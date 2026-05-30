# Design Spec: Riverine process — MeshView compatibility, chunk-safe bridge, correct depth, provider check

**Date:** 2026-05-30 (final design)
**Components:**
- `src/clearwater_modules_v3/processes/riverine.py`, `src/clearwater_modules_v3/model.py` (ClearWater-modules)
- `clearwater_riverine/model.py`, `clearwater_riverine/io/hdf.py`, `clearwater_riverine/utilities.py`, `clearwater_riverine/variables.py` (ClearWater-Riverine, branch `steissberg-riverine-merged`)

**Severity:** Blocking for the `init_from_file` (config-driven) coupling path — the LimnoTech-facing production path.
**Status:** Final. Two-repo change.

## Scope boundary (explicit)

This change stays **inside LimnoTech's coupling design**: shared `VariableRegistry`, `init_from_file`, the `Riverine` process as the coupling point, per-process substepping, chunking. It makes that design **correct and robust**; it does **not** re-architect it. Specifically **out of scope** and deliberately *not* done here:

- No data-source/two-way split, no provider abstraction, no restructuring of how forcings flow (the earlier "(b)/(c)" options).
- No BMI convergence. The BMI interface (`papers/paper2_nutrients/bmi_implementation/`) is a separate track for the Report 2 appendix and stays in its lane.
- No full-NSM-I-state bridge (a separate, deferred design item).

Every item below either fixes a real defect in LimnoTech's design or hardens it against silent failure.

## Defects addressed

1. **Attribute access vs. MeshView.** `init_process` reads `mesh.Ap`; the mesh is now a `MeshView` exposing constituents only by item access. → item access.
2. **Stale-after-chunk-1 bridge.** Constituents bridged once via `copy(deep=False)` aliases; chunked mode re-registers fresh objects per chunk, stranding the aliases. → chunk-safe re-bridge each substep (mirrors the validated manual path, which re-pushes state every step).
3. **Inorganic-P name mismatch.** Bridges `TIP` as `phosphorus_total_inorganic`; v3 Phosphorus reads `tip`. → register as `tip` (algae consumers already prefer it).
4. **Wrong depth.** Registers `depth` as a static alias of `wetted_surface_area` (area, not length). → correct depth, resolved by precedence on the CWR side (below).
5. **Silent missing-input failure.** A `process.variables` entry that no provider supplies surfaces as a runtime `KeyError` mid-substep. → one-shot init-time coverage check (Change C).

## Correct depth: precedence + on-demand computation

**Definition (settled).** The depth NSM needs is the cell mean water-column depth. CWR resolves it by **precedence**, preferring RAS's own output:

1. **RAS "Cell Hydraulic Depth"** (`FACE_HYD_DEPTH`, read from the HDF) when present — authoritative.
2. else **`volume / wetted_surface_area`** (the existing `calculate_average_depth`) — matches the validated manual path.
3. else **`WSE − bed_elevation`** (`maximum_depth`/`calculate_face_hyd_depth`) — last resort, emit a warning (this is *max* depth and overestimates the mean).

This is parity-preserving for the validated Santiam–Salem case (its HDF lacks Cell Hydraulic Depth → falls through to `volume/wsa`, identical to today's manual path).

**On-demand, not default-on.** Depth is computed **only when the coupling asks for it**; standalone transport runs are unaffected. No `average_depth` is forced on for every run.

**Provenance.** `FACE_HYD_DEPTH` currently holds *either* the RAS-read value *or* a synthesized `WSE−bed` (and `__populate_wet_mask` can overwrite it). Presence alone cannot drive precedence. Record RAS Cell Hydraulic Depth availability **at HDF-read time** (a boolean on the model, set where `io/hdf.py` probes `Cell Hydraulic Depth` via `_optional_temporal_variables` + `hdf_path in infile`) and key the resolver off that flag, not off registry presence.

## Change A — ClearWater-Riverine (`steissberg-riverine-merged`)

1. **Revert the default-on `average_depth`** added in the prior iteration. Standalone behavior is unchanged: nothing new is computed unless the coupling enables it. (Keep the `__calculated_variables` None-normalization only if needed to avoid the pre-existing `.items()`-on-None crash; do not force any variable on.)
2. **Record RAS-hydraulic-depth availability.** In the HDF reader, set e.g. `self.__ras_cell_hydraulic_depth_available: bool` when the `Cell Hydraulic Depth` dataset is present in the file. Expose it read-only if needed for tests.
3. **Precedence resolver.** Add a model method (it needs the availability flag, so it cannot use the registry-only `CALCULATED_VARIABLE_MAP` signature) that returns the resolved depth DataArray for the currently loaded window:
   - flag set → the RAS-read `FACE_HYD_DEPTH`;
   - else wsa derivable → `calculate_average_depth(registry)` (volume/wsa);
   - else → `WSE − bed` with `warnings.warn(...)`.
   Guard the wsa branch so a config lacking elevation–volume lookups falls through to the WSE−bed branch rather than raising.
4. **On-demand enablement + per-chunk refresh.** Add public `enable_coupling_depth()`: idempotent; sets an enabled flag, computes the resolved depth for the current window now, and registers it under a new constant **`COUPLING_DEPTH = 'coupling_depth'`** (in `variables.py`). When the flag is set, recompute/re-register `COUPLING_DEPTH` in the per-chunk refresh path (alongside `__update_calculated_variables`/`__load_new_chunk`) so it tracks each chunk's time window. When the flag is unset (standalone), it is never computed or registered.
5. **Public `is_chunked` accessor** (keep from the prior iteration) — used by coupling tests/diagnostics.

## Change B — ClearWater-modules (`processes/riverine.py`)

Adjust the already-implemented chunk-safe bridge so depth comes from the resolved coupling depth:

- In `init_process`, after the existing `water_temperature`/`volume`/`wetted_surface_area` checks, call `self.riverine_instance.enable_coupling_depth()` once (turns depth computation on for this coupled run), then seed via `self._bridge_mesh_to_registry(registry)`.
- In `_bridge_mesh_to_registry`, change the depth source from `mesh["average_depth"]` to **`mesh["coupling_depth"]`**, keeping the fail-loud `KeyError` when it is absent (message: ensure the riverine model is coupling-enabled / depth could not be resolved). Keep `overwrite=True`.
- Constituent bridging, the `_MESH_TO_CANONICAL` map (`TIP → tip`), the per-substep re-bridge in `run()`, and the present-constituent gate are unchanged from the implemented Change B.

**CWR↔modules depth contract (both repos build to this):**
`riverine_instance.enable_coupling_depth()` enables/seed-computes the resolved depth and registers it under `'coupling_depth'`, refreshed per chunk; the modules bridge reads `mesh["coupling_depth"]` and registers it as `depth`.

## Change C — ClearWater-modules (`model.py`): one-shot provider check

Turn the latent runtime `KeyError` into a clear init-time error. In `Model.__init_model`, **after** data sources are loaded (step 1), all `init_process` calls have run (step 3), and hotstart restore (if any) has completed — i.e., once every provider has had its chance to populate the registry — assert coverage:

```
required = {v for p in self.__processes for v in p.variables}
missing  = [v for v in sorted(required) if v not in self.__registry]
if missing:
    raise <clear error naming each missing variable and the process(es) that declare it,
           noting it must be supplied by a data source, the riverine bridge, hotstart,
           or pre-registration>
```

This is provider-agnostic (it checks actual registry state, so it covers data sources, the bridge, hotstart, and defaults uniformly), mirrors the existing `wet_mask_variable` init-time check (`model.py:393`), and adds no provider abstraction. Place it next to that check. Known limitation to note in the message/comment: constituents read via `registry.get_at_time(...)` that are **not** declared in any `process.variables` (e.g. `algae_floating`) are not covered — only declared inputs are. `depth` *is* declared, so the depth wiring is covered.

## Process order (pre-existing caveat, not changed here)

The schedule fires processes in registered (YAML) order; the `upstream_processes` DAG is only validated, and nothing declares `Riverine` upstream. All bridged state is correct only if `riverine` is listed first in `processes:`. Tracked follow-up: add `"Riverine"` to consumers' `upstream_processes`.

## Verification

**ClearWater-Riverine:**
1. Precedence — RAS Cell Hydraulic Depth present → `coupling_depth` equals the RAS-read values; absent but lookups present → equals `volume/wsa`; neither → equals `WSE−bed` and warns.
2. On-demand — without `enable_coupling_depth()`, `coupling_depth` is absent and standalone construction/run is byte-unchanged from baseline (no new computation).
3. Per-chunk — after `enable_coupling_depth()`, a 2-chunk run has correct `coupling_depth` at the first and last timestep of each chunk.
4. `is_chunked` True/False.

**ClearWater-modules:**
1. Bridge unit (stub `riverine_instance` with a `MeshView`, an `enable_coupling_depth()` method, `is_chunked=False`, and a `coupling_depth` mesh entry): `init_process` calls `enable_coupling_depth()`, registers canonical aliases incl. `depth` from `coupling_depth`; `tip` present, `phosphorus_total_inorganic` absent; subset omitting `Ap` bridges only present constituents; omitting `coupling_depth` → `KeyError`.
2. Chunk-safety unit (already present): swap `Ap`/`coupling_depth` for new chunk-2 objects, re-bridge, assert `algae_floating`/`depth` track the new values.
3. **Provider check (new):** a Model whose processes declare an input that no provider supplies raises a clear `__init_model` error naming the variable; a fully-provided Model initializes cleanly; the riverine-coupled case (depth supplied by the bridge) passes the check.
4. No regressions in `tests/v3`.

## Cross-repo coordination

Land **Change A first** so `coupling_depth` and `enable_coupling_depth()` exist before Change B calls them. Changes B and C are same-repo. The Report 2 nutrient re-run proceeds on the manual `build_v3_modules` coupling and depends on none of this.
