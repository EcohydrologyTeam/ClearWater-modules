# Architecture note: Riverine ↔ APL coupling and the shared state contract

**Date:** 2026-05-31
**Status:** Reference. Describes how ClearWater-modules (the aquatic process
library, "APL") couples to ClearWater-Riverine (the HEC-RAS-2D transport
engine) and ClearWater-data (the shared registry), records the state-access
contract, and traces the origin of the constituent-naming convention so
reviewers do not have to reconstruct it.
**Scope:** the `Riverine` process and the boundary between transport state and
kinetics. No code change; documents the as-built design and its provenance.

## Summary

Three repositories cooperate, with a single shared state container as the
contract between them:

- **ClearWater-Riverine (CWR)** reads the HEC-RAS-2D HDF output, runs
  constituent transport on the RAS mesh, and holds the transport state. It is a
  generic engine — it transports whatever constituents its config declares.
- **ClearWater-data** defines the **`VariableRegistry`** — the shared,
  time-indexed, per-cell state container.
- **ClearWater-modules (APL)** runs the reaction kinetics (TSM heat budget,
  NSM1 constituents). Every kinetics process reads and writes the
  `VariableRegistry`. A single process, `Riverine`, is the bridge that drives a
  CWR instance and reconciles its constituent names into the shared registry.

The HEC-RAS HDF is read **only** by CWR; the APL never opens it (no `h5py`/`.hdf`
access exists in `src/clearwater_modules_v3/`). CWR and the APL **share one
`VariableRegistry`**: `config/init.py` constructs a single registry and passes
it into CWR (`cwr.ClearwaterRiverine(variable_registry=…)`, accepted at
`clearwater_riverine/model.py:158`). There is one source of truth, not two.

## ClearWater-modules internal stages (v1 / v2 / v3)

Internal names on the `steissberg-v3-modules` branch; they collapse into `main`
when that branch merges.

- **v1** — the direct Python port of the TSM and NSM Fortran code.
- **v2** — LimnoTech's refactor of TSM plus a partial refactor of NSM1.
- **v3** (`clearwater_modules_v3/`) — completes the refactor: the processes not
  yet refactored were finished here, then LimnoTech's completed v2 NSM1 refactor
  was folded in, yielding a complete v3 process set.

## CWR is constituent-agnostic (config-driven)

CWR does **not** hardcode any APL constituent name. A search of CWR source finds
zero occurrences of `Ap`/`NH4`/`NO3`/`TIP`/`DOX`/`algae`/`ammonium`/`nitrate`.
Instead `clearwater_riverine/constituents.py` defines a generic
`Constituent(constituent_name: str, …, constituent_config: dict)` — one per
constituent the **config** declares — whose values come from a pluggable
`provider` (e.g. `{"provider": "float", …}` for a constant, or a provider that
reads a dataset). The constituent's registry **name** and its value **source**
are independent.

CWR's only hardcoded vocabulary is mesh geometry and hydraulics
(`clearwater_riverine/variables.py`: `node_x`, `face_x`, `volume`, `wet_mask`,
`water_surface_elev`, `edge_velocity`, …) — what it reads from HEC-RAS. So the
layering is clean: **the consumer (APL, and any future consumer) defines the
variables to be transported; CWR transports them.** It was not custom-built for
APL variables.

## CWR's two state-access lines

The same transport state is exposed two ways across CWR history; which one the
bridge is written against fixes the access idiom:

| | streaming line (`ClearWater-Riverine-streaming`) | canonical line (`ClearWater-riverine`, its `main`) |
|---|---|---|
| State container | xarray `Dataset` on `model.mesh` | `VariableRegistry` on `model.registry` |
| Read at *t* | `model.mesh["…"].isel(time=t)` | `model.registry.get_at_time("…", t)` |
| Access idiom | keyed or attribute (`mesh.Ap` / `mesh["Ap"]`) | keyed only |

Canonical ClearWater-riverine (`main`) is registry-centric
(`model.py:158`, `self.registry = … VariableRegistry()`); its
`design/willamette_validation_plan.md` documents the streaming-vs-canonical API
difference. The branch the APL couples against, `steissberg-riverine-merged`,
adopts `model.registry` and additionally exposes `model.mesh` as a `MeshView`
(`clearwater_riverine/fork_compat.py`): a small, audited wrapper that presents
the streaming line's `model.mesh["…"]` keyed-read API over the canonical
`VariableRegistry`, so streaming-era consumers (e.g. the Phase-2 ESM
orchestrator) keep working. The `MeshView` returns the registry's own
`DataArray` objects (no copy).

## Constituent naming — origin and the v3 delta

The descriptive registry names and the short→descriptive mapping at the bridge
**originate in v2 (LimnoTech)**, not in v3. On `main`,
`clearwater_modules_v2/processes/riverine.py` reads CWR's short mesh names and
registers them under descriptive names:

```
mesh.Ap  → register("algae_floating")
mesh.NH4 → register("ammonium")
mesh.NO3 → register("nitrate")
mesh.TIP → register("phosphorus_total_inorganic")
mesh.DOX → register("oxygen_dissolved")
```

So the two-name split (CWR/HEC-RAS short names ↔ APL descriptive registry names)
and the bridge that reconciles them are LimnoTech's v2 design. v3 inherited them
by porting v2's processes (commits `a248d16` "port FloatingAlgae and
BenthicAlgae into v3-native processes" and `ba90f1b` "port Riverine into
v3-native processes"); the v3 `_MESH_TO_CANONICAL` map is the direct descendant
of the v2 mapping.

The **only** v3-era rename is `phosphorus_total_inorganic → tip`, introduced in
commit `101626a` "NSM1 v3 Phases 3-6: port 8 missing constituents + integration
cleanup". The v3 algae and phosphorus processes read `tip` with a
`phosphorus_total_inorganic` fallback for back-compat.

## Data flow

```
HEC-RAS-2D HDF
      │  (read by CWR only: clearwater_riverine/io/hdf.py)
      ▼
ClearWater-Riverine  ── transport state written into the SHARED VariableRegistry
      │                   under CWR/HEC-RAS short names (Ap, NH4, …)
      │  APL Riverine.run(): riverine_instance.update(); reconcile names
      ▼
shared VariableRegistry (ClearWater-data)   ◄──► NSM1 / TSM kinetics
      (descriptive names: algae_floating, ammonium, nitrate, tip,
       oxygen_dissolved, …)
```

The `Riverine` bridge (`processes/riverine.py`):
- builds a `cwr.ClearwaterRiverine` with the shared registry
  (`from_file_path(..., variable_registry=…)`, `config_filepath` points CWR at
  the HEC-RAS HDF);
- each substep calls `riverine_instance.update()`, then aliases each declared
  constituent from its CWR short name to its descriptive name via
  `_MESH_TO_CANONICAL` (16 entries). Two-way constituents are aliased with
  `copy(deep=False)` — a **shared buffer**, so a kinetics write propagates back
  into transport with no data duplication; one-way constituents get an isolated
  per-step snapshot. The aliasing re-runs each substep so it stays correct under
  chunked transport;
- resolves coupling depth on the CWR side (RAS "Cell Hydraulic Depth" when
  present, else `volume / wetted_surface_area`, else `WSE − bed`) and registers
  it as `depth`; fail-loud if depth cannot be resolved.

## State-access boundary (the rule)

- **Only the `Riverine` bridge touches `model.mesh`.** It is the single point
  that reaches CWR's transport state.
- **Every other process — TSM and all NSM1 kinetics — reads and writes the
  shared `VariableRegistry` via `get_at_time`/`set_at_time`, never the CWR
  mesh.** `processes/temperature.py` (the heat budget) has zero mesh/CWR
  references; it is fully insulated from the transport-state representation.

This boundary is why the access-idiom and naming details below are confined to
one file.

## Why the bridge uses keyed access, and the registry sharing

- v2's bridge read `riverine_instance.mesh.Ap` — attribute access on a
  `model.mesh` Dataset (the streaming-line API).
- v3's bridge reads `riverine_instance.mesh["Ap"]` — keyed access, because the
  canonical CWR state is a `VariableRegistry` surfaced through the `MeshView`,
  which is keyed. Keyed access is the robust idiom for both Datasets and
  registries; attribute access is the fragile one (it fails on names that
  collide with `Dataset` methods or are not valid identifiers).
- Because CWR is given the **shared** registry, `mesh["Ap"]` and
  `riverine_instance.registry.get_variable("Ap")` return the same object. The
  bridge's use of the `MeshView` is therefore a stylistic accessor choice, not a
  structural dependency.

## Recommendation / clean end-state (joint decision with LimnoTech)

Because CWR transports whatever its config names, two simplifications are
available — both naming-convention decisions to be made jointly with LimnoTech,
not unilateral changes:

1. **Declare CWR constituents under the descriptive names** (`algae_floating`,
   `ammonium`, …) in the coupling config. CWR would then write them into the
   shared registry under descriptive names directly, the kinetics would read
   them directly, and `_MESH_TO_CANONICAL` would be retired — the `Riverine`
   bridge collapsing to pure orchestration (`update()` + depth resolution). When
   a constituent's *values* come from a HEC-RAS WQ dataset, the HEC-RAS name
   lives in the provider spec, not as the registry key, so descriptive naming is
   still achievable.
2. **Bridge accessor:** have the bridge read CWR state via the shared `registry`
   API rather than the `MeshView`, removing its dependency on the streaming-era
   shim.

Neither is required for correctness; both reduce the seam between the two repos.

## Provider-coverage and wet-mask (registry-side metadata)

Two `Model` mechanisms key off each process's declared `variables` /
`output_variables`, independent of CWR:

- **Provider-coverage check** (`model.py` `__init_model` Step 8): every entry of
  `process.variables` must be present in the registry by init, or init raises a
  clear error. This is what the config-driven `init_from_file` path exercises.
- **Wet-mask** (`model.py.__apply_wet_mask`): NaN-masks dry cells for each
  process's `output_variables` (falling back to `variables` when
  `output_variables` is undeclared). Only output states are masked, never input
  forcings.

## References

- ClearWater-Riverine: `clearwater_riverine/io/hdf.py` (HEC-RAS HDF reader),
  `clearwater_riverine/constituents.py` (generic `Constituent`),
  `clearwater_riverine/variables.py` (mesh/hydraulics names),
  `clearwater_riverine/model.py` (`self.registry`, `model.mesh` → `MeshView`),
  `clearwater_riverine/fork_compat.py` (`MeshView`),
  `design/willamette_validation_plan.md` (streaming vs canonical state APIs).
- ClearWater-data: `clearwater_data/variables.py` (`VariableRegistry`,
  `DataArrayVariable`).
- ClearWater-modules: `processes/riverine.py` (the bridge, `_MESH_TO_CANONICAL`,
  `_bridge_mesh_to_registry`), `processes/temperature.py` (registry-only),
  `model.py` (provider-coverage check; wet-mask), `config/init.py` (single
  shared registry), `clearwater_modules_v2/processes/riverine.py` on `main` (the
  v2 origin of the descriptive names + mapping), commits `a248d16` / `ba90f1b`
  (v2→v3 port) and `101626a` (the `phosphorus_total_inorganic → tip` rename),
  `design/clearwater_modules_v3_riverine_process_meshview_compat.md` and
  `design/clearwater_modules_v3_riverine_full_state_bridge.md`.
