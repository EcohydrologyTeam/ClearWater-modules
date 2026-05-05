# ClearWater Modules v3 — Multi-Agent Code Review Findings

**Status:** Review complete; remediation pending.
**Date:** 2026-05-04
**Scope:** v3 module architecture (overlay + scaffold + Phase 3 Model + config) and v3 TSM (Temperature class).
**Reviewers:** four parallel `water-quality-model-source-code-reviewer` agents covering non-overlapping domains.
**Method:** line-by-line review against v1 source, upstream v2 baseline, governing equations, riverine chunking pattern, and the v3 design specifications.

This document consolidates findings across the four reviews and prioritizes them for triage.

---

## 1. Executive summary

The v3 implementation is structurally sound: the overlay package layout is clean, ProcessFactory registration works correctly across import permutations, the v3 `Temperature` class merges v1 fixes into the v2 framework as designed, and the v3 `Model` is a competent re-implementation that adds the four advertised orchestration capabilities (kernel optimization, wet-mask, hotstart, chunking) without subclass fragility.

**That said, the review uncovered nine CRITICAL defects that block production use, plus eighteen MAJOR defects affecting correctness, robustness, or maintainability.** Several of the CRITICAL defects are inherited from upstream v2 unchanged — they are not v3 regressions, but v3 ships them and so v3 owns them. Two of the CRITICAL defects (sediment-diffusivity unit mismatch, mixing_ratio_air negative-denominator) directly threaten temperature accuracy in shallow riverine flows, the very scenario v3 is designed to handle. Two more (`simulation_directory` str/Path bug, `process.finalize_process` missing) are crash-class defects that prevent any default-config production run from completing.

The single most important architectural action before LimnoTech review is **F4-1 (vendored v2 broken imports)**: as of this branch, `import clearwater_modules_v3` fails at the module level outside of the test conftest shim. This was flagged in Phase 0 but has not yet been resolved. Until it is, no v3 user can use the package.

**Numerical correctness of the energy budget:** all 14 spot-checked physics quantities (Brutsaert e_sat polynomial, Lv polynomial, freshwater density, Richardson formula, sign conventions, constants) match literature within ~10⁻³. The two physics CRITICALS are not in the polynomials themselves — they are in parameter handling and edge-case guards.

---

## 2. Findings index

Severity counts (after deduplication of cross-reviewer overlaps):

| Severity | Count |
|---|---|
| CRITICAL | 10 (10 resolved 2026-05-04 — C1, C2, C3, C4, C5, C6, C7, C8, C9, C10) |
| MAJOR | 18 (2 resolved 2026-05-04 — M6, M9) |
| MINOR | 19 |
| Observations | 6 |

Cross-reviewer convergence (issues found by multiple agents):

| Issue | Agents |
|---|---|
| `mixing_ratio_air` negative-denominator guard | Physics (CRITICAL F2), Stability (MINOR F2) |
| `celsius_to_kelvin` uses 273.16 not 273.15 | Physics (MAJOR F5), Framework (MINOR F12) |
| Per-process `volume > 0` guard interacts with Phase 3 wet-mask | Stability (MINOR F4), Framework (MINOR F11) |
| Stale streaming-local v2 (broken imports) | Phase 0 gap analysis (own finding), Framework (CRITICAL F1) |

---

## 3. CRITICAL findings (must fix before any production run)

### C1. `simulation_directory` default is a str, used as a Path → TypeError on default config

**Source:** Orchestration F1; `model.py:118-120` and `model.py:261`.
**Inherited from:** v2 verbatim.
**Status:** **RESOLVED 2026-05-04** — `__simulation_directory` is now wrapped in `pathlib.Path` at construction time so the `/` operator works regardless of whether the caller passed `None`, a `str`, or a `Path`. Tests `test_c1_*` in `tests/v3/test_model_orchestration_v3.py` pin all four input variations.

### C2. `process.finalize_process` does not exist on Process base class → AttributeError on every chunked run

**Source:** Orchestration F2; `model.py:233-235` and `model.py:490`.
**Inherited from:** v2 verbatim.
**Status:** **RESOLVED 2026-05-04** — `__finalize_model` now uses `getattr(process, "finalize_process", None)` and calls only if callable, mirroring the optional-method pattern v3 already uses for `to_hotstart` / `from_hotstart`. Processes that opt in by defining `finalize_process` are honored; others are silently skipped (no AttributeError). Same fix also resolves **M6** by adding the previously-missing `__finalize_model` call to `__process_loop_full` so non-chunked runs are symmetric. Tests `test_c2_*` in `tests/v3/test_model_orchestration_v3.py`. Upstream v2's `Process` base could still benefit from a no-op default for clarity; flagged for LimnoTech.

### C3. `sediment_diffusivity` docstring/formula unit mismatch — up to 86400× error

**Source:** Physics F1; `temperature.py:79`, `temperature.py:92`, `temperature.py:287-300`.
**Status:** **RESOLVED 2026-05-04** — see resolution at end of this finding.

```python
sediment_diffusivity: float = 0.0061,
# ...
# docstring: "Sediment thermal diffusivity (m^2/s)"
# formula divides by 86400, treating it as m^2/day
```

Two layered defects:
1. **Docstring vs. formula:** docstring says m²/s; formula's `/86400` only produces W/m² if input is m²/day. A user who trusts the docstring and supplies a true m²/s value (~5e-7) gets a flux **86400× too small** — sediment heat exchange effectively disabled.
2. **Default value:** v3's 0.0061 is roughly 7× below v1's 0.0432, with no documented rationale.

Mission-relevant for **shallow rivers** where sediment heat storage damps diel swings — exactly v3's target use case.

**Resolution (2026-05-04):** Authoritative TSM source confirmed at `HEC-RAS-WQ/RAS-1D-WQ/Kinetics Libraries/{TemperatureEnergyBudget,TemperatureEquilibrium}/Source files/`. Both Fortran modules declare:

```fortran
! modGlobal.f90:14
real(R8), allocatable, dimension(:) :: alphas    ! sediment thermal diffusivity (m2/d)
! modTemperature.f90:100
alphas = 0.0432
! modTemperature.f90:248 (energy budget) and modTemperature.f90:104 (equilibrium)
q_sediment = pb(r) * Cps(r) * alphas(r) / 0.5 / h2(r) * (TsedC - TwaterC) / 86400.0
```

- **Authoritative units:** `m²/day` (declared in the Fortran type comment).
- **Authoritative default:** `0.0432`.
- **Authoritative formula:** identical to v3; `/86400` converts m²/day inputs to the W/m² output via the `(T_sed-T_w)` driving force.

**Origin of the v3/v2 defect:** the value `0.0061` was introduced by Paul Tomasula on 2025-07-28 in commit `303d285` ("Draft sediment flux implementation") and inherited by v3. `0.0061` does not correspond to any unit conversion of `0.0432` (`0.0432 / 86400 = 5e-7`, not `0.0061`). The change was a transcription error, not a deliberate physics decision. v1 (`clearwater_modules`) matches the Fortran (`alphas = 0.0432`); only v2 / v3 carry the wrong value.

**Resolution applied:** v3 Temperature constructor changed: `sediment_diffusivity: float = 0.0432`, docstring "m²/day", formula unchanged. Upstream v2 should be patched separately as part of the LimnoTech review handoff.

### C4. `mixing_ratio_air` does not guard `e_air > P_air`

**Source:** Physics F2 (CRITICAL); Stability F2 (MINOR — flagged but lower-rated).
**Location:** `temperature.py:516-521`.
**Status:** **RESOLVED 2026-05-04** — both `xr.where` calls now use `denom <= 0.0` so the guard fires for both the exact-zero case (the previous v2 posture) and the pathological `e_air > P_air` case (the C4 fix). Tests in `tests/v3/test_mixing_ratio_air_v3.py` cover normal/equality/negative-denom/vectorized-mixed scenarios.

### C5. Wet-mask writes NaN into read-only forcing inputs, not just outputs

**Source:** Orchestration F4; `model.py:336-357`, `temperature.py:57-69`.
**Status:** **RESOLVED 2026-05-04** — option (a) applied: v3 ``Temperature`` declares `output_variables = ["water_temperature", "sediment_temperature"]` (the variables it writes), and `Model.__apply_wet_mask` honors `output_variables` if defined (with `getattr` fallback to `variables` for backward compat with processes that haven't migrated). After this fix, dry-cell forcings (wind_speed, air_temperature, solar_radiation, cloudiness, atmospheric_pressure, atmospheric_vapor_pressure, wetted_surface_area, volume, sediment_thickness) are preserved across substeps; only outputs are NaN-masked. Tests in `tests/v3/test_wet_mask_scope_v3.py` pin the scope contract, the backward-compat fallback, and the no-op case.

### C6. `__build_process_schedule` uses `start_time.timestamp()` — timezone-dependent for naive datetimes

**Source:** Orchestration F3; `model.py:298-303`.
**Status:** **RESOLVED 2026-05-04** — schedule firing semantic changed to **"every Nth substep starting at start_time"** (TZ-independent). Implementation: `delta_seconds = i * time_step_seconds` and fire when `delta_seconds % process.time_step_seconds == 0`. Mirrors the C7 chunk_size validation: `process.time_step_seconds` must be an integer multiple of `time_step_seconds`; otherwise raises `ValueError`. Tests in `tests/v3/test_model_orchestration_v3.py` pin TZ-independence (UTC vs Pacific), multi-cadence firing, non-aligned start_time, and the validation error.

### C7. Chunk-end membership test (`current_time in interior_chunk_ends`) is type-mixed and FP-fragile

**Source:** Orchestration F5; `model.py:444-461`.
**Status:** **RESOLVED 2026-05-04** — chunk-boundary detection now uses **integer step-index comparison**. `interior_chunk_step_indices: set[int]` is precomputed in `__process_loop_chunked` from `chunk_size_seconds / time_step_seconds` for each boundary; the hot loop tests `step_index in interior_chunk_step_indices`. Exact-integer; timezone-independent; immune to floating-point drift in `current_time +=` arithmetic. Validation enforces `chunk_size` is an integer multiple of `time_step` (raises `ValueError` if not). Tests `test_c7_*` in `tests/v3/test_model_orchestration_v3.py` pin the boundary indices, the post-loop final-chunk write, and the validation.

### C8. v3 cannot import standalone — vendored streaming-local v2 has broken bare imports

**Source:** Framework F1 (also flagged in Phase 0 gap analysis).
**Location:** `src/clearwater_modules_v2/config/__init__.py:2`, `src/clearwater_modules_v2/config/init.py:1` on the `streaming` branch.
**Status:** **RESOLVED 2026-05-04** — non-sediment v2 files synced to `upstream/memory-refactor-pytestUpdate` via surgical `git checkout`:
- `src/clearwater_modules_v2/__init__.py`
- `src/clearwater_modules_v2/config/init.py`
- `src/clearwater_modules_v2/model.py`
- `src/clearwater_modules_v2/processes/base.py`
- `src/clearwater_modules_v2/processes/{temperature,floating_algae,nitrogen,riverine}.py`

The user's active sediment SSM work in `src/clearwater_modules_v2/processes/sediment/` was preserved. After the sync, `PYTHONPATH=src python -c "import clearwater_modules_v3"` succeeds without the test conftest shim. The streaming-local copies of these files are now byte-identical to upstream, so the upstream-discovered fixes (`dbe0ec7` mixing_ratio_air, `f7b0967` debug-print toggle-off, `209b67f` skip-first-step) are now present in streaming too.

**Note:** the sync brings the same upstream defects v3 already worked around (e.g., `mixing_ratio_air` scalar-only guard). v3's own `Temperature` already has the vectorized form; the v2 module retains upstream's scalar form for fidelity to LimnoTech's baseline.

### C10. Dynamic sediment temperature evolution dropped from all Python ports

**Source:** Author research (2026-05-04) confirming that the canonical Fortran TSM evolves sediment temperature dynamically; v1, v2, and v3 all treat it as a static forcing.
**Status:** **RESOLVED 2026-05-04** in v3 — see resolution at end of this finding.

The Fortran `modTemperature.f90` (both energy-budget and equilibrium variants) computes a per-substep update of `TsedC` whenever the sediment-flux gate is enabled:

```fortran
! modTemperature.f90 (energy budget) line 317; (equilibrium) line 113
if (use_SedTemp) dTsedCdt = alphas(r) / (0.5 * h2(r) * h2(r)) * (TwaterC - TsedC)
```

This drives sediment temperature toward water temperature with a relaxation time constant `τ = 0.5 · h²/α`. With `α = 0.0432 m²/day` and `h₂ = 0.1 m`, `τ ≈ 0.116 day ≈ 2.78 hours` — a sensible heat-storage memory for the active sediment layer.

**Why it matters for energy conservation:** the Fortran's water and sediment temperature updates are paired so that the sediment flux into water `q_sed = ρ_b · C_ps · α / (0.5 · h₂) · (T_s − T_w) / 86400` and the sediment temperature change `ΔT_s = α / (0.5 · h₂²) · (T_w − T_s) · dt / 86400` exchange exactly the same enthalpy:

```
ΔE_water (per m²)    = q_sed · dt = ρ_b·C_ps·α/(0.5·h₂)·(T_s−T_w)·dt/86400
ΔE_sediment (per m²) = ρ_b·C_ps·h₂·ΔT_s = ρ_b·C_ps·α/(0.5·h₂)·(T_w−T_s)·dt/86400
                     = −ΔE_water       ✓ energy-conservative
```

**The Python ports (v1, v2, v3) drop the sediment update.** They compute `q_sed` (gated by `use_SedTemp` / `use_sediment_temperature`) and apply it to `T_water`, but `T_sed` stays at its initial / forcing value forever. This **breaks energy conservation between water and sediment**: heat is added to or removed from the water without the sediment's enthalpy reservoir actually changing. Over multi-day shallow-river runs with strong diel forcing, the sediment storage term effectively becomes an infinite reservoir at the initial T_sed value, biasing the sediment heat-damping term in whichever direction T_sed differs from the water column's daily mean.

**Why this matters for the v3 mission:** v3 is targeted at shallow riverine flows where sediment heat storage is one of the largest damping terms on diel temperature swings. Dropping the sediment T update underestimates damping under sustained warm or cold forcing and produces an asymmetric error (the magnitude depends on the persistent T_water − T_sed gradient). For a 5-day Sumwere-type run this is a small but systematic bias; for sponsor-facing climate-change scenarios spanning weeks, it accumulates.

**Resolution applied (2026-05-04):**

1. v3 `Temperature` adds an `evolve_sediment_temperature: bool = True` constructor flag (defaulting to True so v3 matches the Fortran out of the box).
2. v3 `Temperature` adds a `sediment_temperature_change(...)` method computing `ΔT_s = α / (0.5 · h₂²) · (T_w − T_s) · dt / 86400`.
3. v3 `Temperature.run` writes the updated sediment temperature back to the registry each substep when both `use_sediment_temperature` and `evolve_sediment_temperature` are True.
4. The dry-cell wet-mask is applied to the sediment update too: when `volume == 0`, `ΔT_s = 0` (sediment doesn't relax toward water if there's no water in contact).
5. Tests verify the relaxation time constant, energy conservation between the water and sediment heat reservoirs, and backward-compatible behavior when `evolve_sediment_temperature=False`.

This brings v3 into parity with the Fortran reference for sediment heat exchange. v1 should be similarly patched if it remains in service; v2 should be flagged for LimnoTech.

### C9. `_v2_init_helper` candidate-name dispatch is built on a misunderstanding

**Source:** Framework F2; `config/init.py:108-114, 198-216`.
**Status:** **RESOLVED 2026-05-04** — the multi-candidate `_v2_init_helper` was replaced with a single-name `_resolve_v2_helper(name)` that does direct `getattr(_v2_init, name)` and re-raises `AttributeError` with a descriptive message naming the missing attribute, the v2 module path, and the v3 files to update. Contract test at `tests/v3/test_v2_helper_contract.py` (6 tests) pins existence, callability, and exact parameter-name lists for both `__init_processes` and `__init_model_data` via `inspect.signature`. Upstream renames now surface as CI failures.

---

## 4. MAJOR findings (fix before sponsor-facing v3 1.0 release)

### M1. No validation of `q_net_depth_ramp_ref` and `dTdt_max_per_hour`

**Source:** Stability F1; `temperature.py:71-114`.

Negative `q_net_depth_ramp_ref`: silent disable (the `if > 0` branch is False). Zero `dTdt_max_per_hour`: temperature field freezes. Negative `dTdt_max_per_hour`: every cell is forced to a constant value per substep. NaN `dTdt_max_per_hour`: entire field becomes NaN. v1 has at least a defensive `xr.ufuncs.maximum(q_net_depth_ramp_ref, 1e-30)`; v3 has neither validation nor clamp. **Realistic for ESTCP/SERDP sponsors who edit YAML by hand.**

**Fix:** Add `__init__` validation that rejects negative/NaN/zero values appropriately.

### M2. `flux_sediment` does not guard `sediment_thickness ≤ 0`

**Source:** Physics F3; `temperature.py:298`.

Division produces inf/NaN. Inherited from v1. Wet cells with positive volume but missing sediment_thickness data propagate NaN into water_temperature.

### M3. NaN propagation through `richardson_number` is unbounded

**Source:** Physics F4; `temperature.py:585-593`.

`xr.where(NaN > 2.0, 2.0, NaN)` returns NaN (both branches False); NaN survives clamping, propagates through `richardson_function`, `wind_function`, `flux_sensible`, `flux_latent_heat`. `wind_speed = 0` produces `richardson_number = -inf` plus `RuntimeWarning: divide by zero`. v1 suppressed the warning; v3 does not.

### M4. `celsius_to_kelvin` uses 273.16 instead of 273.15

**Source:** Physics F5, Framework F12; `utils/conversions.py:10`.

Triple-point of water vs. ice-point. 0.01 K bias on every Kelvin-temp formula (longwave T⁴, e_sat polynomial, density_air, density_air_sat). ~0.05 W/m² on radiation, ~10⁻⁴ K bias per substep — small but always positive, accumulates over long runs. v1 has same bug; v3 inherits with comment "for testing consistency with v1." When v1-parity tests are decommissioned, fix.

### M5. `init_process` runs before `from_hotstart` — internal-state ordering hazard

**Source:** Orchestration F8; `model.py:213-228`.

Process author who adds new internal state and forgets to add it to both `init_process` and `from_hotstart` will silently de-sync between hotstart-resume and fresh-start. v3 hotstart roundtrip tests don't exercise this through the Model lifecycle.

### M6. `__finalize_model` only called in chunked mode

**Source:** Orchestration F11; `model.py:408-420` (full) vs `model.py:486-490` (chunked).
**Status:** **RESOLVED 2026-05-04** — `__process_loop_full` now calls `__finalize_model()` after the post-loop save, mirroring the chunked path. Same fix as C2's getattr guard. Test `test_c2_finalize_model_invoked_in_full_mode` covers the symmetry.

### M7. Non-chunked simulation with chunked data source reads only first time-step

**Source:** Orchestration F12; `model.py:201-209`.

`__init_model` reads `[start, start + (chunk_size or time_step)]` from chunked sources. When `chunk_size=None`, this is one substep of data. A user running unchunked Sumwere with HEC-RAS HDF (chunked by design) will load only the first slice and run for 5 days against it.

### M8. `__seed_from_hotstart` falls back to first dim — silently treats space as time

**Source:** Orchestration F10; `model.py:368-383`.

Single-snapshot dataset with only `nface` dim → fallback grabs `nface` as time dim → `isel(nface=-1)` reduces every variable to a single cell.

### M9. `__apply_wet_mask` bare `except Exception` swallows real errors

**Source:** Orchestration F9; `model.py:349-353`.
**Status:** **RESOLVED 2026-05-04** — narrowed to `except KeyError` (the documented intent: "variable not yet in registry"). Resolved alongside C5 in the same patch since both touch `__apply_wet_mask`. Other exception types now propagate, surfacing real bugs.

### M10. `__init_complete` allows `run()` to be called twice silently

**Source:** Orchestration F6; `model.py:131, 145-152`.

Second `run()` re-iterates from `start_time`, advancing the already-advanced registry. Notebook footgun.

### M11. Wet-mask threshold uses `>` strict-inequality with default 0.0

**Source:** Orchestration F13; `model.py:333-334`.

Cells with `1e-300` depth are "wet" by threshold and get processed, producing NaN in the next substep that the wet-mask cannot clean up.

### M12. `_resolve_hotstart` confusing failure on unsupported file suffix

**Source:** Framework F4; `config/init.py:170-174`.

`.h5`/`.netcdf4`/`.cdf` route to `xr.open_dataset` which may need engine packages not installed. Error has no breadcrumb to the YAML key.

### M13. `_resolve_hotstart` does not validate `timestep` parseability

**Source:** Framework F5; `config/init.py:175-176`.

YAML typo `timestep: '2022-13-01 12:00:00'` (invalid month) surfaces deep inside xarray, not at the YAML resolver.

### M14. `_resolve_wet_mask` does not validate variable exists in registry

**Source:** Framework F6; `config/init.py:179-195`.

Typo `wetted_surface_aera` not detected until first substep. Error message won't mention the YAML key.

### M15. `KeyError`-to-`ValueError` wrap loses deep key path

**Source:** Framework F7; `config/init.py:57-83`.

The `try` block spans seven distinct YAML lookups + iteration over `processes` and `data_sources`. `Missing key in config: 'start_dattetime'` doesn't say which top-level block or which iteration index.

### M16. Process integrator-pattern contract documented but not enforced

**Source:** Framework F3; `processes/base.py:1-19`.

Five-step contract documented; no test or runtime check enforces it. v3 `Temperature.run` itself doesn't fully follow the pattern (computes delta, not rate; multiplies by dt inside helper).

### M17. README Status section is outdated

**Source:** Framework F8; `src/clearwater_modules_v3/README.md:21-37`.

Status says "Phase 1 (scaffold)" but Phases 2-4 are complete. Migration table covers only TSM (3 rows); should include all overlay processes and utility modules.

### M18. Schedule uses bare-`==` mod-op on floats

**Source:** Orchestration (re-evaluation of F3 logic); `model.py:301-303`.

`current_seconds % interval == 0.0` is exact for integer-second intervals, but a user passing `time_step=timedelta(seconds=0.1)` would see floating-point drift. No test for non-integer-second intervals.

---

## 5. MINOR findings (track for v3.1+)

| # | Source | Item |
|---|---|---|
| m1 | Physics F6 | `mixing_ratio_air` nested `xr.where` is correct but more complex than necessary |
| m2 | Physics F7 | `flux_atmospheric_longwave` cloud correction wrongly attributed to "Brunt/Kiehl" — should be "Swinbank/Bolz" |
| m3 | Physics F8 | `wind_function` units opaque; `1_000_000` divisor undocumented |
| m4 | Physics F9 | `water_specific_heat` not NaN-safe (silently returns 4178.0) |
| m5 | Physics F10 | Constructor `wind_a/b/c` have no defaults — convenience class method missing |
| m6 | Physics F11 | `flux_sediment` returns Python scalar `0.0` (not DataArray) when disabled |
| m7 | Stability F3 | Inner `xr.where` in `mixing_ratio_air` is correct but more complex than necessary (overlaps m1) |
| m8 | Stability F6 | `q_net_depth_ramp_ref` is scalar-only; v1 supports per-cell |
| m9 | Stability F7 | Depth clamp `xr.where(depth > 0, ...)` is tautology after surface_area guard |
| m10 | Stability F8 | Floating-point ordering differs from v1 (~10⁻¹⁵ relative) |
| m11 | Orchestration F7 | `__step_index` helper is dead code |
| m12 | Orchestration F14 | Empty `Path("")` is falsy and gets replaced with `"."` |
| m13 | Orchestration F15 | Hardcoded `"water_temperature" → "nface"` v2 hack remains |
| m14 | Orchestration F16 | Schedule precomputed `n_steps + 1`; loop runs `n_steps` |
| m15 | Orchestration F17 | `processes` type hint says `tuple[Process, ...]` but constructor coerces from any iterable |
| m16 | Orchestration F19 | `__output_data_store.write_chunk(data=variable)` passes a `Variable`, not a `DataArray` |
| m17 | Framework F9 | `_v2_init_helper` candidate ordering reversed (overlaps C9) |
| m18 | Framework F10 | `RUN_ORDER` is exported but unused |
| m19 | Framework F11 | Per-process `xr.where(volume > 0, ...)` duplicates `Model.__apply_wet_mask` |

---

## 6. Observations

| # | Source | Item |
|---|---|---|
| O1 | Orchestration F18 | Wet-mask compute called per substep when any process fires; profiler note |
| O2 | Orchestration F20 | `__restore_process_hotstart` shares `attrs` across processes; key collisions by convention only |
| O3 | Orchestration F21 | `LOGGER.info("Running timestep: %s", current_time)` is per-substep — should be DEBUG for production |
| O4 | Orchestration F22 | v2 chunked-loop double-write bug genuinely fixed in v3; new bugs (C2, C7) are different from old |
| O5 | Framework O1 | Process registration ordering correct and stable across import permutations |
| O6 | Framework O3 | Backward compatibility verified: v2 YAML configs without new keys run unchanged on v3 |

---

## 7. Recommended action plan

### Phase R-1 (release-blocker, must fix before any production run) — **COMPLETE 2026-05-04**

1. ✓ **C8 — vendored v2 broken imports.** Surgical `git checkout` of non-sediment v2 files from upstream. v3 imports standalone.
2. ✓ **C1 — simulation_directory str/Path bug.** Wrapped in `Path(...)` at construction.
3. ✓ **C2 — finalize_process AttributeError.** `getattr` callable guard. M6 also resolved (full-mode finalize call added).
4. ✓ **C7 — chunk-end membership test.** Refactored to integer step-index comparison; `chunk_size` validated as integer multiple of `time_step`.

### Phase R-2 (correctness, fix before sponsor demos) — **COMPLETE 2026-05-04**

5. ✓ **C3 — sediment_diffusivity unit mismatch.** Aligned to Fortran (0.0432 m²/day).
6. ✓ **C4 — mixing_ratio_air negative-denominator guard.** `denom <= 0` form.
7. ✓ **C5 — wet-mask read-only-input corruption.** `output_variables` opt-in; M9 also resolved.
8. ✓ **C6 — schedule timezone dependence.** Refactored to delta-seconds-from-start; cadence-multiple validation.
9. ✓ **C9 — _v2_init_helper.** Direct `getattr` + 6-test contract test.

### Phase R-3 (robustness, fix before v3 1.0 ship)

10. M1, M2, M3 — parameter validation, sediment_thickness guard, NaN propagation through Richardson
11. M5–M11 — orchestration robustness items
12. M12–M15 — YAML resolver error reporting
13. M17 — README refresh

### Phase R-4 (deferred, track for v3.x)

14. M4 (273.16 → 273.15) — hold until v1-parity tests are decommissioned
15. M16 — integrator-pattern contract: enforce or remove from docstring
16. M18 — non-integer-second timestep test
17. All MINOR items

### Phase R-5 (test coverage, in parallel with R-1 through R-3)

18. Unit-import test (`import clearwater_modules_v3` standalone)
19. v2-helper-contract test
20. Chunked end-to-end test (catches C1, C2, C7)
21. Wet-mask forcing-preservation test (catches C5)
22. Timezone reproducibility test (catches C6)
23. Wet/dry transition regression test (Phase 2 vs Phase 3 semantics)
24. NaN propagation regression test
25. MMS test for energy conservation

---

## 8. Open questions for author input

1. **Sediment diffusivity (C3):** intended units? v1 says m²/day with default 0.0432; v3 docstring says m²/s with default 0.0061. Which is right?
2. **Schedule firing semantic (C6):** "absolute UNIX seconds modulo interval" (current, TZ-dependent), "delta-seconds-from-start modulo interval" (TZ-independent), or "every Nth substep starting at start_time"?
3. **Process variable scope (C5):** does `Process.variables` represent inputs only, outputs only, or both? Today it's read as "everything declared," but the wet-mask needs "outputs only."
4. **`init_process` / `from_hotstart` ordering (M5):** intentional, or should they be merged?
5. **`__finalize_model` mode asymmetry (M6):** should non-chunked runs also call it?
6. **Per-process volume guard (m19):** remove now (Phase 4 cleanup) or keep until orchestration mask is widely deployed?
7. **273.16 vs 273.15 (M4):** correct now and break v1-parity tests, or hold for v3.x?
8. **Process integrator-pattern contract (M16):** enforce via abstract `compute_rate`, or demote to "guideline"?

---

**End of consolidated findings.** Each finding has a verified file:line reference and a recommended fix. The four full reviewer reports are preserved in the conversation transcript that produced this synthesis.
