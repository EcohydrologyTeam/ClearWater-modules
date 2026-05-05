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
| CRITICAL | 10 (2 resolved 2026-05-04) |
| MAJOR | 18 |
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

```python
self.__simulation_directory = (
    simulation_directory if simulation_directory else "."
)
# later:
store_path=self.__simulation_directory / "model_outputs.zarr",
```

`'.' / 'model_outputs.zarr'` raises `TypeError: unsupported operand type(s) for /: 'str' and 'str'`. Default constructor crashes on any run with `output_variables` non-empty. Masked in tests because v3 tests use `output_variables=[]`.

**Fix:**
```python
self.__simulation_directory = Path(simulation_directory) if simulation_directory else Path(".")
```

### C2. `process.finalize_process` does not exist on Process base class → AttributeError on every chunked run

**Source:** Orchestration F2; `model.py:233-235` and `model.py:490`.
**Inherited from:** v2 verbatim.

`grep -rn "def finalize_process" src/` returns nothing. v2's `Process` base defines `init_process` (no-op default) but not `finalize_process`. v2 chunked path was never exercised in production, hiding the bug. Every v3 chunked run will raise `AttributeError` at the end of `__process_loop_chunked` after the final write.

**Fix (option A, recommended):** add a no-op default `finalize_process` to `clearwater_modules_v2/processes/base.py` mirroring `init_process`. Or fix in v3 by reading `getattr(process, "finalize_process", None)` and calling only if callable.

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

The vectorized guard catches only `denominator == 0.0` (measure-zero set). For pathological inputs where `e_air > P_air` (data-entry error, mis-scaled forcing, sensor noise near saturation), the formula returns a **negative mixing ratio**, which propagates into `density_air` via the `(1+r)/(1+1.61·r)` factor; that factor changes sign near `r = -0.621`, producing **negative or sign-flipped air densities**. Negative density poisons every flux that depends on it (`flux_sensible`, Richardson-stability-dependent `flux_latent_heat`).

**Fix:** Replace `denominator == 0.0` with `denominator <= 0.0` so non-positive denominators all return 0.0. Optionally emit a one-time warning on first encounter.

### C5. Wet-mask writes NaN into read-only forcing inputs, not just outputs

**Source:** Orchestration F4; `model.py:336-357`, `temperature.py:57-69`.

```python
for variable_name in getattr(process, "variables", ()) or ():
    # ... writes NaN into ALL declared variables on dry cells
```

`Process.variables` conflates the variable the process **writes** (`water_temperature`) with variables it **reads** (`wind_speed`, `air_temperature`, `solar_radiation`, `cloudiness`, etc.). After `Temperature.run`, the wet-mask code writes NaN into forcing variables on dry cells. Next substep, those forcings are NaN. v1 explicitly avoided this: NaN-fill state on dry cells, NaN-mask only output slots at write time.

**Mission impact:** in any wet/dry-margin run (Sumwere Creek, any natural channel near baseflow), forcing data is silently corrupted on the margins. After a chunk reload the corruption is partially overwritten but persists between chunks for cells that stay dry.

**Fix:** Either (a) add `Process.output_variables: list[str] = []` and have v3 mask only those, or (b) intersect `process.variables` with the registry's writable-state set and mask only that intersection. (a) is cleaner.

### C6. `__build_process_schedule` uses `start_time.timestamp()` — timezone-dependent for naive datetimes

**Source:** Orchestration F3; `model.py:298-303`.

A naive `datetime` in `.timestamp()` is interpreted in the **local** timezone. The same `datetime(2026,1,1,0,0,0)` produces different POSIX seconds on a Pacific-time laptop vs a UTC cluster. The schedule is baked in at init time, so the bug isn't observable per-step, but the firing schedule **differs by timezone** for any process whose `time_step_seconds` doesn't divide 86400.

For TSM at 5-min substep with model time_step at 5 min, the bug is invisible. For any future process with non-divisor cadence (25-min, hourly with non-zero start-minute), it's a reproducibility defect.

**Fix:** Compute schedule offsets in delta-seconds from `start_time`, not absolute POSIX seconds. The semantic question is whether the user wants "fire when wall-clock UNIX seconds is divisible by interval" (current, TZ-dependent) or "fire every Nth substep starting at start_time" (TZ-independent). The TSM design spec is silent; clarify with author.

### C7. Chunk-end membership test (`current_time in interior_chunk_ends`) is type-mixed and FP-fragile

**Source:** Orchestration F5; `model.py:444-461`.

`chunk_ends` is `pd.DatetimeIndex` of `pd.Timestamp` (always nanosecond precision). `current_time` is `datetime` advanced by `+= self.__time_step` each substep. `Timestamp.__hash__` vs `datetime.__hash__` are not guaranteed symmetric across pandas versions when one is tz-aware; may silently miss boundaries that print as equal. Sub-second `time_step` accumulates floating-point error and misses every boundary.

**Consequence:** if a chunk boundary is missed, the chunk is **never written**, the next chunk's data is **never loaded**, the simulation silently produces wrong output for the rest of the run with no warning.

**Fix:** Compare on **chunk-step-index**, not on time identity. Precompute `interior_chunk_step_indices: set[int]` from `(end - start).total_seconds() / time_step_seconds` for each boundary; test `step_index in interior_chunk_step_indices`. Exact-integer; matches the schedule's own indexing.

### C8. v3 cannot import standalone — vendored streaming-local v2 has broken bare imports

**Source:** Framework F1 (also flagged in Phase 0 gap analysis).
**Location:** `src/clearwater_modules_v2/config/__init__.py:2`, `src/clearwater_modules_v2/config/init.py:1` on the `streaming` branch.

```python
# streaming-local v2/config/init.py:1
from model import Model
```

This imports a top-level `model` module that doesn't exist on `sys.path`. Tracing through, `import clearwater_modules_v3` fails with `ModuleNotFoundError: No module named 'model'` because v3's `config/init.py:42` does `from clearwater_modules_v2.config import init as _v2_init`. The test conftest at `tests/v3/conftest.py:26-40` works around this by pre-caching v2 from the sibling modules-repo's editable install.

**Outside the test suite, the package cannot be imported.** This is a release-blocking defect for v3 1.0.

**Fix:** Bring streaming-local v2 up to `upstream/memory-refactor-pytestUpdate`. Five commits behind, including `dbe0ec7` (mixing_ratio_air fix), `f7b0967` (debug-print toggle-off), `209b67f` (skip-first-step), and the absolute-import refactor.

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

The comment claims first candidate exists because "some import paths fold them via name-mangling on intermediate scopes." This is incorrect. Module-level Python name-mangling does not exist; `__name` is mangled only inside class bodies. The first candidate `_init__init_processes` will never match. The wrapper masks real `AttributeError` if v2 ever renames a helper upstream — exactly the silent-break risk the design spec acknowledges.

**Fix:** Replace with direct `getattr(_v2_init, "__init_processes")`. Add a v2-helper-contract test (`tests/v3/test_v2_helper_contract.py`) that asserts presence and signature of `__init_processes` and `__init_model_data` so upstream changes surface as CI failures.

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

Asymmetry between modes. Once C2 is fixed and processes gain real `finalize_process` bodies, non-chunked runs will silently skip finalization.

### M7. Non-chunked simulation with chunked data source reads only first time-step

**Source:** Orchestration F12; `model.py:201-209`.

`__init_model` reads `[start, start + (chunk_size or time_step)]` from chunked sources. When `chunk_size=None`, this is one substep of data. A user running unchunked Sumwere with HEC-RAS HDF (chunked by design) will load only the first slice and run for 5 days against it.

### M8. `__seed_from_hotstart` falls back to first dim — silently treats space as time

**Source:** Orchestration F10; `model.py:368-383`.

Single-snapshot dataset with only `nface` dim → fallback grabs `nface` as time dim → `isel(nface=-1)` reduces every variable to a single cell.

### M9. `__apply_wet_mask` bare `except Exception` swallows real errors

**Source:** Orchestration F9; `model.py:349-353`.

Catches `KeyError` (intent: "variable not yet in registry"), but also `TypeError`, `AttributeError`, `ValueError`. Silent skip on contract regression.

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

### Phase R-1 (release-blocker, must fix before any production run)

1. **C8 — vendored v2 broken imports.** Bring streaming-local v2 up to `upstream/memory-refactor-pytestUpdate`. Without this, v3 cannot be imported standalone.
2. **C1 — simulation_directory str/Path bug.** One-line fix in `model.py:118-120`.
3. **C2 — finalize_process AttributeError.** Add no-op default to v2 `Process` base class, or use `getattr` in v3 Model.
4. **C7 — chunk-end membership test.** Refactor to use chunk-step-index integer comparison.

### Phase R-2 (correctness, fix before sponsor demos)

5. **C3 — sediment_diffusivity unit mismatch.** Decide on units convention; align docstring + default + formula.
6. **C4 — mixing_ratio_air negative-denominator guard.** Replace `denom == 0` with `denom <= 0`.
7. **C5 — wet-mask read-only-input corruption.** Add `Process.output_variables` and mask only those.
8. **C6 — schedule timezone dependence.** Refactor to delta-seconds-from-start.
9. **C9 — _v2_init_helper.** Replace with direct getattr; add contract test.

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
