# ClearWater Modules v3 — Reviewer Materials for LimnoTech

**To:** Paul Tomasula, Anthony Aufdenkampe, Jason Rutyna, Sarah Jordan
**From:** Todd Steissberg (ERDC-EL)
**Re:** v3 convergence of v1 (`clearwater_modules`) and v2 (`clearwater_modules_v2`)
**Date:** 2026-05-04
**Source branch:** `EcohydrologyTeam/ClearWater-modules-streaming` `streaming`
**PR target:** `EcohydrologyTeam/ClearWater-modules` `memory-refactor-pytestUpdate`
**Read time:** 15–20 minutes

This document is the entry point for reviewing the v3 work. It points at the right design specs, summarizes the defects v3 surfaces in upstream v2, and gives a suggested reading order so each of you can land in the file most relevant to your prior work.

---

## 1. What v3 is, in one paragraph

`clearwater_modules_v3` is a thin-overlay package that combines v1's correctness fixes and orchestration optimizations (latent-heat unit fix, thin-water stability guard, kernel optimization, wet-mask gating, hotstart from `xr.Dataset`) with v2's class-based framework (`Process` composition, YAML configuration via `init_from_file`, per-process substepping, chunking). v2's framework is preserved as the architectural baseline; no new framework is introduced. Where v2 and v1 diverge on a specific quantity, v3 picks the v1 value when v1 is correct and v2 carries an inherited error, the v2 value when v2 added a deliberate improvement (e.g., mixing-ratio guard intent), and a synthesized choice when neither is correct (e.g., dynamic sediment-temperature evolution that both v1 and v2 dropped relative to the canonical Fortran).

The v3 package is at `src/clearwater_modules_v3/` on this branch. It is co-located with v1 and v2 and does not modify either of them outside of bringing the streaming branch's vendored v2 up to your `memory-refactor-pytestUpdate` baseline (commit `009fb95`).

## 2. Scope of this PR

**In scope:**
- v3 TSM (`Temperature` class) — full convergence of v1 and v2 with the merged fixes documented below.
- v3 `Model` orchestration — re-implementation of v2's `Model` with kernel optimization, registry-level wet-mask, hotstart, and resolved chunking.
- v3 `init_from_file` — accepts the v2 YAML schema unchanged, with two new optional top-level keys (`hotstart`, `wet_mask`).
- v3 design specs (`design/clearwater_modules_v3_*.md`).
- v3 test suite at `tests/v3/`: 153 passing tests covering TSM regression, latent-heat correctness, thin-water stability, sediment energy conservation, hotstart roundtrip, Model orchestration, wet-mask scope, and Phase R-3/R-4/R-5 robustness.
- Multi-agent code review log (`design/clearwater_modules_v3_review_findings.md`) — 10 CRITICAL + 18 MAJOR + 19 MINOR findings, with disposition for each.
- Streaming-repo pixi env, reproducible from a clean checkout (`pyproject.toml`, `pixi.lock`).

**Not in scope (deferred to follow-up PRs):**
- v3 NSM1 — in progress on a parallel work track per the architecture spec § 1; will land as a separate PR.
- v3 retirement of v1 — v1 stays in the codebase as a deprecated reference for one release cycle per architecture spec § 5.
- v2 retirement — feature-based trigger ("v3 1.0.0 ships"); v2 enters frozen state (bug fixes only) when this PR merges.
- The Sediment Simulation Module (SSM) work in `src/clearwater_modules_v2/processes/sediment/` — separate active work track.
- The `clearwater_riverine` dry→wet cell concentration artifact — documented as a known issue in `ClearWater-Riverine-streaming/design/known_issues.md` and tracked alongside an unrelated `test_mass_end_plan02` mass-balance regression in that repo.

## 3. Defects v3 found in upstream v2 (worth your direct attention)

These are not v3 contributions per se; they are corrections of inherited transcription or design errors in `memory-refactor-pytestUpdate` v2 that v3 has had to address. v3 carries the corrections; v2 still needs them. They are listed in approximate order of impact.

### 3.1 Sediment defaults disagree with the canonical Fortran TSM

| Parameter | v2 default | Fortran default | v1 default | Origin |
|---|---|---|---|---|
| `sediment_diffusivity` | 0.0061 m²/s (docstring) but formula treats as m²/day | 0.0432 m²/day | 0.0432 m²/day | v2 commit `303d285` (2025-07-28; Paul) — transcription error |
| `sediment_specific_heat` | 1000.0 J/kg/C | 1673.0 J/kg/C | 1673.0 J/kg/C | same v2 commit |

The Fortran source (`HEC-RAS-WQ/RAS-1D-WQ/Kinetics Libraries/Temperature*/Source files/modTemperature.f90`) declares `alphas` in m²/day with default `0.0432` and `Cps` with default `1673.0`. v2's defaults disagree with both Fortran and v1. The `sediment_diffusivity` mismatch is layered: the docstring claims m²/s but the formula's `/86400` only produces W/m² if input is m²/day, so a user who trusts the docstring and supplies a true m²/s value gets a flux **86400× too small** — sediment heat exchange effectively disabled in shallow rivers, which is exactly v3's target use case.

v3 corrects both defaults and the docstring to match Fortran/v1. Recommend the same correction in v2 — `clearwater_modules_v2/processes/temperature.py:79` and the surrounding docstring.

### 3.2 Dynamic sediment-temperature evolution is missing in both v1 and v2

The Fortran TSM evolves `T_sed` each substep:

```fortran
if (use_SedTemp) dTsedCdt = alphas(r) / (0.5 * h2(r) * h2(r)) * (TwaterC - TsedC)
```

This is paired with the sediment heat flux `q_sed` so that the water and sediment heat reservoirs exchange identical enthalpy (energy-conservative). Both v1 and v2 dropped this update — they apply `q_sed` to water but keep `T_sed` constant, breaking energy conservation between the reservoirs.

v3 reinstates the update under a new constructor flag `evolve_sediment_temperature: bool = True` (default True). The `tests/v3/test_tsm_sediment_v3.py::test_water_sediment_energy_conservation_per_substep` test verifies the conservation contract holds to floating-point precision. For shallow-river sponsor scenarios with strong diel forcing, neglecting this update biases the sediment damping term in whichever direction `T_water − T_sed` is sustained.

### 3.3 `mixing_ratio_air` guard is scalar-only

v2's guard handles the zero-denominator case correctly for scalar inputs:

```python
if atmospheric_pressure == atmospheric_vapor_pressure:
    return 0.0
```

But the comparison raises `ValueError: The truth value of an array with more than one element is ambiguous` on multi-cell `xr.DataArray` inputs — i.e., on any non-trivial mesh. v3 replaces this with a vectorized `xr.where` form that handles both scalars and arrays and additionally guards `e_air > P_air` (negative denominator → sign-flipped air density via the `(1+r)/(1+1.61r)` factor → poisons every flux dependent on density_air).

`clearwater_modules_v2/processes/temperature.py:511-528` could absorb the same vectorized form. The v3 implementation is at `clearwater_modules_v3/processes/temperature.py:482-505`.

### 3.4 `process.finalize_process` does not exist on the `Process` base class

`v3.Model.__finalize_model` calls `process.finalize_process(model, registry)` for every process. v2's `Process` base class defines `init_process` (no-op default) but not `finalize_process`. Any production-shape v2 chunked run would crash with `AttributeError` after the final chunk write. The defect is hidden because v2's chunked path was never exercised in production.

v3 uses `getattr(process, "finalize_process", None)` and calls only if callable. Recommend adding a no-op default `finalize_process` to `clearwater_modules_v2/processes/base.py` to match the `init_process` pattern.

### 3.5 `simulation_directory="."` (str) used as a `Path / "..."` operator

v2's `Model.__init__` stores `simulation_directory` as a bare string `"."` and later does `self.__simulation_directory / "model_outputs.zarr"`. This raises `TypeError` on the default-constructed Model whenever `output_variables` is non-empty — i.e., the documented default config. Tests that use `output_variables=[]` short-circuit the output store and mask the bug.

v3 wraps in `pathlib.Path` at construction. `clearwater_modules_v2/model.py:118-120` could absorb the same one-line fix.

### 3.6 Schedule firing depends on local timezone

v2's `__process_loop_full` (and v3's optimized fast-path before our fix) computes `current_time.timestamp() % process.time_step_seconds == 0` to decide which processes fire each substep. `naive datetime.timestamp()` is interpreted in the **local** timezone (POSIX rule), so the same model on Pacific-time and UTC-deployed clusters fires processes on different schedules for any process whose `time_step_seconds` does not divide 86400. For TSM at 5-min substep with model time_step at 5 min, the bug is invisible. For any future process with non-divisor cadence (25-min, hourly with non-zero start-minute), it is a reproducibility defect.

v3 refactors to delta-seconds-from-start, with explicit cadence-multiple validation that raises `ValueError` if `process.time_step_seconds` is not an integer multiple of model `time_step_seconds`. Same change is recommended for v2.

### 3.7 Other defects fixed in v3 that v2 should consider

The full list with citations is in `design/clearwater_modules_v3_review_findings.md` Sections 3–5. Beyond the above, the highlights:

- **Wet-mask scope:** v2's per-process `xr.where(volume > 0, ...)` mask is fine for water_temperature but does not generalize. v3 introduces a registry-level wet-mask that masks only `output_variables` (the variables a process writes), preserving forcing inputs.
- **Chunk-end membership test:** v2's `__process_loop_chunked` would silently miss a chunk boundary on tz-aware or sub-second timesteps because `pd.Timestamp.__hash__` and `datetime.__hash__` are not symmetric. v3 uses integer step-index comparison.
- **`__init_complete` non-idempotency:** v2's `Model.run()` can be called twice silently — the second call re-iterates from `start_time` against an already-advanced registry. v3 raises `RuntimeError` on the second call.
- **YAML-loader error reporting:** v2 wraps a 7-key try block with a generic `Missing key in config: ...` message. v3 uses a `_required(d, *path)` helper that names the full key path on failure.

## 4. Numerical correctness

### 4.1 Verified-correct (14 spot-checked items)

The v3 multi-agent review spot-checked these against literature:

- Brutsaert (1982) saturation-vapor-pressure polynomial coefficients (`__A0` … `__A6`) — 6.108/12.27/23.37/31.67/42.43 mb at 0/10/20/25/30 °C, within 10⁻³ of textbook.
- Latent-heat polynomial `2,499,999 − 2385.74·T_C` — 2.500/2.476/2.452/2.440/2.428/2.261 MJ/kg at 0/10/20/25/30/100 °C, within 10⁻³ of textbook (with the v3 K→C correction).
- Freshwater density polynomial — 999.97 / 998.21 / 995.65 kg/m³ at 4 / 20 / 30 °C, within 10⁻³ of UNESCO/IES80.
- Air-density `0.348` constant — `100/R_d = 100/287.058 = 0.34836`, gives 1.224 kg/m³ at standard conditions vs textbook 1.225.
- `(1 + r)/(1 + 1.61·r)` humidity correction — standard form with `M_w/M_d = 1.61` ratio.
- `0.622` molecular-weight ratio (`M_w/M_d = 18.016/28.96`) — used consistently across `mixing_ratio_air`, `flux_latent_heat`, `density_air_sat`.
- Stefan-Boltzmann constant `5.67e-8` — matches the v1-parity comment in v2 (`5.67037442e-8` rounded; ~7×10⁻⁵ relative truncation).
- Water emissivity `0.97` — standard value (Henderson-Sellers, CE-QUAL-W2).
- `(1 + 0.17·C²)` cloud correction — Bolz (1949) form (the v3 docstring corrects v2's mis-attribution to "Brunt"; the correct attribution is Swinbank 1963 for the T² emissivity polynomial and Bolz 1949 for the cloud correction).
- Richardson piecewise stability function — exponents `+0.80` (unstable) and `−0.80` (stable), regime cutoffs at `±0.01`, bounds `[−1, 2]`. Matches v1's `np.select` formulation; v2's chained `xr.where` is logically equivalent though slower.
- `flux_sediment` `/0.5` — sediment active-layer half-thickness convention; matches v1.
- `flux_net` sign convention — `sensible + solar + sediment + atmospheric + upwelling + latent` is mathematically identical to v1's `(qsens + qsol + qsed + qLW_dn − qLW_up − qlat)·86400·dt` because v3 absorbs the sign into `flux_upwelling_longwave` and `flux_latent_heat` returning negative values. The dt scaling moves into `temperature_change`.
- `__skip_first_time_step` v1-coupling-compat — preserved from v2 with hotstart hooks that correctly default the post-hotstart resume to no-skip.

### 4.2 Review-finding triage

`design/clearwater_modules_v3_review_findings.md` is the punch list. Outcome:

| Severity | Total | Resolved | Deferred |
|---|---|---|---|
| **CRITICAL** | 10 | **10** | 0 |
| **MAJOR** | 18 | **17** | 1 (M4) |
| **MINOR** | 19 | 10 | 9 (deferred or out-of-scope with documented rationale) |

**The one deferred MAJOR (M4):** `clearwater_data.utils.conversions.celsius_to_kelvin` uses `+273.16` (the triple-point of water) where the literature value is `+273.15` (the ice-point). The 0.01 K bias accumulates ~0.05 W/m² on radiation and ~10⁻⁴ K per substep — small but always positive. v2's comment says "for testing consistency with v1" — v1 has the same offset. v3 inherits it intentionally so the v1-parity tests in `tests/v3/test_5_tsm_calculations_v3.py` (15 expected values) do not need re-derivation. Worth correcting in a follow-up release alongside v1-parity-test decommissioning.

## 5. Test evidence

### 5.1 v3 suite (153 tests, 0.24 s)

```
tests/v3/test_5_tsm_calculations_v3.py        15  v1-parity TSM calc regression
tests/v3/test_tsm_latent_heat_v3.py           12  Lv unit fix
tests/v3/test_tsm_stability_ramp_v3.py         8  depth ramp + rate cap (T1-T6)
tests/v3/test_tsm_sediment_v3.py              10  Fortran-parity defaults + dynamic T_sed energy conservation
tests/v3/test_hotstart_roundtrip_v3.py        12  to_hotstart / from_hotstart, Model lifecycle
tests/v3/test_mixing_ratio_air_v3.py           4  C4 vectorized guard
tests/v3/test_wet_mask_scope_v3.py             6  C5 output_variables-only mask
tests/v3/test_v2_helper_contract.py            6  v2 helper signatures pinned
tests/v3/test_model_orchestration_v3.py       13  C1, C2, C7, C6 (TZ-independence + multi-cadence)
tests/v3/test_model_robustness_v3.py          11  M5, M7, M8, M10, M11, M14
tests/v3/test_config_init_robustness_v3.py     9  M12, M13, M15 YAML loader hardening
tests/v3/test_tsm_robustness_v3.py            22  M1, M2, M3 stability validation + NaN guards
tests/v3/test_temperature_minor_v3.py          4  m1, m2, m4, m6, m9, m19
tests/v3/test_model_minor_v3.py               10  M16, m11, m12, m15, m18 minor
tests/v3/test_wet_dry_transition_v3.py         4  Phase 2 vs Phase 3 wet/dry semantics
tests/v3/test_nan_propagation_e2e_v3.py        4  M3 NaN propagation E2E
tests/v3/test_schedule_non_integer_v3.py       6  M18 cadence-multiple validation
                                              ---
                                              153  passed
```

### 5.2 v1 + v2 suite (430 + 2 xfail, 97 s)

The pre-existing v1 / v2 suite continues to pass on this branch (430 passed, 2 xfailed) after the upstream-v2 sync (commit `009fb95`). The 2 xfailed cases are documented as the pre-latent-heat-fix expected-value baseline holders.

### 5.3 Sumwere Creek end-to-end run

Coupled TSM+Riverine demo (`ClearWater-modules/examples/03_Example_Coupled_TSM_and _Riverine.ipynb`) executes against v3 in the streaming-repo pixi env with verified output:

```
init_from_file:           6.91 s
model.run():              89.78 s   (v2 baseline ~89 s)
output shape:             (4321 timesteps × 444 cells)
finite-cell coverage:     40 %  (consistent with floodplain mesh)
water_temperature mean:   15.19 C
water_temperature std:    2.11 C
```

The run produces 5 anomalously cold cells (out of 1.92M values; 0.0003 %) due to a `clearwater_riverine` dry→wet cell concentration artifact at the wet/dry margin. The artifact is documented at the source: `ClearWater-Riverine-streaming/design/known_issues.md` (commit `870b486` in that repo). It is not a v3 defect; v3 TSM correctly processes the value the transport solver hands it.

## 6. Suggested reading order

| Reader | Files of primary interest |
|---|---|
| **Paul** | `clearwater_modules_v3/model.py` (chunking refactor + the four `__process_loop_chunked` TODOs you left, all resolved); `design/clearwater_modules_v3_tsm_design_specification.md` § 3.2; `design/clearwater_modules_v3_review_findings.md` §§ C1, C2, C7 (each was a v2-inherited crash-class defect resolved in v3). |
| **Anthony** | `clearwater_modules_v3/processes/temperature.py:482-505` (vectorized `mixing_ratio_air` guard — extends your scalar guard from `dbe0ec7` to multi-cell DataArrays); same file's sediment-flux block (Fortran-parity defaults from your `February 2026` `use_sediment_temperature` work). |
| **Jason** | `clearwater_modules_v3/processes/temperature.py:739-756` (Richardson — your January 2026 `8218962` and `7f4166a` investigation is cited inline; v3 deletes the `−1` TODO per your conclusion); `design/clearwater_modules_v3_review_findings.md` § N5. |
| **Sarah** | `clearwater_modules_v3/__init__.py` and `clearwater_modules_v3/processes/__init__.py` (overlay strategy); `design/clearwater_modules_v3_architecture_specification.md` § 4 (thin-overlay + integrator-pattern decision); `design/clearwater_modules_v3_review_findings.md` § C9 (`_v2_init_helper` simplification + contract test that pins your v2 helper signatures). |

## 7. How to verify locally

```bash
git clone -b streaming git@github.com:EcohydrologyTeam/ClearWater-modules-streaming.git
cd ClearWater-modules-streaming
pixi install -e dev          # ~5 minutes the first time
pixi run -e dev pytest tests/v3/                         # 153 v3 tests
pixi run -e dev pytest tests/ --ignore=tests/v3 \
    --ignore=tests/sediment \
    --ignore=tests/test_5_tsm_calculations_v2.py         # 430 v1 + v2 tests
```

(`tests/test_5_tsm_calculations_v2.py` is a fixture-only file with a missing `import pytest` from a much earlier commit; the gap analysis flagged it as never-run. Worth deleting in this PR but I left it in case you want to keep its fixture.)

End-to-end coupled-run reproduction recipe is in `design/clearwater_modules_v3_review_findings.md` § "Phase 5 verification log" (the Sumwere Creek scenario with v3-modules → refactor-demo riverine). Wall-clock target ~90 s.

## 8. Open questions for LimnoTech

1. **`Process.variables` scope.** Today it conflates inputs and outputs. v3 added an opt-in `output_variables` class attr on `Temperature` so the wet-mask only masks writes. Is there a cleaner v2 contract you would prefer (e.g., split into `inputs` and `outputs` lists with the union derived)?
2. **YAML schema additions.** v3 adds two optional top-level keys (`hotstart`, `wet_mask`). Is the schema you and Anthony envisioned for hotstart aligned with what v3 ships? (`hotstart: {dataset_path, timestep}` per design spec § 3.3.)
3. **`celsius_to_kelvin` 273.16 vs 273.15.** v3 inherits the 273.16 offset for v1-parity. When do you want to retire the v1-parity tests so we can correct the offset?
4. **`_v2_init_helper` direct getattr.** v3 reuses `__init_processes` and `__init_model_data` from `clearwater_modules_v2.config.init` via direct attribute access plus a CI-pinned signature contract. Is this the intended boundary, or would you prefer v3 fork those helpers entirely?

## 9. Risk register for the merge

| Risk | Probability | Impact | Mitigation |
|---|---|---|---|
| LimnoTech disagrees with the thin-overlay strategy | Low (consistent with your own v2-vs-v1 precedent) | High (re-planning) | This document; v3 is reversible — `src/clearwater_modules_v3/` is an isolated subpackage. |
| LimnoTech wants the v2 corrections (§ 3 above) applied directly to v2 in this PR | Medium | Low | Either path works; v3 has the corrections and v2 can adopt them in a follow-up PR or in this one. Tell me which you prefer. |
| The `clearwater_riverine` dry→wet artifact blocks sponsor scenarios with frequent margin transitions | Low (0.0003 % rate observed) | Medium | Documented in the riverine-streaming repo; tracked alongside an unrelated `test_mass_end_plan02` mass-balance regression for a focused future session. |
| The 273.16/273.15 offset surfaces in a sponsor benchmark before retirement | Low | Low (10⁻⁴ K bias) | Documented; correctable when v1-parity is retired. |

## 10. Files changed (high-level)

```
src/clearwater_modules_v3/                             new (Python source + README)
src/clearwater_modules_v2/                             touched only to sync streaming-local v2
                                                        with upstream/memory-refactor-pytestUpdate
                                                        (commit 009fb95)
src/clearwater_modules/{csm,msm,nsm2}/                 deleted (prototype Python ports; user
                                                        retired in commit 8deca2f, see
                                                        Streaming-branch additions in
                                                        the top-level README for rationale)
design/clearwater_modules_v3_architecture_specification.md     new
design/clearwater_modules_v3_tsm_design_specification.md       new
design/clearwater_modules_v3_nsm1_design_specification.md      new (NSM1 v3 spec; implementation in a separate PR)
design/clearwater_modules_v3_tsm_gap_analysis.md               new
design/clearwater_modules_v3_review_findings.md                new
design/clearwater_modules_v3_limnotech_review.md               this document
design/TSM_NSM1_v1_vs_v2_inventory.md                          new (analysis artifact)
design/legacy_modules_validation_status.md                     new (V&V context)
tests/v3/                                              new (153 tests + conftest sys.path shim)
README.md                                              updated (new "v3: convergence" section)
pyproject.toml + pixi.lock                             streaming-repo pixi env recipe
```

Net diff is large but localized. The merge should be reviewable as additions plus the v2 sync.

---

**End of reviewer materials.** I am happy to walk through any of the v3 design choices or the review-finding remediations on a call. Email or Slack works.
