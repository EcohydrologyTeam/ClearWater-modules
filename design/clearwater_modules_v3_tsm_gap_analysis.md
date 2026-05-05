# ClearWater Modules v3 — TSM Gap Analysis (Phase 0 deliverable)

**Status:** Phase 0 deliverable; amended 2026-05-04 (Phase 2) with the N4 ``mixing_ratio_air`` array-guard finding.
**Author:** Todd Steissberg (ERDC), with Claude
**Date:** 2026-05-04 (initial); 2026-05-04 (Phase 2 N4 update)
**Scope:** Row-by-row diff of v1 TSM (with fixes on the streaming branch) vs. v2 TSM (on `upstream/memory-refactor-pytestUpdate`), with each row classified as the disposition for v3.

This file is the work-tracking artifact for v3 TSM Phases 1–5. Companion documents:

- `clearwater_modules_v3_architecture_specification.md` (umbrella architecture)
- `clearwater_modules_v3_tsm_design_specification.md` (TSM-specific design)
- `TSM_NSM1_v1_vs_v2_inventory.md` (longer-form inventory the table below condenses and corrects)

---

## Sources of truth

| Track | Repo / branch | Path |
|---|---|---|
| v1 TSM (with fixes) | `EcohydrologyTeam/ClearWater-modules-streaming` `streaming` | `src/clearwater_modules/tsm/` and `src/clearwater_modules/base.py` |
| v2 TSM (baseline) | `EcohydrologyTeam/ClearWater-modules` `memory-refactor-pytestUpdate` (= `upstream/memory-refactor-pytestUpdate`) | `src/clearwater_modules_v2/processes/temperature.py`, `src/clearwater_modules_v2/model.py`, `src/clearwater_modules_v2/config/init.py` |

**Important caveat (resolved during Phase 0):** the `src/clearwater_modules_v2/` tree on the `streaming` branch is **older** than `upstream/memory-refactor-pytestUpdate`. The streaming-local copy lacks five upstream commits, including `dbe0ec7` (mixing_ratio_air bug fix), `f7b0967` (debug-print toggle-off), and `209b67f` (skip-first-step logic). All v2 references in this gap analysis use the **upstream baseline**, not the streaming-local copy. Before Phase 1, the streaming branch should pull or rebase those commits, or v3 should be branched off `upstream/memory-refactor-pytestUpdate` directly per the architecture spec.

## Disposition legend

| Tag | Meaning |
|---|---|
| **Port v1→v3** | v1 has a fix or feature v2 lacks; v3 takes the v1 logic. |
| **Keep v2 in v3** | v2 has correct/desirable framework or behavior v1 lacks; v3 keeps v2. |
| **Reconcile in v3** | v1 and v2 diverge on convention or implementation; v3 needs an explicit synthesis. |
| **Resolve TODO in v3** | v2 has an unresolved `# TODO`; v3 replaces it with a documented decision. |
| **Resolved upstream** | v2 upstream already addresses the issue; v3 inherits it (no v3-specific work). |

---

## 1. Numerical correctness

| # | Topic | v1 TSM (`clearwater_modules/tsm/...`) | v2 TSM (`upstream/memory-refactor-pytestUpdate`) | v3 disposition |
|---|---|---|---|---|
| N1 | **Latent heat of vaporization** | `mf_latent_heat_vaporization` (`processes.py:213-228`): converts Kelvin → Celsius before applying polynomial `2.5e6 − 2385.74·T_C`. Yields ~2.45 MJ/kg @ 20 °C (correct). Fixed in commit `d9505c6`. | `latent_heat_vaporization` (`temperature.py:455-463`): applies polynomial to Kelvin via `conversions.celsius_to_kelvin(...)`. Yields ~1.80 MJ/kg @ 20 °C (~26 % low → biases evaporative cooling low → simulated water temp warm). | **Port v1→v3.** Replace v2 implementation with v1 logic. Add regression tests from `test_tsm_latent_heat.py`. |
| N2 | **Thin-water stability: depth ramp on `q_net`** | `dTdt_water_c` (`processes.py:456-523`) multiplies `q_net` by `min(1, depth/q_net_depth_ramp_ref)`; default `q_net_depth_ramp_ref = 0.3 m`; set to `0.0` to disable. Constant in `tsm/constants.py::Temperature`. | Not implemented. `temperature_change` (`temperature.py:390-435`) divides directly by `volume`; small volumes → unbounded `dT`. | **Port v1→v3.** Apply ramp inside v3 `temperature_change`. Preserve disable-with-default-of-zero semantics. Port `test_tsm_stability_ramp.py` cases. |
| N3 | **Thin-water stability: rate cap on dT/dt** | `dTdt_water_c` clips `dTdt` to `±dTdt_max_per_hour · dt_hours`; default `5.0 K/hr`; set to `+inf` to disable. | Not implemented. | **Port v1→v3.** Same as N2: insert into v3 `temperature_change`, default-disabled (`+inf`) behavior preserves v2 numerics when not opted in; default-enabled (`5.0 K/hr`) matches v1. |
| N4 | **`mixing_ratio_air` NaN/edge guard** | `mixing_ratio_air` (`processes.py:32-45`): no guard; assumes `pressure_mb > eair_mb`. | `mixing_ratio_air` (`temperature.py:511-528`): `if atmospheric_pressure == atmospheric_vapor_pressure: return 0.0`. Fixed upstream in commit `dbe0ec7`. (Streaming-local copy still has the inverted-tautology bug `eair == eair`; **upstream is correct for scalar inputs**.) | **Reconcile in v3 (Phase 2 deviation, 2026-05-04).** Upstream's scalar guard works for one-cell tests but raises `ValueError: The truth value of an array with more than one element is ambiguous` when ``atmospheric_pressure`` and ``atmospheric_vapor_pressure`` are multi-cell ``xr.DataArray`` inputs (Sumwere Creek and any non-trivial mesh). v3 replaces the scalar ``if`` with a vectorized ``xr.where`` form that produces ``0.0`` for cells where ``pressure == vapor_pressure`` and the standard mixing ratio otherwise. v3 ``mixing_ratio_air`` therefore satisfies the same intent as the upstream fix but works for both scalars and arrays. Verified during Phase 2 smoke testing on a 3-cell DataArray; existing v2 Richardson unit tests still pass against v3 ``Temperature``. |
| N5 | **Richardson `−1` factor** | `ri_number` (`processes.py:150-168`): no leading `−1`. Formula: `g · (ρ_air − ρ_air_sat) · 2 / (ρ_air · u²)`. | `richardson_number` (`temperature.py:579-624`): `−1` is **commented out** with TODO at `temperature.py:607`: `# -1 #TODO: check original equation to see if this multiplication by negative one is needed (not in v1 of code)`. Live formula matches v1. | **Resolve TODO in v3.** Per spec §3.1 and §8 (resolved 2026-05-04 via Jason Rutyna's January 2026 investigation, commits `8218962` and `7f4166a`): delete the `# -1` line and the TODO. v3 formula = v1 formula (no leading `−1`). |
| N6 | **`flux_sediment / 0.5` factor** | `q_sediment` (`processes.py:387-422`): uses `(ρ_b·c_ps·α_s) / (0.5·h_2)` — same `/0.5` factor — gated by `use_sed_temp` flag. Comment explains 86400 conversion only. | `flux_sediment` (`temperature.py:281-300`): same `/ 0.5` factor with `# TODO: determine why we need this 0.5` at line 295. Gated by `self.use_sediment_temperature`. | **Resolve TODO in v3.** Per spec §3.1 (resolved 2026-05-04): replace TODO with one-line docstring documenting `0.5` as the sediment active-layer half-thickness convention. Do not change the formula; matches v1 legacy. |
| N7 | **Energy-balance integration form** | `q_net` (`processes.py:425-453`) returns `(W/m²) · 86400 · dt_days`; `dTdt_water_c` then divides by `V·ρ·c_p`. The 86400·dt factor lives in `q_net`. | `flux_net` returns sum of fluxes in W/m² (no time scaling); `temperature_change` (`temperature.py:390-435`) multiplies `flux_net · surface_area · time_step_seconds / (V·ρ·c_p)`. The dt factor lives in `temperature_change`. | **Reconcile in v3 (Keep v2 framework).** v3 `flux_net` returns W/m² (v2 convention; cleaner unit accounting); the dt factor lives in `temperature_change`, which is also where the depth ramp (N2) and rate cap (N3) are applied. Document the convention shift relative to v1 with a one-line comment so anyone porting v1 tests sees it. |
| N8 | **`xr.where(volume > 0, ...)` post-mask** | Wet-mask gating happens at the orchestration layer (`base.py::increment_timestep(wet_mask=...)`); v1 process bodies do not carry per-cell dry-cell guards. | v2 `Temperature.run` masks via `xr.where(volume > 0, delta_water_temperature, 0)` at `temperature.py:137`. | **Reconcile in v3.** Once v3 Model has registry-level wet-mask (Phase 3), this per-process mask is redundant; remove from v3 `Temperature.run`. v3 architecture spec §4 already specifies this pattern. |

## 2. Process inventory (function-by-function correspondence)

All entries are in v3 unless noted. v3 source location: `src/clearwater_modules_v3/processes/temperature.py` (full file, not overlay).

| # | v1 function | v2 method | v3 disposition |
|---|---|---|---|
| P1 | `air_temp_k`, `water_temp_k` (dynamic vars) | inline `conversions.celsius_to_kelvin(...)` calls | **Keep v2 in v3.** |
| P2 | `mixing_ratio_air` | `mixing_ratio_air` | **Keep v2 in v3** (incl. upstream NaN/edge guard, see N4). |
| P3 | `density_air` | `density_air` | **Keep v2 in v3.** |
| P4 | `mf_density_water` | `water_density` | **Keep v2 in v3** (algebraically identical UNESCO-style anomaly). |
| P5 | `mf_esat_mb` | `saturation_vapor_pressure` | **Keep v2 in v3** (identical Brutsaert 1982 polynomial). |
| P6 | `mf_density_air_sat` | `density_air_sat` | **Keep v2 in v3.** |
| P7 | `ri_number` | `richardson_number` (returns `(rn, rfn)`) | **Keep v2 in v3** + remove `−1` TODO (N5). v2's combined `(rn, rfn)` return is fine. |
| P8 | `ri_function` | inline in `richardson_number` (`xr.where` chain) | **Keep v2 in v3.** Note v2 uses chained `xr.where` (less efficient than v1's `np.select`); not a correctness issue and outside Phase 2 scope. |
| P9 | `mf_latent_heat_vaporization` | `latent_heat_vaporization` | **Port v1→v3** (N1). |
| P10 | `mf_cp_water` | `water_specific_heat` | **Keep v2 in v3** (identical 5-bin lookup). |
| P11 | `emissivity_air` (named dynamic var) | inline at `flux_atmospheric_longwave` (`temperature.py:195`) | **Reconcile in v3 (lean v2).** Inline form is fine; keep as inline unless a separate `emissivity_air` accessor is needed for reuse elsewhere. |
| P12 | `wind_function` | `wind_function` | **Keep v2 in v3.** |
| P13 | `mf_q_longwave_down` | `flux_atmospheric_longwave` | **Keep v2 in v3** (Brunt 1+0.17·C² formulation, identical). |
| P14 | `mf_q_longwave_up` | `flux_upwelling_longwave` | **Keep v2 in v3** (Stefan-Boltzmann with 0.97 emissivity; v1 carries the same factor). |
| P15 | `q_latent` | `flux_latent_heat` | **Keep v2 in v3** + N1 fix flows through via `latent_heat_vaporization`. |
| P16 | `q_sensible` | `flux_sensible` | **Keep v2 in v3.** |
| P17 | `q_sediment` | `flux_sediment` | **Keep v2 in v3** + N6 documentation fix. |
| P18 | `q_net` (W/m² · 86400 · dt) | `flux_net` (W/m²) | **Reconcile in v3 (Keep v2 convention).** See N7. |
| P19 | `dTdt_water_c` (with depth ramp + rate cap) | `temperature_change` (no guards) | **Port v1→v3 guards.** Keep v2 framework shape; add depth-ramp and rate-cap inside v3 `temperature_change` (N2, N3). |
| P20 | `t_water_c` (state update) | inline at `Temperature.run:138` | **Keep v2 in v3** (Forward Euler `T_new = T_old + ΔT` is identical). v3 framework's integrator-pattern contract (architecture spec §4) makes this explicit. |

## 3. Orchestration / framework

| # | Topic | v1 | v2 | v3 disposition |
|---|---|---|---|---|
| O1 | **Process registration** | `@base.register_variable(models=...)` decorator on functions; topological sort via `sorter.py`. | `Process` subclass + `@ProcessFactory.register("name")` on `from_config` staticmethod (`temperature.py:88-91`). | **Keep v2 in v3.** Architecture spec §4 mandates v2 framework as v3 baseline. |
| O2 | **State container** | `xarray.Dataset` with `(time_step, x, y)` dims. | `clearwater_data.VariableRegistry`; `get_at_time` / `set_at_time`. | **Keep v2 in v3.** |
| O3 | **Per-step kernel** | `_iter_computations_fast` (`base.py:534`): cached `(name, callable, [arg_names])` plan, plain `dict[str, np.ndarray]` buffers, direct `DataArray.values[slot] = …` writes. 418× speedup vs. naive `Dataset.__setitem__`. Commit `6daa65e`. | Per-substep Python iteration over `Process.run`; no compute-plan cache. | **Port v1→v3** (architecture spec §3.2). v3 `Model` caches `(process, time_step, variables)` tuples, precomputes next-fire time per process, skips re-checking `current_time_seconds % time_step_seconds`. Variable-level optimization is a deferred follow-up if profiling shows v3 materially slower than v1. |
| O4 | **Wet-mask gating** | `Model.increment_timestep(wet_mask=...)` kwarg; orchestration-level NaN-write for masked cells (commit `3d18965`). | None at orchestration; per-process `xr.where(volume > 0, ...)` only. | **Port v1→v3.** Add a registry-level `wet_mask` variable; v3 `Model.run` honors mask by skipping `set_at_time` for masked cells. Remove redundant per-process masks (architecture spec §4, TSM spec §3.2). |
| O5 | **Hotstart from `xr.Dataset`** | `EnergyBudget(hotstart_dataset=..., hotstart_timestep=...)` (commit `1a226dd`). Loads saved registry state at the saved time. | None. | **Port v1→v3.** Add `hotstart_dataset` and `hotstart_timestep` kwargs to `init_from_file`. Per-process substep state defaults to "fresh start"; processes may opt in via `to_hotstart()` / `from_hotstart()` (TSM spec §3.2 resolution 2026-05-04). v3 `Temperature.from_hotstart` sets `__skip_first_time_step = False`. |
| O6 | **Sub-stepping per process** | One model step per process. | Per-process `time_step` (`timedelta`); `Model` checks `current_time_seconds % process.time_step_seconds == 0` (`model.py:104`). | **Keep v2 in v3.** v3 Phase 3 kernel optimization preserves substepping by precomputing next-fire time per process. |
| O7 | **Chunking execution path** | None in `clearwater_modules`; the streaming/chunking work lives in `clearwater_riverine`. | `__process_loop_chunked` skeleton (Paul, commit `d712c59`); 4 TODOs at `model.py:183, 188, 189, 215`. | **Resolve TODO in v3.** Mirror riverine's chunking conventions per `clearwater_riverine.transport.py` and `clearwater_riverine.constituents.py` (TSM spec §3.2 resolution 2026-05-04). |
| O8 | **Configuration** | Code-only (`EnergyBudget(meteo_parameters=..., temp_parameters=...)`). | YAML via `init_from_file` (`config/init.py`); top-level keys `model.{start_datetime,end_datetime,time_step,root_directory,output_variables}`, `processes`, `data_sources`, `variable_map`. | **Keep v2 in v3.** v3 adds two **optional** top-level keys: `hotstart` and `wet_mask` (TSM spec §3.3). When neither is present, v3 behavior matches v2 exactly (backward compatible). |

## 4. v2-specific behaviors to preserve (Keep v2 in v3)

| # | Topic | Source | v3 disposition |
|---|---|---|---|
| K1 | **`__skip_first_time_step` coupling logic** | `temperature.py:83`, consumed at `temperature.py:98-100`. Comment at line 82: "V1 of the coupling had timestep 1 be skipped and started processing on timestep 2." Added in upstream commit `209b67f`. | **Keep v2 in v3.** Per TSM spec §3.1: carry forward as-is, with a one-line comment documenting the v1-coupling-compat reason. v3 `Temperature.from_hotstart` overrides to `False` when resuming from hotstart (O5). |
| K2 | **Sediment-temperature optional flag** | `use_sediment_temperature` kwarg (`temperature.py:57`) and gate in `flux_sediment` (`temperature.py:288-289`). Anthony, February 2026. | **Keep v2 in v3.** Default `True`. v1 has the same `use_sed_temp` flag; behaviorally aligned. |
| K3 | **Brutsaert-coefficient class-level constants** | `__A0` … `__A6` at `temperature.py:25-31` (private name-mangled). | **Keep v2 in v3.** Equivalent to v1's `tsm/constants.py::Temperature::a0`-`a6`. |
| K4 | **`flux_*` decomposition with sign-bearing returns** | v2 returns negative values from `flux_upwelling_longwave` and `flux_latent_heat` so `flux_net` is a simple sum. | **Keep v2 in v3.** Cleaner than v1's separate sign tracking in `q_net`. |

## 5. Code hygiene (cleanups required in v3)

| # | Item | Location in v2 baseline | v3 disposition |
|---|---|---|---|
| H1 | **Commented-out debug `print(...)` blocks (no information value)** | `temperature.py` lines ~175-188, 227-242, 354-360, 499-504, 613-617, 637-664. | **Delete in v3.** Per TSM spec §3.1, do not just leave them commented. v3 uses Python `logging` at the orchestration layer if log statements are needed. (Streaming-local v2 has these uncommented; that's a streaming-tracking issue, not a v3 issue.) |
| H2 | **`# TODO: We should make the get method handle time selected`** | `temperature.py:116` | **Defer / leave note in v3.** Framework-level concern about `registry.get_at_time`; not TSM-specific. Remove only the TODO comment in v3 and reference the framework-level issue. |
| H3 | **`# TODO: Should this change as a function of temperature?` (emissivity_air)** | `temperature.py:191` | **Resolve TODO in v3 (no behavior change).** Replace with a one-line docstring noting that emissivity is parameterized as `9.37e-6 · T_K²` (Brunt-style) and that the temperature dependence is already present (the TODO question is answered by the formula itself). |
| H4 | **`# TODO: this needs the richardson function`** (above `wind_function`) | `temperature.py:494` | **Delete TODO in v3.** Already resolved — `wind_function` takes `richardson_function` as an arg. |
| H5 | **`# TODO: this needs to be reworked to support array inputs`** (richardson_number) | `temperature.py:609` (deleted in current upstream; verify in v3) | **Delete TODO in v3.** v2 already vectorized via `xr.where`. |
| H6 | **`# TODO: can we find a more efficient way to calculate this?` (richardson_function)** | `temperature.py:626` | **Defer note in v3.** True but not a Phase 2 task; replace with a one-line note pointing to v1's `np.select`-based form as a future optimization. |
| H7 | **`# TODO: We should get Billy to both explain and provide guidance on wind parameters`** | `temperature.py:72` | **Defer TODO in v3.** Wind-parameter physics question is unresolved; preserve TODO but reword to actionable form. |
| H8 | **`# TODO: verify if this equation is correct for both fresh and salt water` (water_density)** | `temperature.py:445` | **Defer note in v3.** Salinity is out of scope for v3 1.0.0; preserve as a note for future v3.x work. |

## 6. Tests

| # | Item | v1 source | v3 target |
|---|---|---|---|
| T1 | TSM calculations regression | `tests/test_5_tsm_calculations.py` (15 tests; rebaselined for latent-heat fix). | **Port** to `tests/v3/test_5_tsm_calculations_v3.py`. |
| T2 | Latent-heat regression | `tests/test_tsm_latent_heat.py` (4 tests). | **Port** to `tests/v3/test_tsm_latent_heat_v3.py`. |
| T3 | Thin-water stability ramp | `tests/test_tsm_stability_ramp.py` (6 tests). | **Port** to `tests/v3/test_tsm_stability_ramp_v3.py`. |
| T4 | Hotstart roundtrip | `tests/test_hotstart_roundtrip.py`. | **Port** to `tests/v3/test_hotstart_roundtrip_v3.py`. |
| T5 | TSM module integration | `tests/test_4_tsm_module.py` (4 tests: init, sort, timestep, sediment-temp toggle). | **Port** the substantive cases; v2 framework changes some assertions. |
| T6 | Existing v2 Richardson unit tests | `tests/unit/temperature/test_richardson.py` (2 tests). | **Keep & extend in v3** under `tests/v3/`. |
| T7 | v2/v3 parity (new) | — | **New**: `tests/v3/test_v2_v3_parity.py`. With v3 corrections disabled (depth ramp `0.0`, rate cap `+inf`, latent-heat fix toggled off via debug switch), v3 must reproduce v2 outputs exactly on Sumwere Creek. |
| T8 | Coupled TSM+Riverine (new) | `examples/03_Example_Coupled_TSM_and _Riverine.ipynb` | **New**: `tests/v3/test_coupled_tsm_riverine_v3.py` programmatic version with wall-time budget assertion. |
| T9 | Empty fixture-only file | `tests/test_5_tsm_calculations_v2.py` (no `def test_…`). | **Delete or replace** in v3; T1 supersedes. |

## 7. Items not in the v3 scope

For completeness, items the architecture spec §3 (umbrella non-goals) and TSM spec §2 (TSM non-goals) place outside v3:

- Ice cover, riparian shading, atmospheric coupling beyond met forcing — neither v1 nor v2 has these; v3 does not add them.
- `clearwater_modules` v1 framework metadata (`Variable` dataclass, topological sort via `sorter.py`) — v3 retains v2's class-based composition and per-process `variables: list[str]` declaration.
- Numba `@njit` kernels in v1 `shared/processes.py` — out of scope for TSM (those kernels are NSM1-relevant).

---

## Phase-1-readiness summary

**To start Phase 1 (v3 scaffold), the following must be in place:**

1. ✓ Phase 0 deliverable (this document) reviewed and accepted.
2. Clear branching choice: branch v3 work off `upstream/memory-refactor-pytestUpdate` (per architecture spec §1) so v3 starts from the **upstream** v2 baseline, not the older streaming-local copy. The streaming-local `src/clearwater_modules_v2/` divergence (5 commits behind) is independent of v3 and should be resolved separately.
3. Pixi env built in the streaming repo (`pixi install -e dev` per architecture spec §1).

**Hot items for Phase 2 (TSM merge), in order of priority:**

1. N1 — latent-heat unit fix (low risk, mechanical port, blocks correct evaporative cooling).
2. N2, N3 — thin-water depth ramp + rate cap (low risk, mechanical port; defaults reproduce v1 hardened behavior).
3. N4 — vectorize ``mixing_ratio_air`` guard (upstream's scalar ``if`` raises on multi-cell DataArrays; replace with ``xr.where``). Required for any multi-cell coupled run including Sumwere Creek.
4. H1 — delete commented debug print blocks.
5. N5 — delete the Richardson `−1` TODO.
6. N6 — replace the `flux_sediment / 0.5` TODO with a one-line docstring.
7. N7, N8 — reconcile `q_net` units convention and remove redundant per-process wet-mask (the latter waits on Phase 3).
8. K1, K2 — preserve `__skip_first_time_step` and `use_sediment_temperature`.

**Hot items for Phase 3 (orchestration), in order of priority:**

1. O3 — kernel optimization adapted to process-level dispatch.
2. O4 — registry-level wet-mask.
3. O5 — hotstart, with optional `to_hotstart()` / `from_hotstart()` per process.
4. O7 — resolve the four chunking TODOs by mirroring riverine.

---

**End of Phase 0 deliverable.** Awaiting review before Phase 1.
