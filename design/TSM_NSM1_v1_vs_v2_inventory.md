# TSM and NSM1: side-by-side inventory of `clearwater_modules` (v1) vs `clearwater_modules_v2` (streaming framework)

Reviewer: water-quality-source-code-reviewer agent
Date: 2026-05-04
Repos under review:

- Repo A (older, mature): `/Users/todd/GitHub/ecohydrology/ClearWater-modules/` (`main`-style branch, last commit `b5f2c8b`, 2026-02-20)
- Repo B (streaming branch, contains both legacy `clearwater_modules` and experimental `clearwater_modules_v2`): `/Users/todd/GitHub/ecohydrology/ClearWater-modules-streaming/` (last commit `5674cbc`, 2026-05-02)

## Headline summary

- **Repo B is two repositories in one.** It carries (1) a hardened, optimized fork of the legacy `clearwater_modules` package with TSM bug fixes, hotstart, wet-mask, and a 418x kernel speedup, AND (2) the experimental `clearwater_modules_v2` streaming framework. The most useful artifacts for production are in (1), not (2).
- **`clearwater_modules_v2` is a working skeleton, not a working model.** TSM is the only v2 process that runs end-to-end (tested and matched to v1 within a few microdegrees). v2 NSM1 covers only nitrogen, floating algae, and benthic algae stubs out of the ~16 NSM1 constituents. Phosphorus, carbon, CBOD, DOX, pathogen, alkalinity, POM, N2/TDG, and sediment diagenesis are absent.
- **The `clearwater_modules_v2.processes.nitrogen` rate update equation is broken** (`new = old + old * rate * dt`, see lines 101 and 115 of `nitrogen.py`); even though sub-rate methods now match v1 within tolerance after the LimnoTech fix, the integrator at line 101/115 multiplies (rather than adds) the rate.
- **The 418x performance kernel (`clearwater_modules-streaming/src/clearwater_modules/base.py`, commit `6daa65e`) lives in the legacy package, not in v2.** v2 has no equivalent vectorized inner loop; it iterates a Python `for` loop over time inside `Model.run` and dispatches each `Process.run` per substep with no compute-plan caching.
- **TSM (v1) has 482 lines of process code and ~150 dynamic variables across NSM1; v2 reimplements TSM in 645 lines as a single class and one nitrogen process class with ~370 lines.** The v2 implementation is methodologically a copy of v1 with refactored composition; it is not a new model.
- **v2's only configuration path is YAML through an `init_from_file` that depends on a `clearwater_data` package not under review here**; the existing example YAML (`clearwater_modules_v2/config/example_config.yml`) drives only the `riverine` and `temperature` processes.
- **No NSM1 v2 integration test exists** that runs more than one timestep of more than one process; the v1-parity tests cover sub-rate methods only.
- **The 418x kernel optimization, hotstart support, wet-mask gating, latent-heat unit fix, and TSM thin-water stability fixes were all added on the streaming branch on top of the legacy package**, not on top of v2.

## 1. Framework architecture comparison

### 1.1 How variables and processes are registered

Repo A (legacy `clearwater_modules`, also present in Repo B's `src/clearwater_modules/` with extensions):

- Variables are dataclasses (`Variable` in `src/clearwater_modules/shared/types.py`) with fields `name`, `long_name`, `units`, `description`, `use ∈ {static, dynamic, state}`, and an optional `process: Callable`.
- Registration is performed at module import time via the `@base.register_variable(models=...)` decorator in `src/clearwater_modules/base.py:537`. Each `Model` subclass owns a class-level `_variables` list.
- Processes are plain Python functions whose argument names must equal registered variable names. The dependency order is resolved automatically by `src/clearwater_modules/sorter.py` from the `__annotations__` of each process function.
- The model loop is `Model.increment_timestep` in `base.py:495` (and the optimized `_iter_computations_fast` at `base.py:534` in Repo B).

Repo B (`clearwater_modules_v2`):

- Processes are `class` objects subclassing `Process` (`src/clearwater_modules_v2/processes/base.py:14`). Each process declares a `variables: list[str]` attribute listing the registry keys it reads, and implements `__init__`, `init_process`, and `run`.
- A `ProcessFactory` decorator (`base.py:61`) registers a class's `from_config` classmethod under a string name (e.g. `@ProcessFactory.register("temperature")`), and `Model.__init__` consumes a tuple of pre-instantiated `Process` objects.
- There is no equivalent of the variable-level `Variable` dataclass; metadata (units, long names) lives only in docstrings.
- Variable storage is delegated to a `VariableRegistry` class imported from a separate `clearwater_data` package (`from clearwater_data.variables import VariableRegistry`). v2 cannot run without the `clearwater_data` package, which is not part of this review.

### 1.2 Model loop and data structures

| Aspect | Repo A / Repo B legacy `clearwater_modules` | Repo B `clearwater_modules_v2` |
|---|---|---|
| State container | `xarray.Dataset` with one DataArray per variable, dims `(time_step, x, y)` | `VariableRegistry` (external); each variable wraps an `xr.DataArray` or scalar (`FloatVariable`) |
| Time iteration | `for ts in range(N): model.increment_timestep()` | `while current_time < end_time: for process in processes: process.run(time, registry); current_time += time_step` (`model.py:96`) |
| Per-step kernel | `_iter_computations_fast` (Repo B): single `for` over cached compute plan, plain numpy dict, direct `DataArray.values[slot] = ...` write at end of step | Each `Process.run` calls `registry.get_at_time` and `registry.set_at_time` per variable; no cached plan |
| Per-step performance (NSM1) | ~3 ms (Repo B, `6daa65e`); ~1292 ms (Repo A, pre-optimization) | Not measured; Python-level loop over processes per substep |
| Sub-stepping per process | One model step per process | Each `Process` carries its own `time_step` (`timedelta`); `Model` checks `current_time_seconds % process.time_step_seconds == 0` (`model.py:104`) to decide whether to run |
| Vectorization across cells | Native via numpy/xarray broadcasting; selected hot kernels are `@numba.njit` decorated in `shared/processes.py` | Native via xarray broadcasting; no numba; no SIMD |
| Parallelism | None beyond NumPy threads and dask | None |

### 1.3 Configuration mechanism

Repo A and Repo B legacy package: pure code-defined. The user instantiates `EnergyBudget(meteo_parameters={...}, temp_parameters={...}, initial_state_values={...})` directly. No YAML.

Repo B `clearwater_modules_v2`: YAML-driven (`config/example_config.yml`, `config/init.py`). The user writes a YAML file enumerating processes, data sources, and a variable map, then calls `init_from_file(yaml_path)`. The current `init_from_config` requires `pandas`, `clearwater_data.io.zarr`, and `clearwater_data.io.csv`, and writes intermediate inputs to a Zarr store. This is the only configuration entry point in v2.

### 1.4 Performance: the 418x kernel optimization

The commit `6daa65e44a74db2341587f7cd470f708faa084cf` ("Optimize Model.increment_timestep kernel to bypass per-variable Dataset.__setitem__", 2026-04-30) lives in the **streaming repo's legacy package** (`src/clearwater_modules/base.py`), not in `clearwater_modules_v2`. Verified by `git show 6daa65e --stat`: the diff modifies only `src/clearwater_modules/base.py` (210 insertions, 52 deletions).

Mechanism: a cached compute plan (`(name, callable, [arg_names])` tuples in dependency order) is built once per `Model` class, stored on the class as `_cached_compute_plan`. The hot loop (`_iter_computations_fast`, `base.py:534`) runs the kernel against a plain `dict[str, np.ndarray]` of buffers and writes results back to the parent dataset via direct `DataArray.values[slot] = ...` assignment, bypassing `Dataset.__setitem__` (which previously triggered `xarray.merge_core -> reindex_all` for every variable assignment, an O(N^2) cost in the variable count).

Reported benchmark (5-cell mesh, 100 timesteps, NSM1 `NutrientBudget`): 1292 ms/step before, 3.1 ms/step after, ~418x. Scaling effectively flat from 1 to 1000 cells. Test suite passes (393 tests, identical to baseline) and runs 3.4x faster.

This optimization is **not present in `clearwater_modules_v2`** and would have to be reinvented if v2 grows to NSM1's variable count.

### 1.5 Test infrastructure

Repo A `tests/`:

- 17 test files, ~830 KB total source.
- Coverage: NSM1 carbon, DOX, alkalinity, CBOD, phosphorus, POM, pathogen (PX), N2, algae, benthic algae, nitrogen; TSM module; TSM calculations; equation sort.
- Style: per-process expected-value tables (`test_*_calculations.py`) imported from manual calculation spreadsheets in `tests/NSM Manual Calcs/`.

Repo B `tests/`:

- All of Repo A's tests plus:
  - `test_5_tsm_calculations.py` (rebaselined for the latent-heat fix; 15 expected values updated)
  - `test_tsm_latent_heat.py` (4 tests)
  - `test_tsm_stability_ramp.py` (6 tests)
  - `test_hotstart_roundtrip.py` (4 cases, NSM1 + TSM, `rtol=1e-12`)
  - `test_5_benthic_algae_calculations_v2.py` (4 tests, v1-parity for v2 `BenthicAlgae` sub-rate methods)
  - `test_5_floating_algae_calculations_v2.py` (6 tests, v1-parity for v2 `FloatingAlgae` sub-rate methods)
  - `test_5_nitrogen_calculations_v2.py` (7 tests, v1-parity for v2 `Nitrogen` sub-rate methods)
  - `test_5_tsm_calculations_v2.py` (0 `def test_` functions; only fixture data)
  - `tests/sediment/` (10 test files, ~190 KB, for the new SSM)
  - `tests/unit/temperature/test_richardson.py` (Richardson number unit tests for v2 `Temperature`)
- Per the streaming branch README, 430 tests pass with 2 expected-failure cases.

### 1.6 Active development status

Repo A activity:

```
b5f2c8b 2026-02-20 Document NSM1 DOX rate bug caused by sentinel values in constants.py
3994c2b 2025-..   Update .gitignore
ff9b488 2025-..   Adding draft Sphinx documentation
79b5673 2025-..   Update README.md
0f867a3 2025-..   Merge pull request #88 from EcohydrologyTeam/KW_Phosphorus
b3ac970 2025-..   Remove kelvin_to_celsius
016b73d 2025-..   Allow Negative Values
fadcbcc 2025-..   Merge branch 'main' into KW_Phosphorus
dc14869 2025-..   Merge pull request #91 from EcohydrologyTeam/performance-profiling
9b66735 2025-..   Update environment.yml
```

Three commits since 2026-04-01 (none after 2026-02-20 in the user-visible history). Repo A is essentially frozen.

Repo B activity since 2026-04-01 (12 commits, all on the streaming branch):

```
5674cbc 2026-05-02 Add SSM (Sediment Simulation Module) — clean-room SEDZLJ port
37e1d18 2026-05-02 Document streaming-branch additions in README
d9505c6 2026-05-01 TSM kernel correctness: thin-water stability + latent-heat unit fix
1a226dd 2026-05-01 Add hotstart_dataset kwargs to NSM1 and TSM kernels
21e28bb 2026-05-01 Add v1-parity tests for v2 BenthicAlgae, FloatingAlgae, Nitrogen
3926a81 2026-05-01 Fix 4 bugs in v2 FloatingAlgae and Nitrogen processes
ae2c5d5 2026-04-30 Add numba to CI test dependencies
f3b4db4 2026-04-30 Add Windows + Linux + macOS CI workflow
4a95da4 2026-04-30 Remove TSM debug print statements that crash on multi-cell arrays
6daa65e 2026-04-30 Optimize Model.increment_timestep kernel to bypass per-variable Dataset.__setitem__
3d18965 2026-04-30 Add wet_mask kwarg to Model.increment_timestep
48cec01 2026-04-30 Merge branch 'main' into memory-refactor
```

All recent development is on the streaming branch. None of these improvements have been backported to Repo A.

### 1.7 Architecture comparison table

| Dimension | Legacy `clearwater_modules` (Repo A and Repo B) | `clearwater_modules_v2` (Repo B only) |
|---|---|---|
| Variable model | `Variable` dataclass + `@register_variable(models=...)` decorator | `Process.variables: list[str]` (no rich metadata) |
| Process registration | Function-style; argument names match variable names; topological sort | `Process` subclass; explicit `run(time, registry)` body |
| Time stepping | Outer Python `for` over `time_steps`; one Model step per call | Outer `while time < end_time`; per-process substep filter on `time_step_seconds` |
| State container | `xarray.Dataset` with `(time_step, x, y)` dims | `VariableRegistry` from `clearwater_data` |
| Configuration | Code-only | YAML via `clearwater_modules_v2/config/init.py` (requires `clearwater_data`) |
| Numba | `@numba.njit` on selected `shared/processes.py` kernels | None |
| Hot loop | `_iter_computations_fast` with cached compute plan (Repo B only) | Per-substep Python iteration; no caching |
| Hotstart from `xr.Dataset` | Yes (Repo B); `EnergyBudget`/`NutrientBudget` accept `hotstart_dataset=` | No |
| Wet-mask gating | Yes (Repo B); `increment_timestep(wet_mask=...)` | No |
| Coupling to ClearWater-Riverine | Direct: `update_state_values=` overrides | Through `Riverine` Process class wrapping `cwr.ClearwaterRiverine` |
| Tests for TSM | 15 calculation tests + latent-heat + stability ramp + hotstart | 1 fixture-only file + 1 Richardson unit-test file (no `def test_` in `test_5_tsm_calculations_v2.py`; the v2 TSM is exercised indirectly through unit tests) |
| Tests for NSM1 | 11 per-constituent calculation files (~830 KB) | 7 sub-rate parity tests for nitrogen only |

## 2. TSM completeness in each repo

### 2.1 State variables, dynamic variables, fluxes

| Item | Legacy TSM (`src/clearwater_modules/tsm/`) | v2 `Temperature` (`src/clearwater_modules_v2/processes/temperature.py`) | Notes |
|---|---|---|---|
| `water_temp_c` (state) | yes (`state_variables.py:12`) | yes (registry key `water_temperature`) | Both |
| `surface_area`, `volume` (state) | yes (`state_variables.py:30`) | yes (registry keys `wetted_surface_area`, `volume`) | Both |
| `sed_temp_c` (boundary) | static (`static_variables.py:152`) | registry key `sediment_temperature` | Both; v2 reads as scalar `FloatVariable` in dev script |
| `air_temp_c`, `q_solar`, `eair_mb`, `pressure_mb`, `cloudiness`, `wind_speed` (forcings) | static (`static_variables.py`) | registry keys (forced via YAML or code) | Both |
| `air_temp_k`, `water_temp_k` (dynamic) | yes (`dynamic_variables.py:12-26`, `processes.py:10-29`) | inline as `conversions.celsius_to_kelvin(...)` | Both compute the same value |
| `mixing_ratio_air` | yes (`processes.py:32`) | yes (`temperature.py:501`) | Both. v2 has a defect: `if atmospheric_vapor_pressure == atmospheric_vapor_pressure: return 0.0` (`temperature.py:511`) -- this branch is `if x == x` which is True for any non-NaN scalar; the function returns 0 unconditionally for non-NaN input. **Suspected bug, needs verification.** |
| `density_air` | yes (`processes.py:48`) | yes (`temperature.py:520`) | Both |
| `density_water` (`mf_density_water`) | yes (`processes.py:67`) | yes (`temperature.py:428`) | Both |
| `esat_mb` (saturation vapor pressure) | yes (`processes.py:92`) | yes (`temperature.py:456`) | Both, identical Brutsaert (1982) coefficients |
| `density_air_sat` | yes (`processes.py:127`) | yes (`temperature.py:541`) | Both |
| `ri_number` (Richardson number) | yes (`processes.py:150`) | yes (`temperature.py:569`) | Both. v2 has a `# TODO` (line 597) noting uncertainty about a sign convention vs v1. |
| `ri_function` (stability function) | yes (`processes.py:171`); uses `np.select` | yes (`temperature.py:618-643`); uses chained `xr.where` | Same numerical formulation. v2 is not vectorized as efficiently. |
| `lv` (latent heat of vaporization) | yes (`processes.py:213`); **fixed in Repo B** to convert K to C before applying polynomial (was a 26-27% bias when input was Kelvin) | yes (`temperature.py:445`); converts to Kelvin before polynomial -- this is the **opposite of the Repo B fix** | **Important discrepancy.** Repo A's pre-fix legacy version applied the polynomial to Kelvin (bug). Repo B's TSM was fixed (`d9505c6`) so that `mf_latent_heat_vaporization` converts Kelvin to Celsius first. v2 (`temperature.py:453`) calls `2499999 - 2385.74 * conversions.celsius_to_kelvin(water_temperature)`, which converts to Kelvin (since the input is already Celsius) -- that is, it reproduces the original v1 bug. **v2 has the latent-heat unit error.** |
| `cp_water` (specific heat of water) | yes (`processes.py:223`); uses `np.select` table | yes (`temperature.py:353`); identical table | Both |
| `emissivity_air` | yes (`processes.py:260`); `0.00000937 * Tk^2` | inline at `temperature.py:186`; `9.37e-6 * Tk^2` | Same |
| `wind_function` | yes (`processes.py:271`) | yes (`temperature.py:486`) | Same |
| `q_longwave_down` (atmospheric LW) | yes (`processes.py:296`) | yes (`flux_atmospheric_longwave`, `temperature.py:153`) | Both implement Brunt-style `(1 + 0.17·C^2)·ε_a·σ·Ta^4` |
| `q_longwave_up` (upwelling LW) | yes (`processes.py:317`) | yes (`flux_upwelling_longwave`, `temperature.py:137`) | Both |
| `q_latent` | yes (`processes.py:331`); `(0.622/p)·ρ_w·Lv·f(U)·(esat - eair)` | yes (`flux_latent_heat`, `temperature.py:198`) | Same form. v2 carries an `embedded print(f"...")` statement (`temperature.py:218-232`) that calls `float(...)` on every term -- this **crashes on multi-cell arrays** (commit `4a95da4` removed equivalent prints from a v1-style codepath but the v2 prints remain in `temperature.py`). |
| `q_sensible` | yes (`processes.py:358`) | yes (`flux_sensible`, `temperature.py:247`) | Same |
| `q_sediment` | yes (`processes.py:387`); guarded by `use_sed_temp` | yes (`flux_sediment`, `temperature.py:271`); guarded by `self.use_sediment_temperature` | Same |
| `q_net` | yes (`processes.py:418`); `(qsens + qsol + qsed + qLW_dn - qLW_up - qlat) * 86400 * dt` | yes (`flux_net`, `temperature.py:292`); `sensible + solar_flux + sediment + atmospheric + upwelling + latent` (signs absorbed into the flux functions) | Same physics but accounting differs. v1 splits sign into `q_net`; v2 puts the sign into each flux (e.g. `flux_upwelling_longwave` returns negative). |
| `dTdt_water_c` | yes (`processes.py:448`); **Repo B adds depth-ramp + rate-cap regularization** | inline in `temperature_change` (`temperature.py:380`); no thin-water guard, no rate cap | **Repo A and v2 are both vulnerable to thin-water blow-up; only Repo B's legacy TSM has the depth ramp and rate cap.** |
| `t_water_c` (state update) | yes (`processes.py:471`); `T_new = T_old + dT` | inline at `temperature.py:130`; `updated_water_temperature = water_temperature + delta_water_temperature`; gated by `xr.where(volume > 0, dT, 0)` | Similar; v2 adds an explicit dry-cell guard. |
| Sediment temperature evolution | Static (held constant) | Static (read once per timestep from registry) | Neither model integrates sediment temperature; both treat it as a forcing. The legacy code has a `use_sed_temp` flag controlling whether the heat flux exchanges with sediment. |
| Ice cover | Not implemented | Not implemented | Neither |
| Riparian shading | Not implemented | Not implemented | Neither |
| Atmospheric coupling beyond met forcing | None | None | Neither |
| Per-segment vs per-cell budget | Per-cell (vectorized over `(x, y)`) | Per-cell (vectorized over registry-array dims) | Both |
| Time stepping | Forward Euler; `dt` from static; no substepping | Forward Euler; `time_step` from `timedelta`; per-substep | Same scheme |
| Stability guard (CFL/thin-water) | Repo B only: depth ramp on `q_net`, rate cap on `dTdt`. Repo A: none. | None | Important gap in v2 |
| Tests | `test_5_tsm_calculations.py` (15 tests, Repo B-rebaselined), `test_tsm_latent_heat.py` (4 tests), `test_tsm_stability_ramp.py` (6 tests), `test_hotstart_roundtrip.py` (TSM included) | `tests/unit/temperature/test_richardson.py` (Richardson number); `test_5_tsm_calculations_v2.py` (no `def test_` functions) | v2 TSM coverage is a small fraction of legacy |

### 2.2 What is missing in v2 TSM relative to legacy TSM

1. **Latent-heat unit fix.** v2 has the original (pre-fix) Kelvin-vs-Celsius bug in `mf_latent_heat_vaporization`; the polynomial is calibrated for Celsius but `temperature.py:453` calls `conversions.celsius_to_kelvin(water_temperature)` before applying it, biasing Lv ~26-27% low and biasing simulated water temperature warm.
2. **Thin-water stability.** Legacy Repo B has `q_net_depth_ramp_ref` (depth ramp on net flux) and `dTdt_max_per_hour` (cadence-invariant rate cap on dT/dt) inside `dTdt_water_c`. v2 has neither. Without these, v2 will diverge at low depths.
3. **Multi-cell `print` debugging.** v2 `temperature.py` lines 167-179, 218-232, 595-607 call `float(...)` on incoming `xr.DataArray` arguments inside `print(f"...")`, which raises if the input is anything other than a scalar.
4. **Hotstart dataset support.** Legacy Repo B accepts `hotstart_dataset=` and `hotstart_timestep=` on `EnergyBudget`. v2 has no equivalent.
5. **Wet-mask gating.** Legacy Repo B accepts `wet_mask=` on `Model.increment_timestep`. v2 has no equivalent and instead masks via `xr.where(volume > 0, dT, 0)` after the fact.
6. **`mixing_ratio_air` always returns 0 for non-NaN input** (line 511: `if atmospheric_vapor_pressure == atmospheric_vapor_pressure: return 0.0`). The intent appears to have been a NaN guard (`if x != x: return 0`), but the comparison was inverted. **Suspected bug, needs verification.**
7. **No `np.select`-based vectorized stability function.** v2 uses chained `xr.where` (lines 618-643), which is correct but slower and less readable.
8. **No numba acceleration.** Legacy `shared/processes.py` decorates several inner kernels with `@numba.njit`. v2 has no equivalent.

## 3. NSM1 completeness in each repo

### 3.1 Legacy NSM1 (Repo A and Repo B `clearwater_modules`)

The legacy NSM1 has 16 state variables in `state_variables.py` and ~150 dynamic variables registered through `dynamic_variables.py` (1400 lines, 158 `Variable(...)` registrations) and `dynamic_variables_global.py`. The process module `processes.py` is 3540 lines with ~290 `def` definitions. By inspection of `state_variables.py:1-157` and `processes.py`, all of the following NSM1 constituents are implemented in legacy:

`Ap, Ab, NH4, NO3, OrgN, N2, TIP, OrgP, POC, DOC, DIC, POM, CBOD, DOX, PX, Alk`.

Sediment diagenesis (Di Toro multi-G) is **not implemented** in legacy NSM1; the `use_SedFlux` flag exists in `constants.py:267` but is `False` by default, and the only sediment terms in `processes.py` are simple parameterized flux constants (`NH4fromBed`, `DIPfromBed`, `NO3_BedDenit`, `DIC_sed_release`). True benthic sediment diagenesis lives in `clearwater_modules/nsm2/sediment_flux/` (NSM2), which is out of scope for this review.

### 3.2 v2 NSM1

`clearwater_modules_v2/processes/`:

- `temperature.py` (TSM, covered above)
- `nitrogen.py` (NH4 and NO3 only; OrgN absent)
- `floating_algae.py` (Ap only)
- `benthic_algae.py` (Ab only)
- `nutrients/__init__.py` (empty file, 0 bytes)
- `sediment/` (the new SSM, sediment transport / SEDZLJ port; out of scope unless coupled to nutrient budget)
- `riverine.py` (transport coupling wrapper)

There is no module for phosphorus, carbon, CBOD, DOX, pathogen, alkalinity, POM, N2, TDG, or sediment diagenesis.

### 3.3 NSM1 coverage matrix

Legend: ✓ full, ◐ partial, ✗ absent

| Constituent / process | Legacy NSM1 (Repo A and Repo B) | v2 (`clearwater_modules_v2`) | Notes |
|---|---|---|---|
| **Ap (floating algae)** state | ✓ `state_variables.py:15` | ✓ registry key `algae_floating` | |
| Ap growth (FT, FN, FP, FL) | ✓ `processes.py:415-606` (mu, mu_max_tc, ApGrowth, FN, FP, FL with multiplicative/limiting/harmonic options) | ◐ `floating_algae.py:169-218` (3 growth options); FT and Arrhenius via `arrhenius_correction` | Numerical formulation differs: v2 returns the rate as `growth_rate * limit_p * limit_n * limit_l * algae`, then runs `algae = old + old * rate * dt * 86400` (line 122), which is `algae(1 + rate*dt*86400)` -- not the linear `algae + rate*dt`. **Likely integrator bug, see section 3.4.** |
| Ap respiration | ✓ `processes.py:621` (ApRespiration) | ◐ `floating_algae.py:231` (`rate_respiration`); does not feed into ammonium production via `floating_algae_process.ammonium_respiration()` (returns 0; `floating_algae.py:401`) | |
| Ap death | ✓ `processes.py:636-650` (ApDeath, ApDeath_OrgN, ApDeath_OrgP) | ◐ `floating_algae.py:220` (`rate_death`); no coupling to OrgN/OrgP (those constituents do not exist in v2) | |
| Ap settling | ✓ `processes.py:650` (ApSettling) | ◐ `floating_algae.py:244` (`rate_settling`) | |
| Ap NH4 preference (`ApUptakeFr_NH4`) | ✓ `processes.py:1206` | ◐ `nitrogen.py:235` (`ammonium_uptake_floating_algae`); contains a defective branch `xr.where(rate == np.nan, ...)` | `x == np.nan` is always False in IEEE 754, so the fallback never fires. **Bug, needs verification with array inputs.** |
| **Ab (benthic algae)** state | ✓ `state_variables.py:24` | ✓ registry key `benthic_algae` | |
| Ab growth (mub, FNb, FPb, FLb, FSb) | ✓ `processes.py:701-1025` | ◐ `benthic_algae.py:90-128` (multiplicative, limiting nutrient; `limit_density` only; `light_limitation_option` 1, 2, 3) | |
| Ab respiration / death | ✓ `processes.py:1040-1054` | ◐ inherited from `FloatingAlgae` | |
| Ab settling/burial | ✓ via Fw/Fb fractions | ✗ | |
| **NH4** state | ✓ | ✓ registry key `ammonium` | |
| Nitrification (with DO inhibition) | ✓ `processes.py:1422-1457` (`NitrificationInhibition` exponential, `NH4_Nitrification`) | ◐ `nitrogen.py:266-283` (`ammonium_nitrification`); inhibition exponential at `nitrogen.py:367` | Both use `1 - exp(-k*DO)`. v2 default `nitrification_oxygen_inhibition_factor = 1.0` is a placeholder, not the typical KsOxna value; **likely missing parameter coupling**. |
| Sediment NH4 release | ✓ `processes.py:1457` (`NH4fromBed`) | ◐ `nitrogen.py:227` (`ammonium_from_bed`) | |
| **NO3** state | ✓ | ✓ registry key `nitrate` | |
| Denitrification (water column, DO-dependent) | ✓ `processes.py:1604-1640` (`NO3_Denit`) | ◐ `nitrogen.py:297` (`nitrate_denitrification`); `half_saturation_oxygen` is hard-coded to `1` at the call site (`nitrogen.py:197`) with `# TODO: need argument` | |
| Sediment denitrification (`vno3_tc`) | ✓ `processes.py:1156-1640` | ◐ `nitrogen.py:321` (`nitrate_bed_denitrification`) | |
| **OrgN** state | ✓ `state_variables.py:51` | ✗ | |
| OrgN hydrolysis to NH4 (`kon_tc`, `OrgN_NH4_Decay`) | ✓ `processes.py:1173-1317` | ✗ | |
| OrgN settling (`vson`, `OrgN_Settling`) | ✓ `processes.py:1333` | ✗ | |
| **TIP** state | ✓ `state_variables.py:69` | ✗ (referenced as `phosphorus_total_inorganic` in `floating_algae.py` but not produced by any v2 process) | v2 reads TIP from the riverine model's mesh (`riverine.py:93`) but does not update it. |
| TIP partitioning (`fdp`) | ✓ `shared/processes.py:257` | ✗ | v2 hard-codes `phosphate_fraction_dissolved=0.5` (`floating_algae.py:113`) with `# TODO`. |
| TIP settling | ✓ `processes.py:1973` (`TIP_Settling`) | ✗ | |
| **OrgP** state and processes | ✓ `state_variables.py:78`, `processes.py:1833-2168` | ✗ | |
| **POC** state and processes | ✓ `state_variables.py:87`, `processes.py:2439-2550` | ✗ | |
| **DOC** state and processes | ✓ `state_variables.py:96`, `processes.py:2565-2671` | ✗ | |
| **DIC** state and processes | ✓ `state_variables.py:105`, `processes.py:2687-2858` (Henry's law, atmospheric reaeration, algal photosynthesis/respiration coupling) | ✗ | |
| **POM** state and processes | ✓ `state_variables.py:114`, `processes.py:2185-2317` | ✗ | |
| **CBOD** state and processes | ✓ `state_variables.py:123`, `processes.py:2334-2422` (multi-group `kbod_tc`, `ksbod_tc`, oxidation, sedimentation) | ✗ | |
| **DOX** state and processes | ✓ `state_variables.py:132`, `processes.py:2876-3123` (saturation `DOX_sat`, atmospheric reaeration, photosynthesis/respiration, nitrification, DOC oxidation, CBOD oxidation, SOD via `shared_processes.SOD_tc`) | ✗ | |
| **PX (pathogen)** state and processes | ✓ `state_variables.py:141`, `processes.py:3141-3227` (decay natural + light + settling) | ✗ | |
| **Alk** state and processes | ✓ `state_variables.py:150`, `processes.py:3246-3435` (denitrification, nitrification, algal growth/respiration, benthic algae) | ✗ | |
| **N2 / TDG** state and processes | ✓ `state_variables.py:60`, `processes.py:3452-3540` (`KHN2_tc`, `N2sat`, `dN2dt`, `TDG`) | ✗ | |
| **Sediment diagenesis (Di Toro multi-G: JNH4, JNO3, JCH4, JSO4, JH2S, JDIC, JDIP, SOD)** | ✗ in NSM1 (only parameterized fluxes); present in `clearwater_modules/nsm2/sediment_flux/` (out of scope) | ✗ | |
| Reaeration menu | ✓ `shared/processes.py:65` (`kah_20`: 9 hydraulic options); `kaw_20`: 13 wind options) | ✗ | v2 has no reaeration menu |
| Light extinction (`L`) | ✓ `shared/processes.py:202` (`L = lambda0 + lambdas*Solid + lambdam*POC/fcom + lambda1*Ap + lambda2*Ap^0.667`) | ◐ `floating_algae.py:306` (`limit_light`) accepts a single `light_attenuation_coefficient`; no Beer-Lambert decomposition over POC, ISS, Ap | |
| Temperature corrections (`theta`) | ✓ `shared/processes.py:16` (`arrhenius_correction`) | ✓ `utils/conversions.py:13` (`arrhenius_correction`) | Identical formula |

### 3.4 v2 NSM1 defects, TODOs, and suspected bugs

The streaming branch already documented the most consequential issues in `tests/test_5_nitrogen_calculations_v2.py:6-12`. Confirmed by inspection:

1. **Nitrogen integrator multiplies instead of adds** (`nitrogen.py:101`):
   ```python
   ammonium = 0 + ammonium * ammonium_rate * self.time_step.total_seconds()
   ```
   The intent is `ammonium_new = ammonium_old + rate * dt`. The current line computes `ammonium_old * rate * dt`, which is dimensionally wrong (rate is [1/s] but the result is [mg-N/L * 1/s * s] = [mg-N/L] only by coincidence) and physically wrong (multiplicative update). Same defect at line 115 for nitrate. **Critical, blocks any production use of v2 nitrogen.**
2. **Floating-algae integrator multiplies instead of adds and includes a stray `* 86400`** (`floating_algae.py:122`):
   ```python
   algae = 0 + algae * rate * self.time_step.total_seconds() * 86400
   ```
   `time_step.total_seconds()` is already in seconds; multiplying by 86400 amounts to multiplying by seconds-per-day a second time. **Critical.**
3. **NaN comparisons that never fire** (`nitrogen.py:256`, `nitrogen.py:319`, `floating_algae.py:274`, `floating_algae.py:300`, `benthic_algae.py:202`):
   ```python
   rate = xr.where(rate == np.nan, 0.0, rate)
   ```
   `x == np.nan` is always False in IEEE 754. The intended behavior is `xr.where(np.isnan(rate), 0.0, rate)` or `rate.fillna(0.0)`. **Bug; the safety net is silently absent.**
4. **`mixing_ratio_air` returns 0 unconditionally for non-NaN input** (`temperature.py:511`):
   ```python
   if atmospheric_vapor_pressure == atmospheric_vapor_pressure:
       return 0.0
   ```
   The intended NaN guard is inverted. Sensible heat flux and downstream evaporation will be wrong if this path is hit. **Bug.**
5. **`half_saturation_oxygen=1` hard-coded** in `Nitrogen.change_nitrate` (`nitrogen.py:197`); flagged with `# TODO: need argument`.
6. **`algea_growth_rate=0` hard-coded** in `Nitrogen.change_nitrate` for `nitrate_uptake_floating_algae` (`nitrogen.py:211`) and `nitrate_uptake_benthic_algae` (`nitrogen.py:217`). With these zero, nitrate uptake by algae is silently absent, even when the user has FloatingAlgae and BenthicAlgae processes registered.
7. **Death rate in `Nitrogen.__init__` is a placeholder** (`nitrogen.py:60`): `# TODO: this should come from floating algae process`.
8. **`limit_phosphorus == 1.0` test treated as harmonic-mean shortcut** (`floating_algae.py:211`): `# TODO: confirm this 1`. The harmonic-mean denominator includes `1/limit_phosphorus`, so a guard of `== 1.0` is unrelated to the actual divide-by-zero risk (which is `limit_phosphorus == 0`).
9. **`-1` factor TODO on Richardson number** (`temperature.py:597`): `# TODO: check original equation to see if this multiplication by negative one is needed (not in v1 of code)`.
10. **`/ 0.5` mystery factor in `flux_sediment`** (`temperature.py:285`): `# TODO: determine why we need this 0.5`. The legacy code at `processes.py:411` carries the same `/ 0.5` with comment "86400 converts the sediment thermal diffusivity from units of m^2/d to m^2/s" -- the `/ 0.5` is unrelated to that conversion and represents the sediment active layer half-thickness factor.
11. **Sentinel-value parameters in legacy `nsm1/constants.py`** that v2 inherits implicitly through the NSM1 reference: `vsop = 999`, `vs = 999`, `SOD_20 = 999`, `SOD_theta = 999`, `kaw_20_user = 999`, `kah_20_user = 999`. These are documented in Repo A's `docs/NSM1_DOX_rate_bug_investigation.md` (commit `b5f2c8b`). `SOD_theta = 999` produces catastrophic Arrhenius blowup (`999^(T-20)`) at T > 20 °C if the user does not override the default. **v2 does not currently expose any DOX or SOD parameters, so the bug cannot be exercised in v2 today, but any port from legacy NSM1 must address this.**

## 4. Shared infrastructure comparison

| Subsystem | Legacy `clearwater_modules` | `clearwater_modules_v2` |
|---|---|---|
| `base.py` (Model framework) | `Model` class with `register_variable`, `_iter_computations_fast`, hotstart, wet-mask, compute-plan caching (Repo B); 850+ lines | `Model` class with `__init__`, `validate`, `run`, `__process_loop_chunked`, `__process_loop_full`; 130 lines, no caching |
| `sorter.py` (dependency resolution) | Topological sort over `Variable.process.__annotations__` (`sorter.py:64`) | Not present; v2 processes declare an explicit `variables: list[str]` and the user provides processes in execution order in the YAML or constructor |
| `utils.py` (xarray utilities) | `validate_arrays`, `_prep_inputs`, `iter_computations` | Not present |
| `shared/processes.py` | Bank of physics primitives: `arrhenius_correction`, `celsius_to_kelvin`, `kah_20` (9 reaeration options), `kaw_20` (13 wind reaeration options), `ka_tc`, `SOD_tc`, `L` (light extinction), `PAR`, `fdp`, plus several numba-decorated kernels | Not present in v2; `utils/conversions.py` carries only `celsius_to_kelvin` and `arrhenius_correction` (35 lines) |
| `shared/variables.py` | `Variable` registrations shared across modules (currently only `arrhenius_correction`) | Not present |
| `shared/types.py` | `Variable` dataclass, `SplitVariablesDict`, `Process`, `InitialVariablesDict`, `VariableTypes` | Not present (`clearwater_data.variables.VariableRegistry` is the substitute) |
| `tsm/constants.py` | `Meteorological`, `Temperature` `TypedDict`s with defaults (Repo B adds `q_net_depth_ramp_ref`, `dTdt_max_per_hour`) | `utils/constants.py` has 4 module-level constants only (`GRAVITY`, `STEFAN_BOLTZMANN`, `EMISSIVITY_WATER`, `AIR_SPECIFIC_HEAT`) |
| `nsm1/constants.py` | 13 `TypedDict`s with defaults for algae, alkalinity, balgae, carbon, CBOD, DOX, nitrogen, POM, N2, phosphorus, pathogen, global parameters, global vars (~250 entries) | None (parameters declared per-process as `__init__` arguments) |
| Configuration validation | None beyond TypedDict | None (raw `dict[str, Any]` from YAML) |
| Hotstart | Yes (Repo B) | No |
| Wet-mask | Yes (Repo B) | No (v2 mask is `xr.where(volume > 0, dT, 0)` after the kinetic step) |

In short: v2 reimplements the model orchestrator at a coarser grain (process-level, not variable-level), foregoing the dependency-sort, the rich variable metadata, the parameter TypedDicts, the constants library, the shared physics primitives, the cached compute plan, the hotstart support, the wet-mask gating, and the numba accelerations. v2 inherits only the `arrhenius_correction` and `celsius_to_kelvin` helpers in `utils/conversions.py`.

## 5. What needs to be improved in `clearwater_modules_v2` to make TSM and NSM1 complete and optimized

This is the punch list. Each item lists the gap, the source files in both repos to consult, the estimated complexity (S/M/L), and whether it is a port-from-A-to-B or a new design.

### 5.1 Bugs that block production use (Critical)

| # | Item | Files | Complexity | Type |
|---|---|---|---|---|
| 1 | Fix `Nitrogen.run` integrator (multiplicative -> additive). | `clearwater_modules_v2/processes/nitrogen.py:101,115` | S | Port (legacy uses additive `Ap = Ap_old + dApdt * dt`) |
| 2 | Fix `FloatingAlgae.run` integrator (drop the stray `* 86400`; switch to additive). | `clearwater_modules_v2/processes/floating_algae.py:122` | S | Same |
| 3 | Replace `xr.where(rate == np.nan, ...)` with `xr.where(np.isnan(rate), ...)` everywhere. | `nitrogen.py:256,319`; `floating_algae.py:274,300`; `benthic_algae.py:202` | S | New |
| 4 | Fix `Temperature.mixing_ratio_air` NaN guard (currently returns 0 for any non-NaN input). | `temperature.py:511` | S | New |
| 5 | Port the latent-heat unit fix from `clearwater_modules-streaming/src/clearwater_modules/tsm/processes.py:213` to `clearwater_modules_v2/processes/temperature.py:445` (do not convert to Kelvin before applying the polynomial). | `temperature.py:445-453` | S | Port |
| 6 | Remove debug `print(float(...))` calls in `Temperature` that crash on multi-cell arrays. | `temperature.py:167-179, 218-232, 595-607` | S | New |
| 7 | Hardcoded `half_saturation_oxygen=1` and `algea_growth_rate=0` in `Nitrogen.change_nitrate` -- replace with parameter or coupled value. | `nitrogen.py:197, 211, 217` | S | New |

### 5.2 Missing TSM features (Major)

| # | Item | Files | Complexity | Type |
|---|---|---|---|---|
| 8 | Port the depth-ramp + rate-cap thin-water stability guard. | Source: `clearwater_modules-streaming/src/clearwater_modules/tsm/processes.py:445-523`; target: `clearwater_modules_v2/processes/temperature.py:380-425` | M | Port |
| 9 | Add hotstart support (load state from `xr.Dataset` checkpoint and resume). | Source: `clearwater_modules-streaming/src/clearwater_modules/tsm/model.py:38-60`; target: v2 has no `Model.from_dataset` | M | Port (but v2's Model contract is different) |
| 10 | Add wet-mask gating to `Process.run` so dry cells produce NaN-honest output. | Source: `clearwater_modules-streaming/src/clearwater_modules/base.py:564-680`; target: `clearwater_modules_v2/processes/base.py` and `temperature.py` | M | Port |
| 11 | Decide on Richardson `-1` factor convention (TODO at `temperature.py:597`) and verify against v1 (`tsm/processes.py:150`). | `temperature.py:569-602`; `tsm/processes.py:150-168` | S | Verification + 1-line fix |
| 12 | Document and unit-test the `/ 0.5` factor in `flux_sediment` (sediment active-layer half-thickness). | `temperature.py:271-289`; `tsm/processes.py:387-415` | S | Documentation |

### 5.3 Missing NSM1 constituents (Critical for production NSM1)

| # | Item | Files | Complexity | Type |
|---|---|---|---|---|
| 13 | OrgN process module (hydrolysis, settling, algal death contribution). | Source: `nsm1/processes.py:1173-1405`; target: new `clearwater_modules_v2/processes/nitrogen.py::OrganicNitrogen` or new file | M | Port |
| 14 | TIP / OrgP / DIP partitioning + settling. | Source: `nsm1/processes.py:1833-2168`; target: new `clearwater_modules_v2/processes/phosphorus.py` | M | Port |
| 15 | POC / DOC / DIC carbon module (hydrolysis, oxidation, atmospheric CO2 reaeration, photosynthesis/respiration coupling, sediment release). | Source: `nsm1/processes.py:2439-2858`; target: new `clearwater_modules_v2/processes/carbon.py` | L | Port |
| 16 | CBOD module (multi-group oxidation and sedimentation). | Source: `nsm1/processes.py:2334-2422`; target: new `clearwater_modules_v2/processes/cbod.py` | M | Port |
| 17 | DOX module (saturation, atmospheric reaeration with 13 wind/9 hydraulic options, photosynthesis/respiration, nitrification, DOC and CBOD oxidation, SOD). | Source: `nsm1/processes.py:2876-3123` plus `shared/processes.py:65-200`; target: new `clearwater_modules_v2/processes/dox.py` | L | Port |
| 18 | POM module (settling, dissolution, benthic algae mortality, burial). | Source: `nsm1/processes.py:2185-2317`; target: new `clearwater_modules_v2/processes/pom.py` | M | Port |
| 19 | Pathogen (PX) decay + settling. | Source: `nsm1/processes.py:3141-3227`; target: new `clearwater_modules_v2/processes/pathogen.py` | S | Port |
| 20 | Alkalinity (Alk) module. | Source: `nsm1/processes.py:3246-3435`; target: new `clearwater_modules_v2/processes/alkalinity.py` | M | Port |
| 21 | N2 / TDG module. | Source: `nsm1/processes.py:3452-3540`; target: new `clearwater_modules_v2/processes/n2_tdg.py` | M | Port |
| 22 | Reaeration menu (`kah_20` 9 options, `kaw_20` 13 options, `ka_tc`). | Source: `shared/processes.py:65-200`; target: new `clearwater_modules_v2/utils/reaeration.py` | M | Port |
| 23 | Light extinction `L` (Beer-Lambert with ISS, POC, Chl-a contributions) and `PAR`. | Source: `shared/processes.py:202-273`; target: shared utility in v2 | S | Port |
| 24 | TIP partitioning `fdp` and proper coupling to floating-algae phosphorus uptake (replace hard-coded 0.5). | Source: `shared/processes.py:257-272`; target: `floating_algae.py:113` | S | Port |

### 5.4 Missing process completeness within constituents v2 already has (Major)

| # | Item | Files | Complexity | Type |
|---|---|---|---|---|
| 25 | Couple algal NH4 preference and N uptake fractions to actual algal growth (currently `algea_growth_rate=0` placeholders). | `nitrogen.py:211, 217`; `floating_algae.py:401-408` | M | New |
| 26 | Implement `floating_algae.ammonium_respiration()` and `ammonium_growth()` (currently return 0). | `floating_algae.py:401-408` | S | Port |
| 27 | Implement `benthic_algae.ammonium_respiration()` and `ammonium_growth()`. | `benthic_algae.py` (not yet defined) | S | Port |
| 28 | Add benthic algae burial / fraction-to-water terms (Fw, Fb). | Source: `nsm1/constants.py:104` + `nsm1/processes.py:1086`; target: new in `benthic_algae.py` | S | Port |
| 29 | Add light limitation harmonic mean for benthic algae (currently only multiplicative + limiting). | `benthic_algae.py:90-128` | S | Port |
| 30 | Sediment diagenesis (Di Toro multi-G) -- not in NSM1 itself but in NSM2; if planned, this is a major new module. | Out of scope for NSM1; reference `clearwater_modules/nsm2/sediment_flux/` | L | New design (not a port) |

### 5.5 Performance and architectural debt (Major)

| # | Item | Files | Complexity | Type |
|---|---|---|---|---|
| 31 | Port the 418x kernel optimization concept to v2 once v2 has enough variables to need it. | Source: `clearwater_modules-streaming/src/clearwater_modules/base.py:480-565` (`_build_compute_plan`, `_iter_computations_fast`); target: v2 `Model.run` | L | Port (concept), but v2's data model differs |
| 32 | Adopt the variable-level dependency sort instead of process-level. v2 currently relies on the user supplying processes in correct order in the YAML; this does not scale beyond a half-dozen processes. | Source: `clearwater_modules/sorter.py`; target: new in `clearwater_modules_v2` | L | Port |
| 33 | Numba-decorate the inner kernels (Arrhenius, Beer-Lambert, reaeration). | `utils/conversions.py`, future `utils/reaeration.py` | S | Port |
| 34 | Validate config against a schema (pydantic, dataclass, or jsonschema) instead of `dict[str, Any]`. | `config/init.py:209-238` | M | New |

### 5.6 Test coverage gaps (Major)

| # | Item | Files | Complexity | Type |
|---|---|---|---|---|
| 35 | Add full-step integration test for v2 Nitrogen (current parity tests cover sub-rates only; the integrator bug at lines 101/115 is not caught by sub-rate tests). | `tests/test_5_nitrogen_calculations_v2.py` | S | New |
| 36 | Add `def test_*` functions to `tests/test_5_tsm_calculations_v2.py` (file currently contains only fixture data, 0 tests). | `tests/test_5_tsm_calculations_v2.py` | S | New |
| 37 | Add v2 latent-heat regression test analogous to legacy `test_tsm_latent_heat.py`. | `tests/test_tsm_latent_heat.py` (exists for legacy); add v2 analogue | S | Port |
| 38 | Add v2 thin-water stability tests analogous to `test_tsm_stability_ramp.py`. | Same | S | Port |
| 39 | Add hotstart roundtrip test for v2 once hotstart is implemented (item 9). | `test_hotstart_roundtrip.py` | M | Port |
| 40 | Add v2 conservation tests (mass conservation across coupled processes). Legacy has none either; this is a new requirement. | New | M | New design |
| 41 | Add MMS or analytical test for TSM thin-water response. | New | M | New design |

### 5.7 Configuration and documentation gaps (Minor to Major)

| # | Item | Files | Complexity | Type |
|---|---|---|---|---|
| 42 | The example YAML (`config/example_config.yml`) only configures `riverine` and `temperature`. Add example blocks for nitrogen, floating algae, benthic algae once stable. | `config/example_config.yml` | S | New |
| 43 | The YAML uses both `start_time`/`end_time` (in `example_config.yml`) and `start_datetime`/`end_datetime` (in `init.py`). Reconcile. | `config/init.py:29-31`; `config/example_config.yml:6-10` | S | New |
| 44 | `nutrients/__init__.py` is a 0-byte placeholder; remove or populate. | `processes/nutrients/__init__.py` | S | New |
| 45 | Generate v2 docstrings into Sphinx and align with legacy module's `docs_sphinx/`. | `docs_sphinx/` | M | New |
| 46 | Document the `/ 0.5` sediment-thickness convention in `flux_sediment`. | `temperature.py:285` | S | New |
| 47 | The v2 `Process.from_config` contract is inconsistent: some classes register with `@ProcessFactory.register("name")` decorating a `staticmethod` (`temperature.py:84`), some bind on `__init__` (`riverine.py:53`). Standardize. | `processes/*.py` | M | New |

### 5.8 Coupling-readiness gaps (Major if v2 is intended for coupled simulation)

| # | Item | Files | Complexity | Type |
|---|---|---|---|---|
| 48 | The v2 `Riverine` process registers state variables it expects from the riverine model (`Ap, NH4, NO3, TIP, DOX`) but does **not** register OrgN, OrgP, POC, DOC, DIC, POM, CBOD, PX, Alk, N2 from the riverine mesh. Once those constituents are ported (items 13-21), the riverine wrapper must be extended. | `processes/riverine.py:79-108` | M | New |
| 49 | v2 has no equivalent of the `update_state_values` override in the legacy `increment_timestep`. If the riverine transport update is supposed to overwrite v2 state, this needs an explicit API. | New | M | New design |
| 50 | The streaming branch's new SSM (sediment transport, SEDZLJ port) lives in `clearwater_modules_v2/processes/sediment/` and has its own bed-state model. If NSM1 in v2 is to use SSM-derived bed concentrations for benthic flux, a coupling contract is needed. | New | L | New design |

## 6. Recommendation

### 6.1 Where to invest

**Recommendation: option (b), bring v2 features into A, with a small carve-out for v2 sediment.** The legacy `clearwater_modules` package on the streaming branch is currently the only viable production path for TSM and NSM1. It carries (i) every NSM1 constituent, (ii) the 418x kernel optimization, (iii) hotstart support, (iv) wet-mask gating, (v) the latent-heat unit fix, (vi) the thin-water stability guard, and (vii) ~430 passing tests. `clearwater_modules_v2`, by contrast, has a working TSM (with several known defects), a partially working nitrogen process (with a critical integrator bug), and stubs for floating and benthic algae. Porting the legacy package's 16 NSM1 constituents into v2 is a multi-month effort that would duplicate ~3500 lines of well-tested code into a framework that is itself missing the dependency sort, the parameter TypedDicts, the cached compute plan, and the rich variable metadata. Conversely, the architectural ideas v2 is exploring (process-level composition, YAML configuration, per-process substepping) can be retrofitted onto the legacy package incrementally without throwing away the working NSM1.

The exception is the SSM (sediment transport) module that landed on the streaming branch in `clearwater_modules_v2/processes/sediment/` (commit `5674cbc`). That work is genuinely new and need not be ported.

### 6.2 Smallest unit of work for an end-to-end TSM+NSM1 in v2

If the v2 path is nonetheless taken, the minimum end-to-end TSM+NSM1 in v2 requires, in order:

1. Items 1, 2, 5, 6, 7 (5 small fixes: nitrogen and algae integrator additivity; latent-heat unit fix; debug-print removal; nitrogen parameter wiring). One day.
2. Items 13, 14, 16, 17, 22, 23, 24 (port OrgN, P module, CBOD, DOX, reaeration menu, light extinction, TIP partitioning). Two to three weeks.
3. Item 25 (couple algal growth to nutrient uptake). Three days.
4. Items 35, 36, 37 (test coverage for the integrator and TSM). One week.

That is roughly 4-6 weeks of focused work to a credible v2 TSM + nitrogen + phosphorus + carbon-demand + DOX system. POC, DOC, DIC, POM, PX, Alk, N2/TDG, and full benthic algae are additional weeks each.

### 6.3 Punch-list priority order

1. **Block 1 (1 day): Critical fixes in v2.** Items 1, 2, 3, 4, 5, 6 (integrator bugs, NaN guards, latent-heat unit fix, debug-print removal). Without these, v2 cannot give correct answers for any case it nominally supports.
2. **Block 2 (1 week): TSM hardening in v2.** Items 8, 11, 37, 38 (thin-water guard, Richardson sign verification, latent-heat regression test, stability tests).
3. **Block 3 (2-3 weeks): Core nutrient kinetics in v2.** Items 13, 14, 24, 25 (OrgN, phosphorus module + partitioning, algal-growth coupling). After this v2 has a credible C/N/P/algae system but no DOX or CBOD.
4. **Block 4 (3-4 weeks): Oxygen and BOD in v2.** Items 16, 17, 22 (CBOD, DOX, reaeration menu). After this v2 can simulate a complete oxygen budget.
5. **Block 5 (2-3 weeks): Carbon and POM in v2.** Items 15, 18, 23 (POC/DOC/DIC, POM, light extinction).
6. **Block 6 (2-3 weeks): Coupling, hotstart, performance.** Items 9, 10, 31, 32, 49 (hotstart, wet-mask, kernel optimization concept, dependency sort, transport override API).
7. **Block 7 (2 weeks): Remaining constituents.** Items 19, 20, 21 (Pathogen, Alkalinity, N2/TDG).

If the choice is option (b) (consolidate on the legacy package), the only essential carry-over from v2 is the YAML configuration system (items 42, 43, 47), which is roughly two weeks of work to bolt onto the legacy package.
