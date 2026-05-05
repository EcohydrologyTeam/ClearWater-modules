[![Tests Status](https://github.com/EcohydrologyTeam/ClearWater-modules/actions/workflows/tests.yml/badge.svg)](https://github.com/EcohydrologyTeam/ClearWater-modules/actions/workflows/tests.yml)
[![Coverage](https://codecov.io/gh/EcohydrologyTeam/ClearWater-modules/graph/badge.svg)](https://codecov.io/gh/EcohydrologyTeam/ClearWater-modules)

# ClearWater Modules in Python

The [ClearWater-modules](https://github.com/EcohydrologyTeam/ClearWater-modules) package is a collection of water quality and vegetation process simulation modules written in modern Python and designed to flexibly couple with a variety of water transport models, such as HEC-RAS-2D, GSSHA, CE-QUAL-W2, [AdH](https://www.erdc.usace.army.mil/Locations/CHL/AdH/), and others. The U.S. Army Engineer Research and Development Center (ERDC), Environmental Laboratory (EL) develops these modules.

- [TSM: Temperature Simulation Module](src/clearwater_modules/tsm) (formerly TEMP)
- [NSM1: Nutrient Simulation Module I](src/clearwater_modules/nsm1)
- [GSM: General Constituent Simulation Module](src/clearwater_modules/gsm)
- SSM: Solids Simulation Module (Fortran only)
- RVSM: Riparian Vegetation Simulation Module (Fortran only)

NSM-II, CSM (Contaminant Simulation Module), and MSM (Mercury Simulation Module) prototypes that previously lived under `src/clearwater_modules/{nsm2,csm,msm}/` were removed on this branch. Future contaminant- and mercury-process work is expected to go through the PhreeqcRM coupling path rather than reviving the prototype Python ports; benthic sediment diagenesis (the substantive content of the former NSM-II) is on the v3.x roadmap as a separate process module.

These water quality modules form the central capabilities of the [ClearWater (Corps Library for Environmental Analysis and Restoration of Watersheds)](https://ui.adsabs.harvard.edu/abs/2023EGUGA..2512470S/abstract) software system. The overall goal of the ClearWater system is to couple these water quality simulation capabilities to state-of-the-art hydrologic and hydraulic modeling tools, such as HEC-RAS-2D, CE-QUAL-W2, and GSSHA, allowing users to leverage existing river, reservoir, and watershed models for water quality studies. The Temperature Simulation Module (TSM) and Nutrient Simulation Module (NSM) have been successfully coupled to HEC-RAS-2D models via the [ClearWater-riverine](https://github.com/EcohydrologyTeam/ClearWater-riverine) package.

A secondary goal is to develop a suite of modern Python tools that build on community-developed scientific workflows, standards, and libraries to automate model setup, prepare input datasets, store output data, and visualize results using Python-based user interfaces such as Jupyter Notebooks.

This Python library is a port and modernization of the algorithms and structures originally written in Fortran 95, released as version 1.0 in 2021, and described in:

- Zhang, Zhonglong and Billy E. Johnson. 2016. Aquatic nutrient simulation modules (NSMs) developed for hydrologic and hydraulic models. Vicksburg, MS: Environmental Laboratory, U. S. Army Engineer Research and Development Center (ERDC). Ecosystem Management and Restoration Research Program (EMRRP). ERDC/EL Technical Report 16-1. https://hdl.handle.net/11681/10112
- Zhang, Zhonglong and Billy E. Johnson. 2016. Aquatic contaminant and mercury simulation modules developed for hydrologic and hydraulic models. Vicksburg, MS: Environmental Laboratory, U. S. Army Engineer Research and Development Center (ERDC). Environmental Quality Technology Research Program (EQTRP). ERDC/EL Technical Report 16-8. https://hdl.handle.net/11681/20249
- Johnson, Billy E. and Zhonglong Zhang. 2016. Testing and Validation Studies of the NSMII-Benthic Sediment Diagenesis Module. Vicksburg, MS: Environmental Laboratory, U. S. Army Engineer Research and Development Center (ERDC). Ecosystem Management and Restoration Research Program (EMRRP). ERDC/EL Technical Report 16-11. https://hdl.handle.net/11681/20343

## Streaming-branch additions

This is the `streaming` branch. The branch carries a set of additions on top of the public `main` that target real-world coupling with two-dimensional transport drivers (in particular ClearWater-Riverine on HEC-RAS-2D) and that improve numerical robustness, hotstart support, and per-timestep performance. The branch hosts three coexisting tracks:

1. **`clearwater_modules` (v1)** — the function-style framework hardened in April–May 2026 with the latent-heat unit fix, thin-water stability guard, hotstart-from-`xr.Dataset`, multi-cell-safe debug-print removal, wet-mask gating, and the 418× kernel optimization. Documented in the sections below.
2. **`clearwater_modules_v2`** — the class-based framework authored by LimnoTech (Paul Tomasula, Anthony Aufdenkampe, Jason Rutyna, Sarah Jordan), mirrored from `upstream/memory-refactor-pytestUpdate`. v2 introduces YAML-driven configuration via `init_from_file`, per-process substepping, and a chunking execution path. It is the architectural baseline for v3.
3. **`clearwater_modules_v3`** — the convergence of v1 and v2 into a single coherent codebase. v3 keeps v2's framework and adds v1's optimization and correctness work. v3 development happens entirely on this branch; the v3 1.0.0 release will go upstream when the work is ready for LimnoTech review. See the dedicated v3 section below.

The remaining sections document the three tracks separately; the rest of the upstream README continues to apply unchanged.

## TSM stability and correctness

The `streaming` branch lands three independent corrections to the TSM kernel under a single commit. The function signatures of `dTdt_water_c` and `mf_latent_heat_vaporization` did not change, so existing callers continue to work without modification.

### Depth-ramp damping of net surface heat flux

In shallow cells, the explicit-Euler temperature update is numerically stiff: a full surface heat flux delivered into a column with vanishing thermal mass produces per-substep temperature changes that exceed physical bounds. The kernel now multiplies `q_net` by `min(1.0, depth / q_net_depth_ramp_ref)`, where `depth = volume / surface_area`. The ramp smoothly attenuates the flux as depth approaches zero. Setting `q_net_depth_ramp_ref = 0.0` disables the ramp and reproduces legacy behavior bit-exactly. The kernel also clamps depth to zero defensively when `surface_area` or `volume` are non-positive.

### Per-substep rate cap on temperature change

The kernel applies a cadence-invariant rate cap on the magnitude of `dTdt_water_c`. The per-substep cap is `dTdt_max_per_hour * dt_hours`, where `dt_hours` is derived from the substep size `dt`. The cap is a belt-and-suspenders guard against rare residual stiffness; under the depth ramp it should not normally activate. Setting `dTdt_max_per_hour = +inf` disables the cap.

### Latent heat of vaporization unit fix

The function `mf_latent_heat_vaporization` previously applied its polynomial coefficients (intercept 2,499,999 J/kg, slope -2385.74 J/kg/K) to a Kelvin temperature, but the polynomial is calibrated for input in degrees Celsius. The mismatch produced a latent heat roughly 26 to 27 percent too low across the typical surface-water range, which underestimated evaporative cooling and biased simulated water temperatures warm. The fix converts Kelvin to Celsius (`water_temp_c = water_temp_k - 273.15`) before applying the polynomial. The function still accepts Kelvin per the kernel's dynamic-variable registration.

### Parameters

| Parameter | Units | Default | Scope | Description |
|-----------|-------|---------|-------|-------------|
| `q_net_depth_ramp_ref` | m | 0.3 | per-cell static | Reference depth for the smooth flux ramp inside `dTdt_water_c`. The kernel multiplies `q_net` by `min(1.0, depth / q_net_depth_ramp_ref)`. Set to 0.0 to disable the ramp (legacy behavior, bit-exact). |
| `dTdt_max_per_hour` | K hr-1 | 5.0 | per-cell static | Maximum magnitude of `dTdt_water_c` expressed as a rate. The kernel applies a per-substep cap of `dTdt_max_per_hour * dt_hours`. Set to `+inf` to disable. |

You can supply both parameters through the `temp_parameters` dictionary when constructing `EnergyBudget`. The defaults live in `clearwater_modules.tsm.constants.DEFAULT_TEMPERATURE`.

## Coupled-simulation features

The base `Model` class on this branch exposes three features that exist primarily to support coupling with an external transport driver such as ClearWater-Riverine. All three are opt-in and preserve legacy behavior when omitted.

### Hotstart dataset round-trip

The `EnergyBudget` (TSM) and `NutrientBudget` (NSM1) constructors accept the keyword arguments `hotstart_dataset` (an `xarray.Dataset`) and `hotstart_timestep` (an integer). When you pass a hotstart dataset, the constructor skips the from-dicts dataset construction path entirely and forwards the dataset and starting timestep to `base.Model._init_from_dataset`. The model resumes from the supplied state without re-broadcasting initial conditions or static parameters. Tests in `tests/test_hotstart_roundtrip.py` cover both modules with `rtol=1e-12`.

```python
import xarray as xr
from clearwater_modules.tsm.model import EnergyBudget

ds = xr.open_dataset("tsm_checkpoint.nc")
model = EnergyBudget(
    time_steps=24,
    hotstart_dataset=ds,
    hotstart_timestep=0,
)
```

### Per-cell wet-mask gating

The method `Model.increment_timestep` accepts an optional `wet_mask` keyword argument (a Boolean `numpy.ndarray` or `xarray.DataArray`). Cells where `wet_mask == True` run kinetics normally. Cells where `wet_mask == False` are still iterated numerically so the array shapes and the dependency-order kernel stay consistent, but the kernel writes NaN into the new-timestep slot for every temporal variable in those cells. To keep kinetics finite during the iteration over dry cells, the kernel fills any inherited NaNs in the previous-step state from the model's `initial_state_values` defaults. Passing `wet_mask=None` (the default) preserves the legacy all-cells-wet behavior. This feature exists so that NaN-honest output is produced for cells that the transport driver reports as dry on a given step.

### `update_state_values` override

The same `Model.increment_timestep` method also accepts an optional `update_state_values` argument: a dict mapping state-variable names to `xarray.DataArray` values. The kernel applies the overrides to the previous-timestep state before kinetics run. Use this entry point when an external coupler (for example, ClearWater-Riverine transport) holds a more authoritative value for a given state variable than the model's own integration. The method validates that override keys are members of the model's state variables or its `updateable_static_variables` list and that override dimensions match the parent dataset's spatial dimensions.

```python
model.increment_timestep(
    update_state_values={"water_temp_c": transported_temp},
    wet_mask=cell_is_wet,
)
```

## Performance

The hot path of `Model.increment_timestep` was reworked on this branch (commit `6daa65e`) for roughly a 418x per-step speedup. Two changes drive the gain. First, the compute plan (kernel name, callable, and ordered argument names per variable) is built once per Model class and cached on the class; subsequent steps reuse the cached plan. Second, the inner loop runs against a plain `dict[str, np.ndarray]` of buffers, and the kernel writes new-timestep results directly into the underlying numpy buffers via `DataArray.values[slot] = ...` rather than going through `Dataset.__setitem__`, which avoids the per-assignment merge_core/align/reindex machinery that previously dominated runtime on NSM1-sized variable graphs (~150 dynamic variables). The result is essentially flat scaling from 1 to 1000 cells, with a per-step time of roughly 3 ms on a 5-cell mesh after the optimization compared to roughly 1292 ms before.

## v3: convergence of v1 and v2

`src/clearwater_modules_v3/` is the convergence track and the active development target on this branch. v3 keeps v2's class-based framework as the architectural baseline and folds in v1's optimization and correctness work, plus the synthesis required to make them coexist. v3 is being co-developed with LimnoTech as the v2 framework authors, and is intended to eventually supersede both v1 and v2.

### Status (as of 2026-05-04)

- Phases 0 through 4 of the original v3 plan are complete.
- A multi-agent code review on 2026-05-04 produced a punch list of 10 CRITICAL, 18 MAJOR, 19 MINOR findings. After five remediation rounds (Phase R-1 through R-5), all 10 CRITICAL findings are closed; 17 of 18 MAJOR are closed (only M4, the 273.16→273.15 ice-point/triple-point offset, is intentionally deferred until v1-parity tests are decommissioned); 10 of 19 MINOR are closed and the remainder are deferred or out-of-scope with documented rationale.
- 153 v3 tests pass in 0.24 s; v3 imports standalone from `PYTHONPATH=src python -c "import clearwater_modules_v3"`.

### What's new in v3 (relative to v2)

The following capabilities exist only in v3. In two cases (sediment-diffusivity Fortran parity, dynamic sediment-temperature evolution), they are not in v1 either.

- **Latent-heat unit fix.** v2 evaluated the `Lv` polynomial on a Kelvin temperature where the polynomial expects Celsius. The error was approximately 26 % in `Lv` at 20 °C, propagating into the latent-heat flux. v3 evaluates the polynomial in Celsius.
- **Thin-water depth ramp + per-hour `dT/dt` cap.** v3 ramps the net heat flux down as wetted depth approaches the configured `q_net_depth_ramp_ref` and caps `dT/dt` at the configured `dTdt_max_per_hour`. Both guards are required for stability in shallow riverine cells where transient depth approaches zero.
- **Vectorized `mixing_ratio_air` guard.** v2's guard against zero or negative `(P_air − e_air)` was scalar-valued and silently broke on multi-cell `xr.DataArray` inputs. v3's guard uses `xr.where` and works for any shape.
- **Sediment-diffusivity Fortran parity.** The Fortran TSM reference (`HEC-RAS-WQ` `modGlobal.f90`) declares `alphas` in m²/day with default `0.0432`. v2 inherited a transcription error introducing `0.0061` with an inconsistent docstring claiming m²/s. v3 restores the Fortran value and units.
- **Dynamic sediment-temperature evolution.** The Fortran TSM evolves `T_sed` on each substep so that water and sediment exchange identical enthalpy (energy-conservative). v1 and v2 dropped this update; v3 reinstates it under `evolve_sediment_temperature=True` (default).
- **Kernel optimization with precomputed schedule.** v3 `Model` precomputes the per-process firing schedule indexed by integer step number. The schedule is timezone-independent and exact under floating-point time arithmetic (no `start_time.timestamp()` modular drift).
- **Registry-level wet-mask gating.** v3 `Model.__apply_wet_mask` honors `process.output_variables` (the variables a process *writes*) rather than `process.variables` (which includes inputs). Dry-cell forcings such as `wind_speed`, `air_temperature`, and `volume` are preserved across substeps; only outputs are NaN-masked.
- **Hotstart from `xr.Dataset`.** v3 `Model` accepts a `hotstart_dataset` and `hotstart_timestep`, seeds the registry from the dataset slice at that timestep, and offers each `Process` an optional pair of opt-in hooks `to_hotstart() -> dict` / `from_hotstart(state: dict)` for substep-internal state.
- **Chunking aligned to integer step indices.** v3 precomputes `interior_chunk_step_indices: set[int]` from `chunk_size_seconds / time_step_seconds`. Boundary detection uses integer comparison, immune to floating-point drift in repeated `current_time +=` arithmetic and immune to timezone effects.
- **YAML schema extensions.** Two new optional top-level keys, `hotstart` and `wet_mask`. v2 YAML configs without these keys run unchanged on v3 (verified backward compatibility).

### Design and review documentation

Under `design/`:

- `clearwater_modules_v3_architecture_specification.md` — umbrella architecture spec for the v3 package.
- `clearwater_modules_v3_tsm_design_specification.md` — TSM-specific design spec.
- `clearwater_modules_v3_nsm1_design_specification.md` — NSM1-specific design spec (deferred to v3.1).
- `clearwater_modules_v3_tsm_gap_analysis.md` — Phase 0 inventory and v1-vs-v2 diff table.
- `clearwater_modules_v3_review_findings.md` — Multi-agent review findings, Phase R-1 through R-5 resolution log.
- `TSM_NSM1_v1_vs_v2_inventory.md` — v1 vs v2 framework inventory.
- `legacy_modules_validation_status.md` — V&V status for legacy Fortran NSM/TSM/MSM/CSM.

The v3 package's own README is at `src/clearwater_modules_v3/README.md` and covers the user-facing migration path, `What's new`, and YAML schema in detail.

## Module maturity at a glance

| Module | Path | Status |
|--------|------|--------|
| v3 (convergence) | `src/clearwater_modules_v3/` | Active development; v3 `Temperature`, `Model`, `init_from_file` are v3-native; v3 1.0.0 targets TSM-complete + multi-agent-review-resolved (10/10 CRITICAL, 17/18 MAJOR, 10/19 MINOR closed); 153 tests passing |
| TSM (v1) | `src/clearwater_modules/tsm/` | Production-ready; latent-heat fix, thin-water stability guard, hotstart, multi-cell debug-print removal |
| NSM1 (v1) | `src/clearwater_modules/nsm1/` | Production-ready; ~150 dynamic variables; hotstart support |
| v2 modules | `src/clearwater_modules_v2/` | LimnoTech framework, mirrored from `upstream/memory-refactor-pytestUpdate`; refactored Temperature, BenthicAlgae, FloatingAlgae, and Nitrogen with v1-parity tests passing; the architectural baseline for v3 |
| GSM (v1) | `src/clearwater_modules/gsm/` | Stub; limited process implementation |

## Variable and process registration

The base `Model` design uses a decorator, `@register_variable(models=...)`, to wire process functions into a `Model` class. Each variable is declared with one of three roles: `static` (a parameter or constant supplied at construction), `state` (a quantity the model integrates over time), or `dynamic` (a per-step diagnostic computed from other variables). The topological sort in `src/clearwater_modules/utils.py` and `src/clearwater_modules/sorter.py` resolves dependency order so that each step's kernel runs each registered process once, in an order that satisfies its argument dependencies. To add a new process, you write a Python function whose argument names match registered variable names, then register the variable and bind the function as its `process` attribute. The compute plan is rebuilt automatically when you register or unregister a variable.

## Tests

### v1 + v2 (legacy and refactored)

The branch ships 430 passing tests with 2 expected-failure cases that hold the pre-latent-heat-fix expected-value baseline. Continuous integration runs on Windows, Linux, and macOS (commit `ae2c5d5`). Notable test modules include:

- `test_5_tsm_calculations.py` (15 rebaselined expected values, energy-balance regression)
- `test_tsm_latent_heat.py` (9 tests, Kelvin-to-Celsius regression guard)
- `test_tsm_stability_ramp.py` (6 tests covering the depth ramp and rate cap)
- `test_hotstart_roundtrip.py` (4 cases, NSM1 + TSM, `rtol=1e-12`)
- `test_5_benthic_algae_calculations_v2.py`, `test_5_floating_algae_calculations_v2.py`, `test_5_nitrogen_calculations_v2.py` (v1-parity for the v2 modules)

### v3

`tests/v3/` ships 153 passing tests in 0.24 s, organized into roughly four groups: regression-suite ports of the v1 TSM tests against the v3 `Temperature` API (calculations, latent-heat, stability-ramp), v3-specific tests for the merge work (sediment-diffusivity Fortran parity, dynamic sediment-T evolution, hotstart roundtrip), Phase R-1 through R-3 robustness tests for the resolved review findings, and Phase R-4/R-5 wet/dry-transition, NaN-propagation, and non-integer-second-time_step tests. Notable test files include:

- `test_5_tsm_calculations_v3.py` (15 v1-parity TSM calculation tests)
- `test_tsm_latent_heat_v3.py` (12 tests pinning the Lv unit fix)
- `test_tsm_stability_ramp_v3.py` (8 tests covering the depth ramp and rate cap)
- `test_tsm_sediment_v3.py` (10 tests covering Fortran-parity sediment defaults and energy-conservative T_sed evolution)
- `test_hotstart_roundtrip_v3.py` (12 tests covering Temperature `to_hotstart` / `from_hotstart` and Model lifecycle hooks)
- `test_model_orchestration_v3.py`, `test_model_robustness_v3.py`, `test_model_minor_v3.py` (Model orchestration, kernel-optimization, wet-mask, hotstart, chunking, idempotency, validation)
- `test_config_init_robustness_v3.py` (YAML loader error reporting; deep-key paths, hotstart and wet_mask validation)
- `test_v2_helper_contract.py` (pins the v2 helper-function signatures so upstream changes surface as CI failures)
- `test_wet_dry_transition_v3.py`, `test_nan_propagation_e2e_v3.py`, `test_schedule_non_integer_v3.py` (R-4/R-5 coverage)

Run the v3 suite with `pixi shell -e dev` then `pytest tests/v3/`.

## Repository Directories

## Getting Started

### Installation

ClearWater-modules was developed against **Python 3.11** (see `pyproject.toml`).

Follow these steps to install.

#### Production

Production installation instructions are forthcoming. For now, follow the [Developer](#developer) instructions.

#### Developer

##### 1. Install Pixi

We recommend installing [pixi](https://pixi.prefix.dev/latest/), a fast, modern, and reproducible package management tool.

##### 2. Clone the ClearWater family of repositories

Three repositories house the dependencies for this project. Navigate to each repository listed below and follow its instructions to clone it to your local machine:

- [clearwater-modules](https://github.com/EcohydrologyTeam/ClearWater-modules)
- [clearwater-riverine](https://github.com/EcohydrologyTeam/ClearWater-riverine)
- [clearwater-data](https://github.com/EcohydrologyTeam/ClearWater-data)

From the GitHub site, click the green "Code" dropdown near the upper right. Choose either "Open with GitHub Desktop" or "Download ZIP." We recommend GitHub Desktop so you can pull updates easily.

Place your copy of each repository in any convenient location on your computer.

##### 3. Create the Python environment using pixi

Navigate to the directory of the `clearwater-modules` repository you cloned. To install the development environment, run:

```bash
pixi install -e dev
```

To activate the newly created environment, run:

```bash
pixi shell -e dev
```

You should now be able to run the example notebooks and create your own.

## Examples

The `examples/` directory contains Jupyter notebooks and supporting data:

- `Implement_TSM_Example.ipynb`: a worked example that builds and runs the TSM `EnergyBudget` model
- `model_architecture.ipynb`: a walkthrough of the variable and process registration design
- `V2/`: example workflows for the experimental v2 modules

## Contributing

We welcome your pull request.

## Acknowledgements

The vision for modernizing this library, including the initial port to Python from Fortran, was developed by:

- Dr. Todd E. Steissberg (ERDC-EL)

The algorithms and structure of this program were adapted from the Fortran 95 version 1.0 of these modules, originally developed by:

- Dr. Billy E. Johnson (ERDC-EL, LimnoTech)
- Dr. Zhonglong Zhang (Portland State University, LimnoTech)
- Mr. Mark Jensen (USACE HEC)
