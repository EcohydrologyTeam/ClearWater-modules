[![Tests Status](https://github.com/EcohydrologyTeam/ClearWater-modules/actions/workflows/tests.yml/badge.svg)](https://github.com/EcohydrologyTeam/ClearWater-modules/actions/workflows/tests.yml)
[![Coverage](https://codecov.io/gh/EcohydrologyTeam/ClearWater-modules/graph/badge.svg)](https://codecov.io/gh/EcohydrologyTeam/ClearWater-modules)

# ClearWater Modules in Python

The [ClearWater-modules](https://github.com/EcohydrologyTeam/ClearWater-modules) package is a collection of water quality and vegetation process simulation modules written in modern Python and designed to flexibly couple with a variety of water transport models, such as HEC-RAS-2D, GSSHA, CE-QUAL-W2, [AdH](https://www.erdc.usace.army.mil/Locations/CHL/AdH/), and others. The U.S. Army Engineer Research and Development Center (ERDC), Environmental Laboratory (EL) develops these modules.

- [TSM: Temperature Simulation Module](src/clearwater_modules/tsm) (formerly TEMP)
- [NSM: Nutrient Simulation Modules](src/clearwater_modules/nsm1) ([NSM-I](src/clearwater_modules/nsm1) and [NSM-II](src/clearwater_modules/nsm2))
- [GSM: General Constituent Simulation Module](src/clearwater_modules/gsm)
- [CSM: Contaminant Simulation Module](src/clearwater_modules/csm)
- [MSM: Mercury Simulation Module](src/clearwater_modules/msm)
- SSM: Solids Simulation Module (Fortran only)
- RVSM: Riparian Vegetation Simulation Module (Fortran only)

These water quality modules form the central capabilities of the [ClearWater (Corps Library for Environmental Analysis and Restoration of Watersheds)](https://ui.adsabs.harvard.edu/abs/2023EGUGA..2512470S/abstract) software system. The overall goal of the ClearWater system is to couple these water quality simulation capabilities to state-of-the-art hydrologic and hydraulic modeling tools, such as HEC-RAS-2D, CE-QUAL-W2, and GSSHA, allowing users to leverage existing river, reservoir, and watershed models for water quality studies. The Temperature Simulation Module (TSM) and Nutrient Simulation Module (NSM) have been successfully coupled to HEC-RAS-2D models via the [ClearWater-riverine](https://github.com/EcohydrologyTeam/ClearWater-riverine) package.

A secondary goal is to develop a suite of modern Python tools that build on community-developed scientific workflows, standards, and libraries to automate model setup, prepare input datasets, store output data, and visualize results using Python-based user interfaces such as Jupyter Notebooks.

This Python library is a port and modernization of the algorithms and structures originally written in Fortran 95, released as version 1.0 in 2021, and described in:

- Zhang, Zhonglong and Billy E. Johnson. 2016. Aquatic nutrient simulation modules (NSMs) developed for hydrologic and hydraulic models. Vicksburg, MS: Environmental Laboratory, U. S. Army Engineer Research and Development Center (ERDC). Ecosystem Management and Restoration Research Program (EMRRP). ERDC/EL Technical Report 16-1. https://hdl.handle.net/11681/10112
- Zhang, Zhonglong and Billy E. Johnson. 2016. Aquatic contaminant and mercury simulation modules developed for hydrologic and hydraulic models. Vicksburg, MS: Environmental Laboratory, U. S. Army Engineer Research and Development Center (ERDC). Environmental Quality Technology Research Program (EQTRP). ERDC/EL Technical Report 16-8. https://hdl.handle.net/11681/20249
- Johnson, Billy E. and Zhonglong Zhang. 2016. Testing and Validation Studies of the NSMII-Benthic Sediment Diagenesis Module. Vicksburg, MS: Environmental Laboratory, U. S. Army Engineer Research and Development Center (ERDC). Ecosystem Management and Restoration Research Program (EMRRP). ERDC/EL Technical Report 16-11. https://hdl.handle.net/11681/20343

## Streaming-branch additions

This is the `streaming` branch. The branch carries a set of additions on top of the public `main` that target real-world coupling with two-dimensional transport drivers (in particular ClearWater-Riverine on HEC-RAS-2D) and that improve numerical robustness, hotstart support, and per-timestep performance. The sections below document those additions; the rest of the README continues to apply unchanged.

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

## Module maturity at a glance

| Module | Path | Status |
|--------|------|--------|
| TSM | `src/clearwater_modules/tsm/` | Production-ready; recent stability and latent-heat fixes |
| NSM1 | `src/clearwater_modules/nsm1/` | Production-ready; ~150 dynamic variables; hotstart support |
| NSM2 | `src/clearwater_modules/nsm2/` | Production; submodules for algae, alkalinity, benthic algae, carbon, CBOD, dissolved oxygen, N2, nitrogen, pathogens, phosphorus, particulate organic matter, and sediment flux |
| v2 modules | `src/clearwater_modules_v2/` | Experimental; refactored temperature, benthic algae, floating algae, and nitrogen with v1-parity tests passing |
| GSM | `src/clearwater_modules/gsm/` | Stub; limited process implementation |
| CSM | `src/clearwater_modules/csm/` | Stub; no dynamic-variable registration |
| MSM | `src/clearwater_modules/msm/` | Stub; no process definitions |

## Variable and process registration

The base `Model` design uses a decorator, `@register_variable(models=...)`, to wire process functions into a `Model` class. Each variable is declared with one of three roles: `static` (a parameter or constant supplied at construction), `state` (a quantity the model integrates over time), or `dynamic` (a per-step diagnostic computed from other variables). The topological sort in `src/clearwater_modules/utils.py` and `src/clearwater_modules/sorter.py` resolves dependency order so that each step's kernel runs each registered process once, in an order that satisfies its argument dependencies. To add a new process, you write a Python function whose argument names match registered variable names, then register the variable and bind the function as its `process` attribute. The compute plan is rebuilt automatically when you register or unregister a variable.

## Tests

The branch ships 430 passing tests with 2 expected-failure cases that hold the pre-latent-heat-fix expected-value baseline. Continuous integration runs on Windows, Linux, and macOS (commit `ae2c5d5`). Notable test modules include:

- `test_5_tsm_calculations.py` (15 rebaselined expected values, energy-balance regression)
- `test_tsm_latent_heat.py` (9 tests, Kelvin-to-Celsius regression guard)
- `test_tsm_stability_ramp.py` (6 tests covering the depth ramp and rate cap)
- `test_hotstart_roundtrip.py` (4 cases, NSM1 + TSM, `rtol=1e-12`)
- `test_5_benthic_algae_calculations_v2.py`, `test_5_floating_algae_calculations_v2.py`, `test_5_nitrogen_calculations_v2.py` (v1-parity for the v2 modules)

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
