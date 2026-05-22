# Changelog

All notable changes to ClearWater-Modules are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/).

---

## [v0.4.0] — 2026-05-22

> 104 commits to `main` by @ptomasula, @jrutyna, @aufdenkampe, and @sjordan29
> spanning 2025-05-14 through 2026-05-05

This release is a ground-up architectural refactor with a custom variable-registry framework, YAML-driven configuration, protocol-decoupled I/O, chunked Zarr output, and tightly coupled integration with [ClearWater-Riverine](https://github.com/EcohydrologyTeam/ClearWater-riverine). A key enabler is the new **[ClearWater-data](https://github.com/EcohydrologyTeam/ClearWater-data)** package that serves as a purpose-built shared data layer (Xarray + Dask + Zarr) that provides zero-copy data exchange across all ClearWater components, such as  ClearWater-Riverine. The refactor is presently introduced with a new package namespace  `clearwater_modules_v2`, which will temporily sit alongside the legacy `clearwater_modules` package for backwards compatibility and testing. In a future release, the legacy code will be fully replaced by the refactored code under a single `clearwater_modules` namespace.

 and ClearWater-data v0.1.0 was released alongside the spring 2026 refactors of both repositories.

---

### ⚠️ Breaking Changes

- **New package namespace: `clearwater_modules_v2`** — the v2 refactor lives in an entirely
  new top-level package. Users of `clearwater_modules` are not affected immediately, but
  migration to `clearwater_modules_v2` is the intended path forward.
- **Process initialization changed to factory pattern** — the previous conditional loop for
  building process instances is replaced by a factory; any code that directly instantiated
  process classes will need to be updated.
- **`time_step` renamed to `time` in process signatures** — process `__call__` signatures
  now use `time` for the timestep argument.
- **YAML config is now the primary entry point** — the model is initialized from a `.yml`
  config file rather than direct constructor arguments.

---

### 🚀 Major New Features

#### `clearwater_modules_v2`: New Architecture
- A **custom variable-registry framework** inspired by xarray-simlab but without the
  variable-intent system; all state is accessed through a shared `VariableRegistry`.
- **Single shared registry instance** between ClearWater-Modules and ClearWater-Riverine,
  enabling tight coupling without data copying.
- **Space dimension support** propagated to all `DataArrayVariable` objects, aligned with
  the riverine mesh.
- New top-level `model.py` orchestrates process execution, chunking, and output.

#### YAML Configuration System
- Full model setup from a `.yml` config file (`d3ed447`, `c679601`).
- `config/` subpackage: `read.py`, `init.py`, `example.py`, and `example_config.yml`.
- CSV data provider reads boundary/forcing data and writes intermediate input files.
- `__main__.py` entry point for running the model from the command line.

#### Protocol-Decoupled I/O via ClearWater-data

- I/O layer refactored around Python `Protocol` classes, decoupling data sources from
  the model (`a3f6b55`).
- Integrated with the new **[ClearWater-data](https://github.com/EcohydrologyTeam/ClearWater-data)**
  package (`e392bc8`) — created specifically to provide **zero-copy data exchange** between
  ClearWater-Modules, ClearWater-Riverine, and future ClearWater components. It combines
  Xarray (array data model), Dask (lazy/parallel computation), and Zarr (cloud-native chunked
  storage) to scale simulations to millions of grid cells and time steps. ClearWater-data was
  first integrated into ClearWater-Riverine, with ClearWater-Modules adopting the same
  protocols as part of this refactor.
- Output loop with **lazy Zarr writes** (`fef304da`, `931e809`), building on the
  ClearWater-data Zarr storage layer.
- `read zarr` added to getting-started example notebook.

#### Chunked Execution
- **Chunked simulation mode** wired up to match ClearWater-Riverine's chunk architecture
  (`d712c59`).
- Chunk config naming convention aligned with riverine (`3f11bf0`).
- Model uses variable space-dimension information for chunk handling (`4de7c01`).
- Logic added to skip the first TSM step to align process timesteps with riverine
  (`4b1a5f9`).
- Iterative script for running coupled models across different chunk sizes and run
  lengths as part of scalability profiling (#98, `b86c16a`).

#### New and Refactored Processes

**Temperature Simulation Module (TSM) v2**
- Fully rewritten `temperature.py` in `clearwater_modules_v2/processes/`.
- **Richardson number** rework: switched from scalar `if` logic to `xarray.where`
  calls to support array inputs (`d12bb09`).
- **Flipped temperature flux equations fixed**: atmospheric and upwelling longwave flux
  equations were swapped; corrected (`5d8321a`).
- **Optional sediment temperature flag** added (`--use_sediment_temperature`); new pytest
  case implemented for this path (`85f24cf`, `c49a108`).
- Fixed `xarray.where` usage in temperature to work with non-Dataset array-like types
  (`60d08b7`).
- **TSM v2 validated against TSM v1**: first 14 pytest cases from v1 now pass in v2
  (`26ab3ba`, `634246d`).

**Floating Algae**
- New `floating_algae.py` process class added to `clearwater_modules_v2/processes/`
  (`2fa8067`).

**Other Processes**
- `nitrogen.py`, `benthic_algae.py`, `nutrients/` subpackage, and `riverine.py` process
  stubs added to the v2 processes directory.
- `base.py` defines the shared `Process` base class.

#### Configurable Timestep Frequency
- Processes can now run on **decoupled timesteps** — a frequency setting controls how
  often each process fires relative to the model's global timestep (`bd0c522`).

#### Coupled Riverine + TSM Simulation
- Full **linked Riverine/TSM v2 simulation** demonstrated end-to-end (`125dc7b`).
- Fine-mesh coupled simulation working (`91b8e5b`, `d1f900d`).
- Riverine **dynamic plot** integrated with linked simulation output (`21a48b9`).
- Logic to skip first TSM step ensures time alignment between the two models (`209b67f`).

---

### ✨ Improvements

#### Relative Pathing & Configuration
- All file paths resolved **relative to the notebook working directory**, fixing coupling
  bugs when `modules.yml` is called from within riverine runs (related to riverine #133).
- `project_path` scoped to the model, not the repo root (`0d51a18`).
- Example data moved into `examples/data/` subtree; configs aligned accordingly
  (related to ClearWater-data #7, `ffdbfe6`).
- `sumwere_creek_coarse` dataset and associated config YAMLs added to `examples/data/`
  for self-contained demos.

#### Riverine Integration
- `Update Riverine class for format dates` — date formatting aligned between modules
  and riverine (`5ac4400`).
- Riverine set to skip for TS to align process timesteps (`4b1a5f9`).

#### Examples & Documentation
- New **V2 example notebooks**: `01_Getting_Started.ipynb`, `02_Configuration_Files.ipynb`,
  `03_Example_TSM.ipynb` in `examples/V2/`.
- Linked Riverine/TSM v2 example notebook (coarse mesh and fine mesh variants) added to
  `examples/dev_sandbox/`.
- Older notebooks moved to `examples/archive/`.
- `README.md` updated.
- Documentation stubs added to v2 source code (`c6d5098`).

#### Profiling Infrastructure
- New `profiling/modules_runtime/` directory with `coupled.py` and `shared.py` scripts
  for measuring runtime performance of the coupled model.

---

### 🧪 Testing & Validation

| Date | Description |
|------|-------------|
| 2026-01-16 | Initial pytest setup for v2 (`ca19a33`) |
| 2026-01-28 | v2 reproduces v1 default pytest case temperature value (`634246d`) |
| 2026-01-29 | pytest v2 running for first 2 TSM v1 cases (`f8ae787`) |
| 2026-01-29 | First 14 TSM v1 pytest cases passing in v2 (`26ab3ba`) |
| 2026-02-03 | `use_sediment_temperature` pytest case implemented (`c49a108`) |
| 2026-02-03 | Runtime warnings reviewed for two pytest cases (issue #97, `326d45c`) |
| 2026-02-04 | `pytest-cov` working; fixed bug in `mixing_ratio_air` (`dbe0ec7`) |
| — | New `tests/unit/temperature/test_richardson.py` for Richardson number unit tests |

---

### 🐛 Bug Fixes

| Date | Description |
|------|-------------|
| 2025-08-04 | Fixed missed argument update to variable registry (`27bdafb`) |
| 2025-11-12 | Fixed indentation in example config file (`e6f99c9`) |
| 2026-01-16 | Fixed `xarray.where` call to work with non-Dataset array types in temperature (`60d08b7`) |
| 2026-01-22 | Fixed flipped atmospheric / upwelling longwave flux equations in temperature (`5d8321a`) |
| 2026-01-29 | Fixed Richardson number to use `xarray.where` arrays instead of scalar `if` logic (`d12bb09`) |
| 2026-01-29 | Added missing lines in temperature process (`c97388f`) |
| 2026-02-04 | Fixed bug in `mixing_ratio_air` calculation (`dbe0ec7`) |
| 2026-02-20 | Toggled off leftover print statements in model (`f7b0967`) |
| 2026-03-27 | Removed print statements that prevented v1 riverine from running (`c5fcacb`) |
| 2026-04-13 | Fixed relative-pathing bug when modules invoked from riverine context (`0d51a18`) |

---

### 🔧 Dependencies & Tooling

- **`pyyaml` added** as an explicit dependency for YAML config parsing (`53b558`).
- Python environment updated (`pyproject.toml`, `environment.yml`, `2df827c`).
- Pixi lock updated; macOS compatibility fixes applied (`7be3a83`).
- `.vscode/` added to `.gitignore` (`46e5c7c`).
- `conda develop` path updated in examples (`8f90081`).
- `pyproject.toml` version set to `0.4.0-alpha`; version is now dynamic, sourced from
  `src/clearwater_modules_v2/__init__.py`.

---

### 👥 Contributors

| Contributor | Role |
|-------------|------|
| **Paul Tomasula** (@ptomasula) | Lead v2 architect — framework, config, I/O, chunking, temperature processes |
| **Jason M Rutyna** (@jrutyna) | TSM v2 validation, pytest updates, coupled examples |
| **Anthony Aufdenkampe** (@aufdenkampe) | xarray/simlab exploration, pathing, release management |
| **Sarah Jordan** (@sjordan29) | Zarr examples, iterative coupled-run profiling script |

---

## [v0.3.1] — 2026-05-05

**Add profiling pre-refactor** — Baseline profiling of the v1 model added before the v2
architectural refactor began. (PR #102, @aufdenkampe)

---

## [v0.3.0] — 2023-09-25

**Nutrient Simulation Module (NSM)** — Introduced the refactored NSM (originally Fortran 95
by Zhang & Johnson); added phosphorus and pathogen modules, improved xarray pre-initialization,
and comprehensive algae/nitrogen test suites. See [GitHub Release](https://github.com/EcohydrologyTeam/ClearWater-modules/releases/tag/v0.3.0).

---

## [v0.2.0] — 2023-12-04

**OOP Refactor — TSM Fully Functional** — Major architectural modernization to
object-oriented design; TSM fully validated with pytest; ~20–35% performance improvements
via Numba. See [GitHub Release](https://github.com/EcohydrologyTeam/ClearWater-modules/releases/tag/v0.2.0).

---

## [v0.1.0] — 2022-10-12

**Initial Release** — Python ports of TSM, NSM-I, and GSM Fortran modules.
See [GitHub Release](https://github.com/EcohydrologyTeam/ClearWater-modules/releases/tag/v0.1.0).
