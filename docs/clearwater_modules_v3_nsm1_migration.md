# v3 NSM1 Migration Notes — v1→v3 and v2→v3

This document covers migration paths to `clearwater_modules_v3` NSM1 1.0.0
from the two predecessor codebases:

- **v1 NSM1** at `src/clearwater_modules/nsm1/` — full 16-constituent
  implementation with the `NutrientBudget` xarray-DAG model.
- **v2 NSM1** at `src/clearwater_modules_v2/processes/{nitrogen,
  floating_algae, benthic_algae}.py` — partial 4-of-16 implementation
  using the v2 `Process`/`Model` framework.

Reference materials:

- `design/clearwater_modules_v3_nsm1_design_specification.md` Section 8
  (Migration Strategy)
- `design/clearwater_modules_v3_nsm1_design_specification.md` Section 6
  (Bug Fix Inventory — 16 items)
- `design/clearwater_modules_v3_nsm1_design_specification.md` Section 7
  (Critical Default-Value Corrections — 7 items)
- `src/clearwater_modules_v3/parameter_defaults_corrections.md` Sections
  1, 2, 3

---

## 1. v2 NSM1 → v3 NSM1

### What changed for v2 users

Existing users of `clearwater_modules_v2.processes.nitrogen.Nitrogen`,
`floating_algae.FloatingAlgae`, `benthic_algae.BenthicAlgae` see the
following behavioral and surface changes when upgrading to v3:

#### Bug fixes (16 items, design spec Section 6)

1. **NH4 multiplicative integrator** (`nitrogen.py:101`):
   `ammonium = 0 + ammonium * rate * dt` → `ammonium + rate * dt`.
2. **NO3 multiplicative integrator** (`nitrogen.py:115`): same fix.
3. **`time_step_frequency` typo** (`nitrogen.py:115`).
4. **Ap multiplicative integrator + stray ×86400**
   (`floating_algae.py:122`): `algae * rate * dt * 86400` →
   `algae + rate * dt`.
5–8. **NaN guards using `rate == np.nan`** (always False per IEEE 754):
   replaced with `.isnull()` at four locations.
9. **Hard-coded `half_saturation_oxygen=1`** (`nitrogen.py:197`): wired to
   parameter library.
10–11. **Hard-coded `algea_growth_rate=0`** for nitrate uptake
   (`nitrogen.py:211, 217`): now reads from `FloatingAlgae` /
   `BenthicAlgae` rate cache. Algal nitrate uptake actually happens.
12. **Death rate placeholder** (`nitrogen.py:60`): wired to floating-algae
   rate cache.
13. **`ammonium_respiration()` returns 0** (`floating_algae.py:401`):
   implemented per v1 line 1273.
14. **`ammonium_growth()` returns 0** (`floating_algae.py:407`):
   implemented per v1 line 1206.
15. **Hard-coded `phosphate_fraction_dissolved=0.5`**
   (`floating_algae.py:113`): replaced with
   `clearwater_modules_v3.utils.partitioning.fdp` shared utility.
16. **`set_at_time` persistence** (state computed but dropped): fixed in
   `Nitrogen.run` and `FloatingAlgae.run` so state actually evolves
   between substeps.

#### New state variable: OrgN

The v3 `Nitrogen` Process owns three state variables: `ammonium`,
`nitrate`, **`organic_nitrogen`** (new in v3). v3 adds OrgN hydrolysis to
NH4 (`kon_tc * OrgN`) and OrgN settling to bed (`vson_tc / depth * OrgN`)
following the v1 kinetics. Sources: algal mortality routed through
`algal_orgn_from_mortality_rate` and
`balgae_orgn_from_mortality_rate`. Optional via `use_OrgN` flag (default
`True`).

#### Behavioral consequences

- **Existing v2 simulations will produce different state values after
  migration.** All four v2 NSM1 processes had at least one
  state-affecting bug; correcting them changes the trajectory.
- **Algal nitrate uptake now actually occurs** (v2 silently set the rate
  to zero). Expect lower steady-state nitrate concentrations whenever
  algae are growing.
- **State actually evolves.** v2's persistence regression dropped
  computed state values on each substep, producing near-constant outputs
  on production paths. v3 writes back via `set_at_time` and state
  trajectories reflect the integrated rates.
- **Phosphorus partitioning is correct.** v2's
  `phosphate_fraction_dissolved=0.5` was a hard-coded placeholder. v3
  computes `fdp` from the configured `kdpo4`, suspended-solids
  concentration, and water density.

#### YAML schema

The v3 YAML schema accepts every v2 NSM1 YAML config that does not
conflict with the new `parameters:` block convention. v3 NSM1 Process
classes share the v2 framework: `Process` subclass, `ProcessFactory.register`,
per-process `time_step` declared in the YAML, registry-driven coupling.

#### Worked example

```yaml
# v2 NSM1 YAML
processes:
  - floating_algae:
      time_step: '30s'
  - benthic_algae:
      time_step: '30s'
  - nitrogen:
      time_step: '30s'
```

```yaml
# v3 NSM1 YAML — additive: parameters now optional, defaults from DEFAULTS
processes:
  - floating_algae:
      time_step: '30s'
      parameters:
        knit_20: 0.10          # optional override; defaults pulled from
                               # parameters/algae.py if omitted
  - benthic_algae:
      time_step: '30s'
  - nitrogen:
      time_step: '30s'
      parameters:
        use_OrgN: true         # optional; default true
  # 8 new Process declarations to exercise the full constituent set:
  - phosphorus: { time_step: '30s' }
  - carbon:     { time_step: '30s' }
  - pom:        { time_step: '30s' }
  - cbod:       { time_step: '30s' }
  - dox:        { time_step: '30s' }
  - pathogen:   { time_step: '30s' }
  - alkalinity: { time_step: '30s' }
  - n2:         { time_step: '30s' }
```

The v2 → v3 migration is import-path-only for the existing
`Nitrogen`/`FloatingAlgae`/`BenthicAlgae` Process classes and additive
for the 8 new constituents.

#### Code-level migration

```python
# v2
from clearwater_modules_v2.processes.nitrogen import Nitrogen
from clearwater_modules_v2.processes.floating_algae import FloatingAlgae
from clearwater_modules_v2.processes.benthic_algae import BenthicAlgae

# v3 — same Process classes, with bugs fixed; same constructor surface
from clearwater_modules_v3.processes import (
    Nitrogen, FloatingAlgae, BenthicAlgae,
    Phosphorus, Carbon, POM, CBOD, DOX, Pathogen, Alkalinity, N2,
)
```

### v2 → v3 API mapping table

| v2 import                                                      | v3 equivalent                                                  |
| ---                                                            | ---                                                            |
| `clearwater_modules_v2.processes.nitrogen.Nitrogen`            | `clearwater_modules_v3.processes.nitrogen.Nitrogen`            |
| `clearwater_modules_v2.processes.floating_algae.FloatingAlgae` | `clearwater_modules_v3.processes.floating_algae.FloatingAlgae` |
| `clearwater_modules_v2.processes.benthic_algae.BenthicAlgae`   | `clearwater_modules_v3.processes.benthic_algae.BenthicAlgae`   |
| (not in v2)                                                    | `clearwater_modules_v3.processes.phosphorus.Phosphorus`        |
| (not in v2)                                                    | `clearwater_modules_v3.processes.carbon.Carbon`                |
| (not in v2)                                                    | `clearwater_modules_v3.processes.pom.POM`                      |
| (not in v2)                                                    | `clearwater_modules_v3.processes.cbod.CBOD`                    |
| (not in v2)                                                    | `clearwater_modules_v3.processes.dox.DOX`                      |
| (not in v2)                                                    | `clearwater_modules_v3.processes.pathogen.Pathogen`            |
| (not in v2)                                                    | `clearwater_modules_v3.processes.alkalinity.Alkalinity`        |
| (not in v2)                                                    | `clearwater_modules_v3.processes.n2.N2`                        |
| `clearwater_modules_v2.config.init_from_file`                  | `clearwater_modules_v3.config.init_from_file`                  |
| `clearwater_modules_v2.Model`                                  | `clearwater_modules_v3.Model`                                  |

---

## 2. v1 NSM1 → v3 NSM1

### What changed for v1 users

Existing users of `clearwater_modules.nsm1.NutrientBudget` (the v1
xarray-DAG model) see the following changes when upgrading to v3:

#### New API surface — Process classes + YAML config

v1's xarray-DAG dependency model is replaced with explicit producer/
consumer rate variables routed through a registry. Each Process owns its
state variables, computes net rates, integrates by Forward Euler, and
writes results back to the registry.

The dispatch order follows v1's `ComputeKinetics` sequence:
`FloatingAlgae → BenthicAlgae → Nitrogen → Phosphorus → Carbon → POM →
CBOD → DOX → Pathogen → Alkalinity → N2`. The order is declared by the
YAML `processes:` block and honored by `clearwater_modules_v3.Model`.

#### Sentinel-`999` defaults corrected (7 items, design spec Section 7)

v1's `nsm1/constants.py` carried sentinel-`999` defaults that produced
catastrophic blowup under any non-trivial run. v3 corrects these at the
port:

| Parameter      | v1 default | v3 default       | Risk                                        |
| ---            | ---        | ---              | ---                                         |
| `vsop`         | 999 m/d    | 0.1 m/d          | Drains all OrgP per timestep                |
| `vs`           | 999 m/d    | 0.1 m/d          | Drains all TIP per timestep                 |
| `SOD_20`       | 999 g-O2/m²/d | 1.0 g-O2/m²/d | Drives DOX immediately negative             |
| `SOD_theta`    | 999        | 1.060            | `999^(T-20)` blows up at T > 20 °C          |
| `kah_20_user`  | 999 1/d    | 0.0 1/d          | Runaway atmospheric flux on user-override   |
| `kaw_20_user`  | 999 m/d    | 0.0 m/d          | Same                                        |
| `pressure_mb`  | 2026.5     | 1013.25 hPa      | ~2× sea-level pressure; biases `O2sat`/`N2sat` |

See `src/clearwater_modules_v3/parameter_defaults_corrections.md` Section
1 for the full rationale per item.

#### Hotstart and wet-mask

v1's hotstart support (added 2026-05-01 via `hotstart_dataset` and
`hotstart_timestep` kwargs on `NutrientBudget`) lifts unchanged into v3.
v3 declares hotstart in the YAML rather than in code:

```yaml
hotstart:
  dataset_path: previous_run.nc
  timestep: '2022-05-13 12:00:00'

wet_mask:
  variable: wetted_surface_area
  threshold: 1.0
```

See the umbrella v3 README (`src/clearwater_modules_v3/README.md`) for
the hotstart and wet-mask schema specification.

#### Worked example

```python
# v1 — direct invocation of NutrientBudget
import xarray as xr
from clearwater_modules.nsm1.model import NutrientBudget
from clearwater_modules.nsm1 import constants

algae_pars = constants.DEFAULT_ALGAE
nitrogen_pars = constants.DEFAULT_NITROGEN
# ... define all 13 TypedDict groups in code ...

initial = xr.Dataset({
    "Ap":  ("cell", [10.0] * 5),
    "NH4": ("cell", [0.1] * 5),
    # ... 14 more state variables ...
})

nsm = NutrientBudget(
    initial_state_values=initial,
    algae_parameters=algae_pars,
    nitrogen_parameters=nitrogen_pars,
    # ... 11 more parameter groups ...
    time_steps=10,
    time_dim="time",
)

for _ in range(10):
    nsm.increment_timestep()

ammonium_trajectory = nsm.dataset["NH4"]
```

```python
# v3 — YAML-driven Model
from clearwater_modules_v3.config import init_from_file

model = init_from_file("nsm1_run.yml")
model.run()

ammonium_trajectory = model.registry.get_history("ammonium")
```

with `nsm1_run.yml`:

```yaml
simulation_directory: ./outputs

time_step: '30s'
duration: '1h'

# YAML-only block; in v1 this was 13 dicts in Python code
processes:
  - floating_algae:
      time_step: '30s'
      parameters:
        AWa: 100.0
        knit_20: 0.10
  - benthic_algae:  { time_step: '30s' }
  - nitrogen:       { time_step: '30s' }
  - phosphorus:     { time_step: '30s' }
  - carbon:         { time_step: '30s' }
  - pom:            { time_step: '30s' }
  - cbod:           { time_step: '30s' }
  - dox:            { time_step: '30s' }
  - pathogen:       { time_step: '30s' }
  - alkalinity:     { time_step: '30s' }
  - n2:             { time_step: '30s' }

initial_state:
  # 16 state variables seeded from a NetCDF or inline
  ...
```

#### v1 expected-value tables

v1's hard-coded expected values in `tests/NSM Manual Calcs/` were
generated under v1's mixed Jacobi/Gauss-Seidel state-read pattern. v3
enforces strict Jacobi-state semantics. For tests where the state-read
pattern affects the result, v3 regenerates expected values rather than
treating the v1 numbers as ground truth. See gap analysis Section 1
finding 3.

### v1 → v3 API mapping table

| v1 reference                                               | v3 equivalent                                                                                                  |
| ---                                                        | ---                                                                                                            |
| `clearwater_modules.nsm1.model.NutrientBudget`             | `clearwater_modules_v3.config.init_from_file(yaml_path)` (returns a `Model`)                                   |
| Code-defined parameters: `algae_parameters={...}, ...`     | YAML `parameters:` blocks per Process                                                                          |
| `nsm.increment_timestep()`                                 | `model.run()`                                                                                                  |
| Hotstart kwargs: `hotstart_dataset=ds, hotstart_timestep=0` | YAML `hotstart:` block with `dataset_path` and `timestep`                                                      |
| State access: `nsm.dataset['NH4']`                         | `model.registry.get('ammonium')` (or `get_history` for trajectories)                                           |
| `clearwater_modules.shared.processes.kah_20`               | `clearwater_modules_v3.utils.reaeration.kah_20`                                                                |
| `clearwater_modules.shared.processes.SOD_tc`               | `clearwater_modules_v3.utils.sediment.SOD_tc`                                                                  |
| `clearwater_modules.shared.processes.fdp`                  | `clearwater_modules_v3.utils.partitioning.fdp`                                                                 |
| `clearwater_modules.shared.processes.PAR`, `L`             | `clearwater_modules_v3.utils.light.PAR`, `L`                                                                   |
| `clearwater_modules.shared.processes.celsius_to_kelvin`    | `clearwater_modules_v3.utils.conversions.celsius_to_kelvin` (273.15 instead of v1's 273.16; see deviations) |

---

## 3. v1↔v3 numerical deviations

For users running A/B comparisons between v1 and v3, the following
numerical deviations are documented in detail in
`src/clearwater_modules_v3/parameter_defaults_corrections.md` Section 3:

- **Carbon POC hydrolysis** — v3 adds DOX-Monod attenuation; v1 has none.
- **DOX SOD** — v3 uses pure-Arrhenius `SOD_tc`; v1 has Monod inline.
- **Alkalinity DOX-attenuation flow** — v3 reads pre-attenuated
  nitrification / denitrification flux from Nitrogen rate cache.
- **Pathogen light decay** — v3 uses `PAR = q_solar * Fr_PAR`; v1 uses
  raw `q_solar`.
- **CBOD sedimentation** — v3 applies `ksbod_tc / depth` (m/d → 1/d); v1
  treats `ksbod_tc` directly as 1/d.
- **Kelvin conversion** — v3 273.15 vs v1 273.16 (0.01 K offset).
- **mb→atm conversion** — v3 `1.0 / 1013.25` vs v1 literal `0.000986923`
  (~7 sig fig agreement).

Per-test docstrings in `tests/test_5_*_calculations_v2.py` document the
empirical verification for each deviation.

---

## 4. Retirement timeline

Consistent with the overall v3 plan (architecture spec Section 5):

- **v3 NSM1 1.0.0 ships by 2026-05-31** (target).
- At that point, v1 NSM1 (`clearwater_modules.nsm1`) and v2 NSM1
  (`clearwater_modules_v2.processes.nitrogen` etc.) are frozen.
- **v3 1.1.0 (post-2026-05-31) removes** v1 and v2 NSM1 source from the
  tree alongside the broader v2 retirement.

Until v3 1.1.0, both v1 and v2 NSM1 remain importable for downstream
projects mid-migration.

---

## 5. References

- `docs/clearwater_modules_v3_nsm1_README.md` — quick-start and
  Process inventory.
- `design/clearwater_modules_v3_nsm1_design_specification.md` Section 6
  (bug list), Section 7 (defaults), Section 8 (migration), Section 11
  (phased plan), Section 14 (resolved decisions).
- `src/clearwater_modules_v3/parameter_defaults_corrections.md` — every
  default-value correction and runtime-numerical deviation.
- `docs/clearwater_modules_v3_nsm1_limnotech_review.md` — packaged
  review materials for external first-pass review.
- `examples/V3/04_Example_NSM1.ipynb` — coupled end-to-end demo.
