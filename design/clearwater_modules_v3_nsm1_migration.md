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
  (Critical Default-Value Corrections)
- `src/clearwater_modules_v3/parameter_defaults_corrections.md` Sections
  1, 2, 3, 4
- `design/clearwater_modules_v3_nsm1_audit_summary.md` (consolidated
  three-way audit)
- `design/clearwater_modules_v3_nsm1_README.md` Section 2.1 (audit
  history)

> **Important.** The audit-driven corrections in Phases 9.A, 9.B, and 9.C
> are behavior-affecting for users who relied on default-instantiated
> Process classes. See Section 3 of this document for the full list and
> the magnitudes of the change.

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
| `SOD_20`       | 999 g-O2/m^2/d | 1.0 g-O2/m^2/d | Drives DOX immediately negative           |
| `SOD_theta`    | 999        | 1.060            | `999^(T-20)` blows up at T > 20 C            |
| `kah_20_user`  | 999 1/d    | 0.0 1/d          | Runaway atmospheric flux on user-override   |
| `kaw_20_user`  | 999 m/d    | 0.0 m/d          | Same                                        |
| `pressure_mb`  | 2026.5     | 1013.25 hPa      | ~2x sea-level pressure; biases `O2sat`/`N2sat` |

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
| `clearwater_modules.shared.processes.celsius_to_kelvin`    | `clearwater_modules_v3.utils.conversions.celsius_to_kelvin` (re-exports v2's +273.16; see deviations)          |

---

## 3. Audit-driven corrections (Phases 9.A, 9.B, 9.C)

This section documents the corrections applied after the line-by-line
three-way audit (Fortran NSM1 vs v1 Python NSM1 vs v3 Python NSM1).
**These corrections affect simulation outputs.** Users who relied on
default-instantiated Process classes pre-Phase-9 will see substantively
different results post-fix; the changes are bug fixes that bring v3 into
agreement with v1/Fortran NSM1.

The full audit record is at `design/clearwater_modules_v3_nsm1_audit_summary.md`;
the canonical corrections record is at
`src/clearwater_modules_v3/parameter_defaults_corrections.md`.

### 3.1 Algae and Nitrogen wiring sweep (Phase 9.A)

**Symptom.** Pre-Phase-9.A, `FloatingAlgae`, `BenthicAlgae`, and `Nitrogen`
maintained two parallel parameter surfaces: a v3-style `DEFAULTS` dict
merged into `self` at construction time, plus legacy v2-style kwargs that
shadowed them with defaults of `0.0` (rates) or `1.0` (thetas, fractions).
The kinetic methods read the legacy kwargs, not the DEFAULTS. A user who
constructed a Process with no kwargs (`FloatingAlgae()`) silently received:

- `mu_max = 1.0/d` (kwarg default) instead of `1.0/d` from DEFAULTS — same
  by coincidence, but `mu_max_theta = 1.0` instead of 1.047 (no temperature
  correction).
- `respiration_rate = 0.0` instead of `krp_20 = 0.2/d` from DEFAULTS.
- `death_rate = 0.0` instead of `kdp_20 = 0.15/d` from DEFAULTS.
- `settling_velocity = 0.0` instead of `vsap = 0.15 m/d` from DEFAULTS.
- For BenthicAlgae: `growth_rate_max = 1.0` instead of `mub_max_20 = 0.4/d`
  (2.5x); `KsN = 0.04` (FloatingAlgae's value, ~6x too small for benthic);
  `KsP = 0.0012` (~100x too small); `Ksb = 1.0` (10x too small).
- For Nitrogen: `nitrification_rate = 1.0` (kwarg default) instead of
  `knit_20 = 0.1/d` from NITROGEN_DEFAULTS (10x off); `kdnit = 1.0` instead
  of `0.002/d` (500x); legacy `rnh4`, `vno3` defaults of 1.0 instead of
  0.0 (sediment-flux fluxes silently enabled at unit rate).

The existing parity tests did not catch this because they all called
`FloatingAlgae(mu_max=...)`, `Nitrogen(nitrification_rate=...)` etc. with
explicit kwargs that matched the legacy kwarg names but used v1/Fortran
values. Default instantiation was never tested.

**Fix.** Phase 9.A rewired all 13 algae kinetic methods and 5 Nitrogen
kinetic methods to read the `DEFAULTS` keys directly (`self.mu_max_20`,
`self.knit_20`, etc.). Legacy kwargs now default to `None` and are deprecated;
passing them still works but emits a deprecation warning. Three algae
formula bugs were also fixed:

- `FloatingAlgae.limit_light` option-1 parenthesization (the `(KL+PAR)`
  numerator and `(KL+PAR*exp(-Ld))` denominator were not both inside
  `np.log(...)`).
- `FloatingAlgae` harmonic-mean growth zero-guard: `where(FP == 1, 0, ...)`
  was inverted to `where(FP == 0, 0, ...)`.
- `BenthicAlgae` Steele light limitation: `x*exp(x-1)` was sign-flipped to
  the correct `x*exp(1-x)`.

**Behavioral consequence.** Default-instantiated Algae and Nitrogen
Processes pre-Phase-9.A produced effectively non-physical kinetics
(zero respiration, zero death, zero settling, runaway nitrification at
10x the canonical rate, runaway denitrification at 500x). Post-Phase-9.A
defaults match v1/Fortran NSM1 exactly. **If your pre-9.A simulation
relied on default instantiation, your results were almost certainly
non-physical; re-run with the post-9.A defaults.**

### 3.2 Phantom NH4 source removed (Phase 9.A)

**Symptom.** Pre-Phase-9.A, `change_ammonium` included a positive
first-order source term `+ ammonium_decay_rate * NH4` with kwarg default
`1.0/d`. Neither v1 nor Fortran NSM1 has any analog. NH4 grew exponentially
without bound at default kwargs (which is why default-instantiated Nitrogen
runs blew up immediately).

**Fix.** Phase 9.A dropped the term from the rate sum. v3's NH4 budget now
matches v1/Fortran exactly: nitrification sink, NH4 algal-uptake sink,
OrgN-hydrolysis source, ammonium-respiration source, sediment NH4
release source.

### 3.3 NO3 algal-uptake split made mass-conservative (Phase 9.A)

**Symptom.** Pre-Phase-9.A, the NH4 algal-uptake sink and NO3 algal-uptake
sink did not sum to `rna * AlgalGrowth`. The NH4 path used a dynamic
fraction `algal_nh4_uptake_fraction` (correct); the NO3 path used a static
parameter `float_algea_faction_uptake_from_nitrate = 1.0` (wrong). At
default `PN`, the two paths sum to ~1.5x of the actual algal-N demand,
violating algal-N mass balance.

**Fix.** Phase 9.A rewired `nitrate_uptake_floating_algae` to use the
dynamic `1 - algal_nh4_uptake_fraction`. The static parameter was retired.
NH4 + NO3 algal-uptake sinks now sum to `rna * AlgalGrowth` to within
floating-point tolerance.

### 3.4 Benthic NO3 uptake reconstructed (Phase 9.A)

**Symptom.** Pre-Phase-9.A, `nitrate_uptake_benthic_algae` divided by
`algal_chlorophyll = AWa = 1000` (the floating-algae denominator) instead
of `BWd = 100` (the canonical benthic-algae denominator), missed the
`/depth` divisor, used a static `fraction_bottom_area = 1.0` instead of
the dynamic `Fb = 0.9`, and used a static
`benthic_algea_faction_uptake_from_nitrate = 0.5` instead of the dynamic
`1 - balgae_nh4_uptake_fraction`. Stoichiometry, units, and dynamics all
disagreed with v1 / Fortran.

**Fix.** Phase 9.A reconstructed the formula as
`(1 - balgae_nh4_uptake_fraction) * (BWn / BWd) * Fb * balgae_growth_rate / depth`.
A `balgae_no3_uptake_fraction` cache was added to BenthicAlgae for clean
inter-Process routing.

### 3.5 Carbon and DOX `rca` / `rcb` stoichiometric ratios (Phase 9.B)

**Symptom.** Pre-Phase-9.B, the v3 `Carbon` and `DOX` Processes used the
raw stoichiometric weights `self.AWc = 40 mg-C/ug-Chla` and
`self.BWc = 40 mg-C/mg-D` directly in the algal coupling terms, where
v1 and Fortran use the *ratios* `rca = AWc/AWa = 40/1000 = 0.04` and
`rcb = BWc/BWd = 40/100 = 0.4`. The raw weights are stoichiometric mass
fractions; the ratios convert algal biomass (in chlorophyll-a or dry
weight units) to mg-C / L.

Affected sites: `dic_algal_resp`, `dic_algal_photo`, `dic_balgae_resp`,
`dic_balgae_photo` in `carbon.py`; `_floating_algae_growth_flux`,
`_floating_algae_respiration_flux`, `_benthic_algae_growth_flux`,
`_benthic_algae_respiration_flux` in `dox.py`. Pre-fix: floating-algae
DIC and O2 algal coupling was 1000x too large; benthic-algae 100x too
large.

The existing parity tests did not catch this because they explicitly
called the v1 helper functions with `rca = AWc = 40` (passing the same
wrong value into both sides of the parity comparison).

**Fix.** Phase 9.B replaced `self.AWc` with `self.AWc / self.AWa` (or
pre-computed `rca` once in `run`) at all four floating-algae call sites
in `carbon.py` and `dox.py`. Same for `self.BWc -> self.BWc / self.BWd`
at all four benthic call sites. The v1 Python helper functions define
the `rca`/`rcb` derivation; Fortran does the same. v3 1.0.0 uses the
derived values.

**Behavioral consequence.** Carbon and DOX algal coupling is now
quantitatively correct. Users running v3 NSM1 with active algae and
either Carbon or DOX enabled will see DIC and DOX trajectories shift by
a factor of 100 (benthic) or 1000 (floating) on the algal-coupling terms.
Other carbon and DOX terms (CBOD oxidation, reaeration, respiration)
are unaffected.

### 3.6 DOX SOD attenuation under hypoxia (Phase 9.B)

**Symptom.** v3's `utils/sediment.SOD_tc` is a pure-Arrhenius helper
(architectural improvement; see Section 4 below). The DOX-Monod
attenuation factor `DOX/(DOX+KsSod)` that Fortran applies inside its
`SOD_tc` was meant to be re-applied at the v3 DOX consumer site. Pre-9.B,
the consumer did not apply it. Under hypoxic conditions, v3 sediment kept
consuming oxygen at the full Arrhenius rate, with the negative-state
clip masking the conservation violation.

**Fix.** Phase 9.B multiplies `_sod_flux` by `DOX / (DOX + KsSod)` at the
consumer site (gated on `use_DOX`). The hypoxic-attenuation behavior now
matches Fortran while preserving the cleaner v3 architectural split.

### 3.7 Carbon DIC budget includes CBOD oxidation source (Phase 9.B)

**Symptom.** Pre-Phase-9.B, `dDIC/dt` omitted the CBOD oxidation source.
Fortran (`modCarbon.f90:262-266`) and v1 (`processes.py:2854`) both add
`DIC_CBOD_oxidation = sum(CBOD_oxidation) / roc / 12000` to dDIC/dt.

**Fix.** Phase 9.B adds `+ cbod_process.cbod_oxidation_rate / self.roc / 12000.0`
to dDIC/dt (gated on `use_cbod`). The CBOD process already exposes
`cbod_oxidation_rate` in its rate cache.

### 3.8 POC hydrolysis no longer DOX-attenuated (Phase 9.B)

**Symptom.** Pre-Phase-9.B, v3 multiplied POC hydrolysis by `DOX /
(KsOxmc + DOX)`. Neither Fortran nor v1 attenuates POC hydrolysis. POC ->
DOC is physical/chemical (cell-wall fragmentation, leaching), not
biochemically O2-limited.

**Fix.** Phase 9.B drops the `dox_attenuation` factor from POC hydrolysis.
v3 now matches v1/Fortran: `kpoc_tc * POC` only.

### 3.9 fdp partitioning unit factor (Phase 9.B)

**Symptom.** v3's `utils/partitioning.fdp` had `fdp = 1 / (1 + kdpo4 *
Solid / 0.000001)`. Fortran has `1 / (1 + kdpo4 * Solid / 1.0E6)`. The
v3 form divides by 1e-6 (multiplies by 1e6); Fortran divides by 1e6
(multiplies by 1e-6). Factor of 1e12 wrong.

This was inherited from v1's `shared/processes.py`. Both forms collapse
to `fdp = 1` at the v3 default `kdpo4 = 0` (latent), but break the
moment a user sets `kdpo4 > 0`.

**Fix.** Phase 9.B changed the literal to `1.0E6`. An MMS-style test now
covers the `kdpo4 > 0` regime against Fortran.

### 3.10 Parameter corrections (Phase 9.C)

Two inherited parameter values were corrected in Phase 9.C:

- **`vson_20` 0.1 → 0.01 (organic-N settling velocity).** The v3 nitrogen
  group's `vson_20` had drifted to 0.1 m/d while v3's own `global_vars.vson`,
  v1 `GlobalVars.vson`, and Fortran `vson` were all 0.01 m/d. Internal
  v3 inconsistency; corrected to 0.01 m/d.
- **`lambdam` 0.0174 → 0.174 (POM contribution to Beer-Lambert light
  extinction).** Likely v1 typo (10x too small); Fortran and QUAL2K
  Table 6 use 0.174. Multiple legacy v1 NSM tests already overrode the
  v1 default to 0.174, indirectly confirming the v1 default was wrong.

See `src/clearwater_modules_v3/parameter_defaults_corrections.md`
Sections 1.8 and 1.9 for the full rationale.

### 3.11 Open items flagged for LimnoTech reconciliation

Four parameter items remain open after Phases 9.A/B/C; see
`src/clearwater_modules_v3/parameter_defaults_corrections.md` Section 4
and the LimnoTech review packet
`design/clearwater_modules_v3_nsm1_limnotech_review.md`. Of particular
note for migration:

- **Nitrogen 4-way theta swap.** `kon_theta`, `kdnit_theta`,
  `rnh4_theta`, `vno3_theta` v3/v1 values appear to swap with Fortran's
  values pairwise (4-way pattern: v3/v1 have `1.074, 1.08, 1.047, 1.045`
  while Fortran has `1.047, 1.045, 1.074, 1.08`). Pending LimnoTech
  reconciliation. v3 inherits v1; users running v3 vs Fortran NSM1
  side-by-side at non-20-C temperatures will see ~10-15% differences in
  Arrhenius-corrected nitrification, denitrification, sediment NH4
  release, and sediment NO3 uptake rates from this swap alone.

---

## 4. v1↔v3 deliberate runtime numerical deviations

For users running A/B comparisons between v1 and v3 NSM1, the following
deliberate numerical differences exist independent of the audit-driven
corrections in Section 3. Each is documented in detail in
`src/clearwater_modules_v3/parameter_defaults_corrections.md` Section 3:

- **DOX SOD architectural split** — v3's `utils/sediment.SOD_tc` is pure
  Arrhenius; v1 has Monod inline (architectural choice). Phase 9.B
  re-applied the Monod factor at the consumer site, so the hypoxic-
  attenuation behavior now matches v1/Fortran.
- **Alkalinity DOX-attenuation flow** — v3 reads pre-attenuated
  nitrification / denitrification flux from the Nitrogen rate cache; v1
  applies the attenuation locally. Same numerical answer.
- **Pathogen light decay** — v3 uses `PAR = q_solar * Fr_PAR`; v1 uses
  raw `q_solar`. Calibration target `apx` absorbs the difference.
- **CBOD sedimentation** — v3 applies `ksbod_tc / depth` (m/d → 1/d); v1
  treats `ksbod_tc` directly as 1/d.
- **Kelvin conversion** — v3 utility re-exports v2's +273.16 form for
  v2-parity test stability; v1 NSM1 `nsm1/processes.py` and Fortran NSM1
  use +273.15.
- **mb→atm conversion** — v3 `1.0 / 1013.25` vs v1 literal `0.000986923`
  (~7 sig fig agreement).

Per-test docstrings in `tests/test_5_*_calculations_v2.py` document the
empirical verification for each deviation. Note that **Carbon POC
hydrolysis** is no longer in this list: Phase 9.B retired the v3-only
DOX-Monod attenuation that had been documented as a deliberate v3
improvement (it was found to disagree with both v1 and Fortran on the
underlying physics; POC hydrolysis is not aerobically limited).

---

## 5. Retirement timeline

Consistent with the overall v3 plan (architecture spec Section 5):

- **v3 NSM1 1.0.0 ships by 2026-05-31** (target).
- At that point, v1 NSM1 (`clearwater_modules.nsm1`) and v2 NSM1
  (`clearwater_modules_v2.processes.nitrogen` etc.) are frozen.
- **v3 1.1.0 (post-2026-05-31) removes** v1 and v2 NSM1 source from the
  tree alongside the broader v2 retirement.

Until v3 1.1.0, both v1 and v2 NSM1 remain importable for downstream
projects mid-migration.

---

## 6. References

- `design/clearwater_modules_v3_nsm1_README.md` — quick-start and
  Process inventory.
- `design/clearwater_modules_v3_nsm1_design_specification.md` Section 6
  (bug list), Section 7 (defaults), Section 8 (migration), Section 11
  (phased plan), Section 14 (resolved decisions).
- `src/clearwater_modules_v3/parameter_defaults_corrections.md` — every
  default-value correction and runtime-numerical deviation; Section 4
  is the LimnoTech reconciliation list.
- `design/clearwater_modules_v3_nsm1_audit_summary.md` — consolidated
  three-way audit (Fortran NSM1 vs v1 Python NSM1 vs v3 Python NSM1).
- `design/clearwater_modules_v3_nsm1_audit_*.md` — the five sub-audits
  (algae, n_p, c_dox, simple_constituents, utilities_params).
- `design/clearwater_modules_v3_nsm1_limnotech_review.md` — packaged
  review materials for external review.
- `examples/V3/04_Example_NSM1.ipynb` — coupled end-to-end demo.
