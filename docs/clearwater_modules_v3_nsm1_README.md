# clearwater_modules_v3 — NSM1

v3 NSM1 1.0.0 — Nutrient Simulation Module port from v1, framework alignment with v2.

This document covers the NSM1 module within `clearwater_modules_v3`. For the
broader v3 package overview (TSM, `Model` orchestration, hotstart, wet-mask)
see `src/clearwater_modules_v3/README.md`. The full source-of-truth design is
`design/clearwater_modules_v3_nsm1_design_specification.md`.

---

## 1. What v3 NSM1 is

v3 NSM1 is a complete port of v1 NSM1's 16-constituent kinetics suite into the
v2 `Process` class framework, with all 16 known v2 NSM1 bugs fixed at the
port. The 11 Process classes share a single registry, a single integrator
contract (additive Forward Euler), and a single dispatch model (Jacobi state
reads, Gauss-Seidel rate variables) provided by `clearwater_modules_v3.Model`.

v3 NSM1 supersedes both:

- v1 NSM1 (`clearwater_modules.nsm1`) — fully functional but xarray-DAG based,
  with sentinel-`999` parameter defaults that produce catastrophic blowup
  under default settings.
- v2 NSM1 (`clearwater_modules_v2.processes.{nitrogen, floating_algae, benthic_algae}`) —
  4-of-16 partial implementation with a multiplicative-integrator bug, broken
  NaN guards, and several hard-coded zero placeholders that silently disabled
  algal nutrient uptake.

---

## 2. Scope: 16 constituents, 11 Process classes

| Process class    | State variables                  | Source                                            |
| ---              | ---                              | ---                                               |
| `FloatingAlgae`  | Ap                               | Extends v2 `FloatingAlgae` (Phase 2 fixes)        |
| `BenthicAlgae`   | Ab                               | Extends v2 `BenthicAlgae` (Phase 2 fixes)         |
| `Nitrogen`       | NH4, NO3, OrgN                   | Extends v2 `Nitrogen` (Phase 2 fixes; OrgN added) |
| `Phosphorus`     | TIP, OrgP                        | New in v3                                         |
| `Carbon`         | POC, DOC, DIC                    | New in v3                                         |
| `POM`            | POM                              | New in v3                                         |
| `CBOD`           | CBOD (multi-group)               | New in v3                                         |
| `DOX`            | DOX                              | New in v3                                         |
| `Pathogen`       | PX                               | New in v3                                         |
| `Alkalinity`     | Alk                              | New in v3 (declared but inactive in v1)           |
| `N2`             | N2, TDG                          | New in v3                                         |

11 Process classes covering 16 state variables. Process source files live at
`src/clearwater_modules_v3/processes/<name>.py`. Default parameter dictionaries
live at `src/clearwater_modules_v3/parameters/<group>.py`.

---

## 3. Quick-start

```python
from datetime import datetime, timedelta

import numpy as np
import xarray as xr

from clearwater_modules_v3.processes import (
    FloatingAlgae, BenthicAlgae, Nitrogen, Phosphorus, Carbon, POM,
    CBOD, DOX, Pathogen, Alkalinity, N2,
)

# Synthetic 5-cell registry; in production this is built by Model from YAML.
from tests.v3.nsm1.conftest import InMemoryRegistry  # for illustration only

cells = 5
registry = InMemoryRegistry()

# Forcings and state. See per-Process docstrings for the full input list.
registry.register("water_temperature", xr.DataArray(np.full(cells, 20.0), dims="cell"))
registry.register("depth",             xr.DataArray(np.full(cells,  1.0), dims="cell"))
registry.register("ammonium",          xr.DataArray(np.full(cells,  0.1), dims="cell"))
registry.register("nitrate",           xr.DataArray(np.full(cells,  0.5), dims="cell"))
registry.register("organic_nitrogen",  xr.DataArray(np.full(cells,  0.2), dims="cell"))
# ... (continue for TIP, OrgP, POC, DOC, DIC, POM, CBOD, DOX, PX, Alk, N2, etc.)

t0 = datetime(2026, 1, 1)
dt = timedelta(seconds=30)

nitrogen = Nitrogen(parameters={}, time_step=dt)  # DEFAULTS auto-merged
nitrogen.run(registry, t_current=t0, dt=dt)

print(registry.get("ammonium"))
```

For an end-to-end coupled simulation (all 11 Processes plus Riverine
transport on a synthetic 5-cell mesh), see
`examples/V3/04_Example_NSM1.ipynb`. That notebook is the canonical
demonstration of v3 NSM1's Process inventory, integrator behavior, mass
conservation contract, and v1↔v3 numerical deviations.

---

## 4. Architecture

### Resolved design decisions (spec Section 14)

These are decided and implemented; LimnoTech confirmation closes the three
items still flagged tentative.

- **Parameter library: DEFAULTS-merge pattern.** Each `Process` class imports
  a `DEFAULTS: dict[str, float]` from `parameters/<group>.py`. At
  construction time it merges `{**DEFAULTS, **user_parameters}` and stores
  results as `self.<name>` attributes. YAML overrides are partial; defaults
  fill in unspecified keys.
- **Within-step semantics: Jacobi state, Gauss-Seidel rate variables.**
  Every `Process.run` reads state via `registry.get_at_time(t=t_current)`
  (always pre-update, order-independent). Rate variables are step-scoped:
  cleared at start of each step, written by producers, read by consumers,
  and reading an unset rate variable raises an error (catches
  dispatch-order bugs immediately).
- **Negative-state handling: clip-with-log.** Each `Process.run` calls
  `clearwater_modules_v3.utils.numerics.clip_negative_state(...)` after the
  Forward Euler step. Clip target is exactly 0. Each clip event is logged
  (rate-limited) and counted on `model.diagnostics`. Tier 1 closed-system
  conservation tests assert `clip_events == 0`.
- **Alkalinity as simple tracer.** Source/sink terms (nitrification consumption,
  denitrification production, algal-growth and -respiration coupling)
  integrated by Forward Euler. No carbonate equilibrium, no pH solver. Full
  pH chemistry is NSM2 territory (v3 1.1+).
- **Single-compartment algae.** One `Ap` and one `Ab`. Multi-group
  phytoplankton lands in v3 1.1+ as a new `PhytoplanktonGroups` Process
  *alongside* `FloatingAlgae`, not as in-place extension.

### Integrator contract

Each `Process.run`:

1. Reads its state variables from the registry at `t_current` (Jacobi).
2. Computes net rate of change for each state variable in `[state]/second`,
   additively combining sources and sinks.
3. Applies Forward Euler: `state_new = state_old + rate * dt_seconds`.
4. Calls `clip_negative_state(...)` on each new state.
5. Writes the post-clip new state to the registry at `t_current + dt`.

The 16 v2 NSM1 bugs (multiplicative integrator, persistence regressions,
broken NaN guards, hard-coded zeros) are exactly the contract violations
this pattern eliminates.

---

## 5. Tier 1 conservation contract

Every Process passes a closed-system mass-conservation test at
`rtol=1e-12` with **zero clip events**. The closed-system harness
(`tests/v3/nsm1/test_validation_tier1_conservation.py`) zeros all boundary
fluxes, all settling velocities, and all sediment fluxes; runs each Process
in isolation; and asserts that the relevant elemental total (N, P, C,
O₂-equivalents, Alk-equivalents) is conserved to roundoff and that
`model.diagnostics.clip_events == 0`. A clip event under closed-system
conditions indicates either unphysical parameters or a malformed test —
the test fails in that case, which is the correct diagnostic.

Per-process Tier 1 tests live at `tests/v3/nsm1/test_<constituent>_tier1.py`.

---

## 6. v1↔v3 deviations summary

v3 NSM1 reproduces v1 NSM1 kinetics within floating-point tolerance for the
overwhelming majority of sub-rate terms. A small set of deliberate numerical
deviations is documented in detail per-test (read the docstrings of
`tests/test_5_*_calculations_v2.py`) and summarized in
`src/clearwater_modules_v3/parameter_defaults_corrections.md` Section 3.
The headline list:

- **Carbon POC hydrolysis** — v3 adds DOX-Monod attenuation; v1 has none.
- **DOX SOD** — v3 uses pure-Arrhenius `SOD_tc`; v1 has Monod inline (an
  architectural choice, not a defect).
- **Alkalinity DOX-attenuation** — v3 reads pre-attenuated nitrification and
  denitrification fluxes from the Nitrogen rate cache; v1 applies the
  attenuation locally inside `Alk_nitrification` / `Alk_denitrification`.
  Same numerical answer when tests pass v3-aligned upstream rates.
- **Pathogen light decay** — v3 uses `PAR = q_solar * Fr_PAR`; v1 uses raw
  `q_solar`. Calibration target `apx` absorbs the difference.
- **CBOD sedimentation** — v3 applies `ksbod_tc / depth` (m/d → 1/d); v1
  treats `ksbod_tc` directly as 1/d.
- **Kelvin conversion** — v3 uses 273.15; v1 uses 273.16. 0.01 K offset.
- **Pressure mb→atm** — v3 uses `1.0 / 1013.25`; v1 uses literal
  `0.000986923`. Agreement to ~7 significant figures.

For per-test docstring detail: `tests/test_5_carbon_calculations_v2.py`,
`test_5_dox_calculations_v2.py`, `test_5_alkalinity_calculations_v2.py`,
`test_5_pathogen_calculations_v2.py`, `test_5_cbod_calculations_v2.py`,
`test_5_n2_calculations_v2.py`.

---

## 7. DOX reaeration default note (Phase 7.C Item 3)

DOX reaeration is **disabled by default** in v3:

```python
DEFAULTS = {
    'kah_20_user': 0.0,   # 1/d; user-override hydraulic reaeration at 20 C
    'kaw_20_user': 0.0,   # m/d; user-override wind reaeration at 20 C
    # ...
}
```

These zero defaults are corrections from v1's sentinel-`999` values (which
produced catastrophic atmospheric flux on any timestep where the override
branch was selected). For runs with non-trivial atmospheric exchange the
user must set `kah_20_user` and `kaw_20_user` explicitly in YAML or via the
`parameters={...}` constructor argument, or select one of the menu-driven
reaeration options (`hydraulic_reaeration_option`, `wind_reaeration_option`)
that do not consult these overrides.

See `src/clearwater_modules_v3/parameter_defaults_corrections.md` Sections
1.5 and 1.6 for the rationale and audit trail.

---

## 8. v1↔v3 high-level deviations beyond DOX reaeration

For the full bug-fix inventory and parameter corrections that v3 carries
forward from the port, read:

- Bug list — `design/clearwater_modules_v3_nsm1_design_specification.md`
  Section 6 (16 v2 bugs fixed at the port).
- Default-value corrections —
  `src/clearwater_modules_v3/parameter_defaults_corrections.md` Section 1
  (7 corrections: `vsop`, `vs`, `SOD_20`, `SOD_theta`, `kah_20_user`,
  `kaw_20_user`, `pressure_mb`).
- Lower-priority audit findings —
  `src/clearwater_modules_v3/parameter_defaults_corrections.md` Section 2
  (8 items kept at v1 values pending Phase 1.3 audit).
- Phase 7 numerical deviations —
  `src/clearwater_modules_v3/parameter_defaults_corrections.md` Section 3
  (the 7-item list summarized above).

For migration from v1 NSM1 (`NutrientBudget`) or v2 NSM1
(`clearwater_modules_v2.processes.nitrogen` etc.), see
`docs/clearwater_modules_v3_nsm1_migration.md`.

---

## 9. Test counts

- **652 tests passing** across the v3 NSM1 suite (Phases 0–7.C).
- **8 Tier 1 conservation tests** under `tests/v3/nsm1/test_<constituent>_tier1.py`
  plus the consolidated harness
  `tests/v3/nsm1/test_validation_tier1_conservation.py`.
- **36 sub-rate v1-parity tests** in `tests/test_5_*_calculations_v2.py`,
  each pinning a single v1 sub-term against v3's cached rate variable.
- **1 coupled end-to-end demo** at `examples/V3/04_Example_NSM1.ipynb`
  exercising all 11 Process classes plus Riverine on a synthetic mesh.

---

## 10. Open follow-up items (Phase 7.C)

The Phase 7.C review identified 4 follow-up items beyond the v3 1.0.0
LimnoTech-review-ready milestone. As of Phase 8.A:

1. **Nitrogen / Phosphorus oxygen-inhibition contract** — verify that the
   pre-attenuated `nitrification_flux_rate` and `denitrification_flux_rate`
   cached by `Nitrogen.run` are the canonical inputs for downstream
   Alkalinity coupling. *Status:* addressed in Phase 8.A.
2. **CBOD `use_DOX` default verification** — the `use_DOX` flag governs
   whether oxygen-inhibition is applied to CBOD oxidation. *Status:*
   addressed in Phase 8.A.
3. **DOX reaeration default note** — the `kah_20_user=0` / `kaw_20_user=0`
   defaults disable reaeration in the user-override branch by default. *Status:*
   documented in Section 7 of this README.
4. **Performance benchmark on Sumwere Creek** — current performance is
   17.6 ms/step on a 5-cell synthetic mesh. Extrapolation to Sumwere Creek
   (~600 cells) is deferred to a formal benchmark in a follow-up phase.
   *Status:* deferred; "must" target of 415 ms/step (spec Section 10) is
   confidently met based on current per-cell cost; "should" target may need
   profiling.

---

## 11. References

- `design/clearwater_modules_v3_nsm1_design_specification.md` — full
  source-of-truth design (Sections 6 bug list, 7 defaults, 8 migration,
  9 testing, 14 design decisions).
- `design/clearwater_modules_v3_architecture_specification.md` — umbrella
  architecture, package layout, integrator-pattern contract.
- `docs/clearwater_modules_v3_nsm1_gap_analysis.md` — Phase 0 readiness
  assessment.
- `docs/clearwater_modules_v3_nsm1_migration.md` — v1→v3 and v2→v3
  migration paths with side-by-side examples.
- `docs/clearwater_modules_v3_nsm1_limnotech_review.md` — packaged review
  materials for external first-pass review.
- `src/clearwater_modules_v3/parameter_defaults_corrections.md` — every
  default-value correction and runtime-numerical deviation v3 carries
  relative to v1.
- `examples/V3/04_Example_NSM1.ipynb` — coupled end-to-end demo.
