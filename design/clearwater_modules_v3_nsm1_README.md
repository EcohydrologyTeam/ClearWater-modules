# clearwater_modules_v3 — NSM1

v3 NSM1 1.0.0 — Nutrient Simulation Module port from v1, framework alignment with v2.

This document covers the NSM1 module within `clearwater_modules_v3`. For the
broader v3 package overview (TSM, `Model` orchestration, hotstart, wet-mask)
see `src/clearwater_modules_v3/README.md`. The full source-of-truth design is
`design/clearwater_modules_v3_nsm1_design_specification.md`.

ClearWater is an ERDC product released under an open-source license; v3 NSM1
is part of the v3 line of the ClearWater modeling framework.

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

## 2. Status (post-audit)

v3 NSM1 1.0.0 is the post-audit deliverable. An earlier Phase 8 README claimed
"review-ready"; that framing was retracted after a line-by-line three-way
audit of every Process and shared utility. The current status is:

- Code in `src/clearwater_modules_v3/processes/` reflects the audit-driven
  corrections from Phases 9.A, 9.B, and 9.C.
- The full test suite reports **705 tests passing, 3 xfailed** as of Phase
  9.C.
- 4 items remain open and are tracked for LimnoTech reconciliation; see
  Section 11 below and `src/clearwater_modules_v3/parameter_defaults_corrections.md`
  Section 4.
- The LimnoTech review packet at
  `design/clearwater_modules_v3_nsm1_limnotech_review.md` lists the focus
  areas for external review.

### 2.1 Audit history

A consolidated three-way audit (Fortran NSM1 in `modGlobalParam.f90` and
sibling files vs. v1 Python NSM1 vs. v3 Python NSM1) was performed across
five Process families and the shared utilities in early May 2026. Findings:

- 64 total findings deduplicated across five sub-audits.
- 22 critical correctness defects identified in v3 (algae wiring, Nitrogen
  wiring, rca/rcb stoichiometric ratios used as raw weights, missing CBOD
  source on DIC budget, missing DOX-Monod attenuation on SOD, phantom
  ammonium decay term, broken NO3 algal-uptake split, fdp partitioning
  unit factor, two parameter typos).
- 8 v3 deliberate improvements vindicated by the Fortran-coded defaults
  (sentinel-999 corrections, SOD_tc / PAR architectural splits, denit-source
  closure on N2 budget, etc.).
- 6 three-way disagreements where v3 inherits v1 values that differ from
  Fortran; flagged for LimnoTech reconciliation.

Phases 9.A, 9.B, and 9.C resolved 22 of the 26 critical items:

- **Phase 9.A** (algae + Nitrogen wiring sweep + 3 algae formula bugs):
  rewired 13 algae kinetic methods and 5 Nitrogen kinetic methods to read
  the `DEFAULTS` dict from `parameters/<group>.py` rather than legacy v2
  kwargs (which defaulted to 0.0 or 1.0, masking the wiring at default
  instantiation); fixed FloatingAlgae light-limit option-1 parenthesization,
  harmonic-mean zero-guard, and BenthicAlgae Steele exponent sign.
- **Phase 9.B** (Carbon + DOX rca/rcb stoichiometry + missing terms +
  fdp): replaced raw-weight references (`self.AWc=40`, `self.BWc=40`)
  with the correctly derived ratios `rca = AWc/AWa = 0.04 mg-C/ug-Chla`
  and `rcb = BWc/BWd = 0.4 mg-C/mg-D` at all eight call sites in
  `carbon.py` and `dox.py` (1000x error for floating algae, 100x for
  benthic before fix); added the missing CBOD oxidation source to the
  DIC budget; restored DOX-Monod attenuation `DOX/(DOX+KsSod)` on the
  SOD sink; dropped the spurious DOX attenuation from POC hydrolysis;
  corrected the `fdp` partitioning denominator from `0.000001` to
  `1.0E6` (latent at default `kdpo4=0`, breaks under sorption).
- **Phase 9.C** (parameter reconciliation + docs cleanup): corrected two
  inherited / drift-related parameter values (`vson_20` 0.1 → 0.01;
  `lambdam` 0.0174 → 0.174); reconciled `utils/reaeration.py` author
  attributions against Fortran source comments; relocated NSM1 docs
  from `docs/` to `design/`.

The audit summary at
`design/clearwater_modules_v3_nsm1_audit_summary.md` is the consolidated
reference; the five sub-audits in
`design/clearwater_modules_v3_nsm1_audit_*.md` carry the per-Process detail.

---

## 3. Scope: 16 constituents, 11 Process classes

| Process class    | State variables                  | Source                                            |
| ---              | ---                              | ---                                               |
| `FloatingAlgae`  | Ap                               | Extends v2 `FloatingAlgae` (Phase 2 + 9.A fixes)  |
| `BenthicAlgae`   | Ab                               | Extends v2 `BenthicAlgae` (Phase 2 + 9.A fixes)   |
| `Nitrogen`       | NH4, NO3, OrgN                   | Extends v2 `Nitrogen` (Phase 2 + 9.A fixes; OrgN added) |
| `Phosphorus`     | TIP, OrgP                        | New in v3                                         |
| `Carbon`         | POC, DOC, DIC                    | New in v3 (Phase 9.B stoichiometry fix)           |
| `POM`            | POM                              | New in v3                                         |
| `CBOD`           | CBOD (multi-group)               | New in v3                                         |
| `DOX`            | DOX                              | New in v3 (Phase 9.B stoichiometry + SOD-Monod fix) |
| `Pathogen`       | PX                               | New in v3                                         |
| `Alkalinity`     | Alk                              | New in v3 (declared but inactive in v1)           |
| `N2`             | N2, TDG                          | New in v3                                         |

11 Process classes covering 16 state variables. Process source files live at
`src/clearwater_modules_v3/processes/<name>.py`. Default parameter dictionaries
live at `src/clearwater_modules_v3/parameters/<group>.py`.

---

## 4. Quick-start

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

## 5. Architecture

### Resolved design decisions (spec Section 14)

These are decided and implemented; LimnoTech confirmation closes the three
items still flagged tentative.

- **Parameter library: DEFAULTS-merge pattern.** Each `Process` class imports
  a `DEFAULTS: dict[str, float]` from `parameters/<group>.py`. At
  construction time it merges `{**DEFAULTS, **user_parameters}` and stores
  results as `self.<name>` attributes. YAML overrides are partial; defaults
  fill in unspecified keys.
  - **Phase 9.A correction:** the kinetic methods now read the DEFAULTS
    keys (e.g., `self.mu_max_20`, `self.knit_20`) instead of the legacy
    v2 kwargs that previously defaulted to 0.0 or 1.0. Pre-Phase-9.A,
    default-instantiated Processes silently produced zero-respiration,
    zero-death, zero-settling algae and 5x-to-500x-inflated nitrogen
    rates; the regression tests didn't catch this because every existing
    parity test passed the legacy kwargs explicitly with v1/Fortran
    values. See migration notes Section 3.1.
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

## 6. Tier 1 conservation contract

Every Process passes a closed-system mass-conservation test at
`rtol=1e-12` with **zero clip events**. The closed-system harness
(`tests/v3/nsm1/test_validation_tier1_conservation.py`) zeros all boundary
fluxes, all settling velocities, and all sediment fluxes; runs each Process
in isolation; and asserts that the relevant elemental total (N, P, C,
O2-equivalents, Alk-equivalents) is conserved to roundoff and that
`model.diagnostics.clip_events == 0`. A clip event under closed-system
conditions indicates either unphysical parameters or a malformed test —
the test fails in that case, which is the correct diagnostic.

Per-process Tier 1 tests live at `tests/v3/nsm1/test_<constituent>_tier1.py`.

---

## 7. v1↔v3 deviations summary

v3 NSM1 reproduces v1 NSM1 kinetics within floating-point tolerance for the
overwhelming majority of sub-rate terms. Two classes of deviation exist:

**Deliberate runtime numerical deviations** (architectural choices and unit
or convention cleanups):

- **Carbon POC hydrolysis** — pre-9.B v3 added a DOX-Monod attenuation that
  v1 and Fortran do not apply. Phase 9.B removed the attenuation; v3 now
  matches v1 / Fortran (`kpoc_tc * POC` only).
- **DOX SOD** — v3's `SOD_tc` utility is pure Arrhenius; v1 has Monod
  inline (architectural choice). Phase 9.B added the `DOX/(DOX+KsSod)`
  Monod factor at the consumer site so the hypoxic-attenuation behavior
  matches Fortran while preserving the cleaner architectural split.
- **Alkalinity DOX-attenuation** — v3 reads pre-attenuated nitrification and
  denitrification fluxes from the Nitrogen rate cache; v1 applies the
  attenuation locally inside `Alk_nitrification` / `Alk_denitrification`.
  Same numerical answer when tests pass v3-aligned upstream rates.
- **Pathogen light decay** — v3 uses `PAR = q_solar * Fr_PAR`; v1 uses raw
  `q_solar`. Calibration target `apx` absorbs the difference.
- **CBOD sedimentation** — v3 applies `ksbod_tc / depth` (m/d → 1/d); v1
  treats `ksbod_tc` directly as 1/d.
- **Kelvin conversion** — v3 utility re-exports v2's `+273.16` form for
  v2-parity test stability; v1 NSM1 `nsm1/processes.py` and Fortran NSM1
  use `+273.15`. The 0.01 K offset propagates weakly into Henry's-law
  saturation.
- **Pressure mb→atm** — v3 uses `1.0 / 1013.25`; v1 uses literal
  `0.000986923`. Agreement to ~7 significant figures.

**Audit-driven default-value corrections (Phase 9.B, 9.C):**

- **Carbon and DOX algal coupling: `rca` and `rcb` stoichiometric ratios.**
  Phase 9.B replaced raw-weight uses (`self.AWc=40`, `self.BWc=40`) with
  the correctly derived ratios `rca = AWc/AWa = 0.04` and `rcb = BWc/BWd
  = 0.4`. Pre-fix DIC and O2 algal coupling was off by ~1000x for floating
  algae and ~100x for benthic algae. The v1 Python source defines the
  helper functions to derive `rca`/`rcb` this way; Fortran derives them
  the same. The v3 1.0.0 release uses the derived values.
- **`vson_20` 0.1 → 0.01 (Phase 9.C).** Internal v3 inconsistency: the
  nitrogen-group default had drifted to 0.1 m/d while v3's own
  `global_vars.vson`, v1 `GlobalVars.vson`, and Fortran `vson` were all
  0.01 m/d.
- **`lambdam` 0.0174 → 0.174 (Phase 9.C).** Likely v1 typo (10x too
  small). Fortran and QUAL2K Table 6 use 0.174; multiple legacy v1 NSM
  test fixtures already overrode the v1 default to 0.174, indirectly
  confirming the v1 default was wrong.

For per-test docstring detail on the deliberate deviations:
`tests/test_5_carbon_calculations_v2.py`,
`test_5_dox_calculations_v2.py`, `test_5_alkalinity_calculations_v2.py`,
`test_5_pathogen_calculations_v2.py`, `test_5_cbod_calculations_v2.py`,
`test_5_n2_calculations_v2.py`. The full corrections record is at
`src/clearwater_modules_v3/parameter_defaults_corrections.md`.

---

## 8. DOX reaeration default note (Phase 7.C Item 3 / Phase 9.C audit)

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

### 8.1 v3 vs Fortran NSM1 default-DOX-recovery divergence (important)

The Phase 9.C three-way audit
(`design/clearwater_modules_v3_nsm1_audit_utilities_params.md`) found that
v3's `kah_20_user=0.0` disagrees with **Fortran NSM1**'s default of
`kah_20_user=1.0` (`modGlobalParam.f90:113`):

* At default `hydraulic_reaeration_option=1` (the user-override branch),
  v3 produces **zero** atmospheric hydraulic reaeration.
* At the same setting, Fortran NSM1 produces **1.0 1/d** hydraulic
  reaeration.
* Side-by-side runs of v3 vs Fortran NSM1 with all-default settings will
  therefore show DOX recovery in Fortran but not in v3.

Users who want non-zero default reaeration matching Fortran's behavior
should either:

  * (a) explicitly set `kah_20_user > 0` (e.g., `kah_20_user=1.0` to mimic
    Fortran's default), **or**
  * (b) select a different `hydraulic_reaeration_option` from the menu
    (options 2-9 use empirical formulas based on velocity, depth, flow,
    topwidth, slope, or shear velocity, and **do not consult**
    `kah_20_user`).

Author attributions for the 9 hydraulic and 13 wind reaeration options
were reconciled in Phase 9.C against the Fortran source comments
(`modGlobalParam.f90:268-414`); see `src/clearwater_modules_v3/utils/reaeration.py`
docstrings for the canonical attributions.

See `src/clearwater_modules_v3/parameter_defaults_corrections.md` Sections
1.5 and 1.6 for the rationale and audit trail.

---

## 9. v1↔v3 high-level deviations beyond DOX reaeration

For the full bug-fix inventory and parameter corrections that v3 carries
forward from the port, read:

- Bug list — `design/clearwater_modules_v3_nsm1_design_specification.md`
  Section 6 (16 v2 bugs fixed at the port).
- Default-value corrections —
  `src/clearwater_modules_v3/parameter_defaults_corrections.md` Section 1
  (9 corrections: `vsop`, `vs`, `SOD_20`, `SOD_theta`, `kah_20_user`,
  `kaw_20_user`, `pressure_mb`, plus Phase 9.C `vson_20` and `lambdam`).
- Lower-priority audit findings —
  `src/clearwater_modules_v3/parameter_defaults_corrections.md` Section 2
  (8 items kept at v1 values pending Phase 1.3 audit).
- Phase 7 numerical deviations —
  `src/clearwater_modules_v3/parameter_defaults_corrections.md` Section 3
  (the 7-item list summarized above).
- Audit-driven corrections —
  `src/clearwater_modules_v3/parameter_defaults_corrections.md` Section 4
  (the LimnoTech reconciliation list).

For migration from v1 NSM1 (`NutrientBudget`) or v2 NSM1
(`clearwater_modules_v2.processes.nitrogen` etc.), see
`design/clearwater_modules_v3_nsm1_migration.md`.

---

## 10. Test counts

- **705 tests passing, 3 xfailed** across the v3 NSM1 suite (Phases 0-9.C).
  Phase 8 reported 652 passing; Phase 9.A/B/C added approximately 35
  Fortran-anchored regression tests covering the audit findings (algae
  default-instantiation regressions, rca/rcb numerical references,
  NH4+NO3 mass-balance closure, hypoxic SOD attenuation, CBOD-DIC source,
  the `kdpo4>0` partitioning regime, and the `vson_20` / `lambdam`
  parameter pinnings).
- **8 Tier 1 conservation tests** under `tests/v3/nsm1/test_<constituent>_tier1.py`
  plus the consolidated harness
  `tests/v3/nsm1/test_validation_tier1_conservation.py`.
- **36 sub-rate v1-parity tests** in `tests/test_5_*_calculations_v2.py`,
  each pinning a single v1 sub-term against v3's cached rate variable.
- **1 coupled end-to-end demo** at `examples/V3/04_Example_NSM1.ipynb`
  exercising all 11 Process classes plus Riverine on a synthetic mesh.

---

## 11. Open follow-up items

The audit identified four items that remain open after Phases 9.A/B/C and
are flagged for LimnoTech reconciliation. Each is documented in detail in
`src/clearwater_modules_v3/parameter_defaults_corrections.md` Section 4.

1. **Nitrogen Arrhenius theta values — possible 4-way swap with Fortran.**
   `kon_theta`, `kdnit_theta`, `rnh4_theta`, `vno3_theta` v3/v1 values
   appear to swap with Fortran's values pairwise (4-way pattern). v3
   inherits v1; flagged with `FIXME(phase9c-audit):` inline comments in
   `parameters/nitrogen.py`. Recommend LimnoTech confirm the canonical
   QUAL2K-Table reference before changing.
2. **`BWa` benthic-algae chlorophyll-a stoichiometry** — Fortran 5000 vs
   v1/v3 3500 (~43% difference). Both within the literature range; flagged
   for confirmation.
3. **`vsop` value** — Fortran 0.01 m/d vs v3 0.1 m/d (10x). v3 chose
   mid-range over Fortran's lower bound; both defensible.
4. **`SOD_20` value** — Fortran 0.2 g-O2/m^2/d vs v3 1.0 g-O2/m^2/d (5x).
   v3 chose conservative midpoint over Fortran's lower bound; both
   defensible from Chapra (1997).

Two further items disclosed but not blocking:

- **`kah_20_user=0.0` (v3) vs `1.0` (Fortran).** Documented in Section 8.1
  above as a deliberate v3 design choice with an explicit behavioral note.
- **DOX salinity correction on `O2sat`.** Fortran applies a salinity
  correction; v1 omits; v3 inherits the v1 omission. Freshwater impact is
  zero; estuarine impact up to ~20%. Defer if v3 1.x targets fresh water
  only.
- **DIC unit reconciliation (mol/L vs mg/L).** Long-standing v1/Fortran
  mismatch in the dDIC/dt unit factors; v3 inherits the v1 form. Escalated
  for LimnoTech review; decision is whether to land a unit cleanup in
  1.0.x or defer to the v3 1.1 carbonate solver.

---

## 12. References

- `design/clearwater_modules_v3_nsm1_design_specification.md` — full
  source-of-truth design (Sections 6 bug list, 7 defaults, 8 migration,
  9 testing, 14 design decisions).
- `design/clearwater_modules_v3_architecture_specification.md` — umbrella
  architecture, package layout, integrator-pattern contract.
- `design/clearwater_modules_v3_nsm1_gap_analysis.md` — Phase 0 readiness
  assessment.
- `design/clearwater_modules_v3_nsm1_migration.md` — v1→v3 and v2→v3
  migration paths with side-by-side examples.
- `design/clearwater_modules_v3_nsm1_limnotech_review.md` — packaged
  review materials for external review.
- `design/clearwater_modules_v3_nsm1_audit_summary.md` — consolidated
  three-way audit (Fortran NSM1 vs v1 Python NSM1 vs v3 Python NSM1).
- `design/clearwater_modules_v3_nsm1_audit_algae.md`,
  `..._n_p.md`, `..._c_dox.md`, `..._simple_constituents.md`,
  `..._utilities_params.md` — the five sub-audits.
- `src/clearwater_modules_v3/parameter_defaults_corrections.md` — every
  default-value correction and runtime-numerical deviation v3 carries
  relative to v1, including the LimnoTech escalation list (Section 4).
- `examples/V3/04_Example_NSM1.ipynb` — coupled end-to-end demo.
