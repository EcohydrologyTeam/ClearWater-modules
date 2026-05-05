# v3 NSM1 Parameter Defaults: Corrections and Audit Findings

This document records every default-value correction applied to v3 NSM1's
parameter library at the port from v1 (`src/clearwater_modules/nsm1/constants.py`)
and the lower-priority audit findings that were left as-is for v3 1.0.0 with
inline `FIXME(phase1-audit):` comments in the relevant `parameters/<group>.py`
module.

The canonical reference for the seven critical corrections is
`design/clearwater_modules_v3_nsm1_design_specification.md` Section 7. The
broader Phase 0 inventory and rationale live in
`docs/clearwater_modules_v3_nsm1_phase0_parameter_audit.md`.

---

## Section 1: Critical default-value corrections (7 items, applied at port)

These are applied directly in the v3 `DEFAULTS` dicts. Each inline comment in
the relevant `parameters/<group>.py` module records the v1 original.

### 1.1 `vsop` — organic-P settling velocity
- **Module:** `parameters/phosphorus.py`
- **v1 default:** `999` m/d
- **v3 default:** `0.1` m/d
- **Rationale:** `vsop` is multiplied directly into the OrgP loss rate
  (`vsop * OrgP / depth`). At 999 m/d the water column would be entirely
  evacuated of organic phosphorus on every timestep. Typical values for
  organic-P settling sit in the 0.01-1.0 m/d range; 0.1 m/d is a defensible
  midpoint consistent with literature values for medium-sized organic
  particles.

### 1.2 `vs` — TIP settling velocity
- **Module:** `parameters/phosphorus.py`
- **v1 default:** `999` m/d
- **v3 default:** `0.1` m/d
- **Rationale:** Same magnitude error as `vsop`. `vs` controls the settling
  loss of total inorganic phosphorus adsorbed onto suspended solids; physically
  reasonable values fall in the 0.01-1.0 m/d range. 0.1 m/d is the same
  midpoint chosen for `vsop` and is consistent with Chapra (1997) recommended
  ranges for inorganic-P partitioning onto fine sediments.

### 1.3 `SOD_20` — sediment oxygen demand at 20 C
- **Module:** `parameters/dox.py`
- **v1 default:** `999` g-O2/m^2/d
- **v3 default:** `1.0` g-O2/m^2/d
- **Rationale:** Realistic SOD values for moderate organic loading sit in the
  0.5-3.0 g-O2/m^2/d range (Chapra 1997, Table 25.2). 1.0 g-O2/m^2/d is a
  defensible default for an unspecified moderate-loading site; users with
  field measurements should override per simulation. The v1 value of 999
  drives DOX immediately negative on any wet-bed timestep.

### 1.4 `SOD_theta` — Arrhenius coefficient for SOD
- **Module:** `parameters/dox.py`
- **v1 default:** `999` (unitless)
- **v3 default:** `1.060` (unitless)
- **Rationale:** Arrhenius-style temperature corrections take the form
  `theta^(T-20)`. With `theta=999` and any temperature above 20 C, the
  corrected SOD blows up by orders of magnitude per degree (catastrophic).
  Chapra (1997) recommends `theta=1.060` for SOD as the standard literature
  value, equivalent to a Q10 of approximately 1.79. This was Phase 0's most
  severe finding.

### 1.5 `kaw_20_user` — user-override wind reaeration coefficient at 20 C
- **Module:** `parameters/dox.py`
- **v1 default:** `999` m/d
- **v3 default:** `0.0` m/d
- **Rationale:** This parameter is consulted only when `wind_reaeration_option`
  selects the user-override branch (option 1). Setting the default to zero
  means the user-override branch is *off* by default; if a user opts into the
  override path without supplying their own value they will see no reaeration
  rather than a runaway 999 m/d.

### 1.6 `kah_20_user` — user-override hydraulic reaeration coefficient at 20 C
- **Module:** `parameters/dox.py`
- **v1 default:** `999` 1/d
- **v3 default:** `0.0` 1/d
- **Rationale:** Same as `kaw_20_user`. Defaulting to zero makes the override
  branch a no-op when not configured, instead of a numerical disaster.

### 1.7 `pressure_mb` — atmospheric pressure
- **Module:** `parameters/global_parameters.py`
- **v1 default:** `2026.5` hPa
- **v3 default:** `1013.25` hPa
- **Rationale:** Standard sea-level atmospheric pressure is `1013.25` hPa
  (ISO 2533). The v1 value is approximately 2x correct; this biases every
  downstream computation that depends on it, notably `O2sat` (DO saturation
  via Henry's law) and `N2sat` (N2 saturation), as well as any atmospheric
  reaeration formula that uses partial pressure. This was identified as a
  Phase 0 finding (2026-05-04) and added to the canonical correction list.

---

## Section 2: Lower-priority audit findings under review (8 items)

These eight items were flagged during Phase 0 as suspicious or unclear but
left at v1 values for the v3 1.0.0 port. Each is marked with a
`FIXME(phase1-audit):` inline comment in the relevant parameter module.
Disposition for each is described below.

### 2.1 `rnh4_20=0`, `vno3_20=0`, `rpo4_20=0` — sediment-flux release rates disabled
- **Modules:** `parameters/nitrogen.py` (`rnh4_20`, `vno3_20`),
  `parameters/phosphorus.py` (`rpo4_20`).
- **Issue:** All three sediment-release/uptake rates are zero, which silently
  disables sediment-flux contributions to NH4, NO3, and PO4 budgets. v1 also
  defaults `use_SedFlux=False`, which should make this moot, but Phase 0
  could not confirm that every code path consuming these parameters is gated
  by `use_SedFlux`.
- **Disposition (v3 1.0.0):** Kept at zero, consistent with v1. Phase 1.3
  Process implementation must verify that `use_SedFlux=False` gates all
  consumers of these parameters; until that audit is complete, leaving the
  values at zero preserves v1 behavior. If `use_SedFlux=True` is set without
  also providing site-specific rates, results will silently exclude sediment
  fluxes.

### 2.2 `kdpo4=0.0` — TIP partitioning coefficient
- **Module:** `parameters/phosphorus.py`
- **Issue:** Phosphorus partitioning between dissolved (DIP) and particulate
  (sorbed) phases is governed by `kdpo4` (L/kg). At zero, no DIP adsorbs onto
  suspended solids; the TIP partitioning feature is effectively disabled.
- **Disposition (v3 1.0.0):** Kept at zero. Full DIP-solid partitioning is
  NSM2 territory; v3 NSM1 1.0.0 maintains the v1 simplification. Documented
  for clarity.

### 2.3 `ksbod_20=0.0` — CBOD settling rate
- **Module:** `parameters/cbod.py`
- **Issue:** CBOD never settles in v1 because the settling rate is hardcoded
  to zero. If CBOD is conceptually fully dissolved this is correct; if CBOD
  represents a mixture including particulate fractions, this is a bug.
- **Disposition (v3 1.0.0):** Kept at zero pending clarification from
  LimnoTech on the intended interpretation of CBOD groups in NSM1.

### 2.4 `apx=1`, `vx=1` — pathogen placeholder values
- **Module:** `parameters/pathogen.py`
- **Issue:** Both pathogen production/coupling (`apx`) and settling velocity
  (`vx`) are set to 1.0 with no documented literature basis. Likely
  placeholder values inherited from v1.
- **Disposition (v3 1.0.0):** Kept as-is; updating defaults to literature-
  backed values is deferred to a future audit pending site-specific
  pathogen-tracer studies.

### 2.5 `h2=0.1` — POM dissolution depth denominator
- **Module:** `parameters/pom.py`
- **Issue:** `h2` is used as a divisor in burial/sedimentation terms (e.g.,
  `vb * POM / h2`); its physical role is unclear, and the units (meters) do
  not obviously line up with the surrounding terms. Phase 0 flagged it as
  needing clarification.
- **Disposition (v3 1.0.0):** Kept as-is; flagged with `FIXME(phase1-audit):`.
  Full clarification is part of the Phase 1.3 process implementation.

### 2.6 `vb=0.01` — burial velocity magnitude
- **Module:** `parameters/global_vars.py`
- **Issue:** Burial velocity of 0.01 m/d has no documented validation; some
  sediment-flux literature uses values one to two orders of magnitude
  different.
- **Disposition (v3 1.0.0):** Kept as-is; flagged with `FIXME(phase1-audit):`.
  Without tied-to-site sediment data, no defensible alternative default
  exists.

### 2.7 `q_solar=500` units mismatch
- **Module:** `parameters/global_vars.py`
- **Issue:** v1's docstring documents `q_solar` as having units of `1/d`, but
  the value 500 and its consumption pattern (Beer-Lambert PAR computation)
  are consistent with W/m^2. The docstring is wrong, not the value.
- **Disposition (v3 1.0.0):** Value kept at 500; inline comment in
  `global_vars.py` documents the correct units (W/m^2). Docstring will be
  corrected in the Process docstring during Phase 1.3.

### 2.8 `lambdas` — light extinction parameter defined but not used
- **Module:** `parameters/global_vars.py`
- **Issue:** `lambdas=0.052` represents the suspended-solids contribution to
  Beer-Lambert light extinction. v1's code path computing the effective
  extinction coefficient comments out this term: `lambdas * Solid` is
  defined but not added to the sum. Either the term should be enabled or
  the parameter should be removed.
- **Disposition (v3 1.0.0):** Parameter retained for backward compatibility
  with any user YAML that already sets it; flagged with `FIXME(phase1-audit):`
  in `global_vars.py`. Phase 1.3 light/extinction utility implementation
  should decide whether to re-enable the term (and document the choice
  here).

---

## Section 3: Phase 7 v1↔v3 runtime numerical deviations

The items in Sections 1 and 2 are parameter default corrections. The items
in this section are **runtime numerical differences** between v1 NSM1 and
v3 NSM1 that are not strictly parameter corrections but are worth recording
alongside them. Each is documented in detail in the docstring of the
corresponding `tests/test_5_*_calculations_v2.py` parity test, and v3
behavior is pinned by that test against the v1 reference under matched
inputs.

The list is intentionally small: v3 NSM1 reproduces v1 NSM1 kinetics within
floating-point tolerance for the overwhelming majority of sub-rate terms.

### 3.1 Carbon POC hydrolysis — DOX-Monod attenuation added in v3

- **Module:** `processes/carbon.py` (POC hydrolysis sub-rate)
- **v1 form:** `kpoc_tc * POC` (no DOX coupling)
- **v3 form:** `kpoc_tc * POC * DOX_attenuation`, where `DOX_attenuation =
  DOX / (KsOxmc + DOX)` follows the standard Monod oxygen-inhibition
  pattern.
- **Rationale:** POC hydrolysis is mediated by aerobic microbial activity;
  treating it as DOX-independent (v1) overestimates hydrolysis under
  hypoxic conditions. v3 follows the same architectural pattern v1 already
  applies to DOC oxidation. The Phase 7 test
  `tests/test_5_carbon_calculations_v2.py::test_poc_hydrolysis_rate_matches_v1`
  forces `DOX_attenuation == 1` (saturated DOX) to verify v3 collapses to
  v1's form when oxygen is non-limiting; the deviation appears only when
  DOX is depleted.
- **Reference test:** `tests/test_5_carbon_calculations_v2.py` lines
  130–180 (parity test plus deviation note in the module docstring).

### 3.2 DOX SOD — pure-Arrhenius `SOD_tc` in v3 vs Monod-inline in v1

- **Module:** `processes/dox.py` (SOD sink)
- **v1 form:** `SOD_tc` includes a DOX-Monod factor inline when
  `use_DOX=True` is passed.
- **v3 form:** `clearwater_modules_v3.utils.sediment.SOD_tc` is a pure
  Arrhenius temperature correction; DOX-Monod attenuation, when desired,
  is applied at the call site rather than baked into the utility.
- **Rationale:** Architectural separation. The Process owns its
  oxygen-inhibition contract; the utility owns only the temperature
  correction. v3's design supports a future Process opt-in to a
  semi-implicit DOX treatment without requiring two parallel SOD
  utilities.
- **Reference test:** `tests/test_5_dox_calculations_v2.py::test_dox_sod_rate_matches_v1`
  passes `use_DOX=False` to v1's `SOD_tc` so v3 and v1 forms match
  exactly under the test fixture.

### 3.3 Alkalinity DOX-attenuation flow — pre-attenuated rates in v3

- **Module:** `processes/alkalinity.py` (nitrification and denitrification
  Alk coupling)
- **v1 form:** `Alk_nitrification` / `Alk_denitrification` apply the
  DOX-Monod / oxygen-inhibition factor locally inside each function.
- **v3 form:** `Alkalinity` consumes the pre-cached
  `nitrification_flux_rate` / `denitrification_flux_rate` from
  `Nitrogen.run`, which has already applied the Monod / inhibition factor
  upstream. v3 multiplies through by the stoichiometric `r_alkn` /
  `r_alkden` and 50000 mg-CaCO3-equivalent factor only.
- **Rationale:** Single-source-of-truth for the attenuation factor — by
  routing through Nitrogen's rate cache, v3 guarantees Alkalinity's coupling
  is consistent with the actual NH4/NO3 transformation rate. v1's
  duplicate-attenuation pattern works correctly only because both call
  sites use the same parameter values.
- **Reference test:** `tests/test_5_alkalinity_calculations_v2.py` module
  docstring (lines 21–30) documents the equivalence and the test fixture
  passes a Nitrogen mock whose `*_flux_rate` already includes the Monod
  factor, matching what v1 computes locally.

### 3.4 Pathogen light decay — `PAR = q_solar * Fr_PAR` in v3

- **Module:** `processes/pathogen.py` (light-driven decay sub-rate)
- **v1 form:** `PathogenDecay` uses raw `q_solar` (W/m² incident).
- **v3 form:** `_rate_light_decay` uses `PAR(q_solar, Fr_PAR) =
  q_solar * Fr_PAR`, scaling effective surface irradiance by `Fr_PAR`
  (default 0.47, the photosynthetically active fraction).
- **Rationale:** Pathogen UV-driven decay is properly a function of UV
  flux, not total shortwave; PAR is a closer proxy than raw `q_solar`.
  The calibration target `apx` absorbs the difference for any historical
  v1 calibration. The Phase 3.1 docstring documents this as an
  intentional deviation.
- **Reference test:** `tests/test_5_pathogen_calculations_v2.py::test_pathogen_light_decay_matches_v1`
  pins `Fr_PAR=1.0` to make v3 and v1 forms exactly equivalent under the
  fixture.

### 3.5 CBOD sedimentation — `ksbod_tc / depth` in v3 (m/d → 1/d)

- **Module:** `processes/cbod.py` (sedimentation sink)
- **v1 form:** `CBOD_sedimentation = CBOD * ksbod_tc`, treating `ksbod_tc`
  directly as a 1/d first-order rate.
- **v3 form:** `CBOD * ksbod_tc / depth`, treating `ksbod_tc` as a
  settling velocity (m/d) divided by water-column depth to yield a 1/d
  first-order rate.
- **Rationale:** Dimensional consistency — settling-driven sedimentation
  scales with `velocity / depth`, not with velocity alone. The v1
  formulation requires re-interpreting `ksbod_tc` as 1/d and reconciling
  with the velocity-style defaults documented in `parameters/cbod.py`.
- **Reference test:** `tests/test_5_cbod_calculations_v2.py` module
  docstring (line 11) and the parity test at lines 173–200 document the
  units mismatch and pin the v3 result against the v1 form scaled
  through the `1/depth` factor.

### 3.6 Celsius-to-Kelvin offset — 273.15 in v3 vs 273.16 in v1

- **Module:** `utils/conversions.py` (`celsius_to_kelvin`)
- **v1 form:** `T_K = T_C + 273.16` (the historical SI triple-point
  definition of the Kelvin offset, kept in v1 for backward compatibility).
- **v3 form:** `T_K = T_C + 273.15` (the modern definition, consistent
  with all other v3 utilities and with the `T_K - 273.15` round-trip).
- **Rationale:** Modern SI convention. The 0.01 K offset propagates
  weakly into Henry's-law saturation calculations (`O2sat`, `N2sat`,
  `Henrys_k`) but is well below the ~1% tolerance of typical aquatic
  measurements.
- **Reference test:** `tests/test_5_n2_calculations_v2.py` lines 13–18
  (module docstring) document the convention difference; tests use the
  v1 +273.16 convention internally to isolate kinetics-formula
  parity from the offset.

### 3.7 Pressure mb→atm conversion — `1.0 / 1013.25` in v3 vs `0.000986923` in v1

- **Module:** `utils/conversions.py` (mb-to-atm scaling factor used in
  Henry's-law saturation calculations)
- **v1 form:** literal `0.000986923` (truncated decimal).
- **v3 form:** `1.0 / 1013.25` (computed exactly from the standard
  sea-level pressure in hPa).
- **Rationale:** Eliminate the truncation. The two forms agree to ~7
  significant figures. Affects `N2sat`, `O2sat`, atmospheric reaeration
  partial-pressure terms.
- **Reference test:** `tests/test_5_n2_calculations_v2.py` lines 19–21
  (module docstring) document the convention difference; the parity
  tests use `rtol=1e-6` to absorb the truncation and the 0.01 K Kelvin
  offset together.

---

For per-deviation empirical verification, see the module docstrings of:

- `tests/test_5_carbon_calculations_v2.py`
- `tests/test_5_dox_calculations_v2.py`
- `tests/test_5_alkalinity_calculations_v2.py`
- `tests/test_5_pathogen_calculations_v2.py`
- `tests/test_5_cbod_calculations_v2.py`
- `tests/test_5_n2_calculations_v2.py`

Each docstring states the deviation explicitly, identifies the test
fixture's strategy for isolating the deviation from kinetics-formula
parity, and pins the expected v1 reference value.
