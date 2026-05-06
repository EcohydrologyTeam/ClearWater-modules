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

## Section 1: Critical default-value corrections (11 items, applied at port)

These are applied directly in the v3 `DEFAULTS` dicts. Each inline comment in
the relevant `parameters/<group>.py` module records the v1 original. Items 1.1
through 1.7 are the original 7 sentinel-999 / 2026.5 corrections applied at the
v1->v3 port. Items 1.8 and 1.9 were added in Phase 9.C after the three-way
v1/v3/Fortran audit (`design/clearwater_modules_v3_nsm1_audit_utilities_params.md`)
identified them as additional v1 flaws or internal inconsistencies.

**Three-way verification summary (v1 / Fortran modGlobalParam.f90 / v3):**

| Section | Parameter | v1 default | Fortran default | v3 default | Verdict |
|---|---|---|---|---|---|
| 1.1 | `vsop` | 999 (sentinel) | 0.01 | 0.1 | v3 chose mid-range over Fortran's lower bound |
| 1.2 | `vs` | 999 (sentinel) | 0.1 | 0.1 | v3 matches Fortran |
| 1.3 | `SOD_20` | 999 (sentinel) | 0.2 | 1.0 | v3 chose conservative midpoint over Fortran's lower bound |
| 1.4 | `SOD_theta` | 999 (sentinel) | 1.06 | 1.060 | v3 matches Fortran (Chapra 1997) |
| 1.5 | `kaw_20_user` | 999 (sentinel) | 0.0 | 0.0 | v3 matches Fortran |
| 1.6 | `kah_20_user` | 999 (sentinel) | 1.0 | 0.0 | v3 disables user-override branch by default; behavioral note below |
| 1.7 | `pressure_mb` | 2026.5 (2x error) | n/a (Fortran uses pressure_atm) | 1013.25 | v3 ISO 2533 standard |
| 1.8 | `vson_20` | 0.01 (in v1 GlobalVars as `vson`) | 0.01 | 0.01 (was 0.1, Phase 9.C fix) | v3 corrected to match Fortran/v1 |
| 1.9 | `lambdam` | 0.0174 (likely v1 typo) | 0.174 | 0.174 (was 0.0174, Phase 9.C fix) | v3 corrected to match Fortran/QUAL2K Table 6 |

For 4 of the original 7 (`vs`, `SOD_theta`, `kaw_20_user`, `pressure_mb`) v3
chose values matching or equivalent to Fortran. For 3 of 7 (`vsop`, `SOD_20`,
`kah_20_user`) v3 chose values differing from Fortran by O(1)-O(10) but
defensible from literature; rationale recorded per-item below. The Phase 9.C
additions (1.8, 1.9) bring v3 into closer Fortran/QUAL2K alignment.

### 1.1 `vsop` — organic-P settling velocity
- **Module:** `parameters/phosphorus.py`
- **v1 default:** `999` m/d (sentinel)
- **Fortran default:** `0.01` m/d (`modGlobalParam.f90:98`)
- **v3 default:** `0.1` m/d
- **Rationale:** `vsop` is multiplied directly into the OrgP loss rate
  (`vsop * OrgP / depth`). At 999 m/d the water column would be entirely
  evacuated of organic phosphorus on every timestep. The literature
  range cited for OrgP settling (Chapra 1997, QUAL2K manual section
  5.5.16) is 0.01-1.0 m/d, but the QUAL2K manual leaves `vop` as a
  calibration parameter without a pinned default.

  **Physical-consistency basis for v3 = 0.1 m/d:** organic-P in NSM1
  originates predominantly from dead-algae detritus (via the algal
  mortality routing `algal_orgp_from_mortality_rate` and the benthic
  `balgae_orgp_from_mortality_rate`). The algae from which OrgP
  derives have a settling velocity of `vsap = 0.15` m/d (Fortran/v1/v3
  agree on this value, the universal NSM1 default). For internal
  consistency, OrgP detritus inherited from the same algal pool
  should settle at a comparable rate — typical 0.1-0.2 m/d. v3's
  `vsop = 0.1` m/d is consistent with this physical-consistency
  argument; Fortran's 0.01 m/d is 15x slower than the algae from
  which the OrgP derives, which is implausible for a representative
  default. v3 deliberately deviates from Fortran here on
  physical-consistency grounds. Users with site-specific data
  (especially fine-detritus or colloidal organic-P regimes) should
  override per project. Regression coverage in
  `tests/test_5_phosphorus_calculations_v2.py::test_phase9e_vsop_consistent_with_vsap`.

### 1.2 `vs` — TIP settling velocity
- **Module:** `parameters/phosphorus.py`
- **v1 default:** `999` m/d (sentinel)
- **Fortran default:** `0.1` m/d (`modGlobalParam.f90:87`)
- **v3 default:** `0.1` m/d
- **Rationale:** Same magnitude sentinel error as `vsop` in v1. `vs` controls
  the settling loss of total inorganic phosphorus adsorbed onto suspended
  solids. v3's 0.1 m/d matches Fortran exactly and falls within the typical
  literature range (Chapra 1997).

### 1.3 `SOD_20` — sediment oxygen demand at 20 C
- **Module:** `parameters/dox.py`
- **v1 default:** `999` g-O2/m^2/d (sentinel)
- **Fortran default:** `0.2` g-O2/m^2/d (`modGlobalParam.f90:122`)
- **v3 default:** `1.0` g-O2/m^2/d
- **Rationale:** Realistic SOD values for moderate organic loading sit in the
  0.5-3.0 g-O2/m^2/d range (Chapra 1997, Table 25.2). v3 chose 1.0 g-O2/m^2/d
  as a defensible midpoint for an unspecified moderate-loading site; this is
  5x Fortran's 0.2 (which sits at the low end of clean-substrate values).
  Users with field measurements should override per simulation. The v1
  sentinel of 999 drives DOX immediately negative on any wet-bed timestep.
  The v3 vs Fortran difference is flagged for LimnoTech reconciliation; both
  values are defensible from literature.

### 1.4 `SOD_theta` — Arrhenius coefficient for SOD
- **Module:** `parameters/dox.py`
- **v1 default:** `999` (sentinel)
- **Fortran default:** `1.06` (`modGlobalParam.f90:122`)
- **v3 default:** `1.060`
- **Rationale:** Arrhenius-style temperature corrections take the form
  `theta^(T-20)`. With `theta=999` and any temperature above 20 C, the
  corrected SOD blows up by orders of magnitude per degree (catastrophic).
  Chapra (1997) recommends `theta=1.060` for SOD as the standard literature
  value, equivalent to a Q10 of approximately 1.79. v3 matches Fortran
  exactly. This was Phase 0's most severe finding.

### 1.5 `kaw_20_user` — user-override wind reaeration coefficient at 20 C
- **Module:** `parameters/dox.py`
- **v1 default:** `999` m/d (sentinel)
- **Fortran default:** `0.0` m/d (`modGlobalParam.f90:117`)
- **v3 default:** `0.0` m/d
- **Rationale:** This parameter is consulted only when `wind_reaeration_option`
  selects the user-override branch (option 1). Setting the default to zero
  means the user-override branch is *off* by default; if a user opts into the
  override path without supplying their own value they will see no reaeration
  rather than a runaway 999 m/d. v3 matches Fortran exactly.

### 1.6 `kah_20_user` — user-override hydraulic reaeration coefficient at 20 C
- **Module:** `parameters/dox.py`
- **v1 default:** `999` 1/d (sentinel)
- **Fortran default:** `1.0` 1/d (`modGlobalParam.f90:113`)
- **v3 default:** `0.0` 1/d
- **Rationale:** v3 chose 0.0 (disabled-by-default) over Fortran's 1.0 to make
  the user-override branch a no-op when not configured, mirroring the
  `kaw_20_user=0` choice. v1's sentinel of 999 was clearly a flaw; the
  subsequent v3 vs Fortran disagreement is a deliberate v3 design choice.

  **Behavioral note (important for v3 vs Fortran side-by-side runs):** at the
  default `hydraulic_reaeration_option=1` (user-defined path), v3's
  `kah_20_user=0.0` produces zero atmospheric hydraulic reaeration. Fortran's
  default of `kah_20_user=1.0` produces 1.0 1/d hydraulic reaeration at the
  same option setting. This means side-by-side runs of v3 vs Fortran NSM1 with
  all-default settings will show DOX recovery in Fortran but not in v3. Users
  who want non-zero default reaeration should either:

  * (a) explicitly set ``kah_20_user > 0`` (e.g., 1.0 to mimic Fortran's
    default), or
  * (b) select a different ``hydraulic_reaeration_option`` from the menu
    (options 2-9 use empirical formulas based on velocity, depth, flow,
    topwidth, slope, or shear velocity, and do not depend on
    ``kah_20_user``).

  This behavioral divergence is mirrored in
  `design/clearwater_modules_v3_nsm1_README.md` Section 7.

### 1.7 `pressure_mb` — atmospheric pressure
- **Module:** `parameters/global_parameters.py`
- **v1 default:** `2026.5` hPa (2x error; possibly intended as Pa-> mb conversion error)
- **Fortran default:** n/a (Fortran NSM1 uses `pressure_atm` in atm, not mb)
- **v3 default:** `1013.25` hPa
- **Rationale:** Standard sea-level atmospheric pressure is `1013.25` hPa
  (ISO 2533). The v1 value is approximately 2x correct; this biases every
  downstream computation that depends on it, notably `O2sat` (DO saturation
  via Henry's law) and `N2sat` (N2 saturation), as well as any atmospheric
  reaeration formula that uses partial pressure. This was identified as a
  Phase 0 finding (2026-05-04) and added to the canonical correction list.

### 1.8 `vson_20` — organic-N settling velocity at 20 C (Phase 9.C fix)
- **Module:** `parameters/nitrogen.py`
- **v1 default:** `0.01` m/d (in `GlobalVars` as `vson`, not in nitrogen group)
- **Fortran default:** `0.01` m/d (`modGlobalParam.f90:92`, `vson` not `vson_20`)
- **v3 default:** `0.01` m/d (corrected in Phase 9.C; was `0.1`)
- **Rationale:** Phase 9.C three-way audit found an internal v3 inconsistency:
  the v3 nitrogen group's `vson_20` was 0.1 m/d while v3 `global_vars.vson`,
  v1 `GlobalVars.vson`, and Fortran `vson` were all 0.01 m/d. The 0.1 in
  the nitrogen group was 10x v1, 10x Fortran, and 10x v3's own `global_vars`
  value, with no documented basis. Phase 9.C corrected v3 to `0.01` m/d to
  match Fortran/v1. Note: `vson` was migrated from `global_vars` to the
  nitrogen group in v3 because it is a nitrogen-specific settling velocity,
  and renamed to `vson_20` for consistency with the other Arrhenius
  rate-base parameters (`knit_20`, `kon_20`, etc.) since v3 added an
  Arrhenius temperature correction (see also `vson_theta` below).

  **v3 addition (`vson_theta=1.024`):** v1 uses `vson` raw (no Arrhenius
  correction); Fortran also has no `vson_theta`. Phase 2.B added the
  Arrhenius correction `vson_tc = arrhenius_correction(T, vson_20,
  vson_theta)` for consistency with the other settling-velocity parameters.
  At T=20 C this collapses exactly to the v1/Fortran behavior. The
  `theta=1.024` value follows the convention used for `kah_theta`,
  `kaw_theta` (other reaeration/settling parameters with mild temperature
  dependence). Documented as a v3 enhancement; flagged for LimnoTech review.

### 1.9 `lambdam` — POM contribution to Beer-Lambert light extinction (Phase 9.C fix)
- **Module:** `parameters/global_vars.py`
- **v1 default:** `0.0174` L/(mg*m) (likely typo)
- **Fortran default:** `0.174` L/(mg*m) (`modGlobalParam.f90:68`)
- **v3 default:** `0.174` L/(mg*m) (corrected in Phase 9.C; was `0.0174`)
- **Rationale:** Phase 9.C three-way audit found a 10x discrepancy: v1
  `GlobalVars` used `0.0174` (10x lower than canonical) and v3 inherited
  the v1 value. Fortran's `0.174` matches QUAL2K Table 6 for the POM
  contribution to Beer-Lambert light extinction and is used throughout the
  legacy v1 NSM test suite (`test_7_nsm_algae_calculations.py:340`,
  `test_10_nsm_carbon_calculations.py:335`, `test_17_nsm_N2_calculations.py:339`,
  etc. all override the v1 default with `lambdam = 0.174` confirming the
  test authors recognized the v1 default as wrong). Phase 9.C corrected v3
  to `0.174` to match Fortran and QUAL2K. The same value is mirrored in
  the inline `_LIGHT_DEFAULTS` dict in `processes/pathogen.py`.

### 1.10 Nitrogen Arrhenius theta transposition (Phase 9.E fix)
- **Module:** `parameters/nitrogen.py`
- **Pre-9.E values:** `kon_theta=1.074`, `rnh4_theta=1.047`, `kdnit_theta=1.08`, `vno3_theta=1.045`
- **Phase 9.E values (matching Fortran):** `kon_theta=1.047`, `rnh4_theta=1.074`, `kdnit_theta=1.045`, `vno3_theta=1.08`
- **Fortran source:** `modNitrogen.f90:82, 89, 95, 100`
- **Rationale:** The four nitrogen Arrhenius theta values were transposed
  in pairs during v1's port from Fortran. v3 inherited the transposition.
  Three independent lines of evidence support the correction:
  (1) Direct Fortran source confirms the canonical values in
  `modNitrogen.f90` initializers. (2) Phosphorus parallel-process check:
  Fortran/v1/v3 all agree on `kop_theta=1.047` (organic-P hydrolysis,
  parallel to `kon`) and `rpo4_theta=1.074` (sediment-P release, parallel
  to `rnh4`); the nitrogen pair should mirror this and does in Fortran but
  did not in v1/v3 pre-9.E. (3) Literature convention (Chapra 1997, QUAL2K
  manual, EPA Bowie et al. 1985): organic-matter hydrolysis universally
  uses `theta=1.047` (matches all other v3 NSM1 organic-matter Arrhenius
  defaults: `mu_max_theta`, `kdp_theta`, `krp_theta`, `kpoc_theta`,
  `kdoc_theta`, `kop_theta`, `kpom_theta`, `kbod_theta`); water-column
  denitrification uses ~1.045; sediment-water exchange velocities use
  steeper temperature dependence (~1.074-1.08).
  Regression coverage in
  `tests/test_5_nitrogen_calculations_v2.py::test_phase9e_*` (5 tests
  pinning each theta value plus the nitrogen-phosphorus parallel-process
  consistency).

### 1.11 DIC unit reconciliation (Phase 9.E fix)
- **Module:** `processes/carbon.py`
- **Issue:** Fortran `modCarbon.f90:268` integrates `dDICdt` in
  mol-C/L/d (every explicit-formula term divides mass by 12000). However
  Fortran `modMain.f90:301` labels DIC as mg-C/L. v1 inherited Fortran's
  formula (with the `/12000` divisions) but stores DIC as mg-C/L
  throughout — the rate is then implicitly added to a state in mg-C/L,
  producing a 12000x scaling error that effectively freezes DIC dynamics.
  v3 inherited v1's mixed convention.
- **Phase 9.E correction:** v3 now integrates `dDICdt` in mg-C/L/d
  throughout. Removed the legacy `/12000.0` divisions from every
  explicit-formula DIC source/sink:
  `dic_algal_resp`, `dic_algal_photo`, `dic_balgae_resp`, `dic_balgae_photo`,
  `dic_sed_release`, `dic_cbod_oxidation`. Converted the Henry's-law
  atmospheric equilibrium term from mol-C/L to mg-C/L by multiplying
  `KH * pCO2 / 1e6` by `MG_C_PER_MOL_C = 12000` (= 12 g-C/mol * 1000 mg-C/g).
  After the fix every term in the dDIC/dt sum is in mg-C/L/d, consistent
  with the mg-C/L DIC state.
- **Magnitude impact:** DIC dynamics are 12000x larger (i.e., physically
  meaningful) under all kinetic conditions where the previously-scaled
  terms were active.
- **v3 deviates from Fortran here.** This is a deliberate v3 correction
  over the legacy reference. Fortran is internally inconsistent (the
  rate is mol-C/L/d but the labeled state is mg-C/L); v3 standardizes on
  the mg-C/L convention used elsewhere in NSM1 (NH4, NO3, OrgN, POC, DOC,
  POM, CBOD, DOX all in mg/L). Regression coverage in
  `tests/test_5_carbon_calculations_v2.py` (the Phase 9.B audit-anchored
  tests now assert mg-C/L/d magnitudes; the
  `test_dic_co2_reaeration_matches_v1` test explicitly documents the
  v3-vs-v1 12000x relation and asserts the corrected v3 form).

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

### 2.8 `lambdas` — light extinction parameter (suspended-solids contribution)
- **Module:** `parameters/global_vars.py`
- **Issue (corrected in Phase 9.C audit):** `lambdas=0.052` represents the
  suspended-solids contribution to Beer-Lambert light extinction. **Earlier
  versions of this section incorrectly claimed v1's `lambdas * Solid` term
  was "commented out / defined but not used."** The Phase 9.C three-way audit
  (`design/clearwater_modules_v3_nsm1_audit_utilities_params.md`) verified
  that v1's `shared/processes.py:232` applies `lambdas * Solid` unconditionally
  in the Beer-Lambert sum, matching Fortran `modGlobalParam.f90:LightExtCoefficient`
  (which loops over `nGS` solid groups summing `lambdas(i) * Solid(i)`). v3's
  `utils/light.py:13-53` also applies the term unconditionally, matching v1
  and Fortran for the single-solid-class case. The previous "commented out"
  claim was a documentation defect.
- **Disposition (v3 1.0.0):** Parameter is fully active in the v3 light
  utility; default value 0.052 matches v1, Fortran, and QUAL2K Table 6. The
  legacy `FIXME(phase1-audit):` inline comment in `global_vars.py` is now
  obsolete and refers to the documentation defect that has been corrected
  here. Multi-solid-class generalization (Fortran's `nGS > 1` loop form) is
  out of scope for v3 NSM1 1.0.0 and would require utility extension.

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

### 3.6 Celsius-to-Kelvin offset — v3 uses 273.15 (SI); v2 uses 273.16 (triple point)

- **Module:** `clearwater_modules_v3/utils/conversions.py` (v3-native;
  no longer a re-export from v2).
- **SI canonical form:** `T_K = T_C + 273.15`. 273.15 K is the absolute
  temperature of 0 deg C; 273.16 K is the triple point of water (a
  separate, slightly higher reference). The 0-deg-C-to-Kelvin offset
  is 273.15.
- **Fortran-A (HEC-RAS-WQ) and Fortran-B (WQM1D):** `T_K = T_C + 273.16`
  in all four sites (audit 2026-05-05 finding — both Fortrans pick the
  triple point as their offset).
- **v1 form:** mixed. `clearwater_modules.shared.processes.celsius_to_kelvin`
  uses 273.16 (for TSM); `clearwater_modules.nsm1.processes.celsius_to_kelvin`
  uses 273.15 (for NSM1).
- **v2 form:** `clearwater_modules_v2.utils.conversions.celsius_to_kelvin`
  uses 273.16, with inline comment "for testing consistency with v1".
- **v3 form (audit 2026-05-05 resolution):** v3's
  `utils/conversions.py` defines `celsius_to_kelvin` locally as
  `T_C + KELVIN_OFFSET` where `KELVIN_OFFSET = 273.15` (in
  `utils/constants.py`). Companion `kelvin_to_celsius` uses the same
  offset. v3 NSM1 modules (`carbon.py`, `dox.py`, `n2.py`) already
  used the literal `+273.15` and are unchanged. v3 TSM
  (`processes/temperature.py`) now picks up 273.15 via the
  `conversions.celsius_to_kelvin` import.
- **Rationale:** the user goal is correct units. 273.15 is the right
  Kelvin offset for 0 deg C; 273.16 is wrong by 0.01 K. The bias was
  small in absolute terms (~3.4e-5 relative at 293 K) but was the
  wrong unit choice for SI temperature physics. Fixed in audit 2026-05-05
  open question 5 resolution.
- **Test impact:** v1-parity tests in `tests/v3/test_5_tsm_calculations_v3.py`
  that were pinned against v1's 273.16-based outputs were re-derived
  at 273.15 (option (b) per user direction). v2-direct-parity tests
  in `tests/test_5_*_calculations_v2.py` are unaffected because they
  exercise v2's `clearwater_modules_v2.utils.conversions.celsius_to_kelvin`
  (still 273.16 there).

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

---

## Section 4: Items flagged for LimnoTech reconciliation

These items were identified by the Phase 9.C three-way audit
(`design/clearwater_modules_v3_nsm1_audit_utilities_params.md`) as
discrepancies between v1 and Fortran that v3 inherited from v1. They were
not corrected in Phase 9.C because the literature basis for the canonical
value is not unambiguous from the codebase alone, or because the change
would alter long-standing v1 behavior that downstream applications may
depend on.

### 4.1 Nitrogen Arrhenius theta values — RESOLVED in Phase 9.E

(Originally flagged here in Phase 9.C as "possible v1 swap with Fortran"
pending reconciliation.) Phase 9.E confirmed the transposition with
three independent lines of evidence (direct Fortran source, phosphorus
parallel-process check, literature convention) and applied the
correction. See Section 1.10 above. This entry remains here as a
historical record of the audit-to-fix path.

### 4.2 `BWa` benthic-algae chlorophyll-a stoichiometry (5000 vs 3500)

- **Fortran default:** `5000.0` g-D / mg-Chla (`modBenthicAlgae.f90:87`)
- **v1 / v3 default:** `3500.0` (`parameters/balgae.py`)
- **Comment:** v3 inherits v1's 3500. Fortran's 5000 differs by ~43%; both
  are within the literature range for benthic algal chlorophyll
  stoichiometry but the canonical QUAL2K value should be confirmed.
  Flagged for LimnoTech review.

### 4.3 `vsop` value — RESOLVED in Phase 9.E

(Originally flagged here in Phase 9.C as "Fortran 0.01 vs v3 0.1
pending reconciliation.") Phase 9.E researched the canonical literature
(Chapra 1997, QUAL2K manual section 5.5.16, EPA Bowie et al. 1985, and
the QUAL2Kw application set) and confirmed that no published source
pins a fixed default — `vop` is a calibration parameter in QUAL2K. The
defensible-default question reduces to a physical-consistency argument:
organic-P in NSM1 originates predominantly from dead-algae detritus,
and the algal settling velocity `vsap = 0.15` m/d is the universal
NSM1 default agreed across Fortran/v1/v3. v3's `vsop = 0.1` m/d is
consistent with the algal-detritus origin; Fortran's 0.01 m/d is 15x
slower than the algae from which the OrgP derives and is implausible
as a representative default. v3 keeps 0.1; rationale and regression
test pinned in Section 1.1 above.

### 4.4 `SOD_20` value (Fortran 0.2 vs v3 1.0)

See Section 1.3 above. v3 chose conservative midpoint (1.0 g-O2/m^2/d)
over Fortran's lower-bound (0.2). Both defensible from Chapra (1997)
Table 25.2; flagged for LimnoTech review.

### 4.5 `kah_20_user` value (Fortran 1.0 vs v3 0.0)

See Section 1.6 above. v3 chose 0.0 (disabled-by-default) over Fortran's
1.0 to make the user-override branch a no-op when not configured. This is
a deliberate v3 design choice; the resulting behavioral divergence is
documented in detail and mirrored in
`design/clearwater_modules_v3_nsm1_README.md` Section 7.

### 4.6 `vson_theta=1.024` is a v3 addition

See Section 1.8 above. v3 added an Arrhenius temperature correction to
the OrgN settling velocity (`vson_tc = arrhenius_correction(T, vson_20,
vson_theta)`); v1 and Fortran apply `vson` raw with no temperature
correction. At T=20 C v3 collapses exactly to v1/Fortran behavior; for
other temperatures v3 differs by `theta^(T-20)`. Flagged for LimnoTech
review; `theta=1.024` is consistent with the values used for related
reaeration/settling parameters but is not directly traceable to a
literature reference for OrgN settling specifically.
