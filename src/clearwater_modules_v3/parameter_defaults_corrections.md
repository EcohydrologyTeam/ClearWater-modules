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

## Section 1: Critical default-value corrections (16 items, applied at port)

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

### 1.6 `kah_20_user` and `hydraulic_reaeration_option` — Phase 9.E correction

This section covers two related parameters: the user-override hydraulic
reaeration coefficient (`kah_20_user`) and the menu selector that
determines which formula computes reaeration (`hydraulic_reaeration_option`).
Phase 9.E researched the cross-model convention and applied a coordinated
correction.

| Parameter | v1 default | Fortran default | v3 (pre-9.E) | v3 (Phase 9.E) |
|---|---|---|---|---|
| `hydraulic_reaeration_option` | 1 (user-supplied) | 1 (user-supplied) | 1 (user-supplied) | **5 (Cover 1976 / Internal)** |
| `kah_20_user` | 999 (sentinel) | 1.0 1/d | 0.0 1/d | 0.0 1/d (unchanged) |

**Cross-model convention research (Phase 9.E):**

| Model | Default reaeration approach | Reference |
|---|---|---|
| **QUAL2K** | "**Internal**" option (Covar 1976 depth-piecewise blend of Owens-Gibbs / O'Connor-Dobbins / Churchill); explicit quote from manual p56: *"if no option is specified, the Internal option is the default."* | Chapra & Pelletier 2008 manual |
| **WASP7** | Empirical formula (O'Connor-Dobbins family) computed from stream hydraulics | EPA WASP documentation |
| **CE-QUAL-W2** | Empirical (REAERC switchable; defaults to formula appropriate to water-body type RIV/LAK) | ERDC/USACE manual |
| **Chapra (1997) textbook** | Empirical hydraulic formulas (velocity, depth based) as the default; "user-supplied constant" only as an advanced calibration option | Surface Water-Quality Modeling, ch. 21 |

The peer water-quality-model standard is to **default to an empirical
hydraulic formula** based on stream hydraulics, with "user-supplied
constant" reserved as an opt-in for advanced calibration. NSM1's
historical default (option 1, user-supplied) plus a sentinel/zero
`kah_20_user` produces a silent constant or zero reaeration that does
not depend on stream hydraulics — inconsistent with peer-model
convention and physically misleading.

**Phase 9.E correction:** v3 default `hydraulic_reaeration_option`
changed from 1 to **5** (Cover 1976 / Internal) to match QUAL2K's
documented default. NSM1 option 5 is the same depth-piecewise blend
that QUAL2K calls "Internal":

* depth < 0.61 m → Owens-Gibbs (`(3.93 v^0.5) / d^1.5`)
* depth > 0.61 m → O'Connor-Dobbins (`(5.32 v^0.67) / d^1.85`)
* depth = 0.61 m → Churchill (`5.026 v / d^1.67`)

Under the new default, v3 produces meaningful empirical reaeration
based on stream hydraulics out of the box. `kah_20_user = 0.0` remains
unchanged but is no longer on the default code path; it is consulted
only when a user explicitly opts into option 1.

This is a deliberate v3 correction over Fortran/v1; v1's sentinel 999
was clearly a flaw, but Fortran's default of "option 1 + 1.0 1/d" was
also non-standard relative to peer-model convention. Phase 9.E
harmonizes v3 with QUAL2K/WASP/CE-QUAL-W2/Chapra-textbook conventions.

**Side-by-side comparison with legacy NSM1:** at default settings,
v3 now produces hydraulic reaeration computed from stream velocity and
depth (Cover 1976 / Internal), whereas Fortran NSM1 produces 1.0 1/d
regardless of hydraulics. The two will diverge for any stream where
the hydraulic-formula reaeration differs from 1.0 1/d (i.e., almost
all real streams). v3's behavior is the physically meaningful one.

**Users who want to mimic Fortran NSM1's behavior:** set
`hydraulic_reaeration_option = 1` and `kah_20_user = 1.0` explicitly.
Users with site-specific calibration data should set their measured
value via the same option-1 + `kah_20_user` path.

Regression coverage in `tests/test_5_dox_calculations_v2.py` (the
existing `test_dox_atmospheric_reaeration_matches_v1_with_user_kah`
test wires the user-supplied path explicitly so the formula change is
isolated to the default behavior); a new
`test_phase9e_default_hydraulic_reaeration_option_is_5` pins the
option-5 default and asserts that default-instantiated DOX produces
non-zero reaeration on a representative stream.

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
  and renamed to `vson_20` (`_20` suffix retained for naming consistency
  with the other parameters in the nitrogen group, even though no
  Arrhenius correction is applied — see Section 1.12 below).

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

### 1.12 `vson_theta` removed in Phase 9.E
- **Module:** `parameters/nitrogen.py`
- **Pre-9.E v3 default:** `1.024` (v3-only addition)
- **Phase 9.E:** removed
- **Fortran:** no `vson_theta`; raw `vson` used (`modNitrogen.f90:233` —
  `OrgN_Settling = vson(r) / depth * OrgN`)
- **v1:** no `vson_theta`; raw `vson` used (`processes.py:1333` —
  `return vson / depth * OrgN`)
- **Rationale:** `vson_theta` was added by Phase 1.2's parameter-library
  agent by analogy with the other nitrogen-group theta values (`knit_theta`,
  `kon_theta`, etc.). Phase 2.B's `organic_nitrogen_settling` then
  consumed it via `arrhenius_correction(T, vson_20, vson_theta)` with a
  docstring claiming "for parity with v1" — but the parity claim was
  false: both v1 and Fortran use raw `vson` without Arrhenius correction.

  Phase 9.E researched Fortran's deliberate type distinction in
  `modGlobalParam.f90`:
  - **Reaction rate constants** (`kon`, `knit`, `kdnit`, `kpoc`, `kdoc`,
    `kop`, `kpom`, `kbod`, `mu_max`, `kdp`, `krp`, ...) are declared as
    `TempCorrectionStruct` with `%rc20` and `%theta` fields — these get
    Arrhenius corrections.
  - **Settling velocities** (`vson`, `vsop`, `vsoc`, `vsap`, `vsbp`,
    `vs`, `vb`) are declared as plain `real(R8)` with no theta — these
    do not get Arrhenius corrections.

  The physical convention behind this distinction: reaction rates
  represent biochemical / microbial activity that scales strongly with
  temperature (Q10 ≈ 2-3, theta ≈ 1.04-1.08). Settling velocities depend
  on particle size/density and water viscosity. Water viscosity decreases
  ~30% over 0-30°C — settling velocity scales correspondingly with
  theta ≈ 1.009, NOT the rate-constant-magnitude theta=1.024 that v3 had
  applied (which overstates the temperature dependence by ~3x).

  Phase 9.E removes `vson_theta` from `parameters/nitrogen.py` and
  changes `Nitrogen.organic_nitrogen_settling` to use raw `self.vson_20`
  directly. This restores parity with v1 and Fortran exactly. The
  `temperature` argument is retained on the method signature for API
  stability with other rate methods on the Process, but is no longer
  consumed. Regression coverage in
  `tests/test_5_nitrogen_calculations_v2.py::test_phase9e_orgn_settling_matches_v1_no_arrhenius`.

### 1.13 `BWa` harmonized to WASP7 canonical (Phase 9.E follow-up)
- **Module:** `parameters/balgae.py`
- **v1 default:** `3500.0` ug-Chla per stoichiometric unit
- **Fortran default:** `5000.0` (`modBenthicAlgae.f90:87`)
- **Pre-9.E v3 default:** `3500.0` (matched v1)
- **Phase 9.E v3 default:** `1000.0`
- **Rationale:** `BWa` is the benthic algae chlorophyll-a stoichiometric
  weight (ug-Chla per stoichiometric unit; v1 docstring at
  ``processes.py:797``). The physically meaningful derived ratio is
  `rab = BWa / BWd` in mg-Chla/g-DW.

  | Source | `BWa` | `rab` (mg-Chla/g-DW) |
  |---|---|---|
  | WASP7 Benthic Algae User's Guide Table 1 | n/a* | **10** |
  | Periphyton literature typical | n/a | 1-15 |
  | NSM1 floating (`AWa`/`AWd` = 1000/100) | 1000 | 10 |
  | Fortran NSM1 benthic | 5000 | 50 (5x WASP7) |
  | v1 / pre-9.E v3 NSM1 benthic | 3500 | 35 (3.5x WASP7) |
  | **Phase 9.E v3 NSM1 benthic** | **1000** | **10 (matches WASP7)** |

  *WASP7 documents Chla:C = 0.025 mg-Chla/mg-C and DW:C = 2.5 mg-DW/mg-C
  -> Chla:DW = 0.025 / 2.5 = 0.01 mg-Chla/mg-DW = **10 mg-Chla/g-DW**.

  Phase 9.E originally kept v1's 3500 with documentation noting it was
  3.5x WASP7's canonical value but had v1-application history. After
  explicit reconciliation, Phase 9.E follow-up harmonizes v3 to the
  WASP7 canonical: `BWa = 1000` gives `rab = 10` mg-Chla/g-DW, which
  (a) matches WASP7's documented benthic stoichiometry, (b) matches
  NSM1's own floating-algae Chla:DW (`AWa/AWd = 10`), bringing benthic
  and floating algae onto the same Chla:DW basis as WASP7 does, and
  (c) sits in the middle of the published periphyton range.

  QUAL2K does not expose a directly comparable parameter; its bottom
  algae kinetics use a cell-quota model (variable Chla per cell driven
  by internal nitrogen and phosphorus quotas) rather than a fixed
  Chla:DW stoichiometry. WASP7 is therefore the authoritative peer-EPA
  reference for this default.

  Calibration impact: any v3 simulation that exercises benthic algae
  with the default `BWa` will see the effective `rab` value drop by
  3.5x. The most visible consequence is that the modeled benthic
  chlorophyll output (computed from biomass via `rab`) will be
  proportionally lower for the same biomass, bringing it into line
  with WASP7 and typical periphyton field measurements. Users with
  v1-calibrated configurations who want to preserve the old behavior
  should explicitly set `BWa = 3500` in their YAML; users with
  Fortran-calibrated configurations should set `BWa = 5000`.

  Regression coverage in
  `tests/test_5_benthic_algae_calculations_v2.py::test_phase9e_bwa_harmonized_to_wasp7_canonical`.
  See also Section 4.2 (audit-history record of how this got from
  "flagged for review" to "RESOLVED in 9.E keep-v1" to the harmonization).

### 1.14 `vb` burial velocity unit-conversion bug (Phase 9.F.A fix)
- **Module:** `parameters/global_vars.py` and inline fallback in `processes/pom.py`
- **v1 default:** `0.01` m/d (likely v1 unit-conversion error)
- **Fortran default:** `0.0025` m/yr ≡ `6.85e-6` m/d (`modGlobalParam.f90:138`)
- **Pre-9.F v3 default:** `0.01` m/d (inherited from v1)
- **Phase 9.F.A v3 default:** `6.85e-6` m/d (= 0.0025 m/yr = 0.25 cm/yr)
- **Rationale:** Phase 9.F.5 research found that v3's `vb = 0.01 m/d` was
  **1460x too high**. Three independent canonical references converge on
  `~6.85e-6 m/d`:
  - **Fortran NSM1** (`modGlobalParam.f90:138`): `vb = 0.0025` m/yr,
    converted at runtime via `/365` to ~6.85e-6 m/d.
  - **WASP7/WASP8** Appendix A parameter table: documented default
    `6.85e-6 m/d` verbatim.
  - **Di Toro (2001)** sediment-flux model and the Chesapeake Bay
    calibration: `~0.25 cm/yr` is the canonical sediment burial rate
    for stream/lake systems. Lake/stream sediment-accumulation
    literature spans `0.05–1.3 cm/yr`.

  Provenance of the v1 bug: v1's `processes.py:2293` removed Fortran's
  runtime `/365` conversion factor with the comment `"removed 365 from
  FORTRAN"`, but did **not** rescale the numerical default's
  magnitude from `m/yr` to `m/d`. `0.0025 m/yr` got silently relabeled
  as `0.01 m/d` (and rounded up). v3 inherited the broken default
  verbatim. This is the same class of bug as `lambdam` (Section 1.9):
  a v1 unit/typo error that v3 carried forward.

  Dimensional smell test: at `vb = 0.01 m/d` the `vb / h2 = 0.1 d⁻¹`
  burial timescale is ~10 days, accidentally matching the typical POM
  dissolution rate `kpom_tc`. Burial cannot physically equal
  mineralization in a well-functioning sediment model. At
  `vb = 6.85e-6 m/d`, `vb / h2 = 6.85e-5 d⁻¹` ≈ 40-year e-folding
  burial timescale, consistent with the Di Toro literature.

  No formula change required — only the numerical default. The Phase 0
  `FIXME(phase1-audit) magnitude not validated` flag was correct;
  Phase 9.F.A research validates that the magnitude was wrong by
  ~1500x and identifies the canonical replacement.

  Calibration impact: at the corrected default, POM burial becomes
  effectively negligible (40-year timescale) at typical NSM1 simulation
  durations (days to years). v1-calibrated simulations that relied on
  the previous fast burial may show POM accumulating in the bed
  compartment instead of being silently buried. Users with
  v1-calibrated configurations who want to preserve the old behavior
  should explicitly set `vb = 0.01` in their YAML.

  Regression coverage in
  `tests/test_5_pom_calculations_v2.py::test_phase9fa_vb_value_pinned`
  and `test_phase9fa_vb_dimensional_smell_test`.

### 1.15 `apx` pathogen sunlight-inactivation efficiency (Phase 9.F.B fix)
- **Module:** `parameters/pathogen.py`
- **v1 default:** `1.0` placeholder (no literature basis)
- **Fortran default:** `1.0` placeholder (same; no literature basis)
- **Pre-9.F.B v3 default:** `1.0` (inherited from v1)
- **Phase 9.F.B v3 default:** `0.017` (W/m^2)^-1 d^-1
- **Rationale:** Phase 9.F research (`docs/clearwater_modules_v3_nsm1_research_2_4_pathogen.md`)
  identified Auer & Niehaus (1993, *Wat. Res.* 27(4):693-701) as the
  canonical literature anchor for pathogen sunlight inactivation:
  alpha = 0.00824 cm^2/cal in cgs units, equivalent to **0.017
  (W/m^2)^-1 d^-1** in SI. This is the value Chapra (1997, *Surface
  Water-Quality Modeling*, McGraw-Hill, Ch. 33) cites in the chapter
  that the legacy NSM1 Fortran source (`modPathogen.f90:90`)
  explicitly references; QUAL2K v2.11b8 §5.5.20.1 inherits the same
  formulation. Mancini (1978) reports a ~5x higher composite value
  (~0.085 (W/m^2)^-1 d^-1) from a multi-study synthesis; both Auer
  /Niehaus and Mancini are within plausible literature scatter.

  The v1 docstring claim that `apx` is "dimensionless" was
  dimensionally incorrect: the rate-balance
  ``[1/d] = apx * q_solar * (dimensionless optical factor)`` requires
  ``apx`` to carry units ``(W/m^2)^-1 d^-1`` because v3 (and v1, and
  Fortran) consume `q_solar` in W/m^2. The placeholder
  ``apx = 1.0`` masked this dimensional inconsistency.

  Coordinated change: Phase 9.F.B also reverts the Phase 3.1
  substitution ``I0 = q_solar * Fr_PAR`` in
  ``processes/pathogen.py:_rate_light_decay`` and uses total broadband
  ``q_solar`` directly. Pathogen inactivation is largely UVA/UVB-
  mediated, not PAR-mediated, so the PAR substitution was a
  v3-introduced deviation from the canonical formulation. With the
  PAR substitution removed, the new default ``apx = 0.017`` ties
  directly to the canonical Auer/Niehaus calibration without any
  ``1/Fr_PAR`` pre-multiplication.

  Calibration impact: at the corrected default, the magnitude of the
  light-induced decay term changes substantially from the placeholder
  baseline (the magnitude of the placeholder rate was effectively
  unbounded by literature; users of v1-calibrated configurations
  likely already overrode ``apx`` to reasonable values for their
  site). v1-calibrated simulations that explicitly set ``apx`` will
  continue to work; users who relied on the unphysical ``apx = 1.0``
  default should expect light-decay rates to drop by ~60x at typical
  surface irradiances (500 W/m^2).

  Regression coverage in
  `tests/test_5_pathogen_calculations_v2.py::test_phase9fb_apx_canonical_value_pinned`
  and `test_phase9fb_light_decay_uses_raw_q_solar`.

### 1.16 `vx` pathogen settling velocity (Phase 9.F.B fix)
- **Module:** `parameters/pathogen.py`
- **v1 default:** `1.0` placeholder (no literature basis); also
  carries an incorrect "(m)" units docstring (should be m/d)
- **Fortran default:** `1.0` placeholder (same)
- **Pre-9.F.B v3 default:** `1.0` (inherited from v1)
- **Phase 9.F.B v3 default:** `1.38` m/d
- **Rationale:** Phase 9.F research identified Auer & Niehaus (1993)
  as the canonical literature anchor: ``vx = 1.38 m/d``, sediment-trap
  measurement of particle-associated fecal coliform in Onondaga Lake.
  This is the value cited by Chapra (1997, Ch. 33), QUAL2K v2.11b8,
  and adopted in subsequent modeling studies (Steets & Holden 2003
  range: 1.0-1.6 m/d). Bowie et al. (1985) compilation reports a
  0.5-2.5 m/d typical range across studies, bracketing the canonical
  1.38 m/d value. Garcia-Armisen & Servais (2009) particle-class
  settling rates (1.17 m/d small, 2.40 m/d large) further bracket the
  composite value.

  No formula change required, only the numerical default and a
  docstring fix (the v3 parameter library already labeled units as
  m/d correctly, but the v1 ``processes.py:3196`` docstring's "(m)"
  typo was carried into early v3 docstring drafts and is now
  corrected).

  Calibration impact: at the corrected default, the settling-loss
  term ``vx / depth * PX`` increases by ~38% relative to the
  placeholder baseline (`1.38 / 1.0`). v1-calibrated simulations that
  explicitly set ``vx`` are unaffected; users who relied on the
  ``vx = 1.0`` default should expect a modest increase in pathogen
  loss at all depths.

  Regression coverage in
  `tests/test_5_pathogen_calculations_v2.py::test_phase9fb_vx_canonical_value_pinned`.

---

## Section 2: Lower-priority audit findings (8 items; 7 RESOLVED, 1 deferred)

These eight items were flagged during Phase 0 as suspicious or unclear but
left at v1 values for the v3 1.0.0 port. Each is marked with a
`FIXME(phase1-audit):` inline comment in the relevant parameter module.
Disposition for each is described below.

### 2.1 `rnh4_20=0`, `vno3_20=0`, `rpo4_20=0` — sediment-flux release rates — RESOLVED in Phase 9.F.C

(Originally flagged here in Phase 0 as "all three sediment-release/uptake
rates are zero, which silently disables sediment-flux contributions to NH4,
NO3, and PO4 budgets; Phase 0 could not confirm that every code path
consuming these parameters is gated by `use_SedFlux`.")

- **Modules:** `parameters/nitrogen.py` (`rnh4_20`, `vno3_20`),
  `parameters/phosphorus.py` (`rpo4_20`),
  `clearwater_modules_v2/processes/nitrogen.py`,
  `clearwater_modules_v3/processes/phosphorus.py`.

- **Phase 9.F.1 finding (audit):** the Phase 0 concern was correct.
  Three v3 1.0.0 bed-flux consumers run **without** `use_SedFlux`
  gating:
  1. `Nitrogen.ammonium_from_bed` — consumes `rnh4_20`, called
     unconditionally from `change_ammonium`.
  2. `Nitrogen.nitrate_bed_denitrification` — consumes `vno3_20`,
     called unconditionally from the NO3 budget.
  3. `Phosphorus.run` `dip_from_bed` term — consumes `rpo4_20`,
     gated by `use_TIP` only (NOT by `use_SedFlux`).
  Only `Carbon` correctly gates its sediment-flux contribution
  (`JDIC`) by `if self.use_SedFlux:`.

- **Resolution (Phase 9.F.C, defensive Option A + Option B):**
  1. **Option A (framing fix).** Documented in the Section 2.1 above
     and in the inline comments of the three parameter modules: the
     zero defaults for `rnh4_20`, `vno3_20`, and `rpo4_20` are the
     **de facto gate** for sediment-flux contributions in v3 1.0.0,
     because the `use_SedFlux` boolean is only consulted in `Carbon`.
     The value-based gate (zero rates) is intentional and correct for
     a v1-parity port; the boolean-gated implementation belongs to
     the future NSM2 diagenesis path.
  2. **Option B (defensive guard).** Added `NotImplementedError`
     guards in `Nitrogen.__init__` (v2 overlay) and
     `Phosphorus.__init__` (v3-native) that fire if a user passes
     `parameters={"use_SedFlux": True, ...}`. The guards explicitly
     point users at the NSM2 path for the full sediment-flux feature
     and at direct `rnh4_20` / `vno3_20` / `rpo4_20` overrides for
     site-specific constant-flux calibration without `use_SedFlux`.
     The guards prevent the historical silent-partial behavior where
     `use_SedFlux=True` would activate Carbon's sediment-flux but
     leave Nitrogen and Phosphorus at value-based zero, producing an
     inconsistent budget. `Carbon` is unaffected because it correctly
     gates by `use_SedFlux`.

  Regression coverage in
  `tests/test_phase9fc_sedflux_guard.py::test_nitrogen_use_sedflux_true_raises_notimplementederror`
  and `test_phosphorus_use_sedflux_true_raises_notimplementederror`.

### 2.2 `kdpo4=0.0` — TIP partitioning coefficient — DEFERRED to NSM2 path

- **Module:** `parameters/phosphorus.py`
- **Issue:** Phosphorus partitioning between dissolved (DIP) and particulate
  (sorbed) phases is governed by `kdpo4` (L/kg). At zero, no DIP adsorbs
  onto suspended solids; the TIP partitioning feature is effectively
  disabled. The v3 partitioning helper `fdp_partition` in
  `utils/phosphorus.py` returns `fdp = 1.0` whenever `kdpo4 * Solid = 0`,
  collapsing the TIP particulate fraction to zero (settling and
  sediment-flux contributions both vanish).
- **Disposition (v3 1.0.0):** Kept at zero by design. Full DIP-solid
  partitioning is NSM2 territory: a proper implementation requires (1) a
  multi-class suspended-solids model for which `Solid` (single-class) is a
  placeholder; (2) coupling to the NSM2 sediment-diagenesis sediment-flux
  model (paired with the Section 2.1 `use_SedFlux` work — both are gated
  by the same NSM2 path). v3 NSM1 1.0.0 maintains the v1 simplification of
  treating TIP as fully dissolved at the default value, matching v1 and
  Fortran behavior exactly. The `FIXME(phase1-audit):` tag in
  `parameters/phosphorus.py:16` and the inline `FIXME(phase1-audit)`
  references in `processes/phosphorus.py:35,38` are retained as cross-
  references to this NSM2-scope item; they intentionally do **not** mark
  this as a defect to be fixed in 1.0.0.

### 2.3 `ksbod_20=0.0` — CBOD settling rate — RESOLVED in Phase 9.F.C

(Originally flagged here in Phase 0 as "CBOD never settles in v1 because
the settling rate is hardcoded to zero. If CBOD is conceptually fully
dissolved this is correct; if CBOD represents a mixture including
particulate fractions, this is a bug.")

- **Module:** `parameters/cbod.py`

- **Phase 9.F research finding** (Phase 9.F.2 research record at
  `docs/clearwater_modules_v3_nsm1_research_2_3_ksbod.md`): the zero
  default is the **intentional, defensible modern-convention value**.
  - **QUAL2K v2.11b8** (Chapra, Pelletier & Tao 2008), QUAL2Kw, WASP7
    EUTRO, and CE-QUAL-W2 all treat CBOD as a **dissolved-only**
    state variable and provide **no CBOD settling parameter at all**;
    particulate organic matter is carried separately (detritus,
    LPOM/RPOM) with its own settling velocity.
  - **QUAL2E** (Brown & Barnwell 1987, EPA/600/3-87/007), the direct
    ancestor of NSM1, defaults `K_3 = 0` for the same reason.
  - **EPA TMDL Technical Guidance Book II** (Sample Calc B-3)
    explicitly assumes `K_s = 0` for treated effluent.
  - **Yamuna River QUAL2E case** (Parmar & Keshari, citing Kazmi &
    Agrawal 2005) calibrated `K_3 = 0.9 /d` uniformly across 16
    reaches in a heavily polluted urban stretch where particulate-
    laden CBOD settling dominated removal — illustrating that
    nonzero values are site-, source-, and treatment-specific
    calibration parameters.

- **Resolution (Phase 9.F.C):** FIXME cleared from
  `parameters/cbod.py`. The zero default is documented inline with
  the QUAL2K / Brown & Barnwell / EPA-TMDL citations. No value
  change.

- **Two related defects flagged for future audit (NOT addressed in
  Phase 9.F.C):**
  1. **Units form mismatch.** v3 `processes/cbod.py:240` divides
     `ksbod_tc` by depth (`ksbod_tc / depth * cbod`), implementing it
     as a settling velocity (m/d). Fortran NSM1
     (`modCBOD.f90:114`) and QUAL2E both treat the parameter as a
     1/d rate constant (no depth division). With `ksbod_20 = 0` the
     form difference is silent; nonzero user values would diverge
     by a factor of `1/depth`.
  2. **Theta mismatch.** v3 `ksbod_theta = 1.047` differs from
     Fortran/QUAL2E `1.024` (the canonical settling-coefficient
     Arrhenius value per Bowie 1985 / QUAL2E).

  Both are flagged for follow-up when CBOD settling becomes
  actively used (e.g., by a user supplying nonzero `ksbod_20`).

### 2.4 `apx=1`, `vx=1` — pathogen placeholder values — RESOLVED in Phase 9.F.B

(Originally flagged here in Phase 0 as "no documented literature basis;
likely placeholder values inherited from v1.")

- **Module:** `parameters/pathogen.py`
- **Resolution (Phase 9.F.B):** Replaced both placeholders with the
  canonical Auer & Niehaus (1993) / Chapra (1997) values cited by
  QUAL2K and the Fortran source. ``apx = 0.017 (W/m^2)^-1 d^-1`` and
  ``vx = 1.38 m/d``. Coordinated with reverting the Phase 3.1
  substitution ``I0 = q_solar * Fr_PAR`` in ``_rate_light_decay`` so
  that ``apx`` ties directly to the canonical calibration on total
  broadband solar radiation. Full rationale in Section 1.15 (apx) and
  Section 1.16 (vx); Phase 9.F research record in
  ``docs/clearwater_modules_v3_nsm1_research_2_4_pathogen.md``.

### 2.5 `h2=0.1` — active sediment layer thickness — RESOLVED in Phase 9.F.C

(Originally flagged here in Phase 0 as "`h2` is used as a divisor in
burial/sedimentation terms (e.g., `vb * POM / h2`); its physical role
is unclear, and the units (meters) do not obviously line up with the
surrounding terms.")

- **Module:** `parameters/pom.py`

- **Phase 9.F.4 research finding** (research record at
  `docs/clearwater_modules_v3_nsm1_research_2_5_pom_h2.md`): the Phase
  0 audit's "unclear physical role" was a **documentation defect**,
  not a substantive issue.
  - The v1 `static_variables.py:921` declaration and Fortran
    `modGlobalParam.f90:38` both define `h2` unambiguously as the
    **active sediment layer thickness (m)**.
  - `h2 = 0.1 m` matches the Di Toro (2001) / QUAL2K v2.11 §5.6
    convention for the lower anaerobic sediment layer thickness
    `H_2` (approx 10 cm; QUAL2K Eq. 214 uses `H_2` as the volumetric
    divisor in the bed-POM mass balance).
  - NSM1's POM state variable represents bed-sediment POM (Fortran
    `POM2` -- the "2" suffix denoting Di Toro layer 2). v1 and v3
    dropped the `2` subscript when porting from Fortran but the
    conceptual identity is preserved via the `h2` divisor: `h2`
    converts areal water-column fluxes (m * mg/L/d) into bed
    volumetric concentration changes (mg/L/d) in the bed layer.
  - All four `h2` consumers in v3 `POM.run` (`rate_burial`,
    `rate_poc_settling`, `rate_algal_settling`,
    `rate_benthic_mortality`) are dimensionally consistent with this
    interpretation.
  - Implementing the full two-layer Di Toro flux model (separate
    `H_1` aerobic layer plus full nutrient-flux equations) is the
    future NSM2 sediment-diagenesis scope; v3 NSM1 1.0.0 carries
    only the `H_2` layer with first-order burial/dissolution kinetics
    matching v1 and Fortran exactly.

- **Resolution (Phase 9.F.C):** FIXME cleared from
  `parameters/pom.py`. The parameter docstring inline-comment now
  records the Di Toro / QUAL2K provenance. The `processes/pom.py`
  module docstring gained a one-paragraph "Conceptual note"
  explaining that POM represents the bed-sediment compartment
  (Fortran `POM2`) and `h2` is the bed-layer thickness (Di Toro
  `H_2`). No value or formula change.

  Regression coverage in
  `tests/test_phase9fc_documentation.py::test_h2_fixme_cleared`.

### 2.6 `vb` burial velocity — RESOLVED in Phase 9.F.A

(Originally flagged here in Phase 0 as "magnitude not validated; some
sediment-flux literature uses values one to two orders of magnitude
different.") Phase 9.F.5 research found the magnitude was wrong by
**1460x** — a v1 unit-conversion error inherited verbatim. Fortran's
0.0025 m/yr (= 6.85e-6 m/d), WASP7/WASP8's documented 6.85e-6 m/d, and
Di Toro 2001's Chesapeake Bay calibration all agree on the canonical
~0.25 cm/yr value. Phase 9.F.A applied the correction; full rationale
in Section 1.14 above.

### 2.7 `q_solar=500` units mismatch — RESOLVED in Phase 9.F
- **Module:** `parameters/global_vars.py`
- **Issue:** v1's docstring documents `q_solar` as having units of `1/d`, but
  the value 500 and its consumption pattern (Beer-Lambert PAR computation)
  are consistent with W/m^2. The docstring is wrong, not the value.
- **Disposition (v3 1.0.0):** Value kept at 500. The inline comment in
  `parameters/global_vars.py` now states the correct units (W/m^2) with an
  explicit note that v1's docstring was incorrect, and the
  `processes/pathogen.py:_rate_light_decay` docstring (where `q_solar` is
  consumed via `utils.light.PAR`) carries a matching units note. The
  `utils/light.py:PAR` Args block already documents `q_solar | W/m^2 |
  total incident solar radiation at the water surface.` These three sites
  are now consistent on W/m^2 — the only sites in v3 NSM1 that reference
  `q_solar` semantics. The v1 docstring error is a v1-side issue and
  does not propagate into v3. The legacy `FIXME(phase1-audit):` tag was
  removed in Phase 9.F; the parameter is no longer flagged for follow-up.
  Regression coverage in
  `tests/test_5_pathogen_calculations_v2.py::test_phase9f_q_solar_units_are_w_per_m2`.

### 2.8 `lambdas` — light extinction parameter (suspended-solids contribution) — RESOLVED in Phase 9.C, FIXME cleared in Phase 9.F
- **Module:** `parameters/global_vars.py`
- **Issue (corrected in Phase 9.C audit):** `lambdas=0.052` represents the
  suspended-solids contribution to Beer-Lambert light extinction. The
  Phase 9.C three-way audit
  (`design/clearwater_modules_v3_nsm1_audit_utilities_params.md`) verified
  that v1's `shared/processes.py:232` applies `lambdas * Solid` unconditionally
  in the Beer-Lambert sum, matching Fortran `modGlobalParam.f90:LightExtCoefficient`
  (which loops over `nGS` solid groups summing `lambdas(i) * Solid(i)`). v3's
  `utils/light.py:13-53` also applies the term unconditionally, matching v1
  and Fortran for the single-solid-class case. The earlier Phase 0 framing
  that v1's `lambdas * Solid` term was "commented out / defined but not
  used" was a documentation defect; the term is and always has been active
  in v1 `shared/processes.py:232`. That false claim has been removed from
  this section.
- **Disposition (v3 1.0.0):** Parameter is fully active in the v3 light
  utility; default value 0.052 matches v1, Fortran, and QUAL2K Table 6. The
  legacy `FIXME(phase1-audit):` inline comment that flagged `lambdas` as
  suspect in `parameters/global_vars.py` was removed in Phase 9.F; the
  inline comment now records that the parameter is active per the Phase
  9.C verification. Multi-solid-class generalization (Fortran's `nGS > 1`
  loop form) is out of scope for v3 NSM1 1.0.0 and would require utility
  extension.

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

### 3.4 Pathogen light decay — `PAR = q_solar * Fr_PAR` substitution — RESOLVED (reverted) in Phase 9.F.B

(Originally documented here as an intentional Phase 3.1 v3 deviation:
v3 `_rate_light_decay` substituted `I0 = q_solar * Fr_PAR` for v1's
raw `q_solar`, on the rationale that PAR is a closer proxy than total
shortwave for UV-mediated pathogen decay. The calibration target
`apx` was supposed to absorb the difference.)

- **Module:** `processes/pathogen.py` (light-driven decay sub-rate)
- **Resolution (Phase 9.F.B):** the substitution was **reverted**. Phase 9.F
  research (`docs/clearwater_modules_v3_nsm1_research_2_4_pathogen.md`) found
  that pathogen inactivation is largely UVA/UVB-mediated, not PAR-mediated,
  and the canonical Auer & Niehaus (1993) / Chapra (1997) / QUAL2K
  formulation operates on **total broadband solar radiation**. Substituting
  PAR was therefore a v3-introduced deviation away from the canonical
  literature, not a port improvement. With Phase 9.F.B's coordinated change
  to `apx = 0.017 (W/m²)⁻¹·d⁻¹` (Auer & Niehaus canonical, calibrated to
  total broadband), the natural choice is to use raw `q_solar` directly in
  the kinetics. v3 now matches v1 `PathogenDecay` exactly at the kinetics
  level. See Section 1.15 for the coordinated `apx` correction.
- **Reference test:** `tests/test_5_pathogen_calculations_v2.py::test_phase9fb_light_decay_uses_raw_q_solar`
  pins the post-9.F.B behavior (rate is insensitive to `Fr_PAR`).

### 3.5 CBOD sedimentation — `ksbod_tc / depth` in v3 (cross-references Section 2.3)

- **Module:** `processes/cbod.py` (sedimentation sink)
- **v1 / Fortran / QUAL2E form:** `CBOD_sedimentation = CBOD * ksbod_tc`,
  treating `ksbod_tc` directly as a 1/d first-order rate (no depth
  division). Confirmed in Phase 9.F.2 research: legacy Fortran
  `modCBOD.f90:114` and QUAL2E (Brown & Barnwell 1987) both define
  `ksbod_20` (their `K_3`) as a 1/d rate constant, with units stated as
  "(1/day) Range {-0.36-0.36}" in the Fortran source comments.
- **v3 form:** `CBOD * ksbod_tc / depth`, treating `ksbod_tc` as a
  settling velocity (m/d) divided by water-column depth to yield a 1/d
  first-order rate. This conflicts with the Fortran/QUAL2E convention
  but matches the velocity-style units label inherited from v1's
  `constants.py`.
- **Status:** the form mismatch is **silent** at the canonical
  `ksbod_20 = 0` default (both forms produce zero) but would diverge
  by a factor of `1/depth` at any nonzero user value. Phase 9.F.C
  (Section 2.3) documented this and the related `ksbod_theta` mismatch
  (v3 `1.047` vs Fortran/QUAL2E `1.024`) as deferred follow-ups; both
  become actionable only if a user activates CBOD settling. See Section
  2.3 above for the Phase 9.F.2 research record and the
  modern-convention rationale (QUAL2K, WASP, CE-QUAL-W2 omit the
  parameter entirely).
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

### 4.2 `BWa` benthic-algae chlorophyll-a stoichiometry — RESOLVED in Phase 9.E

(Originally flagged here in Phase 9.C as "Fortran 5000 vs v1/v3 3500
pending reconciliation.") Phase 9.E researched the canonical literature
and made the keep-v1 choice.

**Parameter:** `BWa` (ug-Chla per stoichiometric unit; v1 docstring at
`processes.py:797` confirms `BWa = Benthic algae chlorophyll-a
(ug-Chla-a)`). Used in v1 helper `rab(BWa, BWd) = BWa / BWd` returning
the benthic algae Chla:DW ratio in ug-Chla/mg-D = mg-Chla/g-DW.

| Source | `BWa` | Chla:DW (mg-Chla/g-DW) | Comment |
|---|---|---|---|
| WASP7 (EPA peer model) | n/a | 10 (= 0.025 mg-Chla/mg-C / 2.5 mg-DW/mg-C) | Documented Table 1 of WASP7 Benthic Algae User's Guide |
| Periphyton literature typical | n/a | 1-15 | Bothwell 1989, Stevenson et al. 1996 |
| NSM1 floating (`AWa=1000`, `AWd=100`) | 1000 | 10 | Matches WASP7 / typical literature |
| Fortran NSM1 benthic | 5000 | 50 | 5x WASP7; high |
| v1 / pre-9.E v3 NSM1 benthic | 3500 | 35 | 3.5x WASP7; high |
| **Phase 9.E v3 NSM1 benthic** | **1000** | **10** | **Matches WASP7 canonical** |

Audit-history record (this section was updated twice during Phase 9.E):

1. **Phase 9.C audit** flagged the BWa value as a v1/v3-vs-Fortran
   disagreement (3500 vs 5000) pending LimnoTech reconciliation.
2. **Phase 9.E first pass** kept v1's 3500 with documentation noting
   it was 3.5x WASP7 canonical but had v1-application history. Both
   v1's 3500 and Fortran's 5000 were noted as HIGH relative to
   WASP7's documented benthic Chla:DW.
3. **Phase 9.E follow-up** (this commit) harmonized v3 to the WASP7
   canonical: `BWa = 1000` gives `rab = 10` mg-Chla/g-DW, matching
   WASP7 explicitly and bringing v3 benthic and floating algae onto
   the same Chla:DW basis as WASP7's convention. Full rationale,
   calibration-impact analysis, and migration guidance for users with
   v1- or Fortran-calibrated configurations are in Section 1.13 above.

This section is RESOLVED. Regression coverage in
`tests/test_5_benthic_algae_calculations_v2.py::test_phase9e_bwa_harmonized_to_wasp7_canonical`.

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

### 4.4 `SOD_20` value — RESOLVED in Phase 9.E

(Originally flagged here in Phase 9.C as "Fortran 0.2 vs v3 1.0
pending reconciliation.") Phase 9.E confirms v3's `SOD_20 = 1.0`
g-O2/m^2/d as the chosen default. Rationale:

- Chapra (1997) Table 25.2 cites SOD literature range of 0.2-3.0
  g-O2/m^2/d for typical surface waters, with values strongly
  site-dependent (low for oligotrophic, high for eutrophic).
- v3's 1.0 g-O2/m^2/d is a conservative midpoint that produces a
  visible non-zero default and surfaces obvious calibration problems
  if the user does not override.
- Fortran's 0.2 g-O2/m^2/d is the lower-bound of the literature
  range and silently understates SOD in moderately-loaded systems
  (the typical NSM1 application target).
- Users with site-specific data should override per project. The
  conservative default keeps DOX dynamics responsive to organic
  loading even at default settings.

Regression coverage in
`tests/test_5_dox_calculations_v2.py::test_phase9e_sod_20_value_pinned`.
See also Section 1.3 for the sentinel-999 correction history.

### 4.5 `kah_20_user` and hydraulic reaeration default — RESOLVED in Phase 9.E

(Originally flagged here in Phase 9.C as "Fortran 1.0 vs v3 0.0
behavioral divergence pending reconciliation.") Phase 9.E researched
the cross-model convention (QUAL2K, WASP7, CE-QUAL-W2, Chapra 1997)
and found that the peer water-quality-model standard is to default to
an **empirical hydraulic formula** (typically Covar 1976 / Internal,
which is QUAL2K's documented default) rather than to a user-supplied
constant. The original "v3 0.0 vs Fortran 1.0" disagreement was not
the right thing to reconcile — both v3 and Fortran were defaulting to
`hydraulic_reaeration_option = 1` (user-supplied path), which is
itself non-standard.

Phase 9.E corrected v3's default `hydraulic_reaeration_option` from 1
to **5** (Cover 1976 / Internal), matching QUAL2K's documented
default. Under the new default, `kah_20_user = 0.0` is no longer on
the default code path (only consulted when the user explicitly opts
into option 1). The behavioral divergence with legacy NSM1 (Fortran
or v3-pre-9.E) is now *intended*: v3 produces empirically-derived
reaeration from stream hydraulics; legacy NSM1 produced a constant
1.0 1/d (Fortran) or 0 (v3-pre-9.E) regardless of hydraulics. Full
rationale and migration guidance in Section 1.6 above.

### 4.6 `vson_theta` — RESOLVED in Phase 9.E (parameter removed)

(Originally flagged here in Phase 9.C as "v3-only addition pending
LimnoTech reconciliation.") Phase 9.E researched the issue and found
that `vson_theta` was an unjustified v3 addition: Phase 1.2's
parameter-library agent added it by analogy with rate-constant theta
values, and Phase 2.B's `organic_nitrogen_settling` consumed it with
a docstring claiming "parity with v1" — but both v1 and Fortran use
raw `vson` without Arrhenius correction. The parameter has been
removed and the Process now uses raw `vson_20` directly, matching v1
and Fortran exactly. Full rationale in Section 1.12 above.
