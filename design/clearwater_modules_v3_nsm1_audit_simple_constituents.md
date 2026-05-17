# v3 NSM1 simple constituents (POM, CBOD, Pathogen, N2, Alk) — Three-way audit

Date: 2026-05-05
Scope: POM, CBOD, Pathogen, N2, Alkalinity Process classes in
`src/clearwater_modules_v3/processes/` and matching parameter defaults in
`src/clearwater_modules_v3/parameters/`.
References: legacy Fortran NSM1 (`/Users/todd/Downloads/NSM_comparison/NSM1/Source Files/`),
v1 Python NSM1 (`src/clearwater_modules/nsm1/processes.py`,
`src/clearwater_modules/nsm1/constants.py`).

## Summary

- 23 blocks audited across 5 constituents.
- Findings: 0 critical, 7 minor (all documented in corrections doc / parity test
  docstrings), 16 matches.
- Top concerns: none of the 7 minor deviations are correctness defects under
  matched inputs; each is either a documented architectural choice
  (Section 3 of `parameter_defaults_corrections.md`), an absorbed unit
  scaling under recalibration (Pathogen `Fr_PAR`, CBOD `1/depth`), or a
  Fortran-only convention (annual `vb`, oxygen-weighted TDG) that v1 and v3
  already do not follow.

## POM

### 1. Algal settling input (floating algae -> POM)

- Fortran (`modPOM.f90:98`): `ApSettling_POM2 = vsap(r) * Ap * rda(r) / h2(r)`
  (sedimentation to bed sediment layer).
- v1 (`processes.py:2200-2218`):
  `POM_algal_settling = vsap * Ap * rda / h2` (gated by `use_Algae`).
- v3 (`pom.py:280-287`): reads cached
  `floating_algae_process.algal_pom_from_settling_rate` (mg/L/d) populated
  by FloatingAlgae (Phase 2.A), gated by `use_Algae` and presence of process.

Match. v3 routes the same `vsap*Ap*rda/h2` term through the FloatingAlgae rate
cache.

### 2. Benthic algae mortality input -> POM

- Fortran (`modPOM.f90:103`):
  `AbDeath_POM2 = AbDeath * Fb(r) * (1 - Fw(r)) / h2(r)`
  (where `AbDeath = kdb_tc * Ab` is computed in `modBenthicAlgae`).
- v1 (`processes.py:2257-2277`):
  `POM_benthic_algae_mortality = Ab * kdb_tc * Fb * (1 - Fw) / h2`.
- v3 (`pom.py:293-300`): reads cached
  `benthic_algae_process.balgae_pom_from_mortality_rate` (mg/L/d).

Match. Same product, routed through the BenthicAlgae rate cache.

### 3. POM dissolution to DOC

- Fortran (`modPOM.f90:113`): `POM2_Dissolution = kpom2_tc * POM2`,
  with `kpom2_tc = Arrhenius_TempCorrection(kpom2(r), TsedC)`
  (sediment temperature, line 87).
- v1 (`processes.py:2222-2233`): `POM_dissolution = POM * kpom_tc`
  with `kpom_tc = arrhenius_correction(TwaterC, kpom_20, kpom_theta)`
  (water temperature; v1 line 2185-2197).
- v3 (`pom.py:253-262`): `kpom_tc = arrhenius_correction(water_temperature,
  kpom_20, kpom_theta)`; `rate_dissolution = kpom_tc * pom`. Cached as
  `pom_hydrolysis_rate` for Carbon.

Minor deviation (Fortran-only). Fortran applies Arrhenius to sediment
temperature (`TsedC`); v1 and v3 both use `TwaterC`. v3 matches v1.

### 4. POM burial

- Fortran (`modPOM.f90:114`): `POM2_Burial = vb(r) / 365.0 * POM2 / h2(r)`.
  Note Fortran's `vb` is stored as m/year and divided by 365 inline; the
  setter at `modGlobalParam.f90:201` also performs the conversion when
  reading user-supplied values.
- v1 (`processes.py:2281-2293`): `POM_burial = vb * POM / h2`
  (inline comment: "note removed 365 from FORTRAN").
- v3 (`pom.py:265`): `rate_burial = self.vb * pom / self.h2`.

Minor deviation. Fortran uses `vb` in m/year; v1 and v3 use `vb` in m/d
directly. Default magnitudes: Fortran 0.0025 m/yr, v1/v3 `vb=0.01` m/d.
The unit convention is consistent v1<->v3; only the Fortran legacy form
differs. Documented inline in v1.

### 5. POC settling -> POM

- Fortran (`modPOM.f90:108`):
  `POCSettling_POM2 = vsoc(r) * POC / focm(r) / h2(r)`.
- v1 (`processes.py:2236-2254`):
  `POM_POC_settling = vsoc * POC / h2 / fcom` (gated by `use_POC`).
- v3 (`pom.py:269-272`): `rate_poc_settling = self.vsoc * poc / self.h2 /
  self.fcom`.

Match. Same algebraic form; v3 reads `poc` defensively (zeros if not in
registry) and uses `_POM_GLOBAL_DEFAULTS` fallback for `vsoc`, `fcom`.

## CBOD

### 6. CBOD oxidation with DOX-Monod attenuation

- Fortran (`modCBOD.f90:108-113`):
  `CBOD_Oxidation = DOX / (DOX + KsOxbod(i,r)) * kbod_tc(i) * CBOD(i)`
  when `use_DOX`, else `kbod_tc * CBOD`. NaN guard inline at line 110.
- v1 (`processes.py:2370-2388`):
  `(DOX / (KsOxbod + DOX)) * kbod_tc * CBOD` when `use_DOX`,
  `kbod_tc * CBOD` otherwise.
- v3 (`cbod.py:228-233`):
  `kbod_tc * dox / (KsOxbod + dox) * cbod` when `use_DOX`,
  `kbod_tc * cbod` otherwise. Cached as `cbod_oxidation_rate` for DOX.

Match. v3 follows the v1/Fortran form exactly.

### 7. CBOD sedimentation

- Fortran (`modCBOD.f90:114`):
  `CBOD_Sediment = ksbod_tc(i) * CBOD(i)` (no depth divide; treats
  `ksbod_tc` directly as 1/d).
- v1 (`processes.py:2392-2404`):
  `CBOD_sedimentation = CBOD * ksbod_tc` (same as Fortran).
- v3 (`cbod.py:240`):
  `settling_rate = ksbod_tc / depth * cbod` (treats `ksbod_tc` as m/d
  settling velocity divided by depth).

Minor deviation. Documented in `parameter_defaults_corrections.md` Section
3.5 and pinned in `tests/test_5_cbod_calculations_v2.py:173-200`. Under the
v1/v3 default `ksbod_20=0`, both forms are identically zero, so no runtime
divergence under defaults. Under user-supplied `ksbod_20 > 0`, the v3
result is `1/depth` times the v1 result; recalibration of `ksbod_20` per
the dimensional-consistency interpretation is required if porting from v1.

### 8. CBOD multi-group support

- Fortran: native multi-group via `do i = 1, nCBOD` loops.
- v1: same multi-group convention.
- v3 (`cbod.py:88-97, 191-196`): single-group (key `cbod`); multi-group
  documented as a future Phase 4+ extension path in the class docstring.

Acknowledged scope deferral, not a deviation. v3 1.0.0 only supports the
single-group fixture path consumed by the Tier 1 conftest.

## Pathogen

### 9. Natural decay

- Fortran (`modPathogen.f90:79, 87`):
  `kdx_tc = Arrhenius_TempCorrection(kdx(r), TwaterC)`;
  `PathogenDeath = kdx_tc * PX`.
- v1 (`processes.py:3141-3170`):
  `kdx_tc = arrhenius_correction(TwaterC, kdx_20, kdx_theta)`;
  `PathogenDeath = kdx_tc * PX`.
- v3 (`pathogen.py:265-275`): same formula, identical Arrhenius temperature
  correction.

Match.

### 10. Light-induced decay

- Fortran (`modPathogen.f90:91`):
  `PathogenDecay = apx(r) * q_solar / (lambda * depth) *
  (1 - exp(-lambda * depth)) * PX`. Source comment "q_solar units is ly/day
  in original formulation (Chapra, 1997)???" indicates uncertainty.
- v1 (`processes.py:3172-3190`):
  `apx * q_solar / (L * depth) * (1 - np.exp(-L * depth)) * PX`
  (raw `q_solar`, no `Fr_PAR` scaling).
- v3 (`pathogen.py:303-329`): replaces `q_solar` with
  `i0 = PAR(q_solar, Fr_PAR) = q_solar * Fr_PAR`; otherwise identical
  Beer-Lambert depth-averaged form. Adds an `xr.where(kd > 0)` guard
  against (KEXT*depth -> 0) NaN.

Minor deviation. Documented in `parameter_defaults_corrections.md` Section
3.4 and pinned in `tests/test_5_pathogen_calculations_v2.py::test_pathogen_light_decay_matches_v1`
with `Fr_PAR=1.0`. The constant 0.47 is absorbable into a recalibrated
`apx`, so this is a calibration-target adjustment, not a correctness
defect. The NaN-guard at small `kd` is a v3 robustness improvement.

### 11. Settling

- Fortran (`modPathogen.f90:95`): `PathogenSettling = vx(r) / depth * PX`.
- v1 (`processes.py:3193-3206`): `vx / depth * PX`.
- v3 (`pathogen.py:331-338`): `self.vx / depth * px`.

Match.

### 12. Pathogen overall sign convention

- Fortran (`modPathogen.f90:98`):
  `dPXdt = -PathogenDeath - PathogenDecay - PathogenSettling`.
- v1 (`processes.py:3209-3224`): same.
- v3 (`pathogen.py:259-263`): `return -(natural + light + settling)`.

Match.

## N2

### 13. Henry's law constant KHN2_tc

- Fortran (`modN2.f90:40`):
  `KHN2_tc = 0.00065 * exp(1300.0 * (1.0 / TwaterK - 1 / 298.15))`,
  with `TwaterK = TwaterC + 273.15` (line 37).
- v1 (`processes.py:3452-3467`):
  `KHN2_tc = 0.00065 * np.exp(1300.0 * (1.0 / TwaterK - 1 / 298.15))`,
  with `TwaterK = celsius_to_kelvin(TwaterC) = TwaterC + 273.15`
  (`processes.py:9-10`). Note: v1's `celsius_to_kelvin` uses 273.15;
  the +273.16 form is in the v2 utilities module (used by v2 parity
  shims), not the v1 reference.
- v3 (`n2.py:79-87`):
  `0.00065 * np.exp(1300.0 * (1.0 / t_water_k - 1.0 / 298.15))`
  with `_kelvin(t_c) = t_c + 273.15` (line 74-76).

Match (Fortran, v1, and v3 all use 273.15 for the Kelvin offset; the
`273.16` reference in `parameter_defaults_corrections.md` Section 3.6
applies to the v2 parity utility, not the v1 NSM1 implementation).

### 14. Water vapor partial pressure pwv

- Fortran (`modN2.f90:43-44`):
  `P_wv = exp(11.8571 - 3840.70 / TwaterK - 216961.0 / TwaterK^2)`.
- v1 (`processes.py:2878-2886`):
  `np.exp(11.8571 - 3840.70 / TwaterK - 216961 / TwaterK ** 2)`.
- v3 (`n2.py:90-95`):
  `np.exp(11.8571 - 3840.70 / t_water_k - 216961.0 / t_water_k**2)`.

Match.

### 15. N2sat formula

- Fortran (`modN2.f90:47`):
  `N2sat = 2.8E+4 * KHN2_tc * 0.79 * (pressure_atm - p_wv)`,
  with negative-clip to zero at line 50.
  `pressure_atm` is the module-level state in atm (set externally).
- v1 (`processes.py:3470-3487`):
  `N2sat = 2.8E+4 * KHN2_tc * 0.79 * (pressure_mb * 0.000986923 - pwv)`,
  with negative-clip to 1e-6 at line 3485.
- v3 (`n2.py:98-113`):
  `2.8e4 * khn2 * 0.79 * (pressure_mb * MB_TO_ATM - pwv_atm)`,
  with `MB_TO_ATM = 1.0 / 1013.25` (line 66) and negative-clip to 1e-6.

Minor deviation. The mb→atm scalar differs: v1 literal `0.000986923` versus
v3 `1.0/1013.25 ≈ 0.0009869232667...`. Documented in
`parameter_defaults_corrections.md` Section 3.7 and absorbed in
`tests/test_5_n2_calculations_v2.py` with `rtol=1e-6`. Agreement to ~7
significant figures; not a correctness concern.

### 16. Atmospheric exchange flux

- Fortran (`modN2.f90:52`):
  `N2_Reaeration = 1.034 * ka_tc * (N2sat - N2)`, where `ka_tc =
  kah_tc + kaw_tc / depth` (`modGlobalParam.f90:247`).
- v1 (`processes.py:3490-3504`):
  `dN2dt = 1.034 * ka_tc * (N2sat - N2)`.
- v3 (`n2.py:332`):
  `atm_exchange = 1.034 * ka_tc_value * (n2_sat - n2_state)`,
  where `ka_tc_value` is computed via `clearwater_modules_v3.utils.reaeration`
  matching the same `kah + kaw/depth` form.

Match.

### 17. Denitrification source

- Fortran: N2 module does not have an explicit denit source term; the
  `dN2dt` budget at `modN2.f90:54` is `dN2dt = N2_Reaeration` only.
  Denitrification mass loss is tracked in modNitrogen but not added back to
  N2.
- v1 (`processes.py:3490-3504`): `dN2dt = 1.034 * ka_tc * (N2sat - N2)`,
  no denit coupling (matches Fortran).
- v3 (`n2.py:340-350`): adds
  `denit_source = nitrogen_process.denitrification_flux_rate` (mg-N/L/d,
  positive magnitude) to the rate when `use_nitrogen` is wired.
  `rate = atm_exchange + denit_source`.

Minor deviation (extension). v3 closes the N mass balance by routing
denitrification into N2 production, which Fortran and v1 do not do. This
is an intentional Phase 3.4 design decision (Item 1 in the task brief).
Under no-coupling (v1 parity test), `denit_source = 0` and the form
collapses to v1.

### 18. TDG derived variable

- Fortran (`modN2.f90:69-72`):
  `TDG = N2 / N2sat` always; if `use_DOX`, overwrites with the
  oxygen-weighted form `(79 * N2/N2sat + 21 * DOX/O2sat)`.
- v1 (`processes.py:3523-3541`):
  `xr.where(use_DOX, (79.0 * N2 / N2sat) + (21.0 * DOX / DOX_sat),
  N2/N2sat)`.
- v3 (`n2.py:375-383`):
  `tdg = n2_new / n2_sat` only (simple form). Documented at
  module docstring lines 22-24 and inline at line 371-374 as a Phase 3
  scope decision; the oxygen-weighted form is deferred until DOX is
  wired up in Phase 5.

Minor deviation (scope deferral). v3 1.0.0 implements only the
non-oxygen-weighted form. Under `use_DOX=False` the v1/Fortran result
matches v3 exactly. Re-enabling the weighted form is a Phase 5 task per
the inline note.

## Alkalinity

### 19. Nitrification consumption

- Fortran (`modAlkalinity.f90:96`):
  `Alk_Nitrification = ralkn * NH4_Nitrification * 50000`,
  where `NH4_Nitrification` (from `modNitrogen`) already includes the
  `(1 - exp(-KNR*DOX))` Monod attenuation.
- v1 (`processes.py:3284-3319`):
  `r_alkn * (1 - np.exp(-KNR * DOX)) * knit_tc * NH4 * 50000`
  (re-applies the Monod factor locally inside `Alk_nitrification`).
- v3 (`alkalinity.py:251-274`):
  `r_alkn * nitrification_flux_rate * 50000`,
  where `nitrification_flux_rate` is read from `Nitrogen` after Nitrogen.run
  (already includes the Monod factor).

Minor deviation (architectural). Documented in
`parameter_defaults_corrections.md` Section 3.3 and
`tests/test_5_alkalinity_calculations_v2.py` module docstring lines 21-30.
Equivalent to v1 under matched parameters and matched DOX-Monod factor;
v3 follows the Fortran single-source-of-truth pattern (Fortran also reads
`NH4_Nitrification` from modNitrogen rather than recomputing the Monod term).

### 20. Denitrification production

- Fortran (`modAlkalinity.f90:103`):
  `Alk_Denit = ralkden * NO3_Denit * 50000`.
- v1 (`processes.py:3246-3281`):
  `r_alkden * (1.0 - DOX/(DOX + KsOxdn)) * kdnit_tc * NO3 * 50000`
  (re-applies oxygen-inhibition factor locally).
- v3 (`alkalinity.py:276-299`):
  `r_alkden * denitrification_flux_rate * 50000`,
  reading the pre-attenuated flux from Nitrogen.

Minor deviation (architectural). Same pattern as #19; equivalent under
matched parameters. v3 matches the Fortran routing.

### 21. Algal photosynthesis (NH4 vs NO3 fractionation)

- Fortran (`modAlkalinity.f90:78`):
  `Alk_ApGrowth = (ralkca * ApUptakeFr_NH4 - ralkcn * (1 - ApUptakeFr_NH4))
  * rca(r) * ApGrowth * 50000`.
- v1 (`processes.py:3322-3342`):
  `(r_alkaa * ApUptakeFr_NH4 - r_alkan * (1 - ApUptakeFr_NH4))
  * ApGrowth * rca * 50000`.
- v3 (`Alkalinity._floating_algae_growth_alk_flux`):
  `(r_alkaa * ap_uptake_fr_nh4 - r_alkan * (1 - ap_uptake_fr_nh4))
  * ap_growth * rca * EQ_TO_MG_CACO3` where `EQ_TO_MG_CACO3 = 50000` and
  `rca = self.AWc / self.AWa`.

**Defect NSM1-CA-1 (CRITICAL) — CORRECTED 2026-05-16 (gold-standard
spec A1).** This entry previously read `rca = self.AWc` and was
incorrectly verdicted "Match": pre-fix v3 bound the *raw* stoichiometric
weight `AWc` (=40) where Fortran (`rca(r)`, `modAlgae`) and v1
(`processes.py:337-347`, `rca = AWc/AWa`) both use the *intensive*
carbon:chlorophyll ratio `AWc/AWa` (=0.04 mg-C/ug-Chla). The raw form
overstated the floating-algae alkalinity flux by `AWa = 1000x`. The fix
binds `rca = self.AWc / self.AWa`, mirroring `carbon.py:495` and v1.

Now a true Match. Stoichiometric ratio names differ (Fortran
`ralkca/ralkcn`, v1/v3 `r_alkaa/r_alkan`) but the values
`14/106/12/1000` and `18/106/12/1000` are identical, and the
carbon:chlorophyll conversion is now the intensive `AWc/AWa` on all
three (Fortran/v1/v3).

### 22. Algal respiration source

- Fortran (`modAlkalinity.f90:79`):
  `Alk_ApRespiration = ralkca * rca(r) * ApRespiration * 50000`.
- v1 (`processes.py:3345-3361`):
  `ApRespiration * r_alkaa * 50000 * rca`.
- v3 (`Alkalinity._floating_algae_respiration_alk_source`):
  `ap_resp * self.r_alkaa * (self.AWc / self.AWa) * EQ_TO_MG_CACO3`.

**Defect NSM1-CA-1 (CRITICAL) — CORRECTED 2026-05-16 (gold-standard
spec A1).** Previously `ap_resp * self.r_alkaa * self.AWc * ...` and
incorrectly verdicted "Match": v3 used the raw weight `AWc` where
Fortran/v1 use the intensive `rca = AWc/AWa`, a 1000x overstatement of
the algal-respiration alkalinity source. Fixed to `self.AWc / self.AWa`.
Now a true Match.

### 23. Benthic algae growth and respiration

- Fortran (`modAlkalinity.f90:87-88`):
  `Alk_AbGrowth = Fb(r) * (ralkca * AbUptakeFr_NH4 - ralkcn *
  (1 - AbUptakeFr_NH4)) * rcb(r) * AbGrowth / depth * 50000`;
  `Alk_AbRespiration = Fb(r) * ralkca * rcb(r) * AbRespiration / depth
  * 50000`.
- v1 (`processes.py:3364-3410`): same form, with `1/depth` factor and
  `Fb` multiplication; uses `r_alkba`/`r_alkbn`.
- v3 (`_benthic_algae_growth_alk_flux` /
  `_benthic_algae_respiration_alk_source`): same form; uses
  `r_alkba`/`r_alkbn`/`rcb = BWc/BWd`/`Fb`. The `1/depth` divider is
  applied explicitly in both helpers.

**Defect NSM1-CA-1 (CRITICAL) — CORRECTED 2026-05-16 (gold-standard
spec A1).** Previously bound raw `BWc` (=40) and was incorrectly
verdicted "Match": Fortran (`rcb(r)`) and v1 use the intensive
carbon:dry-weight ratio `BWc/BWd` (=0.4 mg-C/mg-D); the raw form
overstated the benthic-algae alkalinity terms by `BWd = 100x`. Fixed
to `rcb = self.BWc / self.BWd`. Now a true Match.

### 24. Net dAlk/dt sign convention

- Fortran (`modAlkalinity.f90:109`):
  `dAlkdt = -Alk_ApGrowth + Alk_ApRespiration - Alk_Nitrification
  + Alk_Denit - Alk_AbGrowth + Alk_AbRespiration`.
- v1 (`processes.py:3413-3431`):
  `Alk_denitrification - Alk_nitrification - Alk_algal_growth
  + Alk_algal_respiration - Alk_benthic_algae_growth
  + Alk_benthic_algae_respiration`.
- v3 (`alkalinity.py:461-468`):
  `denit_source - nitr_sink - algal_growth_sink + algal_resp_source
  - balgae_growth_sink + balgae_resp_source`.

Match. Sign convention: growth terms enter as sinks (subtracted); the
NH4-vs-NO3 fractionation inside the growth flux flips sign internally so
NO3-uptake-dominated growth correctly produces alkalinity in the net.

## Parameter defaults audit

For each constituent, comparison of v3 `DEFAULTS` against Fortran
`modGlobalParam.f90` / per-module `Initialize*` defaults and v1
`constants.py`.

### POM

| Parameter | Fortran | v1 | v3 | Status |
| --- | --- | --- | --- | --- |
| `kpom_20` | 0.01 (`modPOM.f90:37`) | 0.1 | 0.1 | v1<->v3 match; Fortran differs by 10x |
| `kpom_theta` | 1.047 | 1.047 | 1.047 | Match |
| `h2` | 0.1 (`modGlobalParam.f90:134`) | 0.1 | 0.1 | Match (FIXME(phase1-audit) noted) |
| `vsoc` | 0.01 (`modGlobalParam.f90:104`) | 0.01 | 0.01 | Match |
| `fcom` | 0.4 (`modGlobalParam.f90:108`) | 0.4 | 0.4 | Match |
| `vb` | 0.0025 m/yr (`modGlobalParam.f90:138`) | 0.01 m/d | 0.01 m/d | v1<->v3 match; Fortran in m/yr (different unit convention) |

Note: Fortran `kpom2_20=0.01` is 10x smaller than v1/v3 `kpom_20=0.1`.
Both are documented in their respective sources; the v1 constants table
overrides Fortran's 0.01. v3 follows v1.

### CBOD

| Parameter | Fortran | v1 | v3 | Status |
| --- | --- | --- | --- | --- |
| `KsOxbod` | 0.5 (`modCBOD.f90:41`) | 0.5 | 0.5 | Match |
| `kbod_20` | 0.12 (`modCBOD.f90:32`) | 0.12 | 0.12 | Match |
| `ksbod_20` | 0.0 (`modCBOD.f90:36`) | 0.0 | 0.0 | Match (FIXME(phase1-audit) noted) |
| `kbod_theta` | 1.047 | 1.047 | 1.047 | Match |
| `ksbod_theta` | 1.024 (`modCBOD.f90:36`) | 1.047 | 1.047 | v1<->v3 match; Fortran 1.024 |

Minor deviation. Fortran sets `ksbod_theta=1.024` while v1 and v3 both use
`1.047`. Under `ksbod_20=0` the Arrhenius correction has no effect, so
this is dormant under defaults. v3 follows v1.

### Pathogen

| Parameter | Fortran | v1 | v3 | Status |
| --- | --- | --- | --- | --- |
| `kdx_20` | 0.8 (`modPathogen.f90:32`) | 0.8 | 0.8 | Match |
| `kdx_theta` | 1.07 | 1.07 | 1.07 | Match |
| `apx` | 1.0 | 1 | 1.0 | Match (FIXME(phase1-audit) noted) |
| `vx` | 1.0 | 1 | 1.0 | Match (FIXME(phase1-audit) noted) |

### N2

v1 / v3 / Fortran all agree: no N2-specific defaults; saturation derives
from `pressure_mb` (1013.25 in v3 after the Section 1.7 correction;
v1 default 2026.5 was the bug). All Henry's-law and pwv constants are
hard-coded at the formula site:

| Constant | Fortran | v1 | v3 | Status |
| --- | --- | --- | --- | --- |
| `KH(298.15K)` | 0.00065 | 0.00065 | 0.00065 | Match |
| dH/R | 1300 K | 1300 K | 1300 K | Match |
| pwv coeffs | 11.8571, 3840.70, 216961.0 | same | same | Match |
| N2sat factor | 2.8e4 | 2.8e4 | 2.8e4 | Match |
| 0.79 (N2 vol fraction) | 0.79 | 0.79 | 0.79 | Match |
| Reaeration weight | 1.034 | 1.034 | 1.034 | Match |
| mb->atm | externally set in atm | 0.000986923 | 1.0/1013.25 | v3<->v1 ~7 sig fig agreement (Section 3.7) |

### Alkalinity

| Parameter | Fortran | v1 | v3 | Status |
| --- | --- | --- | --- | --- |
| `r_alkaa` | 14/106/12/1000 (`modAlkalinity.f90:49`) | same | same | Match |
| `r_alkan` (v3 / v1; Fortran `ralkcn`) | 18/106/12/1000 (line 50) | same | same | Match |
| `r_alkn` | 2/14/1000 (line 53) | same | same | Match |
| `r_alkden` | 4/14/1000 (line 54) | same | same | Match |
| `r_alkba` | 14/106/12/1000 (Fortran reuses ralkca for benthic) | same | same | Match |
| `r_alkbn` | 18/106/12/1000 (Fortran reuses ralkcn for benthic) | same | same | Match |
| `EQ_TO_MG_CACO3` | 50000 | 50000 | 50000 | Match |

Note: Fortran maintains a single `ralkca`/`ralkcn` pair for both algal and
benthic stoichiometry; v1 and v3 carry separate `r_alkaa`/`r_alkan` and
`r_alkba`/`r_alkbn` with identical numerical defaults. No runtime impact;
allows independent tuning if a future calibration distinguishes them.

## Conclusions

### Required actions before LimnoTech review

None. All seven minor deviations are either:

1. Documented in `src/clearwater_modules_v3/parameter_defaults_corrections.md`
   Section 3 (deviations 3.3 Alkalinity routing, 3.4 Pathogen `Fr_PAR`,
   3.5 CBOD `1/depth`, 3.7 mb→atm).
2. Pinned in parity test docstrings
   (`tests/test_5_pathogen_calculations_v2.py`,
   `tests/test_5_n2_calculations_v2.py`,
   `tests/test_5_cbod_calculations_v2.py`,
   `tests/test_5_alkalinity_calculations_v2.py`).
3. Intentional scope deferrals (N2 oxygen-weighted TDG to Phase 5; CBOD
   multi-group to Phase 4+).
4. Equivalence-under-matched-parameters refactors that match the Fortran
   single-source-of-truth pattern more closely than v1 did (Alkalinity
   nitrif/denit routing).

### Acceptable deviations to document for sponsor

The four substantive deviations from v1 (carried over from earlier audits;
listed here for completeness):

1. CBOD sedimentation `ksbod_tc / depth` (v3) versus `ksbod_tc` (v1).
   Dimensionally consistent in v3; under default `ksbod_20=0` there is
   no runtime difference. If a project uses `ksbod_20 > 0` from a v1
   calibration, recalibration is required.
2. Pathogen light decay scaled by `Fr_PAR=0.47`. Absorbed in calibrated
   `apx`; at `Fr_PAR=1.0` matches v1 exactly.
3. Alkalinity routes nitrification / denitrification through the
   pre-attenuated Nitrogen flux cache rather than recomputing the
   Monod factor. Numerically equivalent under matched parameters and
   matches the Fortran routing pattern.
4. N2 budget includes denitrification source (`Nitrogen.denitrification_flux_rate`).
   This is a v3-only completion of the N mass balance; not present in
   v1 or Fortran. Adds correctness (N is conserved across NO3 -> N2 -> N2sat
   exchange) without breaking parity at zero coupling.

### Items to escalate

None for the simple constituents in scope. The `FIXME(phase1-audit)` items
in the v3 parameter modules (`h2`, `vb`, `apx`, `vx`, `ksbod_20=0`,
`q_solar` units docstring, `lambdas` disabled term) are tracked in
`parameter_defaults_corrections.md` Section 2 and remain open for a
future calibration-targeted audit; none of them produces a runtime
correctness defect under matched inputs.
