# v3 NSM1 Carbon + DOX -- Three-way audit (Fortran vs v1 vs v3)

Date: 2026-05-05
Scope: Carbon (POC, DOC, DIC) and DOX kinetic source terms; cross-process
couplings to FloatingAlgae, BenthicAlgae, Nitrogen, CBOD, POM, SedFlux.
Sources audited:
- Fortran: `Source Files/modCarbon.f90`, `Source Files/modDOX.f90`,
  `Source Files/modGlobalParam.f90` (SOD_tc, ka_tc), `Source Files/modGlobal.f90`,
  `Source Files/modAlgae.f90`, `Source Files/modBenthicAlgae.f90`,
  `Source Files/modNitrogen.f90`.
- v1: `src/clearwater_modules/nsm1/processes.py` (Carbon block 2439-2870,
  DOX block 2876-3135), `src/clearwater_modules/nsm1/constants.py`,
  `src/clearwater_modules/shared/processes.py` (ka_tc, SOD_tc, kah/kaw).
- v3: `src/clearwater_modules_v3/processes/carbon.py`,
      `src/clearwater_modules_v3/processes/dox.py`,
      `src/clearwater_modules_v3/parameters/carbon.py`,
      `src/clearwater_modules_v3/parameters/dox.py`,
      `src/clearwater_modules_v3/utils/sediment.py`,
      `src/clearwater_modules_v3/utils/reaeration.py`,
      and the v2 algae caches (`floating_algae.py`, `benthic_algae.py`)
      that v3 carbon/dox consume via getattr.

## Summary

- Counts: 6 critical, 5 minor, 17 matches, 4 observations.
- Top concerns:
  1. **C1 (Critical) -- v3 uses `AWc` and `BWc` raw stoichiometric weights as
     `rca` and `rcb` directly.** Fortran and v1 derive `rca = AWc / AWa`
     (= 40/1000 = 0.04 mg-C/ug-Chla) and `rcb = BWc / BWd` (= 40/100 = 0.4
     mg-C/mg-D). v3 passes `self.AWc` (= 40) and `self.BWc` (= 40) directly
     into the DIC and DOX algal coupling terms, scaling photosynthesis and
     respiration O2/C fluxes by 1000x for floating algae and 100x for
     benthic algae. Affects v3 Carbon `dic_algal_resp`, `dic_algal_photo`,
     `dic_balgae_resp`, `dic_balgae_photo`, and v3 DOX
     `_floating_algae_growth_flux`, `_floating_algae_respiration_flux`,
     `_benthic_algae_growth_flux`, `_benthic_algae_respiration_flux`.
     The cached `algal_*_from_mortality_rate` and `balgae_*_from_mortality_rate`
     fluxes from FloatingAlgae/BenthicAlgae already bake in the correct
     `rca = AWc/AWa` and `rcb = BWc/BWd`, so the *mortality* routing in
     v3 Carbon is correct; only the *growth and respiration* coupling is
     defective.
  2. **C2 (Critical) -- v3 DOX SOD sink omits the DOX-Monod attenuation that
     Fortran applies in `modGlobalParam.f90:254`.** Fortran computes
     `SOD_tc = SOD * theta_corr * DOX/(DOX+KsSod)`. v1 mirrors this in its
     shared `SOD_tc` helper (`shared/processes.py:180-200`). v3 deliberately
     stripped the Monod factor from `utils/sediment.py:SOD_tc` (per Phase 1.1)
     and never reapplies it inside `DOX._sod_flux`. Documented in v3 dox.py
     module docstring (lines 60-63), but the consequence is that under
     hypoxic conditions the v3 SOD sink stays at its full Arrhenius value
     instead of being throttled toward zero.
  3. **C3 (Critical) -- v3 DIC budget omits the CBOD oxidation source.**
     Fortran (`modCarbon.f90:262-266`) and v1 (`processes.py:2854`) both add
     `DIC_CBOD_oxidation` to dDIC/dt. v3 `carbon.py:451-459` has no CBOD
     coupling in its DIC term sum, even though the v3 CBOD process caches a
     `cbod_oxidation_rate` that is already wired into DOX. Net effect: DIC
     evolution understates the C produced from CBOD oxidation.
  4. **C4 (Critical) -- v3 Carbon adds a DOX-Monod attenuation to POC
     hydrolysis that neither Fortran nor v1 apply.** v3 carbon.py:372 writes
     `poc_hydrolysis = kpoc_tc * poc * dox_attenuation`. Fortran
     (`modCarbon.f90:170`) writes `POC_DOC_Hydrolysis = kpoc_tc * POC` with
     no DOX dependence, and v1 (`processes.py:2455-2465`) likewise does not
     attenuate. POC hydrolysis is a physical/chemical process (cell-wall
     fragmentation, leaching) that is not biochemically O2-limited in the
     reference model.
  5. **C5 (Critical) -- v3 atmospheric reaeration short-circuit blocks O2
     transfer when `kaw_20_user == 0`.** v3 dox.py:620-627 sets
     `ka_tc_value = 0.0` when both menu options are 1 and both user values
     are zero. With v3 corrected defaults `kaw_20_user = 0` and
     `kah_20_user = 0` and the menu defaults at option 1 (user-defined),
     a default-instantiated DOX has zero atmospheric reaeration. Fortran
     uses `kah%rc20 = 1.0; kaw%rc20 = 0.0` defaults (`modGlobalParam.f90:113-117`).
     v1 retains the v1 default `kah_20_user = 999`, which itself is invalid
     (the corrections doc Section 1 records this). v3's correction to 0
     means the Phase 5.B DOX runs with no reaeration unless the user
     explicitly opts in. This is a documented design choice (see corrections
     doc) but the audit flags it because the resulting default behavior is
     physically incorrect.
  6. **C6 (Critical) -- v3 omits the salinity correction on O2sat.** Fortran
     (`modDOX.f90:97-99`) applies `O2sat *= exp(-Salinity * (0.017674 -
     10.754/Tk + 2140.7/Tk^2))`. v3 `dox_sat_apha` and v1 `DOX_sat` both
     omit this term. Effectively zero impact for fresh water (Salinity=0
     gives factor 1.0), but coupling v3 to brackish or estuarine
     applications would silently overstate dissolved oxygen saturation.

## Carbon

### 1. POC mineralization (POC -> DOC hydrolysis)

- Fortran (`modCarbon.f90:170`):
  `POC_DOC_Hydrolysis = kpoc_tc * POC` (no DOX dependence).
- v1 (`processes.py:2455-2465` `POC_hydrolysis`): `kpoc_tc * POC`
  (no DOX dependence).
- v3 (`carbon.py:372`):
  `poc_hydrolysis = kpoc_tc_value * poc * dox_attenuation`
  with `dox_attenuation = dox / (KsOxmc + dox)`.

**Finding C4 (Critical, Scientific correctness).** v3 multiplies POC
hydrolysis by `DOX/(KsOxmc+DOX)`, which neither Fortran nor v1 apply.
Under DOX = 1 mg/L, KsOxmc = 1 mg/L, v3 attenuates the POC -> DOC flux to
50% of the Fortran value. The `f_pocp` and `f_pocb` mortality routings are
unaffected (they bypass hydrolysis). Recommendation: remove the
`dox_attenuation` factor from `poc_hydrolysis` and document the divergence
or restore parity with Fortran.

### 2. POC settling

- Fortran (`modCarbon.f90:171`): `POC_Settling = vsoc(r) / depth * POC`.
- v1 (`processes.py:2469-2481`): `vsoc / depth * POC`.
- v3 (`carbon.py:373`): `self.vsoc / depth * poc`.

Match.

### 3. POC from algal mortality (floating)

- Fortran (`modCarbon.f90:174`):
  `ApDeath_POC = rca(r) * ApDeath * Fpocp(r)` with `rca = AWc/AWa = 0.04`.
- v1 (`processes.py:2484-2502`): `f_pocp * kdp_tc * rca * Ap`.
- v3 (`carbon.py:554-568`): reads `algal_poc_from_mortality_rate` from
  the FloatingAlgae sibling. v2 cache
  (`floating_algae.py:363`) stores `f_pocp * rca * ap_death` with the
  correct `rca = AWc / AWa`.

Match (via cache).

### 4. POC from benthic algal mortality

- Fortran (`modCarbon.f90:180`):
  `AbDeath_POC = rcb(r) * AbDeath * Fb(r) * Fw(r) * Fpocb(r) / depth`,
  where `rcb = BWc / BWd = 0.4`.
- v1 (`processes.py:2505-2529`): `(1/depth) * f_pocb * kdb_tc * rcb * Ab * Fb * Fw`.
- v3 (`carbon.py:607-624`): reads `balgae_poc_from_mortality_rate` cache
  populated by v2 BenthicAlgae as
  `f_pocb * fb * fw * rcb * ab_death / depth` with correct `rcb = BWc / BWd`.

Match (via cache).

### 5. dPOC/dt assembly

- Fortran (`modCarbon.f90:185`):
  `dPOCdt = ApDeath_POC + AbDeath_POC - POC_DOC_Hydrolysis - POC_Settling`.
- v1 (`processes.py:2532-2546`): same.
- v3 (`carbon.py:377-382`): same structure
  (`poc_algal_mortality + poc_balgae_mortality - poc_hydrolysis - poc_settling`).

Match (with the per-term defect in POC hydrolysis under finding C4).

### 6. DOC oxidation to DIC

- Fortran (`modCarbon.f90:198`):
  `DOC_DIC_Oxidation = DOX/(DOX+KsOxmc) * kdoc_tc * DOC` when use_DOX,
  else `kdoc_tc * DOC`.
- v1 (`processes.py:2629-2647` `DOC_DIC_oxidation`):
  `xr.where(use_DOX, DOX/(KsOxmc+DOX) * kdoc_tc * DOC, kdoc_tc * DOC)`.
- v3 (`carbon.py:385`):
  `doc_oxidation = kdoc_tc_value * doc * dox_attenuation`
  unconditionally.

**Finding C7 (Minor, Conditional logic).** v3 always applies the Monod
attenuation; Fortran and v1 gate it on `use_DOX`. With `use_DOX = False`
(rare in practice; default is True), v3 silently zeros the DOC oxidation
flux when DOX = 0 instead of using the unattenuated form. Net consequence
under default `use_DOX = True` is zero. Recommendation: gate the
attenuation factor on `self.use_DOX` for full parity, and to keep
behaviour stable when DOX is intentionally disabled.

### 7. DOC from algal mortality

- Fortran (`modCarbon.f90:207, 213`):
  floating: `rca * ApDeath * (1 - Fpocp)`;
  benthic: `rcb * AbDeath * Fb * Fw * (1 - Fpocb) / depth`.
- v1 (`processes.py:2565-2610`): same form with `rca`, `rcb`.
- v3 (`carbon.py:570-585, 626-642`): reads
  `algal_doc_from_mortality_rate` and `balgae_doc_from_mortality_rate`
  caches; both already bake in correct `rca`, `rcb`, `f_pocp`, `f_pocb`.

Match (via cache).

### 8. DOC from POM hydrolysis

- Fortran: not present in `modCarbon.f90` DOC budget. POM is handled
  separately in `modPOM.f90`; the Fortran DOC path does not include POM
  hydrolysis as a DOC source.
- v1 (`processes.py:2651-2667` `dDOCdt`): does not include POM hydrolysis.
- v3 (`carbon.py:391-404`): adds `pom_hydrolysis_rate` from POM sibling.

**Finding C8 (Observation).** v3 wires POM -> DOC explicitly via the
`pom_hydrolysis_rate` cache. Neither Fortran nor v1 has this coupling in
their DOC equation. v3 module docstring (line 24) calls this out as a
design improvement. The Phase 3.2/Items 2-3 trail confirms this is
intentional. Categorise as a v3 enhancement that LimnoTech may want to
review if they expect strict legacy parity. Verify whether the POM
process's `pom_hydrolysis_rate` correctly excludes the fraction Fortran
routes elsewhere (e.g. to NH4 via OrgN hydrolysis).

### 9. DOC oxidation as DOX sink (cached `doc_dic_oxidation_rate`)

- Fortran (`modDOX.f90:124`):
  `O2_DOC_Oxidation = roc * DOC_DIC_Oxidation`.
- v1 (`processes.py:3002-3015` `DOX_DOC_oxidation`):
  `roc * DOC_DIC_oxidation` when `use_DOC`.
- v3 (`carbon.py:469`): caches
  `self.doc_dic_oxidation_rate = doc_oxidation`.
  v3 (`dox.py:524-540` `_doc_oxidation_flux`):
  `roc * carbon.doc_dic_oxidation_rate` when `use_carbon and use_DOC`.

Match.

### 10. dDOC/dt assembly

- Fortran (`modCarbon.f90:218`):
  `dDOCdt = ApDeath_DOC + AbDeath_DOC + POC_DOC_Hydrolysis - DOC_DIC_Oxidation`.
- v1 (`processes.py:2651-2667`): same (no POM term).
- v3 (`carbon.py:398-404`): `poc_hydrolysis + algal_doc_mort + balgae_doc_mort
  + pom_hydrolysis - doc_oxidation`.

Match modulo the POM addition (finding C8) and the POC hydrolysis Monod
defect propagating into the DOC source (finding C4).

### 11. DIC reaeration with atmosphere (Henry's law for CO2)

- Fortran (`modCarbon.f90:236-238`):
  `KH_tc = 10**(2385.73/Tk + 0.0152642*Tk - 14.0184)`;
  `DIC_Reaeration = 0.923 * ka_tc * (KH_tc * pco2 / 1e6 - Fco2 * DIC)`.
- v1 (`processes.py:2687-2714`): same formula.
- v3 (`carbon.py:124-136, 418-421`): same formula in `henrys_k_co2` and
  `co2_reaeration`.

Match.

**Finding C9 (Observation).** Both Fortran and v1 leave the dimensional
inconsistency that `KH_tc * pco2/1e6` is in mol-C/L while `Fco2 * DIC` is
in (mg-C/L * unitless = mg-C/L), giving a difference of mol/L - mg/L.
Fortran resolves this by integrating dDIC/dt in mol/L for the carbon
section (`DOC_DIC_Oxidation/12000` etc.); v1 inherits the formula but
treats DIC as mg/L throughout. v3 retains the v1 numerical form (line
411-417 in carbon.py docstring acknowledges this). Recommend escalating
this as an open question for the LimnoTech review: should v3 1.x land a
proper carbonate-system unit reconciliation, or is it acceptable to
preserve the v1 numerical form?

### 12. DIC from DOC oxidation

- Fortran (`modCarbon.f90:268`):
  `dDICdt = DOC_DIC_Oxidation / 12000.0 + ...` (DOC oxidation contributes
  to DIC, in mol-C/L since /12000 converts mg-C to mol-C).
- v1 (`processes.py:2834-2854` `dDICdt`): does **not** include
  `DOC_DIC_oxidation` in the DIC budget. This is a v1 omission relative
  to Fortran.
- v3 (`carbon.py:451-459`): includes `+ doc_oxidation` (in mg-C/L/d, no
  /12000 conversion).

**Finding C10 (Minor, Documentation-to-code fidelity).** v3 partially
restores Fortran's DOC -> DIC coupling that v1 dropped, but does so in
mg-C/L/d rather than mol-C/L/d. Combined with the DIC reaeration
mol/L-vs-mg/L mismatch (finding C9), the DIC budget is internally
inconsistent in units. Recommendation: either (a) match v1's omission and
defer DIC -> DOC accounting to a v3 1.x carbonate solver, or (b) follow
through on the unit reconciliation across the entire DIC budget.

### 13. DIC from algal respiration / sink from algal growth

- Fortran (`modCarbon.f90:247-248`):
  `ApRespiration_DIC = rca * ApRespiration / 12000.0` (rca = AWc/AWa = 0.04).
  `DIC_ApGrowth = rca * ApGrowth / 12000.0`.
- v1 (`processes.py:2717-2748`):
  `ApRespiration * rca / 12000`, `ApGrowth * rca / 12000`.
  Caller is responsible for passing `rca = AWc/AWa`.
- v3 (`carbon.py:429-430`):
  `dic_algal_resp = algae_respiration * self.AWc / 12000.0`,
  `dic_algal_photo = algae_growth * self.AWc / 12000.0`.

**Finding C1 (Critical, Scientific correctness).** v3 uses `self.AWc`
(default 40, raw stoichiometric weight) where Fortran uses
`rca = AWc / AWa = 0.04`. Result: v3 DIC algal terms are 1000x larger
than they should be. With default `AWc = 40`, `AWa = 1000`, and an
algal respiration rate of 0.5 ug-Chla/L/d, Fortran/v1 yield
`0.5 * 0.04 / 12000 ≈ 1.7e-6 mg-C/L/d`; v3 yields
`0.5 * 40 / 12000 ≈ 1.7e-3 mg-C/L/d`. The existing parity test
(`tests/test_5_carbon_calculations_v2.py:343-412`) passes only because it
calls `v1.DIC_algal_respiration(rca=AWc=40)`, which is already
mis-parameterising v1 with the same wrong value.
**Recommendation**: change v3 carbon.py:429-430 to
`algae_respiration * self.AWc / self.AWa / 12000.0` (or compute
`rca = self.AWc / self.AWa` once at the top of run, mirroring the
v2 floating_algae.py:358 derivation).

### 14. DIC from benthic algal respiration / sink from benthic algal growth

- Fortran (`modCarbon.f90:255-256`):
  `AbRespiration_DIC = rcb * AbRespiration * Fb / depth / 12000.0`,
  `DIC_AbGrowth = rcb * AbGrowth * Fb / depth / 12000.0`,
  where `rcb = BWc / BWd = 40 / 100 = 0.4`.
- v1 (`processes.py:2751-2789`): same with `rcb` parameter.
- v3 (`carbon.py:436-441`):
  `dic_balgae_resp = balgae_respiration * self.BWc * self.Fb / depth / 12000.0`.

**Finding C1 (Critical, same root cause as floating algae).** v3 uses
`self.BWc` (= 40) instead of `rcb = BWc / BWd` (= 0.4). v3 DIC benthic
algal terms are 100x too large. Recommendation: derive
`rcb = self.BWc / self.BWd` once and use it instead of raw `self.BWc`.

### 15. DIC sediment release (`JDIC`)

- Fortran (`modCarbon.f90:240-244`):
  `if use_SedFlux: DICfromBed = JDIC / depth / 12000.0`
  `else: DICfromBed = SOD_tc / roc / depth / 12000.0`.
- v1 (`processes.py:2817-2830` `DIC_sed_release`):
  `SOD_tc / roc / depth / 12000` unconditionally.
- v3 (`carbon.py:446-449`):
  `if use_SedFlux: dic_sed_release = JDIC / depth / 12000.0`
  `else: dic_sed_release = 0.0`.

**Finding C11 (Minor, Scope/parity).** v3 currently only supports the
SedFlux branch (with `JDIC` user-supplied). Fortran's non-SedFlux fallback
(use SOD-derived DIC release) is not implemented. v3 docstring lines
443-449 documents this as Phase 5.A scope. With default v3 `JDIC = 0.0`
and `use_SedFlux = False`, the v3 DIC sediment release is identically zero,
whereas Fortran/v1 release `SOD_tc / roc / depth / 12000`. Recommendation:
add the SOD-derived fallback in Phase 5.5 or document the deviation in
the corrections doc.

### 16. DIC from CBOD oxidation

- Fortran (`modCarbon.f90:262-266`): sums `CBOD_Oxidation(i)` across
  groups, then `CBOD_DIC_Oxidation = sum / roc / 12000.0`.
- v1 (`processes.py:2793-2814` `DIC_CBOD_oxidation`): per-group
  `(1/roc) * (DOX/(KsOxbod+DOX)) * kbod_tc * CBOD / 12000`.
- v3 (`carbon.py:451-459`): **no CBOD term in dDIC/dt**.

**Finding C3 (Critical, Missing source term).** v3 DIC budget omits the
CBOD oxidation source. Fortran and v1 both include it. The v3 CBOD
process (`processes/cbod.py`) caches `cbod_oxidation_rate` (mg-O2/L/d);
to match Fortran semantics, v3 carbon.py should add
`+ self.cbod_process.cbod_oxidation_rate / self.roc / 12000.0` (or, in
mg-C/L/d if the unit reconciliation lands, `cbod_oxidation_rate / roc`).
Recommendation: wire CBOD into the Carbon DIC budget the same way it is
wired into DOX (`carbon_process.cbod_oxidation_rate`), gated on
`self.use_cbod` if a CBOD process is registered.

### 17. dDIC/dt assembly

- Fortran (`modCarbon.f90:268-269`):
  `DOC_DIC_Oxidation/12000 + DIC_Reaeration + DICfromBed
  + ApRespiration_DIC - DIC_ApGrowth + AbRespiration_DIC - DIC_AbGrowth
  + CBOD_DIC_Oxidation`.
- v1 (`processes.py:2854`):
  `Atm_CO2_reaeration + DIC_algal_respiration - DIC_algal_photosynthesis
  + DIC_benthic_algae_respiration - DIC_benthic_algae_photosynthesis
  + DIC_CBOD_oxidation + DIC_sed_release` (no DOC oxidation term).
- v3 (`carbon.py:451-459`):
  `doc_oxidation + co2_reaeration + dic_algal_resp - dic_algal_photo
  + dic_balgae_resp - dic_balgae_photo + dic_sed_release` (no CBOD term).

Affected by findings C1, C3, C9, C10, C11.

## DOX

### 1. O2 saturation

- Fortran (`modDOX.f90:78-99`): four-coefficient log polynomial in 1/Tk
  giving Benson-Krause O2sat, then pressure correction
  `O2sat *= P_atm * (1 - pwv/P_atm) * (1 - alpha*P_atm) /
  ((1 - pwv) * (1 - alpha))`, then salinity correction
  `O2sat *= exp(-Salinity * (0.017674 - 10.754/Tk + 2140.7/Tk^2))`.
- v1 (`processes.py:2901-2923` `DOX_sat`): four-coefficient log polynomial,
  pressure correction. **No salinity correction.**
- v3 (`dox.py:150-190` `dox_sat_apha`): four-coefficient log polynomial
  matching Fortran, pressure correction with `pressure_atm = pressure_mb * 0.000986923`.
  **No salinity correction.**

Match against v1.

**Finding C6 (Critical, Missing term for non-fresh water).** v3 (and v1)
omit the salinity-based reduction of O2sat. For freshwater (Salinity = 0)
the omitted factor evaluates to 1.0, so freshwater simulations are
unaffected. For brackish or estuarine applications v3 will overstate
O2sat by up to ~20% at typical seawater salinity (35 ppt, ~25 deg C).
Recommendation: add a `salinity_psu` input variable (defaulting to 0.0)
and apply the Fortran-form correction; document the addition in the
corrections doc.

### 2. Atmospheric reaeration

- Fortran (`modDOX.f90:110`): `O2_Reaeration = ka_tc * (O2sat - DOX)`.
- v1 (`processes.py:2927-2939` `Atm_O2_reaeration`):
  `ka_tc * (DOX_sat - DOX)`.
- v3 (`dox.py:391-398, 654`):
  `_atm_reaeration_flux(dox, dox_sat, ka_tc) = ka_tc * (dox_sat - dox)`.

Match in form.

**Finding C5 (Critical, Default behaviour).** v3 short-circuits
`ka_tc = 0` when both `kah_20_user == 0` and `kaw_20_user == 0` and both
menu options are 1 (`dox.py:620-627`). With v3's corrected defaults
(`kah_20_user = 0`, `kaw_20_user = 0`, both menus = 1) the
default-instantiated DOX has zero atmospheric reaeration. This is
arguably worse than v1's `999` placeholder (which Phase 0.2 audit flagged
as invalid), because the v3 short-circuit silently produces a physically
wrong simulation rather than blowing up. Fortran defaults to
`kah%rc20 = 1.0; kaw%rc20 = 0.0` (`modGlobalParam.f90:113-117`).
Recommendation: set v3 `kah_20_user = 1.0` (matching Fortran) or change
the default `hydraulic_reaeration_option` to `2` (Covar 1976) so that
flow-driven reaeration is computed unless the user opts into a different
formula.

### 3. Photosynthesis O2 source from algae (floating)

- Fortran (`modDOX.f90:135`):
  `O2_ApGrowth = (138/106 - 32/106 * ApUptakeFr_NH4) * roc * rca * ApGrowth`,
  with `rca = AWc/AWa = 0.04`.
- v1 (`processes.py:2942-2959` `DOX_ApGrowth`):
  `ApGrowth * rca * roc * (138/106 - 32 * ApUptakeFr_NH4 / 106)`.
- v3 (`dox.py:400-428` `_floating_algae_growth_flux`):
  `ap_growth * self.AWc * self.roc * (138/106 - 32/106 * ap_uptake_fr_nh4)`.

**Finding C1 (Critical, same defect as Carbon).** v3 uses `self.AWc` =
40 in place of `rca = AWc/AWa` = 0.04. v3 photosynthesis O2 source is
1000x larger than Fortran. Recommendation: replace `rca = self.AWc` with
`rca = self.AWc / self.AWa` (line 422 in dox.py).

The Redfield stoichiometric factor (138/106 - 32/106 * NH4_fraction)
matches Fortran and v1 exactly.

### 4. Algal respiration O2 sink (floating)

- Fortran (`modDOX.f90:136`):
  `O2_ApRespiration = roc * rca * ApRespiration`, rca = 0.04.
- v1 (`processes.py:2962-2977`): `ApRespiration * rca * roc`.
- v3 (`dox.py:430-443`): `ap_resp * self.AWc * self.roc`.

**Finding C1 (Critical).** Same 1000x error from `self.AWc`-as-rca.

### 5. Photosynthesis O2 source from benthic algae

- Fortran (`modDOX.f90:143`):
  `O2_AbGrowth = (138/106 - 32/106 * AbUptakeFr_NH4) * roc * rcb * AbGrowth * Fb / depth`,
  rcb = BWc/BWd = 0.4.
- v1 (`processes.py:3032-3054` `DOX_AbGrowth`):
  `(138/106 - 32/106 * AbUptakeFr_NH4) * roc * rcb * AbGrowth * Fb / depth`.
- v3 (`dox.py:445-475` `_benthic_algae_growth_flux`):
  `(138/106 - 32/106 * ab_uptake_fr_nh4) * self.roc * self.BWc * ab_growth * self.Fb / depth`.

**Finding C1 (Critical).** v3 uses `self.BWc` = 40 in place of
`rcb = BWc/BWd` = 0.4. 100x error.

### 6. Benthic algae respiration O2 sink

- Fortran (`modDOX.f90:144`):
  `O2_AbRespiration = roc * rcb * AbRespiration * Fb / depth`.
- v1 (`processes.py:3057-3078`):
  `roc * rcb * AbRespiration * Fb / depth`.
- v3 (`dox.py:477-489`): `self.roc * self.BWc * ab_resp * self.Fb / depth`.

**Finding C1 (Critical).** Same 100x error.

### 7. Nitrification O2 sink

- Fortran (`modDOX.f90:117-121`):
  `O2_Nitrification = ron * NH4_Nitrification` when use_NH4,
  where `NH4_Nitrification = NitrificationInhibition * knit_tc * NH4`
  and `NitrificationInhibition = 1 - exp(-KNR * DOX)` (when use_DOX).
- v1 (`processes.py:2980-2998` `DOX_Nitrification`):
  `(1 - exp(-KNR * DOX)) * ron * knit_tc * NH4` when use_NH4.
- v3 (`dox.py:491-522` `_nitrification_flux`):
  `self.ron * self.nitrogen_process.nitrification_flux_rate`,
  where the cache is computed by v2 Nitrogen as
  `NH4 * knit_tc * (1 - exp(-KNR * DOX))` (`nitrogen.py:493-509`).

Match. v3's design (read pre-attenuated flux from Nitrogen) is the
"registry rate-variable convention" called out in the spec resolved Q10.
The v3 dox.py:519-520 also includes a None-guard which is sensible.

### 8. DOC oxidation O2 sink

- Fortran (`modDOX.f90:124`): `O2_DOC_Oxidation = roc * DOC_DIC_Oxidation`.
- v1 (`processes.py:3002-3015`): `roc * DOC_DIC_oxidation` when use_DOC.
- v3 (`dox.py:524-540`): `self.roc * carbon.doc_dic_oxidation_rate` when
  `use_carbon and use_DOC`.

Match. The cached `doc_dic_oxidation_rate` in v3 carbon.py:469 is
populated each step before DOX runs (per the Phase 5.A/5.B sequencing).

### 9. CBOD oxidation O2 sink

- Fortran (`modDOX.f90:129-132`):
  `O2_CBOD_Oxidation = sum_i CBOD_Oxidation(i)` (already in mg-O2/L/d).
- v1 (`processes.py:3019-3029`): `CBOD_oxidation` (no roc multiplication).
- v3 (`dox.py:542-552`): `cbod_process.cbod_oxidation_rate`.

Match.

### 10. SOD O2 sink

- Fortran (`modGlobalParam.f90:250-256`):
  `SOD_tc = Arrhenius(SOD, TwaterC); SOD_tc *= DOX/(DOX+KsSod)` (when use_DOX).
- Fortran (`modDOX.f90:111-115`):
  `if use_SedFlux: O2_SOD = SOD_Bed / depth`
  `else: O2_SOD = SOD_tc / depth`.
- v1 (`shared/processes.py:180-200` `SOD_tc`):
  Arrhenius then `xr.where(use_DOX, SOD_tc * DOX/(DOX+KsSOD), SOD_tc)`.
  v1 dox.py: `DOX_SOD = SOD_tc / depth`.
- v3 (`utils/sediment.py:16-31` `SOD_tc`):
  pure Arrhenius, no DOX-Monod.
- v3 (`dox.py:554-573` `_sod_flux`):
  `sod = sod_tc_util(SOD_20, SOD_theta, t_water_c); return sod / depth`.
  No use_SedFlux branch, no DOX attenuation, no `KsSOD` use.

**Finding C2 (Critical, Conservation/correctness).** v3 SOD sink stays
at the Arrhenius value regardless of dissolved oxygen. Under hypoxic
conditions (DOX -> 0), Fortran's SOD attenuates to zero (the sediment
cannot deplete oxygen that is not there); v3 keeps consuming oxygen at
the full Arrhenius rate, producing negative DOX after clipping. The
clip-with-log behavior masks the conservation violation. v3 dox.py
docstring (lines 60-63) explicitly notes this is a Phase 5.5 deferred
item; the audit confirms this is a real defect with calibration impact
(SOD calibration in Fortran is implicitly DOX-coupled).

**Finding C12 (Minor, Missing branch).** v3 also lacks the
`use_SedFlux` branch that swaps `SOD_tc` for `SOD_Bed`. Documented in
v3 dox.py docstring; flag for Phase 5.5.

### 11. Cached rates contract for DOX consumption

- v3 reads `nitrification_flux_rate` from Nitrogen (positive magnitude,
  mg-N/L/d).
- v3 reads `doc_dic_oxidation_rate` from Carbon (mg-C/L/d).
- v3 reads `cbod_oxidation_rate` from CBOD (mg-O2/L/d).

Cache contract is internally consistent and matches the spec (registry
rate-variable convention, resolved Q10). Match.

### 12. dDOX/dt assembly

- Fortran (`modDOX.f90:150-151`):
  `dDOXdt = O2_Reaeration - O2_Nitrification - O2_DOC_Oxidation
  - O2_CBOD_Oxidation - O2_SOD + O2_ApGrowth - O2_ApRespiration
  + O2_AbGrowth - O2_AbRespiration`.
- v1 (`processes.py:3095-3119` `dDOXdt`): same sign convention.
- v3 (`dox.py:665-675`): same sign convention
  (`atm_reaer + algal_grow - algal_resp + balgae_grow - balgae_resp
  - nitr_sink - doc_sink - cbod_sink - sod_sink`).

Match in form, defective in magnitude due to findings C1, C2.

## Parameter defaults audit

### Carbon parameters (`parameters/carbon.py`)

| Key          | v3 default | v1 default | Fortran default | Verdict |
|--------------|-----------:|-----------:|----------------:|---------|
| `f_pocp`     | 0.9        | 0.9        | 0.9             | match   |
| `kdoc_20`    | 0.01       | 0.01       | 0.01            | match   |
| `kdoc_theta` | 1.047      | 1.047      | 1.047           | match   |
| `f_pocb`     | 0.9        | 0.9        | 0.9             | match   |
| `kpoc_20`    | 0.005      | 0.005      | 0.005           | match   |
| `kpoc_theta` | 1.047      | 1.047      | 1.047           | match   |
| `KsOxmc`     | 1.0        | 1.0        | 1.0             | match   |
| `pCO2`       | 383.0      | 383.0      | 383.0           | match   |
| `FCO2`       | 0.2        | 0.2        | 0.2             | match   |
| `roc`        | 32/12      | 32/12      | 32/12           | match   |

### DOX parameters (`parameters/dox.py`)

| Key                            | v3 default | v1 default | Fortran default | Verdict |
|-------------------------------|-----------:|-----------:|----------------:|---------|
| `ron`                          | 32/14*2 = 4.5714 | 4.5714 | 4.5714 | match |
| `KsSOD`                        | 1.0        | 1.0        | 1.0 (`KsSod`)   | match   |
| `SOD_20`                       | 1.0        | 999 (invalid) | 0.2          | **C13 (Minor)** |
| `SOD_theta`                    | 1.060      | 999 (invalid) | 1.060        | match (vs Fortran) |
| `kaw_20_user`                  | 0.0        | 999 (invalid) | 0.0           | match (vs Fortran) |
| `kah_20_user`                  | 0.0        | 999 (invalid) | 1.0           | **C5 (Critical)** |
| `kaw_theta`                    | 1.024      | 1.024      | 1.024            | match   |
| `kah_theta`                    | 1.024      | 1.024      | 1.024            | match   |
| `hydraulic_reaeration_option`  | 1          | 1          | 1                | match   |
| `wind_reaeration_option`       | 1          | 1          | 1                | match   |

**Finding C13 (Minor, Default deviates from Fortran).** v3 corrected
`SOD_20 = 1.0` g-O2/m^2/d to replace the invalid v1 sentinel `999`, but
Fortran initialises `SOD%rc20 = 0.2` (`modGlobalParam.f90:122`). The v3
default is 5x larger than Fortran's. Recommendation: align with Fortran
(`SOD_20 = 0.2`) or document the rationale in the corrections doc.

### Algal stoichiometry composed by Carbon and DOX

| Key   | v3 default | v1 default | Fortran default | Verdict |
|-------|-----------:|-----------:|----------------:|---------|
| `AWc` | 40.0       | 40.0       | 40.0            | match (raw value) |
| `AWa` | 1000.0     | 1000.0     | 1000.0          | match (raw value) |
| `BWc` | 40.0       | 40.0       | 40.0            | match (raw value) |
| `BWd` | 100.0      | 100.0      | 100.0           | match (raw value) |
| `BWa` | 3500.0     | 5000.0     | 5000.0          | **C14 (Minor)** |
| `Fb`  | 0.9        | 0.9        | 0.9             | match   |
| `Fw`  | 0.9        | 0.9        | 0.9             | match   |

**Finding C14 (Minor, Default deviation).** v3 BWa = 3500 but Fortran and
v1 use BWa = 5000. The benthic-algae chlorophyll-a stoichiometry
ratio enters the v3 DIC budget only via `rab` derivations that are not
exercised in carbon/dox; cross-check whether v3 BenthicAlgae uses BWa
elsewhere.

### Henry's law and CO2 reaeration constants

The constant 0.923 on the CO2 reaeration term, the formula
`KH = 10^(2385.73/Tk + 0.0152642*Tk - 14.0184)`, and the reference
pressure 1e6 (ppm -> atm) all match across Fortran, v1, and v3.

## Conclusions

### Required actions before LimnoTech review

1. **Fix C1** (rca/rcb derivation in v3 carbon.py and dox.py). Without
   this fix, the v3 algal photosynthesis and respiration coupling to DIC
   and DOX is off by 100--1000x. The existing parity test
   (`test_5_carbon_calculations_v2.py:343`) is not sensitive to this
   defect because it is a same-error parity comparison; replace it (or
   add a sibling test) with a Fortran-anchored value:
   `AWc=40, AWa=1000, ApRespiration=0.5 ug-Chla/L/d -> 1.667e-6 mg-C/L/d`.
2. **Fix C3** (missing CBOD oxidation source in DIC budget). The v3 CBOD
   process already exposes `cbod_oxidation_rate`; wire it into Carbon's
   DIC term sum.
3. **Fix C4** (drop the DOX-Monod factor on POC hydrolysis). Both
   Fortran and v1 treat POC hydrolysis as a non-O2-limited physical
   process. Either remove the attenuation or add a clear opt-in flag.
4. **Decide on C5 (default reaeration short-circuit).** Either set
   `kah_20_user = 1.0` (Fortran default) or change
   `hydraulic_reaeration_option` to a non-trivial Covar/Owens-Gibbs
   formula so that default DOX runs include atmospheric reaeration.

### Acceptable deviations to document

5. **C2 (SOD-DOX Monod attenuation deferred).** Already documented in
   v3 dox.py docstring as Phase 5.5 work. Add a regression test that
   exercises low-DOX SOD attenuation to ensure Phase 5.5 catches it.
6. **C8 (POM -> DOC coupling)**: v3 enhancement beyond Fortran/v1. Add to
   the corrections doc as an intentional v3 addition.
7. **C9, C10 (DIC unit reconciliation)**: long-standing v1/Fortran
   numerical inconsistency. v3 partially restores the missing DOC -> DIC
   term but does not fix the units. Defer to v3 1.x carbonate solver.
8. **C11 (use_SedFlux=False fallback for DIC sediment release)**:
   Phase 5.A scope; document and defer to Phase 5.5.
9. **C12 (use_SedFlux branch for SOD)**: documented as Phase 5.5; defer.
10. **C13 (SOD_20=1.0 vs Fortran 0.2)**: small kinetic-rate calibration
    issue; record in the corrections doc and consider aligning to 0.2 in
    a future release.
11. **C14 (BWa=3500 vs Fortran 5000)**: not used in carbon/dox; verify
    at the BenthicAlgae interface and either align or document.

### Items to escalate

12. **C6 (salinity correction on O2sat)**: significance depends on
    intended deployment scope; if v3 1.x targets brackish or estuarine
    applications, this must be addressed.
13. **C9 (DIC unit consistency)**: open question for the LimnoTech
    review team. The v1 model has an established calibration practice
    around the mol/L vs mg/L mismatch in the DIC budget; v3 should not
    silently change it without team input.

## Cross-references

- `src/clearwater_modules_v3/parameter_defaults_corrections.md` Section 1
  (records the SOD_20, SOD_theta, kah_20_user, kaw_20_user, pressure_mb
  corrections; does not yet record AWc/BWc/AWa/BWd handling).
- `tests/test_5_carbon_calculations_v2.py` line 348-349 docstring
  explicitly states "rca = AWc / AWa = AWc in v3's per-Chla convention",
  but the v3 default `AWc = 40, AWa = 1000` does not satisfy this
  identity. The test passes because it uses the same wrong AWc value on
  both sides.
- `tests/test_5_dox_calculations_v2.py` line 328-330 docstring documents
  the same AWc-as-rca convention; same review applies.
- v3 `parameters/carbon.py` (CARBON_DEFAULTS), `parameters/dox.py`
  (DOX_DEFAULTS), `utils/sediment.py` (Phase 1.1 pure-Arrhenius
  decision), `utils/reaeration.py` (kah_20, kaw_20, ka_tc).
- Phase 5.A and 5.B agent reports in conversation history.
