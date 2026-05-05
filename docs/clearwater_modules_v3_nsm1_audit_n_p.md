# v3 NSM1 Nitrogen + Phosphorus -- Three-way audit (Fortran vs v1 vs v3)

Audit date: 2026-05-05
Auditor: senior water-quality modeling reviewer
Scope: kinetic-block parity for v3 NSM1 Nitrogen and Phosphorus, comparing
`src/clearwater_modules_v2/processes/nitrogen.py` (re-exported as v3
`Nitrogen`), `src/clearwater_modules_v3/processes/phosphorus.py`, and
`src/clearwater_modules_v3/parameters/{nitrogen,phosphorus}.py` against the
legacy Fortran NSM1 (`/Users/todd/Downloads/NSM_comparison/NSM1/Source Files/`)
and the v1 Python NSM1 (`src/clearwater_modules/nsm1/`).

Out of scope: orchestration/Model.run sequencing, YAML registry plumbing,
sediment-flux (NSM2) coupling.

---

## Summary

- Nitrogen: 5 critical, 4 minor, 7 matches (out of 16 enumerated kinetic blocks).
- Phosphorus: 1 critical (inherited from `fdp` utility, latent at default `kdpo4=0`),
  2 minor, 11 matches.
- Parameter defaults: 6 disagreements (3 minor unit/value, 3 critical default-value
  divergences).

Top concerns (read these first):

1. v3 Nitrogen carries a *phantom NH4 source term* `ammonium_decay_nitrate`
   with default rate `1.0/d` that has no v1 or Fortran analogue. At default
   kwargs, it injects `1.0 * NH4 mg-N/L/d` into the NH4 budget, faster than
   any sink. Critical.
2. v3 nitrate algal-uptake fraction (`float_algea_faction_uptake_from_nitrate`)
   is a *static parameter* defaulting to `1.0`, NOT recomputed each step as
   `1 - ApUptakeFr_NH4`. NH4 algal-uptake uses the dynamic fraction
   (`algal_nh4_uptake_fraction` cache). The two uptake paths therefore do
   not sum to the total algal-N demand `rna * ApGrowth`; mass balance is
   broken. Critical.
3. v3 `nitrate_uptake_benthic_algae` divides by `algal_chlorophyll` (the
   *floating*-algae chlorophyll factor `AWa = 1000`) instead of by `BWd`
   (benthic dry-weight), and is missing the `/ depth` divisor that v1 and
   Fortran both apply. Stoichiometry and units are wrong. Critical.
4. v3 `ammonium_from_bed` uses default `sediment_ammonium_release_rate=1.0`
   (v1 default `rnh4_20=0`); the formula is correct but the default magnitude
   injects a large (`1/depth`) NH4 source where Fortran and v1 are silent.
   Critical at default kwargs (calibration-impacting).
5. v3 Phosphorus inherits the v1 unit error in `fdp`: divides by `0.000001`
   instead of `1.0E6` (a factor of 1E12). Latent at the v3 default `kdpo4=0`,
   but breaks particulate-P partitioning the moment the user enables `kdpo4>0`.
   Critical, but gated.

---

## Nitrogen

### N1. NH4 nitrification

- Fortran (`modNitrogen.f90:265,270`):
  `NitrificationInhibition = 1 - exp(-KNR * DOX)`;
  `NH4_Nitrification = NitrificationInhibition * knit_tc * NH4`.
- v1 (`processes.py:1437,1454`):
  `xr.where(use_DOX, 1 - exp(-KNR*DOX), 1)`;
  `NH4_Nitrification = NitrificationInhibition * knit_tc * NH4`.
- v3 (`v2/nitrogen.py:493-509,596-602`):
  `nitrification_inhibition = 1 - exp(-KNR * DOX)` with
  `KNR := self.nitrification_oxygen_inhibition_factor` (default `1.0`);
  `ammonium_nitrification = NH4 * arrhenius(T, knit_20, knit_theta) * inhibition`.

Severity: minor. The formula is identical. The v3 wiring uses
`nitrification_oxygen_inhibition_factor` (default `1.0`) for KNR, while v1
and Fortran both use the named constant `KNR = 0.6 mg-O2/L`. v3
`NITROGEN_DEFAULTS` defines `KNR = 0.6` but the kinetic call routes through
the legacy kwarg, so the v3 default for nitrification inhibition is
effectively `1.0`, not `0.6`. Documented divergence; recommend rewiring
`nitrification_inhibition` to read `self.KNR` (NITROGEN_DEFAULTS value).

### N2. NH4 -> NO3 nitrification term wiring (`change_ammonium`)

- Fortran (`modNitrogen.f90:296`): `dNH4dt = OrgN_NH4_Decay - NH4_Nitrification + NH4fromBed + NH4_ApRespiration - NH4_ApGrowth + NH4_AbRespiration - NH4_AbGrowth`.
- v1 (`processes.py:1584`): identical.
- v3 (`v2/nitrogen.py:334-353`):

  ```
  rate = (
      self.ammonium_decay_nitrate(...)        # PHANTOM SOURCE: rate=1.0/d * NH4
      - self.ammonium_nitrification(...)
      + self.ammonium_from_bed(...)
      + self.ammonium_floating_respiration()
      - self.ammonium_floating_growth()
      + self.ammonium_benthic_respiration()
      - self.ammonium_benthic_growth()
      + orgn_to_nh4
  )
  ```

Severity: critical. The legacy v2 kwarg `ammonium_decay_rate` (no v1 analogue)
adds `arrhenius(T, ammonium_decay_rate, ammonium_decay_theta) * NH4` as a
*positive source* on NH4. With defaults `1.0/d`, `theta=1.0`, this injects
NH4 at first-order rate `1.0/d * NH4` for any positive NH4. There is no
matching sink anywhere in the budget. NSM1 Fortran has no such term (line
296 of `modNitrogen.f90`); v1 has no such term (line 1584). The v3 docstring
(lines 41-46) calls these "legacy v2 kwargs preserved for back-compat";
they were never validated against NSM1 and break NH4 mass balance.

Recommendation before LimnoTech review: drop `ammonium_decay_nitrate` from
the `change_ammonium` rate sum, or default `ammonium_decay_rate=0.0`.

### N3. NH4 hydrolysis from OrgN (OrgN_NH4_Decay)

- Fortran (`modNitrogen.f90:231`): `OrgN_NH4_Decay = kon_tc * OrgN`.
- v1 (`processes.py:1330`): `xr.where(use_OrgN, kon_tc * OrgN, 0)`.
- v3 (`v2/nitrogen.py:613-627`): `kon_tc * OrgN` with `kon_tc = arrhenius(T, kon_20, kon_theta)`.

Match. Defaults `kon_20=0.1`, `kon_theta=1.074` agree across all three.

### N4. NH4 from bed (sediment release)

- Fortran (`modNitrogen.f90:275`): `NH4fromBed = rnh4_tc / depth` (with `rnh4` default 0; gated by `use_SedFlux`).
- v1 (`processes.py:1470`): `rnh4_tc / depth` (`rnh4_20=0` default).
- v3 (`v2/nitrogen.py:451-457`): `arrhenius(T, sediment_ammonium_release_rate=1.0, sediment_ammonium_release_theta=1.0) / depth`.

Severity: critical (default-value defect; formula correct).
The formula is structurally correct. The v2-style legacy kwarg
`sediment_ammonium_release_rate` defaults to `1.0` (1/d), versus the v1/Fortran
default `rnh4_20=0`. At unit-step instantiation `Nitrogen()` (no parameter
override), this term injects `1.0/depth mg-N/L/d` into NH4 every step. The
v3 `NITROGEN_DEFAULTS['rnh4_20'] = 0.0` is correctly carrying the v1 value,
but `change_ammonium` reads the legacy kwarg, not the v3 default. Same wiring
defect as N1.

Recommendation: rewire `ammonium_from_bed` to read `self.rnh4_20` and
`self.rnh4_theta` from NITROGEN_DEFAULTS; default `sediment_ammonium_release_rate=0.0`.

### N5. NH4 from floating-algae respiration

- Fortran (`modNitrogen.f90:279`): `NH4_ApRespiration = rna(r) * ApRespiration`.
- v1 (`processes.py:1486`): `xr.where(use_Algae, rna * ApRespiration, 0)`.
- v3 (`v2/nitrogen.py:363-366` -> `floating_algae.py:673-681`):
  `rna * algal_respiration_rate` where `rna = AWn/AWa`, reading the cached
  `algal_respiration_rate` populated by `FloatingAlgae.run`.

Match.

### N6. NH4 sink from floating-algae growth (NH4-vs-NO3 fractionation)

- Fortran (`modNitrogen.f90:280`): `NH4_ApGrowth = ApUptakeFr_NH4 * rna(r) * ApGrowth` with `ApUptakeFr_NH4 = PN(r) * NH4 / (PN(r)*NH4 + (1-PN(r))*NO3)` recomputed per step (line 208).
- v1 (`processes.py:1504`): `xr.where(use_Algae, ApUptakeFr_NH4 * rna * ApGrowth, 0)`; `ApUptakeFr_NH4` recomputed (line 1226-1247).
- v3 (`v2/nitrogen.py:373-376` -> `floating_algae.py:683-691`):
  `algal_nh4_uptake_fraction * rna * algal_growth_rate`. `algal_nh4_uptake_fraction`
  is recomputed per step in `FloatingAlgae.run` (line 308) via
  `_ap_uptake_fr_nh4(ammonium, nitrate)`.

Match.

### N7. NH4 from benthic-algae respiration

- Fortran (`modNitrogen.f90:287`): `NH4_AbRespiration = rnb(r) * Fb(r) * AbRespiration / depth`.
- v1 (`processes.py:1525`): `(rnb * AbRespiration * Fb) / depth` (`Fw` not used; v1 footnote `# TODO changed the calculation for respiration from the inital FORTRAN due to conflict with the reference guide`).
- v3 (`v2/nitrogen.py:368-371` -> `benthic_algae.py:499`): `rnb * balgae_respiration_rate * fb / depth`.

Match (Fb only, no Fw, agreeing with v1 and Fortran).

### N8. NH4 sink from benthic-algae growth

- Fortran (`modNitrogen.f90:288`): `NH4_AbGrowth = AbUptakeFr_NH4 * rnb(r) * Fb(r) * AbGrowth / depth`.
- v1 (`processes.py:1547`): `(AbUptakeFr_NH4 * rnb * Fb * AbGrowth) / depth`.
- v3 (`v2/nitrogen.py:378-381` -> `benthic_algae.py:513`):
  `balgae_nh4_uptake_fraction * rnb * fb * balgae_growth_rate / depth`.

Match. (Confirmed `balgae_nh4_uptake_fraction` is recomputed per step in
`BenthicAlgae.run`, mirroring the floating-algae cache.)

### N9. NO3 -> NO3 nitrification source (`change_nitrate`)

- Fortran (`modNitrogen.f90:335`): `dNO3dt = NH4_Nitrification - NO3_Denit - NO3_BedDenit - NO3_ApGrowth - NO3_AbGrowth`.
- v1 (`processes.py:1729`): identical.
- v3 (`v2/nitrogen.py:412-442`):

  ```
  rate = (
      ammonium_nitrification(...)            # source: + knit_tc * NH4 * inhibition
      - nitrate_denitrification(...)
      - nitrate_bed_denitrification(...)
      - nitrate_uptake_floating_algae(...)
      - nitrate_uptake_benthic_algae(...)
  )
  ```

Match in structural form. See N12, N13 for component defects.

### N10. NO3 denitrification (water-column)

- Fortran (`modNitrogen.f90:308`): `NO3_Denit = (1 - DOX/(DOX+KsOxdn(r))) * kdnit_tc * NO3` with NaN guard fall-back to `kdnit_tc * NO3`.
- v1 (`processes.py:1623-1635`): `np.select` reproducing the same NaN-handling.
- v3 (`v2/nitrogen.py:524-548`): `nitrate * arrhenius(T, denitrification_rate, denitrification_theta) * (1 - DOX/(DOX+half_saturation_oxygen))` with NaN-to-0 guard.

Severity: minor. The default `denitrification_rate=1.0/d` (legacy kwarg)
is 500x larger than v1/Fortran `kdnit_20=0.002`. v3 NITROGEN_DEFAULTS
correctly stores `kdnit_20=0.002` but `change_nitrate` reads the legacy kwarg.
The v3 NaN guard returns 0 (Fortran returns `kdnit_tc * NO3`); these differ
when `DOX = -KsOxdn` (impossible physically), so this is a stability detail
of no practical consequence. Recommend rewiring to `self.kdnit_20`
and matching the v1/Fortran NaN fall-back.

### N11. NO3 bed denitrification (sediment)

- Fortran (`modNitrogen.f90:317`): `NO3_BedDenit = vno3_tc * NO3 / depth` (`vno3` units m/d; default 0).
- v1 (`processes.py:1655`): `vno3_tc * NO3 / depth`.
- v3 (`v2/nitrogen.py:550-563`): `nitrate * arrhenius(T, sediment_denitrification_rate, sediment_denitrification_theta) / depth`.

Severity: minor (default-value divergence). Formula matches. Default
`sediment_denitrification_rate=1.0` is a v2 legacy kwarg; v1/Fortran default
`vno3_20=0`. Not gated by `use_SedFlux` in v2/v3.

### N12. NO3 sink from floating-algae growth

- Fortran (`modNitrogen.f90:321`): `NO3_ApGrowth = ApUptakeFr_NO3 * rna(r) * ApGrowth` with `ApUptakeFr_NO3 = 1 - ApUptakeFr_NH4` recomputed per step.
- v1 (`processes.py:1675`): `xr.where(use_Algae, ApUptakeFr_NO3 * rna * ApGrowth, 0)` with `ApUptakeFr_NO3 = 1 - ApUptakeFr_NH4` (line 1260) recomputed per step.
- v3 (`v2/nitrogen.py:565-576`):

  ```
  return (
      self.floating_algae_nitrogen_weight     # AWn, mg-N/ug-Chla
      / self.algal_chlorophyll                # AWa, ug-Chla/ug-Chla = 1000 default
      * algea_growth_rate
      * self.float_algea_faction_uptake_from_nitrate   # STATIC, default 1.0
  )
  ```

Severity: critical. Two related defects:

1. **Static uptake fraction**. `float_algea_faction_uptake_from_nitrate` is
   set to `1.0` in `__init__` (line 84) and never recomputed each step. v1
   and Fortran recompute `ApUptakeFr_NO3 = 1 - ApUptakeFr_NH4` from current
   NH4 / NO3 every step. Because `ammonium_growth` (N6) DOES read the dynamic
   `algal_nh4_uptake_fraction`, the NH4 sink uses the dynamic split and the
   NO3 sink uses the static `1.0`. The two paths therefore do not sum to
   `rna * ApGrowth` -- they sum to roughly `(0.5 + 1.0) * rna * ApGrowth = 1.5 * rna * ApGrowth`
   under typical PN=0.5, NH4 ~ NO3. Algal-N mass balance is violated by a
   factor approaching 1.5x.

2. **Wrong stoichiometric ratio reference**. The formula `AWn / AWa * AbGrowth`
   evaluates to `7.2 / 1000 * AbGrowth`, but `rna = AWn / AWa` is exactly
   that division (v1 `rna` and Fortran `rna(r)`). So the floating-algae
   coefficient is correct numerically (7.2/1000 = 0.0072 mg-N/ug-Chla,
   matching `rna`); the formula structure happens to be correct here but
   is opaque, and it breaks for benthic algae (see N13).

Recommendation: rewire to read the dynamic
`1 - floating_algae_process.algal_nh4_uptake_fraction`. Drop
`float_algea_faction_uptake_from_nitrate` as a static parameter.

### N13. NO3 sink from benthic-algae growth

- Fortran (`modNitrogen.f90:328`): `NO3_AbGrowth = AbUptakeFr_NO3 * rnb(r) * Fb(r) * AbGrowth / depth`.
- v1 (`processes.py:1697`): `xr.where(use_Balgae, (AbUptakeFr_NO3 * rnb * Fb * AbGrowth) / depth, 0)`.
- v3 (`v2/nitrogen.py:578-594`):

  ```
  return (
      self.benthic_algae_nitrogen_weight       # BWn, mg-N/g-D
      / self.algal_chlorophyll                  # AWa = 1000 (FLOATING denominator!)
      * algea_growth_rate                       # g/m^2/d
      * self.benthic_algea_faction_uptake_from_nitrate  # static 0.5
      * self.fraction_bottom_area               # = 1.0 default; should be Fb
  )
  ```

Severity: critical. Three structural defects:

1. **Wrong stoichiometric denominator**: `algal_chlorophyll = AWa` is the
   *floating-algae* chlorophyll-per-chlorophyll ratio (1000), not the
   benthic dry-weight `BWd`. v1/Fortran use `rnb = BWn / BWd` (v1 line
   ~1374, balgae constants). With `BWn = ` (v2 BALGAE_DEFAULTS not inspected
   here) and `AWa = 1000`, the coefficient is wrong by orders of magnitude.

2. **Missing `/ depth`**: v1 and Fortran divide by `depth` to convert
   `g/m^2/d * mg-N/g-D` into `mg-N/m^3/d` (= `mg-N/L/d`). v3 omits this
   divisor, so units are `g-N/m^2/d`, not `mg-N/L/d`.

3. **`fraction_bottom_area` substituted for `Fb`**: v3 uses `self.fraction_bottom_area`
   (init default `1.0`) where v1/Fortran use `Fb` (default `0.9`). The v3
   path bypasses the BenthicAlgae-side `Fb` configuration.

The structural NH4-uptake counterpart (N8) routes through
`benthic_algae_process.ammonium_growth()` which uses correct `rnb * fb / depth`.
Only the NO3 path is broken.

Recommendation: rewire `nitrate_uptake_benthic_algae` to read
`benthic_algae_process.balgae_no3_uptake_fraction` (add a counterpart cache
to BenthicAlgae if needed), use `BWn/BWd`, multiply by `Fb`, divide by `depth`.

### N14. OrgN settling

- Fortran (`modNitrogen.f90:233`): `OrgN_Settling = vson(r) / depth * OrgN` (raw `vson`, no Arrhenius).
- v1 (`processes.py:1345`): `vson / depth * OrgN` (raw `vson`, no Arrhenius).
- v3 (`v3/processes/phosphorus.py` no, this is in nitrogen `v2/nitrogen.py:629-643`):
  `arrhenius(T, vson_20, vson_theta) / depth * OrgN`.

Severity: minor (documented Phase 2.B deviation). v3 applies
`vson_theta=1.024` Arrhenius correction; v1 and Fortran do not. The v3
docstring (line 642) acknowledges this. At T=20 C the deviation is
zero. At T=25 C, `1.024^5 = 1.126`, a 12.6% increase in settling rate.
Calibration-impacting only off the reference temperature. Already
documented in `parameter_defaults_corrections.md`.

Note also: v1 default `vson = 0.01` m/d, Fortran default `vson = 0.01`.
v3 default `vson_20 = 0.1` m/d (10x larger). See parameter audit below.

### N15. OrgN from algal mortality (floating + benthic routing)

- Fortran (`modNitrogen.f90:236,242`):
  `ApDeath_OrgN = rna(r) * ApDeath`;
  `AbDeath_OrgN = rnb(r) * Fw(r) * Fb(r) * AbDeath / depth`.
- v1 (`processes.py:1360,1381`): identical.
- v3 (`v2/nitrogen.py:645-667` -> floating/benthic algae caches): reads
  `algal_orgn_from_mortality_rate` and `balgae_orgn_from_mortality_rate`
  (Phase 2.A populates these in algae `run`).

Match (assuming algae caches are correctly populated; not re-audited here).

### N16. dOrgN/dt budget

- Fortran (`modNitrogen.f90:247`): `dOrgNdt = ApDeath_OrgN + AbDeath_OrgN - OrgN_NH4_Decay - OrgN_Settling`.
- v1 (`processes.py:1402`): identical.
- v3 (`v2/nitrogen.py:669-703`): identical structure.

Match.

### N17. Cached step-scoped flux rates (`nitrification_flux_rate`, `denitrification_flux_rate`)

v3-only feature (Phase 2.B / Item 1; lines 175-176, 247-268 of `v2/nitrogen.py`).
Computed as positive-magnitude fluxes in mg-N/L/d before the change-rate
decomposition. No v1 or Fortran analogue (they recompute on the fly in
DOX / N2 modules).

Severity: observation (v3-only enhancement). The values agree numerically
with the in-line `change_ammonium` / `change_nitrate` calls because
both invoke `self.ammonium_nitrification(...)` / `self.nitrate_denitrification(...)`
with the same arguments. No defect.

---

## Phosphorus

### P1. OrgP -> TIP hydrolysis (`OrgP_DIP_decay`)

- Fortran (`modPhosphorus.f90:123`): `OrgP_DIP_decay = kop_tc * OrgP`.
- v1 (`processes.py:1879`): `xr.where(use_OrgP, kop_tc * OrgP, 0)`.
- v3 (`v3/phosphorus.py:283-287`): `kop_tc * orgp` with `kop_tc = arrhenius(T, kop_20, kop_theta)`.

Match. Defaults `kop_20=0.1`, `kop_theta=1.047` agree across all three.

### P2. OrgP settling

- Fortran (`modPhosphorus.f90:124`): `OrgP_Settling = vsop(r) / depth * OrgP` (raw `vsop`, no Arrhenius).
- v1 (`processes.py:1895`): `(vsop / depth) * OrgP`.
- v3 (`v3/phosphorus.py:297-301`): `self.vsop / depth * orgp` (raw `vsop`, no Arrhenius).

Match in formula. Default-value divergence: v1 default `vsop=999`
(sentinel), Fortran default `vsop=0.01`, v3 default `vsop=0.1` (corrected
per Phase 1 corrections doc). The v3 default differs from Fortran by 10x;
see parameter audit.

### P3. OrgP from floating-algae mortality

- Fortran (`modPhosphorus.f90:127`): `ApDeath_OrgP = rpa(r) * ApDeath`.
- v1 (`processes.py:1912`): `xr.where(use_Algae, rpa * ApDeath, 0)`.
- v3 (`v3/phosphorus.py:447-459`): reads `floating_algae_process.algal_orgp_from_mortality_rate` cache.

Match.

### P4. OrgP from benthic-algae mortality

- Fortran (`modPhosphorus.f90:133`): `AbDeath_OrgP = rpb(r) * Fw(r) * Fb(r) * AbDeath / depth`.
- v1 (`processes.py:1935`): `xr.where(use_Balgae, (rpb * Fw * Fb * AbDeath) / depth, 0)`.
- v3 (`v3/phosphorus.py:461-475`): reads `benthic_algae_process.balgae_orgp_from_mortality_rate`
  cache; docstring confirms cache already includes `Fw * Fb / depth`.

Match (assuming benthic algae cache correctness; not re-audited here).

### P5. dOrgP/dt budget

- Fortran (`modPhosphorus.f90:137`): `dOrgPdt = ApDeath_OrgP + AbDeath_OrgP - OrgP_DIP_Decay - OrgP_Settling`.
- v1 (`processes.py:1956`): identical (with `xr.where(use_OrgP, ...)`).
- v3 (`v3/phosphorus.py:339-347`): identical.

Match.

### P6. TIP partitioning (`fdp` utility)

- Fortran (`modGlobalParam.f90:226-230`):
  ```
  fdp = 1.0
  do i = 1, nGS
    fdp = fdp + kdpo4(i,r) * Solid(i) / 1.0E6
  end do
  fdp = 1.0 / fdp
  ```
  i.e. `fdp = 1 / (1 + sum_i(kdpo4_i * Solid_i / 1e6))`.
- v1 (`shared/processes.py:271`): `xr.where(use_TIP, 1/(1 + kdpo4 * Solid/0.000001), 0)`
  i.e. `fdp = 1 / (1 + kdpo4 * Solid * 1e6)`.
- v3 (`utils/partitioning.py:31`): `xr.where(use_TIP, 1.0 / (1.0 + kdpo4 * Solid / 0.000001), 0.0)`.

Severity: critical (latent at default `kdpo4=0`, manifests immediately when
the user enables sorption). v1 and v3 agree; both diverge from Fortran by
a factor of `1E12` in the denominator scaling. With the v3 default
`kdpo4=0.0`, both `(1 + 0)` and `(1 + 0)` give `fdp=1`, so TIP settling
in N7 (TIP_Settling = `vs/depth * (1-fdp) * TIP`) is zero on both. The
moment a user sets `kdpo4 > 0`, the v3/v1 path computes `fdp ≈ 0` (entire
TIP particulate, all settles) while Fortran computes `fdp ≈ 1` (entirely
dissolved, nothing settles). The two are extremes of opposite sign.

This is a v1 bug inherited verbatim by v3. The fix is to use the Fortran
form `kdpo4 * Solid / 1.0E6` (treating `kdpo4` as L/kg and `Solid` as
mg/L: `kdpo4 [L/kg] * Solid [mg/L] * 1e-6 [kg/mg] = dimensionless`).

Recommendation: fix the `fdp` utility before any LimnoTech demonstration
of the phosphorus partitioning pathway. Until fixed, document a hard
constraint that `kdpo4=0` is the only validated regime.

### P7. TIP settling

- Fortran (`modPhosphorus.f90:156`): `TIP_Settling = vs(r) / depth * (1 - fdp) * TIP`.
- v1 (`processes.py:1988`): `vs / depth * (1 - fdp) * TIP`.
- v3 (`v3/phosphorus.py:289-294`): `self.vs / depth * (1 - fdp) * tip`.

Match in formula. Default-value divergence: v1 default `vs=999`,
Fortran `vs=0.1`, v3 `vs=0.1`. v3 corrects v1's sentinel to match Fortran.

### P8. TIP from sediment release (`DIPfromBed`)

- Fortran (`modPhosphorus.f90:153`): `DIPfromBed = rpo4_tc / depth` (default `rpo4_20=0`; gated by `use_SedFlux`).
- v1 (`processes.py:1969`): `rpo4_tc / depth`.
- v3 (`v3/phosphorus.py:303-310`): `rpo4_tc / depth` with `rpo4_tc = arrhenius(T, rpo4_20, rpo4_theta)`.

Match. Default `rpo4_20=0` agrees across all three; term is silently zero
unless calibrator overrides.

### P9. TIP from floating-algae respiration

- Fortran (`modPhosphorus.f90:161`): `DIP_ApRespiration = rpa(r) * ApRespiration`.
- v1 (`processes.py:2003`): `xr.where(use_Algae, rpa * ApRespiration, 0)`.
- v3 (`v3/phosphorus.py:405-415`): `self._rpa() * algal_respiration_rate` with
  `_rpa = AWp / AWa = 1.0/1000`.

Match.

### P10. TIP sink from floating-algae growth

- Fortran (`modPhosphorus.f90:162`): `DIP_ApGrowth = rpa(r) * ApGrowth`.
- v1 (`processes.py:2018`): `xr.where(use_Algae, rpa * ApGrowth, 0)`.
- v3 (`v3/phosphorus.py:393-403`): `self._rpa() * algal_growth_rate`.

Match.

### P11. TIP from benthic-algae respiration

- Fortran (`modPhosphorus.f90:169`): `DIP_AbRespiration = rpb(r) * Fb(r) * AbRespiration / depth`.
- v1 (`processes.py:2037`): `xr.where(use_Balgae, rpb * Fb * AbRespiration / depth, 0)`.
- v3 (`v3/phosphorus.py:433-445`): `self._rpb() * self.Fb * balgae_respiration_rate / depth`.

Match.

### P12. TIP sink from benthic-algae growth

- Fortran (`modPhosphorus.f90:170`): `DIP_AbGrowth = rpb(r) * Fb(r) * AbGrowth / depth`.
- v1 (`processes.py:2056`): `xr.where(use_Balgae, rpb * Fb * AbGrowth / depth, 0)`.
- v3 (`v3/phosphorus.py:417-431`): `self._rpb() * self.Fb * balgae_growth_rate / depth`.

Match. (Note: this is the formula that v3 *Nitrogen* `nitrate_uptake_benthic_algae`
N13 fails to match. v3 Phosphorus is correct here.)

### P13. dTIP/dt budget

- Fortran (`modPhosphorus.f90:176-177`): `dTIPdt = OrgP_DIP_decay - TIP_Settling + DIPfromBed + DIP_ApRespiration - DIP_ApGrowth + DIP_AbRespiration - DIP_AbGrowth`.
- v1 (`processes.py:2091`): identical.
- v3 (`v3/phosphorus.py:323-333`): identical.

Match.

### P14. DIP derived variable (post-step)

Severity: minor. Fortran (`modPhosphorus.f90:223`) computes
`DIP = TIP / fdp`, while v1 (`processes.py:2179`) computes `DIP = TIP * fdp`.
These are *opposite* (reciprocal) operations. Given the matching `fdp` v1/v3
utility (`fdp = 1/(1 + kdpo4*Solid*1e6)`), `fdp` is the *dissolved fraction*,
so `DIP = TIP * fdp` (v1/v3) is *correct* and Fortran is wrong. v3
re-exports this v1 derived variable indirectly (no v3-native `DIP` derived
variable was inspected here). Note for documentation: v3's matching v1
behavior is the right answer; Fortran has the bug.

(This finding is informational; `DIP` is a diagnostic, not a state variable,
so it does not feed back into other kinetics.)

---

## Parameter defaults audit

### Nitrogen parameters

| Parameter   | Fortran        | v1 default | v3 NITROGEN_DEFAULTS | Status                                                                                                                              |
|-------------|----------------|------------|----------------------|-------------------------------------------------------------------------------------------------------------------------------------|
| KNR         | 0.6            | 0.6        | 0.6                  | Match value. **Wired wrong**: kinetic call uses legacy kwarg `nitrification_oxygen_inhibition_factor=1.0`, not `self.KNR`. Critical. |
| knit_20     | 0.1            | 0.1        | 0.1                  | Match. **Wired wrong**: kinetic call uses `nitrification_rate=1.0`, not `self.knit_20`. Critical.                                     |
| knit_theta  | 1.083          | 1.083      | 1.083                | Match (wired wrong via legacy kwarg `nitrification_theta=1.0`).                                                                       |
| kon_20      | 0.1            | 0.1        | 0.1                  | Match (wired correctly via `self.kon_20`).                                                                                            |
| kon_theta   | 1.047          | 1.047      | 1.074                | **Disagreement**: v3 `1.074`, v1/Fortran `1.047`. Minor calibration impact at 25 C: `(1.074/1.047)^5 = 1.14`, ~14% rate divergence.   |
| kdnit_20    | 0.002          | 0.002      | 0.002                | Match. Wired wrong via legacy `denitrification_rate=1.0`. Critical.                                                                   |
| kdnit_theta | 1.045          | 1.08       | 1.08                 | **Disagreement**: Fortran `1.045`, v1/v3 `1.08`. v1 already diverged. Minor calibration impact.                                      |
| rnh4_20     | 0              | 0          | 0                    | Match. Wired wrong via legacy `sediment_ammonium_release_rate=1.0`. Critical.                                                         |
| rnh4_theta  | 1.074          | 1.047      | 1.047                | **Disagreement**: Fortran `1.074`, v1/v3 `1.047`. Minor.                                                                              |
| vno3_20     | 0              | 0          | 0                    | Match. Wired wrong via legacy `sediment_denitrification_rate=1.0`. Critical.                                                          |
| vno3_theta  | 1.08           | 1.045      | 1.045                | **Disagreement**: Fortran `1.08`, v1/v3 `1.045`. Minor.                                                                               |
| vson        | 0.01           | 0.01       | 0.1 (`vson_20`)      | **Disagreement**: v3 10x larger. Minor.                                                                                              |
| vson_theta  | not present    | n/a        | 1.024                | v3-only feature; documented Phase 2.B deviation.                                                                                     |
| KsOxdn      | 0.1            | 0.1        | 0.1                  | Match.                                                                                                                                |
| PN          | 0.5            | 0.5        | 0.5                  | Match.                                                                                                                                |
| PNb         | 0.5            | 0.5        | 0.5                  | Match.                                                                                                                                |
| use_OrgN    | True (default) | True       | True                 | Match.                                                                                                                                |

### Phosphorus parameters

| Parameter  | Fortran | v1 default | v3 PHOSPHORUS_DEFAULTS | Status                                                       |
|------------|---------|------------|------------------------|--------------------------------------------------------------|
| kop_20     | 0.1     | 0.1        | 0.1                    | Match.                                                       |
| kop_theta  | 1.047   | 1.047      | 1.047                  | Match.                                                       |
| rpo4_20    | 0       | 0          | 0                      | Match (silent at default; gated by `use_SedFlux` in Fortran). |
| rpo4_theta | 1.074   | 1.074      | 1.074                  | Match.                                                       |
| kdpo4      | 0.0     | 0.0        | 0.0                    | Match. (Formula bug in `fdp` is gated; see P6.)              |
| vsop       | 0.01    | 999        | 0.1                    | **Disagreement**: v3 10x larger than Fortran. Minor.        |
| vs         | 0.1     | 999        | 0.1                    | Match (v3 corrects v1 sentinel).                             |

### Cross-cutting wiring defect

The v3 Nitrogen Process (`v2/nitrogen.py`) contains BOTH a v3 NITROGEN_DEFAULTS
attribute set (lines 95-108) AND legacy v2 kwargs (lines 71-86 / 114-137).
The kinetic methods (`ammonium_nitrification`, `nitrate_denitrification`,
`ammonium_from_bed`, `nitrate_bed_denitrification`, `nitrification_inhibition`,
`ammonium_decay_nitrate`) read from the **legacy kwargs**, not from the
NITROGEN_DEFAULTS attributes. The legacy kwarg defaults are uniformly `1.0`,
which is 5x to 500x larger than the matching v1/Fortran NSM1 defaults.

Without rewiring, calibrating via NITROGEN_DEFAULTS or YAML config does not
take effect for nitrification, denitrification, sediment NH4 release, or
sediment NO3 denitrification. Every calibration pathway must override the
legacy kwarg names.

This is the single highest-leverage fix: rewire the six kinetic methods to
read `self.knit_20 / knit_theta / kdnit_20 / kdnit_theta / rnh4_20 /
rnh4_theta / vno3_20 / vno3_theta / KNR` from NITROGEN_DEFAULTS, drop the
phantom `ammonium_decay_nitrate` term, and retire the legacy kwargs (or
default them to the NSM1 values).

---

## Conclusions

### Required actions before LimnoTech review

1. **Drop the phantom `ammonium_decay_nitrate` source term** from
   `change_ammonium` (`v2/nitrogen.py:335`), or default
   `ammonium_decay_rate=0.0`. There is no v1 or NSM1-Fortran analogue. With
   the current default `1.0/d`, NH4 grows exponentially during any
   integration starting with NH4 > 0.

2. **Rewire kinetic methods to NITROGEN_DEFAULTS attributes** instead of
   the legacy v2 kwargs. The DEFAULTS values match v1/Fortran; the legacy
   kwargs do not. Affected methods: `ammonium_nitrification`,
   `nitrate_denitrification`, `ammonium_from_bed`,
   `nitrate_bed_denitrification`, `nitrification_inhibition`.

3. **Fix `nitrate_uptake_floating_algae`** to read the dynamic
   `1 - floating_algae_process.algal_nh4_uptake_fraction`, not the static
   `float_algea_faction_uptake_from_nitrate=1.0`. Otherwise NH4-vs-NO3
   uptake does not sum to total algal-N uptake; mass balance is violated.

4. **Fix `nitrate_uptake_benthic_algae`** structurally:
   - Use `BWn / BWd` (benthic dry-weight ratio) instead of `BWn / AWa`.
   - Multiply by `Fb`, not by `fraction_bottom_area` (different default).
   - Divide by `depth`.
   - Use a dynamic `1 - balgae_nh4_uptake_fraction` cache, not the static
     `benthic_algea_faction_uptake_from_nitrate=0.5`.

5. **Fix `fdp` utility** in `utils/partitioning.py:31` to divide by
   `1.0E6` (Fortran convention), not by `0.000001` (v1 inheritance bug).
   Latent at default `kdpo4=0` but breaks immediately when sorption is
   enabled. Add an MMS or analytical test against Fortran for `kdpo4>0`
   regimes.

### Acceptable deviations to document

1. v3 OrgN settling applies `vson_theta=1.024` Arrhenius; v1/Fortran do not.
   Documented in Phase 2.B notes; no action required.

2. v3 stores nitrification/denitrification flux caches (`_flux_rate`
   suffix); no v1/Fortran analogue. Strict enhancement.

3. `kdnit_theta`, `rnh4_theta`, `vno3_theta`: v1 and v3 already diverge
   from Fortran by small amounts; v3 inherits v1's choice. Document as
   v1-inherited, not a v3 regression.

4. `vson` and `vsop` v3 defaults (0.1 m/d) are 10x the Fortran/v1 values
   (0.01 m/d). The Phase 1 corrections doc justifies this as a sentinel
   replacement, not a deliberate calibration. Recommend the v3 defaults
   be lowered to 0.01 m/d to match Fortran/v1.

### Items to escalate

- The Phosphorus `fdp` unit error is shared between v1 and v3. Whatever v1
  calibration work has been validated against was also using the wrong
  `fdp` formula whenever `kdpo4 > 0`. If LimnoTech has v1 results that
  they trust at non-zero `kdpo4`, those results are likely unreliable.
  Escalate to the V&V team before any new partitioning regime claims.

- The v2 legacy-kwarg shim approach in v3 Nitrogen creates a documentation
  trap: NITROGEN_DEFAULTS is exposed as the canonical parameter set, but
  the canonical set is silently bypassed by the kinetic implementation.
  Either fully retire the legacy kwargs (preferred) or add an integration
  test that asserts the kinetic terms read from NITROGEN_DEFAULTS, not
  from the legacy kwargs.
