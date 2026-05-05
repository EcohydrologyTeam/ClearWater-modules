# v3 NSM1 Algae — Three-way audit (Fortran vs v1 vs v3)

Date: 2026-05-05
Scope: floating algae and benthic algae kinetics (constituent math only).
Sources audited:
- Fortran: `Source Files/modAlgae.f90`, `Source Files/modBenthicAlgae.f90`
- v1: `src/clearwater_modules/nsm1/processes.py`, `src/clearwater_modules/nsm1/constants.py`
- v3: `src/clearwater_modules_v2/processes/floating_algae.py`,
      `src/clearwater_modules_v2/processes/benthic_algae.py`,
      `src/clearwater_modules_v3/parameters/algae.py`,
      `src/clearwater_modules_v3/parameters/balgae.py`

## Summary

- Counts: 5 critical, 7 minor, 13 matches, 4 observations.
- Top concerns:
  1. v3 **DEFAULTS dict is silently shadowed** by legacy v2 kwargs (`growth_rate_max`,
     `death_rate`, `repiration_rate`, `light_attenuation_coefficient`,
     `density_michaelis_menton_constant`). With no kwargs supplied, the user gets
     death=respiration=settling=0 and (for benthic) Ksb=1.0 instead of 10.0.
     ALGAE_DEFAULTS keys `mu_max_20`, `kdp_20`, `krp_20`, `vsap`, `KL`, `KsN`,
     `KsP` are assigned to `self` but never read by the rate methods.
  2. **FloatingAlgae `limit_light` option 1 has a misplaced parenthesis**: the
     `np.log(...)` argument lacks the inner ratio. Numerator and denominator are
     split across the `*` and `/` operators rather than wrapped inside `np.log`.
     This breaks the half-saturation depth-averaged light formulation.
  3. **BenthicAlgae `limit_light` option 3 (Steele) has a sign error in the
     exponent**: v3 implements `(x) / exp(1-x)` which equals `x*exp(x-1)`,
     whereas v1 and Fortran implement `x*exp(1-x)`.
  4. **Floating-algae `growth_rate_option == 3` (harmonic mean) zero-guard
     uses the wrong condition** (`limit_phosphorus == 1.0` instead of
     `limit_phosphorus == 0.0`).
  5. **PAR vs solar-radiation coupling**: v3 passes `solar_radiation` directly
     into `limit_light` without the `Fr_PAR=0.47` scaling that converts
     shortwave to PAR; Fortran/v1 use `PAR` as the input to FL/FLb.

## FloatingAlgae

### F1. Maximum growth rate (mu_max_tc)
- Fortran (`modAlgae.f90:240`): `mu_max_tc = Arrhenius(mu_max%rc20, TwaterC)`,
  with `rc20=1.0`, `theta=1.047`.
- v1 (`processes.py:365-378`): `arrhenius_correction(TwaterC, mu_max_20, mu_max_theta)`,
  defaults `mu_max_20=1.0`, `mu_max_theta=1.047`.
- v3 (`floating_algae.py:419-423`):
  `arrhenius_correction(TwaterC, self.growth_rate_max, self.growth_rate_correction)`.
- **Critical**: `growth_rate_max` is the *kwarg* (default 1.0) and
  `growth_rate_correction` is the *kwarg* (default 1.0, not 1.047). The
  ALGAE_DEFAULTS keys `mu_max_20` and `mu_max_theta` are merged into `self`
  (line 175-177) but never consulted by `rate_growth`. Default-instantiated
  `FloatingAlgae()` therefore uses theta=1.0, not 1.047.
- Severity: critical
- Category: v3-deviation (wiring bug)

### F2. Respiration rate (krp_tc) and rate_respiration
- Fortran (`modAlgae.f90:241, 342`): `krp_tc = Arrhenius(krp%rc20=0.2, theta=1.047)`,
  `ApRespiration = krp_tc * Ap`.
- v1 (`processes.py:382-395, 621-632`): identical algebra; defaults `krp_20=0.2`,
  `krp_theta=1.047`.
- v3 (`floating_algae.py:470-482`):
  `arrhenius_correction(T, self.repiration_rate, self.repiration_rate_correction_factor)`,
  then `algae * corrected_rate`.
- **Critical**: `self.repiration_rate` is the kwarg (default 0.0). ALGAE_DEFAULTS
  `krp_20=0.2` is assigned but unread. With no kwargs, respiration is identically
  zero. The unit test (`test_5_floating_algae_calculations_v2.py:65-66`) supplies
  `repiration_rate=0.2` explicitly and therefore passes; the regression hides the
  default-instantiation bug.
- Severity: critical
- Category: v3-deviation (wiring bug)

### F3. Death rate (kdp_tc) and rate_death
- Fortran (`modAlgae.f90:242, 345`): `kdp_tc = Arrhenius(kdp%rc20=0.15, theta=1.047)`,
  `ApDeath = kdp_tc * Ap`.
- v1 (`processes.py:399-412, 636-646`): identical; defaults `kdp_20=0.15`,
  `kdp_theta=1.047`.
- v3 (`floating_algae.py:461-468`):
  `arrhenius_correction(T, self.death_rate, self.death_rate_correction_factor)`,
  then `algae * corrected_rate`.
- **Critical**: same wiring bug as F1, F2. ALGAE_DEFAULTS `kdp_20=0.15` and
  `kdp_theta=1.047` are unread.
- Severity: critical
- Category: v3-deviation (wiring bug)

### F4. Settling rate (ApSettling)
- Fortran (`modAlgae.f90:348`): `ApSettling = vsap/depth * Ap`, default `vsap=0.15`.
- v1 (`processes.py:650-662`): `vsap / depth * Ap`, default `vsap=0.15`.
- v3 (`floating_algae.py:484-488`): `algae / depth * self.settling_velocity`.
- **Critical**: `self.settling_velocity` is the kwarg (default 0.0). ALGAE_DEFAULTS
  `vsap=0.15` is assigned but unread. With no kwargs, ApSettling is identically zero.
- Severity: critical
- Category: v3-deviation (wiring bug)

### F5. Light limitation (FL) — option 1, half-saturation
- Fortran (`modAlgae.f90:267-269`):
  `FL = (1/KEXT) * log( (KL + PAR) / (KL + PAR*exp(-KEXT)) )` with `KEXT = lambda*depth`.
- v1 (`processes.py:447`): `(1/(L*depth)) * log( (KL+PAR) / (KL + PAR*exp(-(L*depth))) )`.
- v3 (`floating_algae.py:556-564`):
  ```
  raw_rate = ((1.0 / (L * depth))
              * np.log(self.light_limitation_constant + surface_light_intensity)
              / (self.light_limitation_constant
                 + surface_light_intensity * np.exp(-(L * depth))))
  ```
- **Critical**: by Python operator precedence (left-to-right `*` and `/`):
  `result = ((1/(L*d)) * log(KL+PAR)) / (KL + PAR*exp(-Ld))`, i.e. the
  `np.log` argument is only `(KL + PAR)` and the `(KL + PAR*exp(-Ld))`
  factor divides the entire expression rather than appearing inside the
  logarithm as the denominator of the ratio.
  The v1/Fortran form is `(1/(L*d)) * log( (KL+PAR) / (KL+PAR*exp(-Ld)) )`.
  These are mathematically distinct under all non-trivial inputs.
- Severity: critical
- Category: actual-bug

### F6. Light limitation (FL) — option 2, Smith
- Fortran (`modAlgae.f90:271-279`): with abs(KL)<1e-10 fallback `FL=1.0`,
  else `FL = (1/KEXT) * log( (PAR/KL + sqrt1) / (PAR*exp(-KEXT)/KL + sqrt2) )`,
  where `sqrt_i` are `(1 + (PAR.../KL)^2)^0.5`.
- v1 (`processes.py:448-449`): identical.
- v3 (`floating_algae.py:566-608`): same formula, with the abs(KL)<1e-10 short
  circuit returning 1 (matching Fortran but disagreeing with v1, which returns 1).
  Note: v1's `np.select` returns 1 in the small-KL Smith branch (line 448);
  v3 also returns 1. Match.
- Severity: match

### F7. Light limitation (FL) — option 3, Steele
- Fortran (`modAlgae.f90:281-287`): with abs(KL)<1e-10 returns `FL=0`, else
  `FL = (2.718/KEXT) * (exp(-PAR/KL * exp(-KEXT)) - exp(-PAR/KL))`.
- v1 (`processes.py:450-451`): identical.
- v3 (`floating_algae.py:610-627`):
  `(2.718 / (L*depth)) * (exp(-PAR/KL * exp(-Ld)) - exp(-PAR/KL))`.
- Match.
- Severity: match

### F8. FL guard for Ap<=0, KEXT<=0, PAR<=0
- Fortran (`modAlgae.f90:264-266`): unified guard before any computation.
- v1 (`processes.py:436-437`): unified guard via `np.select`.
- v3 (`floating_algae.py:633-636`): `algae <= 0 -> 0`; `L*depth <= 0 -> 0`; no PAR
  guard; clamp `raw_rate > 1` to 1.
- **Minor**: the `PAR <= 0` (after-sunset) guard is missing. If `surface_light_intensity`
  is zero and option 1 is selected (with the bug in F5 unfixed), `np.log(KL+0)`
  yields `log(KL)` and the result is finite but unphysical. With the F5 bug
  fixed, log(KL/KL)=0 and the result naturally reduces to 0; even so, the
  explicit guard is a defensive measure that v1/Fortran provide.
- Severity: minor
- Category: v3-deviation

### F9. PAR vs solar_radiation coupling
- Fortran/v1: `PAR` is supplied as input to FL (W/m^2 of photosynthetically active
  radiation). v1 user-facing path computes `PAR = q_solar * Fr_PAR` upstream
  (Fr_PAR=0.47 default; see `parameter_defaults_corrections.md` 3.4).
- v3 (`floating_algae.py:277, 393`): reads `solar_radiation` from registry and
  passes it as `surface_light_intensity` (the FL input slot).
- **Minor / observation**: there is no Fr_PAR scaling in the FloatingAlgae or
  BenthicAlgae path. Either `solar_radiation` is conventionally already PAR
  (in which case the registry contract should say so), or FL/FLb is being
  driven by total shortwave, overstating effective irradiance by ~2x. Needs
  verification against the registry-side variable contract.
- Severity: minor (needs verification)
- Category: v3-deviation

### F10. Nitrogen limitation (FN)
- Fortran (`modAlgae.f90:297-303`): `FN = (NH4+NO3)/(KsN+NH4+NO3)` when
  use_NH4 OR use_NO3; else 1. NaN→0; clamp to 1.
- v1 (`processes.py:474-527`): four explicit branches (both, NH4-only, NO3-only,
  neither) yielding `(N_active)/(KsN+N_active)`. NaN→0; clamp to 1.
- v3 (`floating_algae.py:520-544`):
  `n = nitrate if use_nitrate else 0; n += ammonium if use_ammonium else 0;
   rate = n/(KsN+n)`. NaN→0; clamp to 1.
- Algebra is equivalent in all four sub-cases.
- Severity: match

### F11. Phosphorus limitation (FP)
- Fortran (`modAlgae.f90:308-314`): `FP = fdp*TIP/(KsP+fdp*TIP)` when use_TIP;
  else 1. NaN→0; clamp to 1.
- v1 (`processes.py:530-561`): identical.
- v3 (`floating_algae.py:490-518`): identical formula. Note the gating flag is
  `self.use_phosphate` not `self.use_TIP`.
- **Minor**: the gating flag in `FloatingAlgae.limit_phosphorus` is
  `use_phosphate`; the `fdp` partitioning utility is gated on `self.use_TIP`
  (line 283-285). Both default to `True` so masked. Two flag names for the
  same physical decision is a maintainability hazard.
- Severity: minor
- Category: v3-deviation

### F12. Growth rate combination — option 1 (multiplicative)
- Fortran (`modAlgae.f90:321`): `mu = mu_max_tc * FL * FP * FN`.
- v1 (`processes.py:592`): identical.
- v3 (`floating_algae.py:427`): `growth_rate * limit_phosphorus * limit_nitrogen * limit_light`.
- Match (modulo the F1 wiring bug for `growth_rate`).
- Severity: match

### F13. Growth rate combination — option 2 (limiting nutrient, min)
- Fortran (`modAlgae.f90:326`): `mu = mu_max_tc * FL * min(FP, FN)`.
- v1 (`processes.py:593-594`): two-branch `where` equivalent to min.
- v3 (`floating_algae.py:430-434`):
  `where(FP > FN, growth*FN*FL, growth*FP*FL)`. Equivalent to min.
- Match.
- Severity: match

### F14. Growth rate combination — option 3 (harmonic mean)
- Fortran (`modAlgae.f90:328-336`): if FN==0 OR FP==0, `mu=0`; else
  `mu = mu_max_tc * FL * 2 / (1/FN + 1/FP)`.
- v1 (`processes.py:587-588, 595-596`): same algebra.
- v3 (`floating_algae.py:436-452`):
  ```
  rate_raw = growth * FL * 2 / (1/FN + 1/FP)
  rate = where(FN == 0, 0, rate_raw)
  rate = where(FP == 1, 0, rate)        # <-- wrong condition
  ```
- **Critical**: the second guard zeros out the rate when `FP == 1` (i.e. when
  phosphorus is fully non-limiting). This is the opposite of the intended
  guard, which should zero the rate when `FP == 0` to avoid division by zero
  in `1/FP`. As written, growth shuts down precisely when P is saturating,
  which is the most common steady-state P-replete condition. The TODO at
  line 449 acknowledges the constant is suspect.
- Severity: critical
- Category: actual-bug

### F15. Algal biomass integration (dApdt and Ap update)
- Fortran (`modAlgae.f90:339, 342, 345, 348, 353`):
  `dApdt = ApGrowth - ApRespiration - ApDeath - ApSettling` (per-day).
- v1 (`processes.py:666-697`): same `dApdt`; `Ap = Ap + dApdt * dt` (dt in days).
- v3 (`floating_algae.py:290-326`):
  `rate = growth - death - respiration - settling`; `algae_new = algae + rate * dt_days`
  with `dt_days = time_step.total_seconds() / 86400`. Then `clip_negative_state`
  applied. `set_at_time` persists the new state.
- Match. The Phase 2.A bug-#4/#16 fixes are present and correct.
- Severity: match

### F16. NH4-uptake fraction (ApUptakeFr_NH4)
- Fortran: not implemented in modAlgae (uptake fractionation is in NSM1's
  Nitrogen module).
- v1 (`processes.py:1206-1247`):
  - use_NH4 only → 1.0
  - use_NO3 only → 0.0
  - neither → 0.5
  - both → `PN*NH4 / (PN*NH4 + (1-PN)*NO3)`; NaN→PN.
- v3 (`floating_algae.py:643-671`): identical; uses `self.PN` (default 0.5) when
  the parameter exists, else falls back to 0.5. Both NH4/NO3 branch exactly
  matches v1's algebra; `denom > 0` guard with fallback to PN matches v1's
  NaN-fallback.
- Match.
- Severity: match

### F17. NH4 source from algal respiration (NH4_ApRespiration)
- Fortran: not implemented in modAlgae; surfaced via Nitrogen module.
- v1 (`processes.py:1472-1486`): `rna * ApRespiration` with use_Algae gating.
- v3 (`floating_algae.py:673-681`): `(AWn/AWa) * algal_respiration_rate`.
- Match. Bug #13 (was returning 0) is fixed.
- Severity: match

### F18. NH4 sink from algal growth (NH4_ApGrowth)
- v1 (`processes.py:1488-1504`): `ApUptakeFr_NH4 * rna * ApGrowth`.
- v3 (`floating_algae.py:683-691`): same algebra via cached
  `algal_growth_rate` and `algal_nh4_uptake_fraction`.
- Match. Bug #14 (was returning 0) is fixed.
- Severity: match

### F19. Algal mortality routing (OrgN, OrgP, POC, DOC)
- v1 (`processes.py:1347-1360, plus carbon module`):
  `ApDeath_OrgN = rna * ApDeath`, `ApDeath_OrgP = rpa * ApDeath`,
  `POC_algal_mortality = f_pocp * rca * ApDeath`,
  `DOC_algal_mortality = (1 - f_pocp) * rca * ApDeath`.
- v3 (`floating_algae.py:328-370`): identical algebra. `f_pocp` defaults to 0.5
  (see Observation O1 below); v1 default `f_pocp=0.9`.
- Severity: match (formula); see O1 for default-value note.

### F20. Algal POM source from settling
- v1: `POM_algal_settling = vsap * Ap * (AWd/AWa) / h2`.
- v3 (`floating_algae.py:368-370`): `vsap * algae * (AWd/AWa) / h2`. `h2` is
  composed from `POM_DEFAULTS`.
- Match.
- Severity: match

## BenthicAlgae

### B1. Maximum growth rate (mub_max_tc)
- Fortran (`modBenthicAlgae.f90:257`): `mub_max_tc = Arrhenius(mub_max%rc20=0.4, theta=1.047)`.
- v1 (`processes.py:701-713`): same; defaults `mub_max_20=0.4`, `mub_max_theta=1.047`.
- v3 (`benthic_algae.py:344-359, via FloatingAlgae kwargs`):
  `arrhenius_correction(T, self.growth_rate_max, self.growth_rate_correction)`.
- **Critical**: same wiring bug as F1. BALGAE_DEFAULTS `mub_max_20=0.4` and
  `mub_max_theta=1.047` are merged into `self` but the rate methods read
  `self.growth_rate_max` (kwarg, default 1.0) and `self.growth_rate_correction`
  (kwarg, default 1.0). Default-instantiated `BenthicAlgae()` uses the
  *floating-algae* growth-rate kwarg of 1.0/d rather than the benthic
  reference 0.4/d.
- Severity: critical
- Category: v3-deviation (wiring bug)

### B2. Respiration rate (krb_tc)
- Fortran (`modBenthicAlgae.f90:258`): `krb_tc = Arrhenius(krb%rc20=0.2, theta=1.06)`.
- v1 (`processes.py:717-729`): same; defaults `krb_20=0.2`, `krb_theta=1.06`.
- v3 (`benthic_algae.py:385-401`):
  `arrhenius_correction(T, self.repiration_rate, self.repiration_rate_correction_factor)`.
- **Critical**: same wiring bug. BALGAE_DEFAULTS `krb_20=0.2` and `krb_theta=1.06`
  are unread. Note: `krb_theta=1.06` differs from the typical 1.047; this is
  Fortran-correct but the wiring bug means it never reaches the rate
  computation under default instantiation.
- Severity: critical
- Category: v3-deviation (wiring bug)

### B3. Death rate (kdb_tc)
- Fortran (`modBenthicAlgae.f90:259`): `kdb_tc = Arrhenius(kdb%rc20=0.3, theta=1.047)`.
- v1 (`processes.py:733-745`): same; defaults `kdb_20=0.3`, `kdb_theta=1.047`.
- v3 (`benthic_algae.py:461-468 inherited from FloatingAlgae.rate_death`):
  `arrhenius_correction(T, self.death_rate, self.death_rate_correction_factor)`.
- **Critical**: same wiring bug. BALGAE_DEFAULTS `kdb_20=0.3` and `kdb_theta=1.047`
  are unread under default instantiation.
- Severity: critical
- Category: v3-deviation (wiring bug)

### B4. Light limitation (FLb) — option 1, half-saturation
- Fortran (`modBenthicAlgae.f90:269, 284`): `KEXT = exp(-lambda*depth)`;
  `FLb = PAR*KEXT / (KLb + PAR*KEXT)`.
- v1 (`processes.py:825, 838`): identical.
- v3 (`benthic_algae.py:411-422`):
  `light_at_depth_coefficent = exp(-L*depth);
   raw_rate = surface * coef / (KLb + surface * coef)`.
- Match.
- Severity: match

### B5. Light limitation (FLb) — option 2, Smith
- Fortran (`modBenthicAlgae.f90:287`): `FLb = PAR*KEXT / ((KLb^2 + (PAR*KEXT)^2)^0.5)`.
- v1 (`processes.py:839`): identical.
- v3 (`benthic_algae.py:424-435`): same algebra.
- Match.
- Severity: match

### B6. Light limitation (FLb) — option 3, Steele
- Fortran (`modBenthicAlgae.f90:293`):
  `FLb = PAR*KEXT/KLb * exp(1 - PAR*KEXT/KLb)`. With `x = PAR*KEXT/KLb`,
  this equals `x * exp(1 - x)`.
- v1 (`processes.py:841`):
  `PAR*KEXT/KLb * np.exp(1.0 - PAR*KEXT/KLb)`. Same as Fortran.
- v3 (`benthic_algae.py:438-454`):
  ```
  raw_rate = surface * coef / (KLb * exp(1 - surface*coef/KLb))
  ```
  By precedence, this is `(surface*coef) / (KLb * exp(1-x))` = `x / exp(1-x)`
  = `x * exp(x - 1)`.
- **Critical**: v3 implements `x * exp(x - 1)` instead of `x * exp(1 - x)`.
  At x=1 (light = KLb) both forms give 1; for x<1 v3 underestimates and for
  x>1 v3 overestimates the limitation factor. The error is the sign of the
  exponent argument.
- Severity: critical
- Category: actual-bug

### B7. Nitrogen limitation (FNb)
- Fortran (`modBenthicAlgae.f90:302-308`): same form as FN with KsNb instead
  of KsN; default `KsNb=0.25`.
- v1 (`processes.py:864-916`): four-branch form, same algebra.
- v3: BenthicAlgae inherits `limit_nitrogen` from FloatingAlgae.
- **Critical**: the inherited method uses `self.nitrogen_michaelis_menton_constant`
  (kwarg, default 0.04 = `KsN`) rather than `self.KsNb=0.25` from BALGAE_DEFAULTS.
  Under default instantiation, BenthicAlgae sees the *floating-algae*
  half-saturation constant of 0.04 mg-N/L instead of the benthic 0.25 mg-N/L,
  understating the half-saturation by ~6x and consequently overstating FNb
  at low N concentrations.
- Severity: critical
- Category: v3-deviation (wiring bug)

### B8. Phosphorus limitation (FPb)
- Fortran (`modBenthicAlgae.f90:311-317`): `FPb = fdp*TIP/(KsPb+fdp*TIP)`,
  default `KsPb=0.125`.
- v1 (`processes.py:919-950`): identical.
- v3: BenthicAlgae inherits `limit_phosphorus` from FloatingAlgae which uses
  `self.phosphorus_michaelis_menton_constant` (kwarg, default 0.0012 = `KsP`).
- **Critical**: same class of wiring bug as B7. BenthicAlgae uses the
  floating-algae KsP=0.0012 instead of the benthic KsPb=0.125, understating
  the half-saturation by ~100x.
- Severity: critical
- Category: v3-deviation (wiring bug)

### B9. Density limitation (FSb)
- Fortran (`modBenthicAlgae.f90:320`): `FSb = 1 - Ab/(Ab + KSb)`, default
  `KSb=10.0` g-D/m^2; clamp NaN→1, FSb>1→1.
- v1 (`processes.py:953-982`): identical.
- v3 (`benthic_algae.py:467-480`):
  `1 - algae/(algae + self.density_michaelis_menton_constant)`,
  with `density_michaelis_menton_constant` from kwarg (default 1.0).
- **Critical**: BALGAE_DEFAULTS provides `Ksb=10.0` but `rate_growth` reads
  the kwarg `density_michaelis_menton_constant` (default 1.0), 10x smaller.
  At default Ab values this overstates the density limitation. Note also
  the v1 NaN→0 vs v3 NaN→0 (both v1 line 974 and v3 line 478 use 0); v1
  comments suggest 1, but the implementation uses 0.
- Severity: critical
- Category: v3-deviation (wiring bug)

### B10. mub combination — option 1 (multiplicative)
- Fortran (`modBenthicAlgae.f90:330`): `mub = mub_max_tc * FLb * FPb * FNb * FSb`.
- v1 (`processes.py:1014`): identical.
- v3 (`benthic_algae.py:362-369`):
  `growth * limit_phosphorus * limit_nitrogen * limit_light * limit_density`.
- Match (modulo wiring bugs feeding into each factor).
- Severity: match

### B11. mub combination — option 2 (limiting nutrient)
- Fortran (`modBenthicAlgae.f90:335`):
  `mub = mub_max_tc * FLb * FSb * min(FPb, FNb)`.
- v1 (`processes.py:1015-1016`): two-branch `where` equivalent to min.
- v3 (`benthic_algae.py:371-376`):
  `where(FP > FN, growth*FN*FL*FSb, growth*FP*FL*FSb)`. Equivalent to min.
- Match.
- Severity: match

### B12. mub combination — option 3 (harmonic mean)
- Fortran: not implemented (only options 1 and 2).
- v1: not implemented (only options 1 and 2).
- v3 (`benthic_algae.py:377-378`): raises `ValueError("Invalid growth rate option")`
  for option != 1, 2.
- Match (deliberate v3 reject of unsupported option).
- Severity: match

### B13. Benthic biomass integration (dAbdt and Ab update)
- Fortran (`modBenthicAlgae.f90:354`):
  `dAbdt = AbGrowth - AbRespiration - AbDeath` (no settling term, since
  benthic algae are by definition attached).
- v1 (`processes.py:1069-1101`): same.
- v3 (`benthic_algae.py:331-342`):
  `rate = growth - death - respiration` (no settling). Forward-Euler
  integrator with dt_days. Match.
- Severity: match

### B14. NH4 coupling — NH4_AbRespiration and NH4_AbGrowth
- v1 (`processes.py:1506-1547`):
  `NH4_AbRespiration = rnb * AbRespiration * Fb / depth`,
  `NH4_AbGrowth = AbUptakeFr_NH4 * rnb * Fb * AbGrowth / depth`.
- v3 (`benthic_algae.py:486-513`): identical algebra; uses
  `self._cached_depth` set during `run`.
- Match.
- Severity: match

### B15. Mortality routing (AbDeath_OrgN, AbDeath_OrgP, POC, DOC, POM)
- v1 (`processes.py:1362-1381` plus carbon/POM modules):
  - `AbDeath_OrgN = rnb * Fw * Fb * AbDeath / depth`
  - `AbDeath_OrgP = rpb * Fw * Fb * AbDeath / depth`
  - `POC_balgae_mortality = (1/depth) * f_pocb * Fb * Fw * rcb * AbDeath`
  - `DOC_balgae_mortality = (1/depth) * (1 - f_pocb) * Fb * Fw * rcb * AbDeath`
  - `POM_balgae_mortality = AbDeath * Fb * (1 - Fw) / h2`
- v3 (`benthic_algae.py:262-284`): identical algebra; pulls `Fw`, `Fb`, `BWn/BWd`,
  `BWp/BWd`, `BWc/BWd`, `f_pocb` from BALGAE/inline defaults.
- Match.
- Severity: match

### B16. Chla derived variable (Chlb)
- Fortran (`modBenthicAlgae.f90:383`): `Chlb = rab * Ab` with
  `rab = BWa/BWd = 5000/100 = 50` µg-Chla/mg-D.
- v1 (`processes.py:1104-1116`): `Chlb = rab * Ab` with default `BWa=3500`,
  `BWd=100`, so rab=35.
- v3: not exposed in `BenthicAlgae` Process; the registry-level Chlb
  derivation is out of scope for this audit but the BWa default disagreement
  carries forward.
- Severity: minor (BWa default disagreement)
- Category: v1-bug-carry-forward (or deliberate v1 deviation from Fortran)

## Parameter defaults audit

For each parameter the table lists Fortran default, v1 default, and v3 default.
Disagreements are flagged.

### Floating-algae parameters

| Param | Fortran (`modAlgae.f90`) | v1 (`constants.py`) | v3 (`parameters/algae.py`) | Status |
|---|---|---|---|---|
| AWd | 100.0 (line 69) | 100 (line 27) | 100.0 | match |
| AWc | 40.0 (line 73) | 40 (line 28) | 40.0 | match |
| AWn | 7.2 (line 77) | 7.2 (line 29) | 7.2 | match |
| AWp | 1.0 (line 81) | 1 (line 30) | 1.0 | match |
| AWa | 1000.0 (line 85) | 1000 (line 31) | 1000.0 | match |
| KL | 10.0 (line 110) | 10 (line 32) | 10.0 | match |
| KsN | 0.04 (line 137) | 0.04 (line 33) | 0.04 | match |
| KsP | 0.0012 (line 143) | 0.0012 (line 34) | 0.0012 | match |
| mu_max_20 | 1.0 (line 106) | 1 (line 35) | 1.0 | match |
| kdp_20 | 0.15 (line 123) | 0.15 (line 36) | 0.15 | match |
| krp_20 | 0.2 (line 128) | 0.2 (line 37) | 0.2 | match |
| mu_max_theta | 1.047 (line 106) | 1.047 (line 38) | 1.047 | match |
| kdp_theta | 1.047 (line 123) | 1.047 (line 39) | 1.047 | match |
| krp_theta | 1.047 (line 128) | 1.047 (line 40) | 1.047 | match |
| vsap | 0.15 (line 132) | 0.15 (line 41) | 0.15 | match |
| PN | 0.5 (line 149) | 0.5 (Nitrogen TypedDict, line 139) | not in algae DEFAULTS | observation |
| Fpocp | 0.9 (line 155) | 0.9 (carbon, `f_pocp`, line 157) | not in algae DEFAULTS | observation |
| growth_rate_option | 1 (line 114) | 1 (line 42) | 1 | match |
| light_limitation_option | 1 (line 118) | 1 (line 43) | 1 | match |

### Benthic-algae parameters

| Param | Fortran (`modBenthicAlgae.f90`) | v1 (`constants.py`) | v3 (`parameters/balgae.py`) | Status |
|---|---|---|---|---|
| BWd | 100.0 (line 71) | 100 (line 87) | 100.0 | match |
| BWc | 40.0 (line 75) | 40 (line 88) | 40.0 | match |
| BWn | 7.2 (line 79) | 7.2 (line 89) | 7.2 | match |
| BWp | 1.0 (line 83) | 1 (line 90) | 1.0 | match |
| **BWa** | **5000.0 (line 87)** | **3500 (line 91)** | **3500.0** | **disagreement: v1+v3 differ from Fortran** |
| KLb | 10.0 (line 113) | 10 (line 93) | 10.0 | match |
| KsNb | 0.25 (line 142) | 0.25 (line 94) | 0.25 | match |
| KsPb | 0.125 (line 148) | 0.125 (line 95) | 0.125 | match |
| Ksb / KSb | 10.0 (line 117) | 10 (line 96) | 10.0 | match |
| mub_max_20 | 0.4 (line 109) | 0.4 (line 97) | 0.4 | match |
| krb_20 | 0.2 (line 122) | 0.2 (line 98) | 0.2 | match |
| kdb_20 | 0.3 (line 127) | 0.3 (line 99) | 0.3 | match |
| mub_max_theta | 1.047 (line 109) | 1.047 (line 100) | 1.047 | match |
| krb_theta | 1.06 (line 122) | 1.06 (line 101) | 1.06 | match |
| kdb_theta | 1.047 (line 127) | 1.047 (line 102) | 1.047 | match |
| Fw | 0.9 (line 136) | 0.9 (line 105) | 0.9 | match |
| Fb | 0.9 (line 131) | 0.9 (line 106) | 0.9 | match |
| Fpocb | 0.9 (line 160) | 0.9 (carbon, `f_pocb`) | not in balgae DEFAULTS | observation |
| PNb | 0.5 (line 154) | 0.5 (Nitrogen) | not in balgae DEFAULTS | observation |
| b_growth_rate_option | 1 (line 166) | 1 (line 103) | 1 | match |
| b_light_limitation_option | 1 (line 170) | 1 (line 104) | 1 | match |

### Observations on parameter defaults

- **O1. PN, PNb, f_pocp, f_pocb live outside the (b)algae DEFAULTS dicts**
  in v3, but are read by the algae Processes via inline fallback dicts
  (`_FDP_DEFAULTS["f_pocp"] = 0.5`, `_BENTHIC_FDP_DEFAULTS["f_pocb"] = 0.5`).
  Fortran/v1 defaults for these are **0.9**. The v3 inline fallback of 0.5
  silently changes the POC vs DOC mortality split from 90/10 to 50/50.
  Severity: minor (consequence: shifts ~40% of mortality C from POC into DOC).
  Category: v3-deviation.

- **O2. BWa default of 3500 in v1/v3 vs 5000 in Fortran.** This propagates to
  `rab = BWa/BWd` which feeds the benthic chlorophyll-a derived variable
  (`Chlb = rab * Ab`). v1 and v3 produce Chlb = 35*Ab; Fortran produces
  Chlb = 50*Ab. Severity: minor. Category: v1-bug-carry-forward (or v1
  deliberate deviation, see also `clearwater_modules_v3_nsm1_phase0_parameter_audit.md`).

- **O3. PN and PNb live in Fortran modAlgae/modBenthicAlgae** but in v1 they
  live in `NitrogenStaticVariables`. v3 places them with neither — they
  must be passed by user or fall back to 0.5 inside the methods. The 0.5
  fallback matches the v1/Fortran default but the routing is unusual and
  should be documented.

- **O4. `light_attenuation_coefficient` default (kwarg = 1.0)** is roughly
  consistent with `lambda` ≈ 1.0 1/m for moderate-turbidity rivers but
  Fortran `lambda` is computed from a sum of contributions in
  modGlobalParam (`lambda0+lambda1*Chla+lambda2*POM*+lambdam*POM`); v3's
  scalar default short-circuits this calculation. Out of scope for this
  algae audit but needs cross-check in the global-vars audit.

## Conclusions

### Required actions before LimnoTech review

The five critical findings (F1, F2, F3, F4, F5, F14, B1, B2, B3, B6, B7,
B8, B9) collectively render default-instantiated v3 floating-algae and
benthic-algae kinetics incorrect in measurable ways. They fall into three
classes:

1. **Wiring**: the v3 ALGAE_DEFAULTS / BALGAE_DEFAULTS keys (`mu_max_20`,
   `kdp_20`, `krp_20`, `vsap`, `KsN`, `KsP`, `mub_max_20`, `krb_20`,
   `kdb_20`, `KsNb`, `KsPb`, `Ksb`, `mu_max_theta`, etc.) are merged into
   `self` but the rate methods never read them; instead they read the
   shadowing v2 kwargs (`growth_rate_max`, `death_rate`, `repiration_rate`,
   `settling_velocity`, `nitrogen_michaelis_menton_constant`,
   `phosphorus_michaelis_menton_constant`, `density_michaelis_menton_constant`,
   `growth_rate_correction`, `death_rate_correction_factor`,
   `repiration_rate_correction_factor`, `light_attenuation_coefficient`,
   `light_limitation_constant`). The kwargs default to values that disagree
   with v1/Fortran (notably 0.0 for the rates and 1.0 for the corrections).
   Either the rate methods must be rewritten to read the DEFAULTS keys, or
   the DEFAULTS-merge must overwrite the kwargs.
2. **Formula errors**: F5 (FL option-1 misplaced parenthesis), F14
   (harmonic-mean wrong zero-guard `FP==1`), B6 (FLb option-3 sign of
   exponent argument).
3. **Inheritance leakage**: B7, B8, B9 — BenthicAlgae inherits FloatingAlgae's
   `limit_nitrogen`/`limit_phosphorus`/`limit_density` and reads the
   floating-algae KsN/KsP/Ksb constants. BenthicAlgae must expose its own
   half-saturation constants (KsNb, KsPb, Ksb) to the limit methods, or
   override them.

The minor findings (F8, F9, F11, O1, O2) are documentation- or routing-level
issues that should be resolved or explicitly documented in the LimnoTech
review packet.

### Acceptable deviations to document

- BWa = 3500 (v1/v3) vs 5000 (Fortran). v1-aligned; document as v1-deviation.
- PN/PNb/f_pocp/f_pocb default values: v3 uses inline fallbacks of 0.5; v1
  uses 0.5 (PN/PNb) and 0.9 (f_pocp/f_pocb). The 0.5 PN matches; the 0.5
  f_pocp/f_pocb is a v3 deviation that should be raised to 0.9 to match
  v1/Fortran, or added to the DEFAULTS dicts and documented.
- Floating-algae FL option 2 (Smith) abs(KL)<1e-10 fallback returns 1.0;
  acceptable.

### Items to escalate

- The systematic shadowing of v3 DEFAULTS by v2 kwargs is the single largest
  source of risk and should be the first item on the Phase 2.A.1 follow-up
  list. It implies that every kinetic test passing today does so only
  because the test fixture explicitly supplies the kwargs that match v1
  defaults. Default-instantiated regression simulations would diverge from
  v1 immediately.
- The three formula bugs (F5, F14, B6) are independent of the wiring issue
  and persist even when the kwargs are explicitly supplied.
- The PAR vs solar-radiation question (F9) requires a cross-check with the
  variable-registry contract for `solar_radiation`. If `solar_radiation` is
  conventionally already PAR, no fix is needed; if it is total shortwave,
  Fr_PAR=0.47 scaling must be applied before passing to `limit_light`.
