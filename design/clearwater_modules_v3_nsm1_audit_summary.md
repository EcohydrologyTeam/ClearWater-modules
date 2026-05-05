# v3 NSM1 — Consolidated Three-way Audit Summary

**Date:** 2026-05-05
**Scope:** Synthesis of five Process-level three-way audits (Fortran NSM1 vs v1 Python NSM1 vs v3) covering FloatingAlgae, BenthicAlgae, Nitrogen, Phosphorus, Carbon, DOX, POM, CBOD, Pathogen, N2, Alkalinity, plus shared utilities and the parameter library.
**Inputs:**
- `docs/clearwater_modules_v3_nsm1_audit_algae.md` (Audit 1)
- `docs/clearwater_modules_v3_nsm1_audit_n_p.md` (Audit 2)
- `docs/clearwater_modules_v3_nsm1_audit_c_dox.md` (Audit 3)
- `docs/clearwater_modules_v3_nsm1_audit_simple_constituents.md` (Audit 4)
- `docs/clearwater_modules_v3_nsm1_audit_utilities_params.md` (Audit 5)

**Framing.** The legacy Fortran NSM1 was not adequately validated and is known to contain flaws. v3 deviation from Fortran is therefore not, by itself, evidence of a v3 defect; in many cases it is the corrective action. Each finding below is categorized by direction (Fortran-flaw / v1-flaw / v3-deliberate / v3-bug) rather than by Fortran-distance.

---

## Executive summary

- **Total findings**: 64 (across 5 audits, deduplicated).
  - Real v3 bugs requiring code fixes: **22** (algae wiring 11; algae formulas 3; Nitrogen 4; Carbon 3; DOX 1).
  - v3 deliberate improvements (already correct, vindicated by Fortran): **8**.
  - Three-way disagreements (need reconciliation): **6**.
  - Likely v1/Fortran flaws bypassed by v3 (no v3 action): **8**.
  - Documentation defects (cheap fixes): **9**.
  - Matches: **~120 blocks across all five audits (not enumerated here).**

- **Headline:** the v3 NSM1 1.0.0 LimnoTech review packet claimed "review-ready" in Phase 8 is **not** review-ready. The audit found a coherent class of correctness defects (algae wiring, rca/rcb stoichiometry as raw weights, missing DIC/CBOD source, missing DOX/SOD attenuation, phantom NH4 source, broken NO3 algal-uptake split) that the existing parity tests do not exercise because they call default-instantiated processes with explicit kwargs that mask the wiring bugs and pass the same wrong stoichiometric ratios into both sides. Estimated time-to-review-ready: **2-3 weeks of focused work** (code fixes ~3-5 days; new Fortran-anchored regression tests ~5-7 days; doc updates and reconciliation discussions ~3-5 days).

- **Wins worth stating up front.** The seven sentinel-999 corrections (`vsop`, `vs`, `SOD_20`, `SOD_theta`, `kaw_20_user`, `kah_20_user`, `pressure_mb`) are all vindicated by the Fortran-coded defaults (see Audit 5). The `SOD_tc` pure-Arrhenius split, the `PAR` toggle inversion (avoiding v1's two-arg `xr.where` NaN bug), and the resolved-Q architectural decisions (clip_negative_state, registry rate-variable convention, fixed-point/forward-Euler integrator) are sound. The Alkalinity, N2, Pathogen, POM, and CBOD Processes match Fortran/v1 to within documented minor deviations.

---

## Categorized findings

### Category 1: Real v3 bugs (must fix before review)

| # | Process | Location | Issue | Severity | Fortran says | v1 says | Recommended fix |
|---|---|---|---|---|---|---|---|
| 1 | FloatingAlgae | `floating_algae.py:419-423` | `rate_growth` reads kwarg `growth_rate_max` (default 1.0) and `growth_rate_correction` (default 1.0); ALGAE_DEFAULTS `mu_max_20=1.0`, `mu_max_theta=1.047` are merged into `self` but never read. Default theta=1.0, not 1.047. | Critical | mu_max=1.0, theta=1.047 | mu_max=1.0, theta=1.047 | Rewire `rate_growth` to read `self.mu_max_20`, `self.mu_max_theta` from ALGAE_DEFAULTS; retire the kwargs or make them read DEFAULTS. |
| 2 | FloatingAlgae | `floating_algae.py:470-482` | `rate_respiration` uses kwarg `repiration_rate` (default 0.0); ALGAE_DEFAULTS `krp_20=0.2` unread. Default-instantiated respiration is identically zero. | Critical | krp=0.2/d | krp=0.2/d | Rewire to `self.krp_20`. |
| 3 | FloatingAlgae | `floating_algae.py:461-468` | `rate_death` uses kwarg `death_rate` (default 0.0); ALGAE_DEFAULTS `kdp_20=0.15` unread. Default-instantiated death is identically zero. | Critical | kdp=0.15/d | kdp=0.15/d | Rewire to `self.kdp_20`. |
| 4 | FloatingAlgae | `floating_algae.py:484-488` | `rate_settling` uses kwarg `settling_velocity` (default 0.0); ALGAE_DEFAULTS `vsap=0.15` unread. Default-instantiated settling is identically zero. | Critical | vsap=0.15 m/d | vsap=0.15 m/d | Rewire to `self.vsap`. |
| 5 | FloatingAlgae | `floating_algae.py:556-564` | `limit_light` option 1: `np.log` is called on `(KL+PAR)` only; the denominator `(KL+PAR*exp(-Ld))` divides the entire expression instead of being inside the log. Operator-precedence bug. | Critical | (1/Ld)*log((KL+PAR)/(KL+PAR*exp(-Ld))) | same as Fortran | Fix parenthesization: pass the ratio to `np.log`. |
| 6 | FloatingAlgae | `floating_algae.py:436-452` | Harmonic-mean growth (option 3) zero-guard: `where(FP == 1, 0, ...)` instead of `where(FP == 0, 0, ...)`. Shuts down growth when P is fully non-limiting. | Critical | guard on FP==0 | guard on FP==0 | Change condition to `FP == 0`. |
| 7 | BenthicAlgae | `benthic_algae.py:344-359` | Same wiring bug as #1: `rate_growth` reads `growth_rate_max`/`growth_rate_correction` instead of `mub_max_20=0.4`, `mub_max_theta=1.047` from BALGAE_DEFAULTS. | Critical | 0.4/d, 1.047 | 0.4/d, 1.047 | Rewire to BALGAE_DEFAULTS keys. |
| 8 | BenthicAlgae | `benthic_algae.py:385-401` | Same wiring bug: respiration reads kwarg (0.0) instead of `krb_20=0.2`, `krb_theta=1.06`. | Critical | 0.2, 1.06 | 0.2, 1.06 | Rewire. |
| 9 | BenthicAlgae | inherited `rate_death` | Same wiring bug: death reads kwarg (0.0) instead of `kdb_20=0.3`. | Critical | 0.3, 1.047 | 0.3, 1.047 | Rewire (or override `rate_death` in BenthicAlgae). |
| 10 | BenthicAlgae | `benthic_algae.py:438-454` | Steele light limitation (option 3): implements `x*exp(x-1)` instead of `x*exp(1-x)`; sign of exponent argument reversed. | Critical | x*exp(1-x) | x*exp(1-x) | Fix exponent sign. |
| 11 | BenthicAlgae | inherits FloatingAlgae.limit_nitrogen | Reads floating-algae KsN=0.04 instead of KsNb=0.25 (~6x too small). | Critical | KsNb=0.25 | KsNb=0.25 | Override in BenthicAlgae or parameterize the half-saturation constant. |
| 12 | BenthicAlgae | inherits FloatingAlgae.limit_phosphorus | Reads floating-algae KsP=0.0012 instead of KsPb=0.125 (~100x too small). | Critical | KsPb=0.125 | KsPb=0.125 | Override or parameterize. |
| 13 | BenthicAlgae | `benthic_algae.py:467-480` | `limit_density` reads kwarg `density_michaelis_menton_constant=1.0` instead of BALGAE_DEFAULTS `Ksb=10.0` (10x too small). | Critical | Ksb=10.0 | Ksb=10.0 | Rewire to `self.Ksb`. |
| 14 | Nitrogen | `v2/nitrogen.py:334-353` | Phantom `ammonium_decay_nitrate` source term: positive first-order source `1.0/d * NH4` with no v1 or Fortran analogue. Causes exponential NH4 growth at default kwargs. | Critical | absent | absent | Drop the term from the rate sum, or default `ammonium_decay_rate=0.0`. |
| 15 | Nitrogen | `v2/nitrogen.py:451-457`; also `nitrification_inhibition`, `nitrate_denitrification`, `nitrate_bed_denitrification`, `ammonium_nitrification` | Kinetic methods read legacy v2 kwargs (uniformly `1.0`) rather than NITROGEN_DEFAULTS (matching v1/Fortran). 5x to 500x rate inflation: `knit` 10x, `kdnit` 500x, `rnh4` infinite-vs-zero, `vno3` infinite-vs-zero, `KNR` ~1.7x. | Critical | matches v1 | matches Fortran | Rewire all six methods to read `self.knit_20/knit_theta/kdnit_20/kdnit_theta/rnh4_20/rnh4_theta/vno3_20/vno3_theta/KNR` from NITROGEN_DEFAULTS; retire or zero the legacy kwargs. |
| 16 | Nitrogen | `v2/nitrogen.py:565-576` | `nitrate_uptake_floating_algae` uses static `float_algea_faction_uptake_from_nitrate=1.0` instead of dynamic `1 - algal_nh4_uptake_fraction`. NH4 sink uses dynamic fraction; NO3 sink uses 1.0. The two paths do not sum to `rna * ApGrowth`; algal-N mass balance violated by ~1.5x at default PN. | Critical | dynamic 1-PN-fraction | dynamic | Rewire to read `1 - floating_algae_process.algal_nh4_uptake_fraction`; retire the static parameter. |
| 17 | Nitrogen | `v2/nitrogen.py:578-594` | `nitrate_uptake_benthic_algae`: divides by `algal_chlorophyll=AWa=1000` (floating denom) instead of `BWd=100`; missing `/depth` divisor; uses `fraction_bottom_area=1.0` instead of `Fb=0.9`; uses static `benthic_algea_faction_uptake_from_nitrate=0.5` instead of dynamic `1 - balgae_nh4_uptake_fraction`. Stoichiometry, units, and dynamics all wrong. | Critical | rnb=BWn/BWd, *Fb/depth, dynamic | rnb=BWn/BWd, *Fb/depth, dynamic | Reconstruct: `(1 - balgae_nh4_uptake_fraction) * (BWn/BWd) * Fb * balgae_growth_rate / depth`. Add `balgae_no3_uptake_fraction` cache to BenthicAlgae if needed. |
| 18 | Phosphorus utilities | `utils/partitioning.py:31` | `fdp = 1/(1 + kdpo4 * Solid / 0.000001)` divides by 1e-6 (i.e. multiplies by 1e6); Fortran divides by 1e6 (i.e. multiplies by 1e-6). Factor of 1e12 wrong. Latent at default `kdpo4=0` (both forms give fdp=1); breaks the moment a user sets `kdpo4>0`. v1 has the same bug; v3 inherits it. | Critical (gated) | `kdpo4 * Solid / 1.0E6` | `kdpo4 * Solid / 0.000001` (v1 bug) | Change literal to `1.0E6` (or multiply Solid by 1e-6). Add an MMS test covering `kdpo4>0`. |
| 19 | Carbon | `carbon.py:429-430, 436-441` | DIC algal coupling: uses `self.AWc=40` and `self.BWc=40` directly instead of `rca=AWc/AWa=0.04` and `rcb=BWc/BWd=0.4`. Floating algae DIC terms 1000x too large; benthic 100x too large. Same defect in DOX (#21). | Critical | rca=AWc/AWa, rcb=BWc/BWd | rca=AWc/AWa, rcb=BWc/BWd | Replace `self.AWc` with `self.AWc/self.AWa` (or pre-compute `rca`); same for `rcb`. |
| 20 | Carbon | `carbon.py:451-459` | dDIC/dt omits the CBOD oxidation source. Fortran adds `CBOD_DIC_Oxidation = sum(CBOD_Oxidation)/roc/12000`; v1 includes it; v3 silently drops it. CBOD process already exposes `cbod_oxidation_rate`. | Critical | included | included | Add `+ cbod_process.cbod_oxidation_rate / self.roc / 12000.0` to dDIC/dt (gated on `use_cbod`). |
| 21 | Carbon | `carbon.py:372` | `poc_hydrolysis = kpoc_tc * poc * dox_attenuation` — v3 multiplies POC hydrolysis by `DOX/(KsOxmc+DOX)`; neither Fortran nor v1 attenuates POC hydrolysis. POC->DOC is physical/chemical (cell-wall fragmentation, leaching), not biochemically O2-limited. | Critical | `kpoc_tc * POC` only | `kpoc_tc * POC` only | Drop `dox_attenuation` factor from POC hydrolysis. |
| 22 | DOX | `dox.py:400-489` | Same rca/rcb defect as #19 in `_floating_algae_growth_flux`, `_floating_algae_respiration_flux`, `_benthic_algae_growth_flux`, `_benthic_algae_respiration_flux`. 1000x error for floating O2 photosynthesis/respiration, 100x for benthic. | Critical | rca=AWc/AWa, rcb=BWc/BWd | same | Same fix as #19. |
| 23 | DOX | `dox.py:554-573` and `utils/sediment.py:16-31` | SOD sink omits `DOX/(DOX+KsSod)` Monod attenuation that Fortran applies in `modGlobalParam.f90:254`. Under hypoxia, v3 sediment keeps consuming O2 at full Arrhenius rate; clip-to-zero masks the conservation violation. v3 doc acknowledges as "Phase 5.5 deferred". | Critical | DOX-attenuated SOD | DOX-attenuated SOD | Multiply `_sod_flux` by `dox/(dox+KsSod)` (gated on `use_DOX`). |

(Findings #1-#13 are the 13 algae-Process bugs from Audit 1; #14-#17 are Nitrogen from Audit 2; #18 is Phosphorus from Audit 2; #19-#21 are Carbon from Audit 3; #22-#23 are DOX from Audit 3.)

### Category 2: v3 deliberate improvements (already correct; document rationale)

| # | Process / Utility | Improvement | Source of rationale |
|---|---|---|---|
| 1 | DOX defaults | `pressure_mb=1013.25` replaces v1's `2026.5` (2x error). Restored to ISO 2533 standard atm. | `parameter_defaults_corrections.md` Section 1.7; Audit 5 confirms. |
| 2 | DOX defaults | `SOD_theta=1.060` replaces v1's sentinel `999`. Matches Fortran `modGlobalParam.f90:122` exactly. | Section 1.4; Audit 5. |
| 3 | DOX defaults | `vs=0.1` replaces v1's `999`. Matches Fortran exactly. | Section 1.2; Audit 5. |
| 4 | DOX defaults | `kaw_20_user=0.0` replaces v1's `999`. Matches Fortran. | Section 1.5; Audit 5. |
| 5 | Phosphorus defaults | `vsop=0.1` replaces v1's `999`. (Differs from Fortran's 0.01 by 10x; literature-defensible — see Cat 3.) | Section 1.1. |
| 6 | DOX defaults | `SOD_20=1.0` replaces v1's `999`. (Differs from Fortran's 0.2 by 5x — see Cat 3.) | Section 1.3. |
| 7 | Utilities | `SOD_tc` pure-Arrhenius split: DOX-Monod factor moved out of the shared utility into the DOX consumer. Cleaner architecture. (Currently *not* re-applied at the consumer — see Cat 1 #23.) | Section 3.2; Phase 1.1. |
| 8 | Utilities | `PAR` toggle inversion: v3 returns `q_solar * Fr_PAR` unconditionally; the `use_Algae|use_Balgae` toggle moved to the consumer Process. Avoids v1's latent two-arg `xr.where` NaN-on-disable bug. | Section 3.4; Phase 1.1. |
| 9 | Utilities | `clip_negative_state` + `Diagnostics`: v3-only safety net. Logic correct (clip target = 0; per-call rate-limited logging; aggregate suppression stub). | Q7 design spec. |
| 10 | Utilities | `np.select` dim-stripping fix in `kah_20`/`kaw_20` (Phase 5.5). Match-preserving. | Phase 5.5 work item. |
| 11 | Nitrogen | Cached `nitrification_flux_rate` and `denitrification_flux_rate` (positive-magnitude flux registry convention). Used by Alkalinity and N2 to avoid recomputation. | Resolved Q10. |
| 12 | N2 | Adds `denit_source = nitrogen_process.denitrification_flux_rate` to dN2/dt — closes the N mass balance (NO3 -> N2 -> N2sat exchange) that Fortran and v1 silently break. Collapses to v1 form when `use_nitrogen=False`. | Phase 3.4 task brief Item 1. |
| 13 | Phosphorus | DIP derived variable: v3 uses `DIP = TIP * fdp` (matching v1, which is correct); Fortran's `DIP = TIP / fdp` is the sign-inverted bug. v3 is right. | Audit 2 finding P14. |
| 14 | FloatingAlgae | Bug #13 (NH4_ApRespiration was returning 0) and bug #14 (NH4_ApGrowth was returning 0) fixed in Phase 2.A. | Phase 2.A. |
| 15 | FloatingAlgae | Forward-Euler integrator with `dt_days = total_seconds()/86400` (Phase 2.A bug #4/#16 fixes). | Phase 2.A. |
| 16 | DOX | Salinity correction on O2sat: v3 omits, but so does v1 (Fortran has it). v3 inherits v1; for freshwater this evaluates to factor 1.0 and is harmless. Document as deferred (see Cat 4). | Audit 3 finding C6. |

### Category 3: Three-way disagreements (need reconciliation)

| # | Process | Parameter / formula | Fortran | v1 | v3 | Recommended resolution |
|---|---|---|---|---|---|---|
| 1 | global_vars | `lambdam` (POM Beer-Lambert coefficient) | 0.174 | 0.0174 | 0.0174 | **Likely v1 typo propagated to v3.** Fortran 0.174 consistent with QUAL2K Table 6. Reconcile with LimnoTech; recommend correcting v3 to 0.174. |
| 2 | DOX | `SOD_20` | 0.2 | 999 (invalid) | 1.0 | v3's 1.0 is 5x Fortran's 0.2. Pick a value with literature backing (Chapra 1997 supports values from 0.2 to ~5 depending on sediment type) and document. |
| 3 | DOX | `kah_20_user` | 1.0 | 999 (invalid) | 0.0 | v3 zero blocks default reaeration; Fortran's 1.0 produces meaningful default reaeration. Recommend either set v3 `kah_20_user=1.0` (Fortran-aligned) or change `hydraulic_reaeration_option` default to a non-trivial Covar/Owens-Gibbs formula. |
| 4 | Phosphorus | `vsop` (OrgP settling velocity) | 0.01 | 999 | 0.1 | v3's 0.1 is 10x Fortran's 0.01. Sentinel rescue, but without explicit literature anchor. Recommend lowering v3 to 0.01 (Fortran-aligned) or document the choice. |
| 5 | Nitrogen | `vson` value & location | 0.01 (in modGlobalParam) | 0.01 (in GlobalVars) | **0.1** in `parameters/nitrogen.py:vson_20`, **0.01** in `parameters/global_vars.py:vson` | Two contradictory values inside v3 itself. Consolidate to 0.01 (Fortran/v1-aligned) and remove the duplicate. |
| 6 | Nitrogen / Carbon | minor theta values: `kon_theta` (v3 1.074 vs v1/Fortran 1.047), `kdnit_theta` (v3/v1 1.08 vs Fortran 1.045), `rnh4_theta` (v3/v1 1.047 vs Fortran 1.074), `vno3_theta` (v3/v1 1.045 vs Fortran 1.08), `ksbod_theta` (v3/v1 1.047 vs Fortran 1.024) | various | various | various | Most of these are v1-inherited choices; v1 already diverged. Document as v1-inherited and assess calibration impact (10-15% off-reference temperature). |

### Category 4: Likely v1/Fortran flaws (no v3 action needed; document for posterity)

| # | Reference | What it does | How v3 handles it |
|---|---|---|---|
| 1 | v1 (`shared/processes.py` `PAR`) | Two-arg `xr.where(cond, value)` returns NaN in the false branch when both algae modules disabled. Latent NaN-propagation bug. | v3 returns `q_solar*Fr_PAR` unconditionally; toggle moved to consumer. Bypassed cleanly. |
| 2 | v1 sentinels (999) | `vsop=999`, `vs=999`, `SOD_20=999`, `SOD_theta=999`, `kaw_20_user=999`, `kah_20_user=999`. Crashes or absurd kinetics if used as-is. | v3 corrected all 7 (see Category 2 #1-#6). Vindicated. |
| 3 | v1 `pressure_mb=2026.5` | 2x ISO standard atm. Doubles all DOX saturation calculations. | v3 corrected to 1013.25. Vindicated. |
| 4 | Fortran `DIP = TIP/fdp` (modPhosphorus.f90:223) | Sign-inverted: should be `DIP = TIP * fdp` since `fdp` is the dissolved fraction. | v1 and v3 use the correct `TIP * fdp`. v3 is right; Fortran is wrong. |
| 5 | Fortran SOD-derived DIC sediment release fallback (`use_SedFlux=False`) | Fortran computes `DICfromBed = SOD_tc/roc/depth/12000` when `use_SedFlux=False`; v1 always uses this form (Fortran-flaw-as-fallback). | v3 returns 0 when `use_SedFlux=False`; documented as Phase 5.A scope. Not technically a v3 bug, but a deferred completion. |
| 6 | Fortran `SOD_tc` DOX-Monod conflation in shared utility | The DOX-attenuation factor lives inside the shared `SOD_tc` helper, coupling sediment temperature dynamics to dissolved-O2 state. Architectural smell. | v3 split it cleanly (utility = pure Arrhenius; consumer applies Monod). The split is correct; the consumer just hasn't yet *re-applied* the Monod factor (see Cat 1 #23). |
| 7 | v1 OrgN/OrgP settling | v1 omits Arrhenius correction on `vson`/`vsop` (matches Fortran). | v3 applies `vson_theta=1.024` Arrhenius (Phase 2.B documented divergence). At T=20 C zero impact; at T=25 C ~12.6% rate increase. Calibration-level deviation. |
| 8 | Fortran `kpom2_20=0.01` vs v1/v3 `kpom_20=0.1` | Fortran sets POM dissolution rate 10x lower than v1/v3. v1's value is documented but Fortran is the literature reference. | v3 inherits v1. Document as v1-inherited; calibrate accordingly. |

### Category 5: Documentation defects (cheap fixes)

| # | Doc location | Current claim | Correct claim |
|---|---|---|---|
| 1 | `parameter_defaults_corrections.md` Section 2.8 | v1's `lambdas * Solid` term is "commented out / defined but not used" | v1 source `shared/processes.py:232` applies `lambdas * Solid` unconditionally; v3 reproduces this. Retract the claim. |
| 2 | `parameter_defaults_corrections.md` Phase 1.1 | "suspicious unit factor `1/(1 + kdpo4 * Solid / 0.000001)`" | The formula is dimensionally consistent with `(L/kg)(mg/L)(1 kg / 1e6 mg)`. But the v1/v3 `/0.000001` is mathematically the *inverse* of Fortran's `/1.0E6` — it's a real bug, not a unit formality. Update to record this as a v1-inherited bug, latent at `kdpo4=0`. |
| 3 | `parameter_defaults_corrections.md` Section 3.6 | Kelvin offset is 273.15 | `utils/conversions.py` re-exports v2's `celsius_to_kelvin` which returns `T_C + 273.16` (with comment "for testing consistency with v1"). Update Section 3.6 to reflect actual v3 behavior. |
| 4 | `parameter_defaults_corrections.md` Section 1.1 | `vsop` correction rationale | Add: "Fortran value is 0.01 m/d; v3 chose 0.1 m/d to match v3's `vs=0.1` and the literature midrange. Acknowledge the 10x divergence from Fortran." |
| 5 | `parameter_defaults_corrections.md` Section 1.3 | `SOD_20` correction rationale | Add: "Fortran value is 0.2 g-O2/m^2/d; v3 chose 1.0 representing a moderate-eutrophic baseline. Acknowledge the 5x divergence from Fortran." |
| 6 | `parameter_defaults_corrections.md` Section 1.6 | `kah_20_user=0.0` | Add: "Fortran value is 1.0; v3 chose 0.0 to make the user-override explicit. **Behavioral consequence:** at default `hydraulic_reaeration_option=1`, v3 default reaeration = 0 1/d while Fortran default = 1.0 1/d. Document prominently in migration notes." |
| 7 | `utils/reaeration.py` docstrings | Author attributions for kah_20 options (Covar/Owens-Gibbs/Churchill/Tsivoglou-Wallace/Padden-Gloyna/etc) | Disagree with Fortran source comments (Owens, O'Connor, Churchill, Cover, Melching-Flores, Tsivoglou-Neal, Thackson-Dawson). Reconcile to Fortran's attributions. |
| 8 | `utils/reaeration.py` docstrings | Author attributions for kaw_20 options (Banks 1975, etc) | Disagree with Fortran (Broecker 1978, etc). Reconcile. |
| 9 | LimnoTech review packet (Phase 8 deliverable) | "v3 NSM1 1.0.0 is review-ready" | Update to reflect this audit's findings: ~22 critical fixes pending, estimated 2-3 weeks to actual review-ready. |

### Category 6: Matches (counts only)

| Process | Match-block count |
|---|---|
| FloatingAlgae | 13 (F6, F7, F10, F12, F13, F15, F16, F17, F18, F19, F20, plus 17 of 19 parameter defaults) |
| BenthicAlgae | 6 (B4, B5, B10, B11, B12, B13, B14, B15) plus 18 of 19 parameter defaults |
| Nitrogen | 7 (N3, N5, N6, N7, N8, N9 structure, N15, N16) plus most theta-corrected blocks |
| Phosphorus | 11 (P1-P5 budgets, P7-P13 structure) plus 4 of 7 parameter defaults |
| Carbon | 8 (POC settling, POC mortality, DOC mortality, DIC reaeration form, DOC-as-DOX-sink, dPOC/dt structure, dDOC/dt structure modulo POM) plus 10 of 10 parameter defaults |
| DOX | 6 (atm reaeration form, nitrification flux structure, DOC oxidation sink, CBOD oxidation sink, dDOX/dt structure, cached-rates contract) plus 7 of 10 parameter defaults |
| POM | 5 (algae settling, balgae mortality, POC settling, sign convention, structure) plus 6 of 6 parameter defaults |
| CBOD | 2 (oxidation Monod) plus 3 of 5 parameter defaults |
| Pathogen | 4 (natural decay, settling, sign convention, plus parameter defaults) |
| N2 | 5 (Henry's law, pwv, N2sat, atm exchange, structure) plus all hard-coded constants |
| Alkalinity | 6 (#19-#24 structure, sign convention) plus all 6 ratios |
| Utilities | `kah_20` (9 options), `kaw_20` (13 options), `ka_tc`, `L`, `fdp` (single-solid form) — all match Fortran exactly |

---

## Top systemic patterns identified

### Pattern A: Wiring defect (algae, Nitrogen)

- **Description.** v3 maintains a v3-style DEFAULTS attribute set that is merged into `self.__dict__` at construction, alongside legacy v2-style kwargs that shadow them. The kinetic methods read the v2 kwargs, not the DEFAULTS. The kwargs default to values that disagree with the v1/Fortran NSM1 reference (typically 0.0 for rates, 1.0 for thetas, 1.0 for sediment release). The DEFAULTS values are correct but never reach the rate methods unless the user manually overrides via the kwarg names.
- **Affected Processes.** FloatingAlgae (`mu_max`, `kdp`, `krp`, `vsap` all wiring-broken), BenthicAlgae (`mub_max`, `krb`, `kdb`, `KsNb`, `KsPb`, `Ksb` all wiring-broken via inheritance and kwarg override), Nitrogen (`KNR`, `knit_20`, `kdnit_20`, `rnh4_20`, `vno3_20` plus their thetas all wiring-broken).
- **Root cause.** The Phase 1-2 migration left both surfaces in place "for back-compat". The kinetic call sites were never updated to read DEFAULTS. The legacy kwargs are documented as v2 back-compat shims but operate as the *only* live wiring.
- **Why tests didn't catch it.** Every existing kinetic parity test calls `FloatingAlgae(mu_max=..., krp_20=..., kdp_20=..., vsap=...)` with explicit kwargs *that match the legacy kwarg names but use v1/Fortran values*. The kwarg overrides hide the wiring defect. A `FloatingAlgae()` call (default args) produces algae with zero respiration, zero death, zero settling, theta=1.0, and growth rate 1.0/d — none of which matches v1.
- **Fix plan.**
  - Sweep all six v3 NSM1 Processes and rewire each kinetic method to read the DEFAULTS keys (`self.mu_max_20`, `self.knit_20`, etc.) instead of the legacy v2 kwarg.
  - Retire the legacy kwargs (preferred), or if back-compat is required, default them to the DEFAULTS values rather than to `0.0/1.0`.
  - Add an integration test that asserts `Process()` (no kwargs) reads from DEFAULTS, not from the legacy kwarg defaults.

### Pattern B: Wrong stoichiometric ratios used as rates (Carbon, DOX)

- **Description.** v3 Carbon and DOX algal coupling terms use `self.AWc=40` and `self.BWc=40` directly where Fortran and v1 use `rca=AWc/AWa=0.04` and `rcb=BWc/BWd=0.4`. AWc and BWc are *raw stoichiometric weights* (mg-C per ug-Chla and per mg-D respectively); the *ratio* converts algal biomass (in chlorophyll or dry weight) to mg-C. v3 omits the division.
- **Affected blocks.**
  - Carbon: `dic_algal_resp`, `dic_algal_photo`, `dic_balgae_resp`, `dic_balgae_photo`.
  - DOX: `_floating_algae_growth_flux`, `_floating_algae_respiration_flux`, `_benthic_algae_growth_flux`, `_benthic_algae_respiration_flux`.
- **Error magnitudes.** Floating algae 1000x (AWa=1000), benthic algae 100x (BWd=100). v3 algal photosynthesis O2 source is 1000x Fortran's; v3 algal respiration O2 sink is 1000x Fortran's; v3 DIC algal terms are 1000x/100x off.
- **Why tests didn't catch it.** The Carbon and DOX parity tests (`test_5_carbon_calculations_v2.py:343-412`, `test_5_dox_calculations_v2.py:328-330`) explicitly call `v1.DIC_algal_respiration(rca=AWc=40)` — passing the same wrong value as `rca`. Same-error parity. Test docstrings even document the (incorrect) "rca = AWc / AWa = AWc in v3's per-Chla convention" identity, which simple arithmetic disproves.
- **Fix plan.**
  - Replace `self.AWc` with `self.AWc / self.AWa` (or pre-compute `rca` once at the top of `run`) at all four floating-algae sites in carbon.py and dox.py.
  - Same for `self.BWc -> self.BWc / self.BWd` at all four benthic sites.
  - Add Fortran-anchored numerical regression tests with explicit reference values (e.g., `AWc=40, AWa=1000, ApRespiration=0.5 ug-Chla/L/d -> 1.667e-6 mg-C/L/d`).
  - Update test docstrings to remove the incorrect identity claim.

### Pattern C: Missing terms (Carbon DIC budget, DOX SOD attenuation, others)

- **C-DIC: Missing CBOD oxidation source in dDIC/dt.** Fortran (`modCarbon.f90:262-266`) and v1 (`processes.py:2854`) both add `DIC_CBOD_oxidation` to dDIC/dt. v3 omits it. Silent missing source term; v3 underestimates DIC production from CBOD oxidation. Fix: add `+ cbod_process.cbod_oxidation_rate / self.roc / 12000.0` to dDIC/dt.
- **C-DIC: Missing SOD-derived DIC release fallback under `use_SedFlux=False`.** Fortran/v1 use `DICfromBed = SOD_tc/roc/depth/12000` when `use_SedFlux=False`; v3 returns 0. Documented as Phase 5.A scope, but flagged because at default `use_SedFlux=False, JDIC=0` v3 silently produces zero sediment DIC release.
- **DOX-SOD: Missing DOX-Monod attenuation.** Fortran's `SOD_tc` is multiplied by `DOX/(DOX+KsSod)`. v3's pure-Arrhenius split (a deliberate architectural improvement) requires the consumer to re-apply the Monod factor. The DOX consumer currently doesn't. Under hypoxia, v3 sediment keeps consuming O2 at full rate. Fix: multiply `_sod_flux` by `dox/(dox+KsSod)` (gated on `use_DOX`).
- **DOX-SedFlux: Missing `use_SedFlux` branch.** Fortran swaps `SOD_tc/depth` for `SOD_Bed/depth` when `use_SedFlux=True`; v3 always uses `SOD_tc/depth`. Phase 5.5 scope.
- **DOX: Missing salinity correction on O2sat.** Fortran applies `O2sat *= exp(-Salinity * (0.017674 - 10.754/Tk + 2140.7/Tk^2))`. v1 omits; v3 inherits the omission. Freshwater impact = 0; brackish/estuarine impact up to ~20%. Defer if v3 1.x targets fresh water only; flag if estuarine deployment is on the roadmap.

### Pattern D: Test suite blind spots

- **Why parity tests didn't catch the bugs.**
  - Algae wiring defects (Pattern A) are masked because every test passes the legacy kwargs explicitly with v1/Fortran values. Default instantiation is never tested.
  - rca/rcb errors (Pattern B) are masked because both sides of the parity comparison receive the same wrong stoichiometric ratio (`rca=AWc=40`).
  - Phantom NH4 source (#14) and broken NO3-uptake split (#16, #17) are masked because the Nitrogen kinetic tests fix `NH4`, `NO3`, `ApGrowth` and check derivative magnitude, not mass balance closure.
  - Missing CBOD-DIC source (#20) and missing SOD-DOX attenuation (#23) are masked because the Carbon and DOX parity tests don't exercise CBOD or hypoxia regimes.
- **Recommended Fortran-anchored regression tests.**
  - Default-instantiation regression: assert `FloatingAlgae()` (no kwargs) produces v1/Fortran-matching `dApdt` for a fixed input fixture. Repeat for BenthicAlgae, Nitrogen, etc.
  - Stoichiometric-ratio regression: assert `Carbon().dic_algal_resp(...)` matches a Fortran-anchored value (numerical, not parity).
  - Mass-balance closure: assert `NH4_ApGrowth + NO3_ApGrowth = rna * ApGrowth` (within tolerance) over a 100-step integration with non-trivial NH4/NO3.
  - CBOD-DIC source: assert dDIC/dt under `use_cbod=True, CBOD>0` includes the CBOD oxidation contribution.
  - Hypoxic SOD: assert SOD sink shrinks toward zero as DOX -> 0 (not just clip).
  - `kdpo4>0` partitioning: MMS or analytical test against Fortran for fdp under non-zero sorption.

---

## Required actions before LimnoTech review

Numbered with effort estimates (1 dev-day = 6 hours focused work). Dependencies marked.

### Code fixes (estimated total: 3-5 days)

1. **Algae wiring sweep (FloatingAlgae + BenthicAlgae).** Rewire 13 kinetic methods to read DEFAULTS keys instead of v2 kwargs. Includes `mu_max_20`, `kdp_20`, `krp_20`, `vsap`, `mub_max_20`, `kdb_20`, `krb_20`, `KsN`/`KsNb`, `KsP`/`KsPb`, `Ksb`. **Effort: 0.5-1.0 day.** Tests: 1 day.

2. **Algae formula bug fixes.** F5 (FL option-1 parenthesization), F14 (harmonic-mean zero-guard), B6 (Steele exponent sign). **Effort: 2 hours code, 0.5 day tests.**

3. **Nitrogen wiring sweep.** Drop phantom `ammonium_decay_nitrate` (#14); rewire 5 kinetic methods to NITROGEN_DEFAULTS (#15); fix `nitrate_uptake_floating_algae` to use dynamic NO3 fraction (#16); reconstruct `nitrate_uptake_benthic_algae` (#17). **Effort: 1.0 day.** Tests: 1 day. Depends on: BenthicAlgae exposing `balgae_no3_uptake_fraction` cache (~2 hours).

4. **Carbon DIC fixes.** Add CBOD oxidation source (#20); fix rca/rcb in 4 sites (#19); drop DOX-Monod from POC hydrolysis (#21). **Effort: 0.5 day.** Tests: 0.5-1 day.

5. **DOX fixes.** Fix rca/rcb in 4 sites (#22); add DOX-Monod attenuation to SOD sink (#23). **Effort: 0.5 day.** Tests: 0.5 day.

6. **fdp partitioning fix.** Change `0.000001` to `1.0E6` in `utils/partitioning.py:31` (#18). **Effort: 5 minutes code, 0.5 day for MMS test.**

### New Fortran-anchored regression tests (estimated total: 5-7 days)

7. Default-instantiation regression suite: `Process()` (no kwargs) vs Fortran-anchored numerical reference. One test class per Process (FloatingAlgae, BenthicAlgae, Nitrogen, Phosphorus, Carbon, DOX). **Effort: 3 days.**

8. Stoichiometric-ratio regression: Carbon DIC algal terms and DOX algal flux terms with explicit numerical references. **Effort: 0.5 day.**

9. Mass-balance closure tests: NH4+NO3 algal-uptake sum equals total N demand; OrgN/POC budgets close over 100-step integrations. **Effort: 1 day.**

10. Regime-coverage tests: hypoxic SOD attenuation; `kdpo4>0` partitioning; CBOD-on DIC source; default reaeration (kah_20_user choice). **Effort: 1-2 days.**

### Documentation updates (estimated total: 0.5-1 day)

11. Update `parameter_defaults_corrections.md` Sections 1.1, 1.3, 1.6, 2.8, 3.6 per Cat 5 #1-#6. **Effort: 1 hour.**

12. Add corrections doc entries for the rca/rcb fix, the algae wiring sweep, the Nitrogen wiring sweep, the fdp fix, the CBOD-DIC source addition, the SOD-DOX attenuation addition. **Effort: 1 hour.**

13. Reconcile `utils/reaeration.py` docstring author attributions with Fortran source comments. **Effort: 1 hour.**

14. Update README and migration notes to flag the v3-vs-Fortran default reaeration divergence (`kah_20_user`), the OrgN/OrgP `vson`/`vsop` defaults, the `vson_20` internal duplication. **Effort: 1-2 hours.**

15. Update LimnoTech review packet (Phase 8 deliverable) to reflect post-audit state. **Effort: 2-3 hours.**

### Reconciliation discussions with LimnoTech (estimated total: ~1 day spread over a week)

16. **lambdam=0.174 vs 0.0174.** Fortran has 0.174 (matches QUAL2K Table 6); v1/v3 have 0.0174. Likely v1 typo. Recommend 0.174 unless LimnoTech has independent calibration evidence.

17. **DIC unit reconciliation.** Long-standing v1/Fortran mol/L vs mg/L mismatch in dDIC/dt. v3 inherits the v1 form. Decide: defer to v3 1.x carbonate solver, or land a unit cleanup in 1.0.0.

18. **fdp partitioning bug** (latent at default). Confirm whether any v1 calibration LimnoTech relies on used `kdpo4>0`; if so, those results are unreliable and need re-running.

19. **Default reaeration choice.** v3 `kah_20_user=0` blocks default reaeration; Fortran `kah_20_user=1.0` produces meaningful default. Decide whether to align with Fortran or document the divergence.

20. **vsop=0.1 vs Fortran 0.01, SOD_20=1.0 vs Fortran 0.2.** Pick literature-anchored values and document.

### Total effort estimate

- **Critical-path code+tests:** 8-12 dev-days.
- **Documentation+reconciliation:** 2-3 dev-days.
- **Round-trip with reviewer/sponsor:** add 3-5 calendar days for clarifications.
- **Realistic timeline to actual review-ready: 2-3 weeks.**

---

## Recommendation

**The v3 NSM1 1.0.0 LimnoTech review packet I claimed in Phase 8 is not review-ready.** The audit found ~22 critical correctness defects across multiple Processes (algae wiring, rca/rcb stoichiometric ratios, Nitrogen wiring, phantom NH4 source, broken NO3-uptake split, missing CBOD-DIC source, missing SOD-DOX attenuation, fdp inherited bug). The existing parity tests pass because they call default-instantiated Processes with explicit kwargs that mask the wiring bugs and pass the same wrong stoichiometric ratios into both sides of the comparison. Default-instantiated v3 NSM1 simulations would diverge from v1 immediately on the first step — algae would have zero respiration/death/settling, NH4 would grow exponentially from a phantom source, and DIC/O2 algal coupling would be off by 100-1000x.

Independent of the bugs, the 7 sentinel-999 corrections, the SOD_tc/PAR architectural refactors, the resolved-Q decisions, and the clip_negative_state safety net are all sound and vindicated by the Fortran-coded defaults. The audit also surfaces several documentation defects in `parameter_defaults_corrections.md` that are cheap to fix.

**Realistic timeline to actual review-ready: 2-3 weeks.** Critical-path code fixes are 3-5 days; new Fortran-anchored regression tests (which the existing test suite needs anyway) are 5-7 days; doc updates and LimnoTech reconciliation discussions add 1-2 days plus calendar time for the reconciliation round-trip.

Recommend pausing the LimnoTech review hand-off, executing the action list above, then re-running this audit synthesis against the post-fix state before re-claiming review-ready.

---

## Appendix: Per-Process finding tally

| Process | Critical | Minor | Match | Notes |
|---|---|---|---|---|
| FloatingAlgae | 6 | 4 | 13 | F1-F4 (4 wiring bugs), F5 (parens), F14 (harmonic guard); F8/F9/F11 minor; F6/F7/F10/F12/F13/F15-F20 match. |
| BenthicAlgae | 7 | 1 | 8 | B1-B3 (3 wiring bugs), B6 (Steele sign), B7-B9 (3 inheritance leakage); B16 minor; B4/B5/B10-B15 match. |
| Nitrogen | 5 | 4 | 7 | N2 (phantom source), N4 (default-value defect), N12 (static fraction), N13 (3 structural defects), wiring sweep (5 methods); N1/N10/N11/N14 minor; N3/N5-N9/N15-N17 match. |
| Phosphorus | 1 | 2 | 11 | P6 (fdp bug, gated); P2 default-value, P14 informational; P1/P3-P5/P7-P13 match. |
| Carbon | 4 | 5 | 8 | C1 (rca/rcb 4 sites), C3 (missing CBOD), C4 (POC Monod), and the C1 root cause; C7-C11 minor; POC settling, POC mortality, DOC mortality, DOC oxidation, DOC-as-DOX-sink, DIC reaeration form, dPOC/dDOC structure match. |
| DOX | 4 | 4 | 6 | C1 (rca/rcb 4 sites), C2 (SOD Monod), C5 (default reaeration), C6 (salinity); C7/C12-C14 minor; atm reaeration form, nitrification flux, DOC sink, CBOD sink, dDOX/dt, cached-rates contract match. |
| POM | 0 | 2 | 5 | Fortran `TsedC` vs water-temperature deviation; Fortran `vb` annual unit; rest match. |
| CBOD | 0 | 2 | 3 | Sedimentation `1/depth` form; `ksbod_theta`; rest match. |
| Pathogen | 0 | 1 | 4 | Light decay `Fr_PAR` scaling absorbable; rest match. |
| N2 | 0 | 3 | 5 | N2sat mb→atm scalar; v3 denit-source extension; oxygen-weighted TDG deferred; rest match. |
| Alkalinity | 0 | 2 | 6 | Architectural routing (nitrif/denit through Nitrogen flux cache) v3-equivalent; rest match. |
| Utilities | 0 | 0 | 5 | All five utilities match Fortran (`kah_20`, `kaw_20`, `ka_tc`, `L`, `fdp` single-solid form); `SOD_tc`/`PAR` deliberate v3 refactors. |
| Parameter library | 0 | 5 | ~135 | `lambdam` (likely v1 flaw), `vson_20` (internal v3 inconsistency), `vsop`/`SOD_20`/`kah_20_user` (v3 vs Fortran value disagreements). 7 sentinel corrections all vindicated. |

**Total critical: 27** (some single defects span multiple sites; deduplicated to 22 unique fixes in Category 1).
**Total minor: 35.**
**Total match: ~120 enumerated blocks plus ~135 parameter entries.**
