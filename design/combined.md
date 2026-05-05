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
- v3 (`alkalinity.py:301-343`):
  `(r_alkaa * ap_uptake_fr_nh4 - r_alkan * (1 - ap_uptake_fr_nh4))
  * ap_growth * rca * EQ_TO_MG_CACO3` where `EQ_TO_MG_CACO3 = 50000` and
  `rca = self.AWc`.

Match. Stoichiometric ratio names differ (Fortran `ralkca/ralkcn`,
v1/v3 `r_alkaa/r_alkan`) but the values
`14/106/12/1000` and `18/106/12/1000` are identical.

### 22. Algal respiration source

- Fortran (`modAlkalinity.f90:79`):
  `Alk_ApRespiration = ralkca * rca(r) * ApRespiration * 50000`.
- v1 (`processes.py:3345-3361`):
  `ApRespiration * r_alkaa * 50000 * rca`.
- v3 (`alkalinity.py:345-360`):
  `ap_resp * self.r_alkaa * self.AWc * EQ_TO_MG_CACO3`.

Match.

### 23. Benthic algae growth and respiration

- Fortran (`modAlkalinity.f90:87-88`):
  `Alk_AbGrowth = Fb(r) * (ralkca * AbUptakeFr_NH4 - ralkcn *
  (1 - AbUptakeFr_NH4)) * rcb(r) * AbGrowth / depth * 50000`;
  `Alk_AbRespiration = Fb(r) * ralkca * rcb(r) * AbRespiration / depth
  * 50000`.
- v1 (`processes.py:3364-3410`): same form, with `1/depth` factor and
  `Fb` multiplication; uses `r_alkba`/`r_alkbn`.
- v3 (`alkalinity.py:362-418`): same form; uses
  `r_alkba`/`r_alkbn`/`BWc`/`Fb`. The `1/depth` divider is applied
  explicitly in `_benthic_algae_growth_alk_flux` and
  `_benthic_algae_respiration_alk_source`.

Match.

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
# v3 NSM1 Utilities + Parameter Library — Three-way audit

**Date:** 2026-05-05
**References:** v3 (streaming branch), v1 (`src/clearwater_modules`), Fortran (`/Users/todd/Downloads/NSM_comparison/NSM1/Source Files`).

## Summary

1. Utilities: 6 functions audited. **5 of 6 match Fortran exactly** (`kah_20`, `kaw_20`, `ka_tc`, `L`, `fdp`). `SOD_tc` is a deliberate v3 architectural refactor (pure Arrhenius). `PAR` is a deliberate v3 refactor (toggle moved to consumer).
2. Parameters: ~145 distinct parameter entries audited across 13 v3 groups. **All 7 critical sentinel-999 corrections are confirmed genuinely needed** — Fortran has the correct physical defaults (e.g., 0.1 m/d for `vsop`, 1013.25 hPa is standard atm), so v1 alone introduced the 999/2026.5 sentinels; v3 restored Fortran-aligned values.
3. **Likely v1 flaw not corrected in v3**: `lambdam=0.0174` in v1/v3 vs Fortran's `0.174` — 10x discrepancy in POM contribution to Beer-Lambert light extinction. Needs reconciliation.
4. v3 deliberate improvements verified: `pressure_mb=1013.25`, sentinel rescue, `PAR` toggle inversion, `SOD_tc` pure-Arrhenius split. All consistent with `parameter_defaults_corrections.md` Sections 1-3.
5. Internal v3 inconsistency: `vson_20=0.1` in `parameters/nitrogen.py` vs `vson=0.01` in `parameters/global_vars.py` (Fortran uses 0.01). Either consolidate or document.

---

## Part 1 — Utility modules

### `reaeration.py`

#### `kah_20` (hydraulic, 9 options)

Fortran source: `modGlobalParam.f90:268-339`, subroutine `O2Reaeration`. v1 source: `shared/processes.py:65-98`. v3 source: `utils/reaeration.py:26-119`.

| Option | Fortran formula | v1 formula | v3 formula | Status |
|---|---|---|---|---|
| 1 | user `kah%rc20` | `kah_20_user` | `kah_20_user` | match |
| 2 | `(3.93 v^0.5)/h^1.5` (line 275) | identical | identical | match |
| 3 | `(5.32 v^0.67)/h^1.85` (line 281) | identical | identical | match |
| 4 | `5.026 v / h^1.67` (line 287) | identical | identical | match |
| 5 | depth-piecewise: `<0.61` Owens, `>0.61` O'Connor, `=0.61` Churchill (lines 289-306) | identical with `==0.61` Churchill branch | identical | match |
| 6 | flow-piecewise (lines 308-315), `<0.556` → `517 (vS)^0.524 Q^-0.242`, else `596 (vS)^0.528 Q^-0.136` | identical | identical | match |
| 7 | flow-piecewise (lines 317-324), `<0.556` → `88 (vS)^0.313 h^-0.353`, else `142 (vS)^0.333 h^-0.66 W^-0.243` | identical | identical | match |
| 8 | flow-piecewise (lines 326-333), `<0.425` → `31183 vS`, else `15308 vS` | identical | identical | match |
| 9 | Froude form: `2.16 (1 + 9 Fr^0.25) u*/h` with `Fr = v/sqrt(g h)` (lines 334-338) | identical | identical | match |

**Finding:** match across all three references. The docstring author attributions in v3 (Covar/Owens-Gibbs/Churchill/Tsivoglou-Wallace/Padden-Gloyna/USGS pool-and-riffle/Thackston-Krenkel/Langbien-Durum) disagree with Fortran's inline comments (Owens, O'Connor, Churchill, Cover, Melching-Flores, Tsivoglou-Neal, Thackson-Dawson). The author attributions differ but the formulas are identical. **Severity:** minor documentation inconsistency.

**v3 implementation detail (Phase 5.5):** `np.select` dim-stripping fix at `utils/reaeration.py:107-119` — match-preserving.

#### `kaw_20` (wind, 13 options)

Fortran source: `modGlobalParam.f90:341-414`. v1 source: `shared/processes.py:117-146`. v3 source: `utils/reaeration.py:122-204`.

Fortran labels its 13 options as: 1=user, 2=Broecker, 3=Gelda, 4=Banks-Herrera, 5=Wanninkhof, 6=Cole-Buchak, 7=Banks, 8=Smith, 9=Liss, 10=Downing-Truesdale, 11=Kanwisher, 12=Yu, 13=Weiler.

v1 and v3 match byte-for-byte numerically:

| Option | Fortran formula | v1 / v3 formula | Status |
|---|---|---|---|
| 1 | user `kaw%rc20` | `kaw_20_user` | match |
| 2 | `0.864 Uw10` | identical | match |
| 3 | `Uw10<=3.5: 0.2 Uw10; else 0.057 Uw10^2` | identical | match |
| 4 | `0.728 Uw10^0.5 - 0.317 Uw10 + 0.0372 Uw10^2` | identical | match |
| 5 | `0.0986 Uw10^1.64` | identical | match |
| 6 | `0.5 + 0.05 Uw10^2` | identical | match |
| 7 | `Uw10<=5.5: 0.362 sqrt(Uw10); else 0.0277 Uw10^2` | identical | match |
| 8 | `0.64 + 0.128 Uw10^2` | identical | match |
| 9 | `Uw10<=4.1: 0.156 Uw10^0.63; else 0.0269 Uw10^1.9` | identical | match |
| 10 | `0.0276 Uw10^2` | identical | match |
| 11 | `0.0432 Uw10^2` | identical | match |
| 12 | `0.319 Uw10` | identical | match |
| 13 | `Uw10<1.6: 0.398; else 0.155 Uw10^2` | identical | match |

`Uw10 = wind_speed * (10/2)^0.143` matches in all three.

**Finding:** match across all three on numerics. v3 docstring author attributions disagree with Fortran source comments (e.g., v3 calls option 2 "Banks 1975" while Fortran credits Broecker 1978). Recommend reconciling author attributions before LimnoTech review. **Severity:** observation.

#### `ka_tc` (combined, temperature-corrected)

Fortran: `modGlobalParam.f90:245-247`, `kah_tc + kaw_tc/depth`. v1: `shared/processes.py:165-178`. v3: `utils/reaeration.py:207-236`.

**Finding:** match. All three apply `Arrhenius_TempCorrection` to each component before summing. **Severity:** match.

### `sediment.py`

#### `SOD_tc`

Fortran: `modGlobalParam.f90:250-256`. v1: `shared/processes.py:180-200`. v3: `utils/sediment.py:16-31`.

- Fortran: `SOD_tc = Arrhenius(SOD%rc20, theta, T)`; if `use_DOX`, multiplies by `DOX/(DOX + KsSod)`.
- v1: same: `arrhenius_correction(...) * xr.where(use_DOX, DOX/(DOX+KsSOD), 1)`.
- v3: pure Arrhenius only. The DOX-Monod factor moved to the DOX Process call site.

**Finding:** v3 deliberate architectural refactor, documented in `parameter_defaults_corrections.md` Section 3.2. Numerically equivalent under matched fixtures (`tests/test_5_dox_calculations_v2.py::test_dox_sod_rate_matches_v1` passes `use_DOX=False` to v1 to obtain parity). **Severity:** intended improvement, not a deviation.

### `light.py`

#### `L` (Beer-Lambert extinction)

Fortran source: `modGlobalParam.f90:420-428`, `LightExtCoefficient`:

```
lambda = lambda0
do i=1,nGS:  lambda += lambdas * Solid(i)
if use_POC:  lambda += lambdam * POC / focm
if use_Algae: lambda += lambda1 * Ap + lambda2 * Ap^0.66667
```

v1 (`shared/processes.py:202-237`): identical, applies `lambdas * Solid` always (matching Fortran), then adds POC term inside `xr.where(use_POC,...)`, then adds algae term inside `xr.where(use_Algae,...)`.

v3 (`utils/light.py:13-53`): identical to v1.

**Finding:** match across all three. **Note:** `parameter_defaults_corrections.md` Section 2.8 incorrectly claims v1's `lambdas * Solid` is "commented out / defined but not used"; in fact `shared/processes.py:232` applies it unconditionally. v3 reproduces this. Recommend correcting the corrections doc Section 2.8.

Multi-solid difference: Fortran loops over `nGS` solid groups summing `lambdas * Solid(i)`. v3/v1 take a single scalar Solid concentration. For single-class use the formulas agree.

#### `PAR`

Fortran (`modGlobalParam.f90:222,234`):
```
real(R8) :: Fr_PAR = 0.47
if (use_Algae .or. use_BAlgae) PAR = q_solar * Fr_PAR
```

v1 (`shared/processes.py:240-253`): wraps with `xr.where(use_Algae or use_Balgae, q_solar*Fr_PAR)`. **Latent v1 bug:** `xr.where` with only two args returns NaN in the false branch.

v3 (`utils/light.py:56-70`): returns `q_solar * Fr_PAR` unconditionally; the `use_Algae/use_BAlgae` toggle moved to the consumer Process per Phase 1.1.

**Finding:** v3 deliberate refactor. Numerically equivalent inside the `use_Algae|use_BAlgae` branch. Avoids latent v1 NaN-propagation bug. **Severity:** v3 deliberate improvement.

### `partitioning.py`

#### `fdp`

Fortran (`modGlobalParam.f90:225-231`):
```
fdp = 1.0
do i=1,nGS:  fdp = fdp + kdpo4(i,r) * Solid(i) / 1.0E6
fdp = 1.0 / fdp
```

v1 (`shared/processes.py:256-271`): `xr.where(use_TIP, 1/(1 + kdpo4*Solid/0.000001), 0)` — single solid class.

v3 (`utils/partitioning.py:12-31`): identical to v1, with `0.000001` literal.

**Finding:** match (under single solid class). The `1.0E6` factor in Fortran (and `1e-6` denominator in v1/v3) is the unit conversion `(L/kg)(mg/L)(1 kg / 1e6 mg) = dimensionless`. The Phase 1.1 audit comment characterizing this as "suspicious" was incorrect — the formula is dimensionally consistent. All three references agree. **Severity:** match. Recommend retracting the Phase 1.1 "suspicious unit factor" flag.

Multi-solid difference: Fortran sums over `nGS` solid groups with possibly different `kdpo4(i,r)` per class. v3/v1 collapse to single scalar. If `nGS > 1` is ever activated, v3 would need a sum-over-classes form to match Fortran.

### `numerics.py` (v3-only)

`Diagnostics` dataclass + `clip_negative_state` at `utils/numerics.py:28-118`. No Fortran or v1 counterpart. Per Q7 in the design spec (Section 14), this is v3-only architectural infrastructure.

**Finding:** v3-only by design. Logic correct: clip target is exactly 0 (line 76), matching the Q7 contract. Detail-limit-per-call (default 10) properly rate-limits log records (line 84). Aggregate suppressed-count stub appended when n_clipped > limit (line 100). Returns DataArray with preserved coords/dims/attrs (line 112). **Severity:** intended addition, no issues.

---

## Part 2 — Parameter library

### Group: `algae` (`parameters/algae.py`, 17 entries)

All 17 v3 algae defaults match v1 exactly (`AWd=100`, `AWc=40`, `AWn=7.2`, `AWp=1`, `AWa=1000`, `KL=10`, `KsN=0.04`, `KsP=0.0012`, `mu_max_20=1`, `kdp_20=0.15`, `krp_20=0.2`, `mu_max_theta=1.047`, `kdp_theta=1.047`, `krp_theta=1.047`, `vsap=0.15`, `growth_rate_option=1`, `light_limitation_option=1`).

### Group: `balgae` (`parameters/balgae.py`, 19 entries)

All 19 entries in v3 match v1 byte-for-byte (`BWd=100`, `BWc=40`, `BWn=7.2`, `BWp=1`, `BWa=3500`, `KLb=10`, `KsNb=0.25`, `KsPb=0.125`, `Ksb=10`, `mub_max_20=0.4`, `krb_20=0.2`, `kdb_20=0.3`, `mub_max_theta=1.047`, `krb_theta=1.06`, `kdb_theta=1.047`, `b_growth_rate_option=1`, `b_light_limitation_option=1`, `Fw=0.9`, `Fb=0.9`).

### Group: `nitrogen` (`parameters/nitrogen.py`, 16 entries)

| Parameter | Fortran | v1 | v3 | Status |
|---|---|---|---|---|
| KNR | (modNitrogen) | 0.6 | 0.6 | match v1 |
| knit_20 | 0.1 | 0.1 | 0.1 | match v1 |
| kon_20 | 0.1 | 0.1 | 0.1 | match v1 |
| kdnit_20 | 0.002 | 0.002 | 0.002 | match v1 |
| rnh4_20 | 0 | 0 | 0.0 | match v1 — FIXME-flagged |
| vno3_20 | 0 | 0 | 0.0 | match v1 — FIXME-flagged |
| **vson_20** | **0.01** (in modGlobalParam.f90:92) | not in nitrogen group; v1 GlobalVars has `vson=0.01` | **0.1** | **structural relocation + 10x value change; undocumented** |
| knit_theta | 1.083 | 1.083 | 1.083 | match v1 |
| kon_theta | 1.074 | 1.074 | 1.074 | match v1 |
| kdnit_theta | 1.08 | 1.08 | 1.08 | match v1 |
| rnh4_theta | 1.047 | 1.047 | 1.047 | match v1 |
| vno3_theta | 1.045 | 1.045 | 1.045 | match v1 |
| **vson_theta** | n/a | not in v1 | **1.024** | **v3 addition, undocumented** |
| KsOxdn | 0.1 | 0.1 | 0.1 | match v1 |
| PN | 0.5 | 0.5 | 0.5 | match v1 |
| PNb | 0.5 | 0.5 | 0.5 | match v1 |
| use_OrgN | True | True (in GlobalParameters) | True | match (relocated) |

**Findings:**
- `vson_20=0.1` in v3 vs `vson=0.01` in v1 GlobalVars and `vson=0.01` in v3's own global_vars (`global_vars.py:26`). The v3 nitrogen group's `vson_20=0.1` is **10x v1's value, 10x v3's own global_vars value, and 10x Fortran's value**. Categorize as **undocumented v3 deviation needing review**.
- `vson_theta=1.024`: also new in v3, not in v1.
- Fortran `vson` (modGlobalParam.f90:92) initializes to `0.01`, matching v1's GlobalVars value.

### Group: `phosphorus` (`parameters/phosphorus.py`, 7 entries)

| Parameter | Fortran modGlobalParam | v1 | v3 | Status |
|---|---|---|---|---|
| kop_20 | (modPhosphorus) | 0.1 | 0.1 | match v1 |
| rpo4_20 | 0 | 0 | 0.0 | match v1 — FIXME-flagged |
| kop_theta | 1.047 | 1.047 | 1.047 | match v1 |
| rpo4_theta | 1.074 | 1.074 | 1.074 | match v1 |
| kdpo4 | line 82: `kdpo4 = 0.0` | 0.0 | 0.0 | match all three — FIXME-flagged |
| **vsop** | line 98: **0.01** | **999** | **0.1** | **likely v1 flaw; v3 corrects but to 10x Fortran's 0.01** |
| **vs** | line 87: **0.1** | **999** | **0.1** | **likely v1 flaw; v3 matches Fortran exactly** |

### Group: `carbon` (`parameters/carbon.py`, 10 entries)

All 10 v3 carbon parameters match v1 exactly: `f_pocp=0.9`, `kdoc_20=0.01`, `kdoc_theta=1.047`, `f_pocb=0.9`, `kpoc_20=0.005`, `kpoc_theta=1.047`, `KsOxmc=1.0`, `pCO2=383.0`, `FCO2=0.2`, `roc=32/12`.

### Group: `cbod` (`parameters/cbod.py`, 5 entries)

5 of 5 match v1 (`KsOxbod=0.5`, `kbod_20=0.12`, `ksbod_20=0.0`, `kbod_theta=1.047`, `ksbod_theta=1.047`).

### Group: `dox` (`parameters/dox.py`, 10 entries)

| Parameter | Fortran modGlobalParam | v1 GlobalVars | v3 dox | Status |
|---|---|---|---|---|
| ron | n/a | 2*32/14 = 4.5714 | 2*32/14 | match v1 |
| KsSOD | line 127: `KsSod = 1.0` | 1 | 1.0 | match all three |
| **SOD_20** | line 122: **0.2** | **999** | **1.0** | **all-three disagreement: Fortran=0.2, v1=999, v3=1.0** |
| **SOD_theta** | line 122: **1.06** | **999** | **1.060** | **likely v1 flaw; v3 matches Fortran** |
| **kaw_20_user** | line 117: **0.0** | **999** | **0.0** | **likely v1 flaw; v3 matches Fortran** |
| **kah_20_user** | line 113: **1.0** | **999** | **0.0** | **all-three disagreement: Fortran=1.0, v1=999, v3=0.0** |
| kaw_theta | line 117: 1.024 | 1.024 | 1.024 | match all three |
| kah_theta | line 113: 1.024 | 1.024 | 1.024 | match all three |
| hydraulic_reaeration_option | line 143: 1 | 1 | 1 | match all three |
| wind_reaeration_option | line 147: 1 | 1 | 1 | match all three |

**Critical findings:**
- **SOD_theta**: v3's correction to 1.060 vindicated by Fortran (Chapra 1997 standard). **Likely v1 flaw**.
- **SOD_20**: Fortran 0.2, v3 chose 1.0 (5x Fortran's value). Both defensible from literature; v3 deviates from Fortran by 5x. Recommend documenting.
- **kaw_20_user**: v3 matches Fortran exactly. **Likely v1 flaw, v3 correctly restored.**
- **kah_20_user**: Fortran 1.0, v3 0.0. **Behavioral consequence:** at default `hydraulic_reaeration_option=1`, v3 reaeration = 0, Fortran reaeration = 1.0 1/d. v3 makes user-override path explicit; Fortran has a non-zero hidden default. v3 + Fortran will produce different DOX trajectories at default settings.

### Group: `pathogen` (`parameters/pathogen.py`, 4 entries)

`kdx_20=0.8`, `kdx_theta=1.07`, `apx=1`, `vx=1`. All 4 match v1. `apx` and `vx` FIXME-flagged for unknown literature basis.

### Group: `alkalinity` (`parameters/alkalinity.py`, 6 entries)

All 6 stoichiometric ratios match v1 exactly.

### Group: `n2` (`parameters/n2.py`, 0 entries)

Empty in v3 and v1. Match.

### Group: `pom` (`parameters/pom.py`, 3 entries)

| Parameter | Fortran modGlobalParam | v1 | v3 | Status |
|---|---|---|---|---|
| kpom_20 | n/a | 0.1 | 0.1 | match v1 |
| h2 | line 134: `h2 = 0.1` | 0.1 | 0.1 | match all three |
| kpom_theta | n/a | 1.047 | 1.047 | match v1 |

### Group: `global_parameters` (`parameters/global_parameters.py`, 17 entries)

All 16 boolean `use_*` flags match v1 (use_NH4, use_NO3, use_OrgN, use_OrgP, use_TIP, use_SedFlux=False, use_POC, use_DOC, use_DOX, use_DIC, use_Algae, use_Balgae, use_N2, use_Pathogen, use_Alk, use_POM all True except SedFlux).

Fortran defaults differ: `use_BAlgae=.false.`, `use_POC=.false.`, `use_DOC=.false.`, `use_DIC=.false.`, `use_N2=.false.`, `use_Pathogen=.false.`, `use_Alk=.false.`, `use_POM2=.false.`. v3 enables all by default (matching v1's "all on" stance), differing from Fortran's selective-enable approach. **Inherited v1 convention**, not a deviation in v3 from v1.

| Parameter | Fortran | v1 GlobalVars | v3 | Status |
|---|---|---|---|---|
| **pressure_mb** | n/a in NSM1 (Fortran uses pressure_atm) | **2026.5** | **1013.25** | **likely v1 flaw, v3 correctly restored to ISO standard** |

### Group: `global_vars` (`parameters/global_vars.py`, 21 entries)

| Parameter | Fortran modGlobalParam | v1 GlobalVars | v3 | Status |
|---|---|---|---|---|
| vson | line 92: 0.01 | 0.01 | 0.01 | match all three |
| vsoc | line 104: 0.01 | 0.01 | 0.01 | match all three |
| theta | n/a | 1.047 | 1.047 | match v1 |
| **vb** | line 138: **0.0025 m/yr** (line 201 divides by 365) | **0.01 m/d** | **0.01 m/d** | **Fortran/v1 unit reconciliation issue**; FIXME-flagged in v3 |
| fcom (Fortran focm) | line 108: `focm = 0.4` | 0.4 | 0.4 | match all three |
| dt | n/a | 1 | 1.0 | match v1 |
| depth | n/a (runtime) | 1.5 | 1.5 | match v1 |
| TwaterC | n/a | 20 | 20.0 | match v1 |
| velocity | n/a | 1 | 1.0 | match v1 |
| flow | n/a | 2 | 2.0 | match v1 |
| topwidth | n/a | 1 | 1.0 | match v1 |
| slope | n/a | 2 | 2.0 | match v1 |
| shear_velocity | n/a | 4 | 4.0 | match v1 |
| wind_speed | n/a | 4 | 4.0 | match v1 |
| q_solar | n/a | 500 | 500.0 | match v1; FIXME for unit-doc mismatch |
| Solid | n/a | 1 | 1 | match v1 |
| lambda0 | line 60: 0.02 | 0.02 | 0.02 | match all three |
| lambda1 | line 72: 0.0088 | 0.0088 | 0.0088 | match all three |
| lambda2 | line 76: 0.054 | 0.054 | 0.054 | match all three |
| lambdas | line 64: 0.052 | 0.052 | 0.052 | match all three |
| **lambdam** | line 68: **0.174** | **0.0174** | **0.0174** | **likely v1 flaw — Fortran has 0.174, v1/v3 have 0.0174 (10x lower)** |
| Fr_PAR | line 222: 0.47 | 0.47 | 0.47 | match all three |

**Critical finding (lambdam):** Fortran `modGlobalParam.f90:68` initializes `lambdam = 0.174` while v1 and v3 use `0.0174` — 10x discrepancy in POM contribution to Beer-Lambert light extinction. Fortran's 0.174 is consistent with QUAL2K Table 6 references; v1's 0.0174 may be a typo. **Likely v1 flaw, propagated to v3. Requires reconciliation with LimnoTech before review.**

---

## Critical-correction verification

For each of the 7 corrections in `parameter_defaults_corrections.md` Section 1:

| Correction | v1 sentinel | v3 fix | Fortran value | Verdict |
|---|---|---|---|---|
| 1.1 vsop | 999 | 0.1 | **0.01** | v1 flaw confirmed; v3 corrects but to 10x Fortran's |
| 1.2 vs | 999 | 0.1 | **0.1** | v1 flaw confirmed; v3 matches Fortran exactly |
| 1.3 SOD_20 | 999 | 1.0 | **0.2** | v1 flaw confirmed; v3 corrects but to 5x Fortran's |
| 1.4 SOD_theta | 999 | 1.060 | **1.060** | v1 flaw confirmed; v3 matches Fortran exactly |
| 1.5 kaw_20_user | 999 | 0.0 | **0.0** | v1 flaw confirmed; v3 matches Fortran exactly |
| 1.6 kah_20_user | 999 | 0.0 | **1.0** | v1 flaw confirmed; v3 zeroes user-override branch but disagrees with Fortran's 1.0 default |
| 1.7 pressure_mb | 2026.5 | 1013.25 | **n/a** (Fortran uses pressure_atm) | v1 flaw confirmed (2x); v3 restored to ISO 2533 standard |

**Summary:** All 7 v3 corrections are vindicated as genuine v1 flaws. For 4 of 7 (`vs`, `SOD_theta`, `kaw_20_user`, `pressure_mb`) v3 chose a value matching or equivalent to Fortran. For 3 of 7 (`vsop`, `SOD_20`, `kah_20_user`) v3 chose a value differing from Fortran by O(1)-O(10) but defensible from literature.

---

## Conclusions

### v3 deliberate improvements (confirmed correct)

1. `pressure_mb=1013.25` (Section 1.7): v1 had 2026.5 (2x error). v3 restored to ISO standard.
2. `SOD_theta=1.060` (Section 1.4): v1 had sentinel 999. v3 matches Fortran's 1.06 and Chapra (1997).
3. `vs=0.1`, `kaw_20_user=0.0` (Sections 1.2, 1.5): v1 had 999. v3 restored Fortran-aligned defaults.
4. `vsop=0.1` (Section 1.1): v1 had 999. v3 chose 0.1 (10x Fortran's 0.01); literature-defensible.
5. `SOD_tc` pure Arrhenius split (Section 3.2): cleaner architecture; DOX-Monod moved to consumer Process.
6. `PAR` toggle inversion (Section 3.4 / Phase 1.1): avoids latent v1 NaN-on-disable bug from `xr.where(cond, value)` two-arg form.
7. `clip_negative_state` + `Diagnostics` (Q7): v3-only safety net, not a deviation.
8. `np.select` dim-stripping fix in `kah_20`/`kaw_20`: v3 implementation detail correctly preserves xarray broadcasting.

### Likely v1 flaws (v3 correctly bypassed; not v3 issues)

1. The 7 sentinel-999 / 2026.5 defaults — all confirmed by Fortran-coded defaults to be flaws in v1.
2. v1 `PAR` two-arg `xr.where` returns NaN when both algae modules disabled (Fortran/v3 avoid).

### Likely v1 flaw NOT corrected in v3 (needs action)

1. **`lambdam=0.0174`** in v1/v3 vs Fortran's `0.174` — 10x discrepancy in POM contribution to Beer-Lambert light extinction. Fortran's value consistent with QUAL2K Table 6; v1 likely has a typo that propagated to v3. **Reconcile with LimnoTech before review.**

### Undocumented v3 deviations (flag for review)

1. **`vson_20=0.1`** in `parameters/nitrogen.py:16`: v1 GlobalVars has `vson=0.01`, v3 global_vars also has `vson=0.01`, Fortran has 0.01. The new `vson_20=0.1` in nitrogen group is 10x v1's value and 10x v3's own global_vars value. **Internal v3 inconsistency.** Recommend either consolidate to 0.01 (matching Fortran/v1) or document the 0.1 choice with rationale.
2. **`vson_theta=1.024`** in `parameters/nitrogen.py:22`: not present in v1. New v3 parameter; document.
3. **`vsop=0.1` (Section 1.1)**: v3 chose 0.1 over Fortran's 0.01. Acknowledged in corrections doc but Fortran value not noted.
4. **`SOD_20=1.0` (Section 1.3)**: v3 chose 1.0 over Fortran's 0.2. Document.
5. **`kah_20_user=0.0` (Section 1.6)**: v3 chose 0.0 over Fortran's 1.0. Behavioral consequence: at default `hydraulic_reaeration_option=1`, v3 reaeration = 0, Fortran reaeration = 1.0 1/d. **Document this divergence prominently** — could surface as "v3 has no DOX recovery" puzzle in side-by-side runs.

### Documentation defects in `parameter_defaults_corrections.md`

1. Section 2.8 claims v1's `lambdas * Solid` term is "commented out / defined but not used"; the actual v1 source (`shared/processes.py:232`) applies it unconditionally. Recommend correcting.
2. Phase 1.1's "suspicious unit factor `1/(1 + kdpo4 * Solid / 0.000001)`" comment for `fdp`: the formula is dimensionally consistent. Recommend retracting the suspicion.
3. Section 3.6/3.7 (Kelvin offset, mb-to-atm scaling): v3's `utils/conversions.py` re-exports from v2, and v2's `celsius_to_kelvin` returns `T_C + 273.16` (with comment "for testing consistency with v1"), **not** 273.15 as Section 3.6 claims. Recommend updating Section 3.6 to reflect actual v3 behavior.

### Required actions before LimnoTech review

1. Resolve the **lambdam 0.174 vs 0.0174** discrepancy (Fortran vs v1/v3 10x mismatch).
2. Resolve the **vson_20 nitrogen group 0.1 vs global_vars 0.01** internal v3 inconsistency.
3. Update `parameter_defaults_corrections.md` Sections 1.1, 1.3, 1.6 to record the Fortran-coded default and the rationale for v3's chosen value where it differs from Fortran.
4. Correct `parameter_defaults_corrections.md` Section 2.8 (lambdas not commented out) and Section 3.6 (Kelvin offset).
5. Reconcile docstring author attributions in `utils/reaeration.py` (option labels disagree with Fortran source comments).

---

## Finding count

- **Utility findings:** 5 matches, 2 deliberate v3 improvements (PAR, SOD_tc), 1 v3-only addition (numerics), 0 undocumented deviations. Documentation issues: 3.
- **Parameter findings across ~145 entries:** ~135 match v1, 7 critical sentinel corrections (all vindicated), 1 likely v1 flaw not corrected (`lambdam`), 5 undocumented v3 deviations or value choices needing rationale (`vson_20`, `vson_theta`, `vsop` value vs Fortran, `SOD_20` value vs Fortran, `kah_20_user` value vs Fortran).
