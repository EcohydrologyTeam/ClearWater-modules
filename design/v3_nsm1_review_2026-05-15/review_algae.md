# v3 NSM1 Algae Code Review (Floating + Benthic)

Review date: 2026-05-15
Reviewer: water-quality-model source-code and science-correctness reviewer
Branch: `streaming`
Repo: `/Users/todd/GitHub/ecohydrology/ClearWater-modules-streaming`

Scope (v3 files reviewed line by line):

- `src/clearwater_modules_v3/processes/floating_algae.py`
- `src/clearwater_modules_v3/processes/benthic_algae.py`
- `src/clearwater_modules_v3/parameters/algae.py`
- `src/clearwater_modules_v3/parameters/balgae.py`

Cross-read for parity:

- `src/clearwater_modules/nsm1/processes.py` (v1 function-style reference: `rna`/`rpa`/`rca`/`rda`, `FL`/`FN`/`FP`, `mu`, `ApGrowth`/`ApRespiration`/`ApDeath`/`ApSettling`, `dApdt`, `FLb`/`FNb`/`FPb`/`FSb`, `mub`, `AbGrowth`/`AbRespiration`/`AbDeath`, `dAbdt`, `ApUptakeFr_NH4`/`AbUptakeFr_NH4`, `NH4_ApRespiration`/`NH4_ApGrowth`/`NH4_AbRespiration`/`NH4_AbGrowth`, `ApDeath_OrgN`/`AbDeath_OrgN`/`ApDeath_OrgP`/`AbDeath_OrgP`, `POC_algal_mortality`/`DOC_algal_mortality`/`POC_benthic_algae_mortality`/`DOC_benthic_algae_mortality`, `POM_algal_settling`/`POM_benthic_algae_mortality`)
- `src/clearwater_modules/nsm1/constants.py` (v1 `DEFAULT_ALGAE`, `DEFAULT_BALGAE`, `f_pocp`/`f_pocb`/`PN`/`PNb`)
- `src/clearwater_modules_v3/utils/conversions.py` (`arrhenius_correction`)
- `src/clearwater_modules_v3/utils/partitioning.py` (`fdp`)
- `src/clearwater_modules_v3/utils/numerics.py` (`clip_negative_state`, `sanitize_rate`)
- `src/clearwater_modules_v3/processes/base.py` (`Process`, `time_step`)

Audit references consulted:

- `design/clearwater_modules_v3_nsm1_audit_algae.md` (dated 2026-05-05; audited the pre-fix `clearwater_modules_v2` tree)
- `design/clearwater_modules_v3_nsm1_audit_summary.md`
- `design/clearwater_modules_v3_nsm1_design_specification.md` (Section 6 bug list)
- `src/clearwater_modules_v3/parameter_defaults_corrections.md`

Tests located and inspected (not re-run; benchmarks out of scope per instruction):

- `tests/v3/nsm1/test_phase9a1_algae_wiring.py` (default-instantiation wiring regression)
- `tests/v3/nsm1/test_floating_algae_v1_parity_v3.py`, `test_benthic_algae_v1_parity_v3.py`

---

## 1. Summary verdict

The v3 Algae code is in good shape on the dimensions this review targeted. Every bug on the original "Bug #N" list (Bug #4 broken multiplicative integrator, Bug #13 `ammonium_respiration` returning 0, Bug #14 `ammonium_growth` returning 0, Bug #15 hard-coded `phosphate_fraction_dissolved=0.5`, Bug #16 missing `set_at_time` persistence) is verifiably fixed in current code, and so are the three formula bugs the algae audit flagged independently (F5 `limit_light` option-1 parenthesization, F14 harmonic-mean zero-guard, B6 Steele-exponent sign). The systematic "wiring" defect class identified in the 2026-05-05 audit (F1-F4, B1-B3, B7-B9 -- DEFAULTS keys silently shadowed by legacy kwargs) has been resolved by the Phase 9.A.1 legacy-kwarg-to-DEFAULTS bridge in both `__init__` methods, and a dedicated regression test (`test_phase9a1_algae_wiring.py`) asserts default-instantiated processes read the DEFAULTS values.

The xarray refactoring is substantially complete. The kinetic helpers are written in vectorized `xr.where` form and the IEEE-754 `== np.nan` defects are corrected to `.isnull()` / `np.isnan`. Three residual array-truthiness / scalar-broadcast patterns remain (the `if not self.use_phosphate: return 1.0` early returns, the `if isinstance(...)` scalar/array branching in `_ap_uptake_fr_nh4`, and a scalar `return 0.5` / `return 1.0` from limitation methods); these are correctness-safe under current single-flag usage but are latent multi-cell broadcasting hazards and are recorded as findings.

The most material remaining issues are documentation, not computation. The `floating_algae.py` module header (lines 1-46) and the `benthic_algae.py` module header (lines 1-20) are written as forward-looking work plans ("apply the FloatingAlgae bug fixes", "implement `ammonium_respiration` (was returning 0)", "replace the broken multiplicative integrator"). Because the work is done, these headers are stale and actively misleading to a future reader, who would reasonably conclude the bugs are still open. Separately, the `init_process` stub comment `# TODO: implement` at `floating_algae.py:376` is stale (the line it annotates is implemented), and the v3 inline default `f_pocp = 0.5` / `f_pocb = 0.5` silently changes the POC-vs-DOC algal-mortality carbon split from v1's documented 90/10 to 50/50 without a corrections-doc entry.

Findings by severity: 0 CRITICAL, 2 MAJOR, 7 MINOR, 6 OBSERVATION.

Overall confidence: high for the bug-fix verification and algorithm-parity conclusions (every claim is traced to specific lines in both v3 and v1); medium for the multi-cell broadcasting findings (classified MINOR/OBSERVATION because they are not reproduced end-to-end -- benchmarks were out of scope -- and depend on registry-supplied array shapes the review did not exercise).

## 2. Findings table

| ID | Severity | File:line | Category | Description | Recommended fix |
|---|---|---|---|---|---|
| A1 | MAJOR | `floating_algae.py:1-46` | documentation / stale-comment | Module docstring is written as a forward-looking work plan ("apply the FloatingAlgae bug fixes", "Bug #13: implement `ammonium_respiration` (was returning 0)", "Bug #4: replace the broken multiplicative integrator"). All listed fixes are present in current code (verified line by line, see Section 3). A reader of the header would conclude the module is still broken. | Rewrite the header in past tense / "as-implemented" form: state that Bug #4/#13/#14/#15/#16 and the NaN-guard defects are corrected, with the current line numbers, and move the historical narrative to a CHANGELOG or the closeout docs. |
| A2 | MAJOR | `floating_algae.py:111-117`, `benthic_algae.py:78-84` | algorithm-parity / open-issue | `_FDP_DEFAULTS["f_pocp"] = 0.5` and `_BENTHIC_FDP_DEFAULTS["f_pocb"] = 0.5`. v1 `constants.py:157,159` set `f_pocp = 0.9`, `f_pocb = 0.9`. This silently shifts ~40% of algal/benthic mortality carbon from POC to DOC versus v1, and is not recorded in `parameter_defaults_corrections.md` (grep confirms no `f_pocp`/`f_pocb` entry). Either an undocumented numeric deviation or a defect. | Either set the inline fallback to `0.9` to match v1/Fortran, or add a numbered entry to `parameter_defaults_corrections.md` justifying `0.5` with literature, mirroring the `mu_max_20`/`kdp_20`/`krp_20` precedent in `algae.py`. |
| A3 | MINOR | `floating_algae.py:370-377` | stale-comment | `# TODO: implement` at line 376 sits directly above `self.use_phosphate = True`, which is implemented (the flag is set and consumed by `limit_phosphorus`). The "is there a phosphorus process" sibling-discovery is hardcoded `True`, matching the pattern; the TODO no longer describes an unimplemented line. Contrast `benthic_algae.py:253-276`, which implements the same concept via `model.has_process`. | Remove the `# TODO: implement` comment, or replace it with the same sibling-discovery pattern used in `BenthicAlgae.init_process` (`model.has_process("Phosphorus")`) and an accurate comment. |
| A4 | MINOR | `floating_algae.py:735-736`, `757-764`, `benthic_algae.py:565-570` | xarray | `limit_phosphorus`/`limit_nitrogen`/`_ap_uptake_fr_nh4`/`_ab_uptake_fr_nh4` early-return Python scalars (`return 1.0`, `return 0.5`, `return 0.0`) when use-flags are off. v1 returns a per-cell array of the same shape (`np.select`/`xr.where`). A scalar in a downstream `xr.where`/array sum broadcasts and is correct numerically, but loses the DataArray coords/dims and can surface as a shape or alignment surprise when a sibling Process expects a coord-bearing array. Latent, not currently triggered. | Return `xr.ones_like(concentration)` / `xr.full_like(..., 0.5)` (guarding the scalar-input unit-test path) so the limitation factors are always the same container type and shape as the inputs. |
| A5 | MINOR | `floating_algae.py:902-907`, `benthic_algae.py:565-570` | xarray | `_ap_uptake_fr_nh4` / `_ab_uptake_fr_nh4` branch on `isinstance(ammonium, xr.DataArray)` to decide between `xr.ones_like(...)`/`1.0`. The FloatingAlgae path handles the array case for the NH4-only branch but `_ab_uptake_fr_nh4` returns bare `1.0`/`0.0`/`0.5` for all special cases regardless of input type. Inconsistent container handling between the two sibling methods; the benthic path is the weaker one. | Make `_ab_uptake_fr_nh4` mirror `_ap_uptake_fr_nh4`'s container-aware returns (or, preferably, both return `xr.where`-shaped arrays unconditionally). |
| A6 | MINOR | `floating_algae.py:62-85`, `723` | numerical-robustness | `_sanitize_cache` replaces only NaN (`isnull` / `np.isnan` / `x != x`), not `inf`. `rate_settling` computes `algae / depth * self.settling_velocity`; at `depth == 0` this yields `+inf`, not NaN, which `_sanitize_cache` passes through. The project's own `sanitize_rate` (utils/numerics.py:210) deliberately catches both NaN and inf for exactly this `x / depth` case. The cached `algal_settling_rate` could therefore deliver `inf` to a consumer. The `Model` wet-mask is described as the primary dry-cell defense, so this is defense-in-depth, not a live failure. | Extend `_sanitize_cache` to also map `inf`/`-inf` to 0 (`xr.where(value.isnull() | np.isinf(value), 0, value)`), matching `sanitize_rate`'s contract, or call `sanitize_rate` instead of the local helper. |
| A7 | MINOR | `floating_algae.py:674-685` | algorithm-parity | Harmonic-mean (`growth_rate_option == 3`) computes `rate_raw = growth * FL * 2 / (1/FN + 1/FP)` *before* the `xr.where((FN==0)|(FP==0), 0, rate_raw)` guard. The division `1.0 / limit_nitrogen` is evaluated at every cell including `FN==0`/`FP==0` cells, producing a transient `inf`/`NaN` inside `rate_raw` that `xr.where` then discards. Result is numerically correct (v1 `np.select` has the same eager-evaluation property), but emits divide-by-zero RuntimeWarnings that v1 explicitly suppresses with `warnings.filterwarnings`. Parity in value, not in warning hygiene. | Optional: compute the reciprocal on a safe denominator (`xr.where(FN==0, 1, FN)`) before the guard, or wrap in `np.errstate(divide="ignore", invalid="ignore")`, to match v1's warning suppression. Not a correctness defect. |
| A8 | MINOR | `floating_algae.py:3-46`, `benthic_algae.py:1-20` | documentation | Header docstrings open with `"""v2 NSM1 FloatingAlgae Process.` and `"""v2 NSM1 BenthicAlgae Process.` though these are the v3 modules under `clearwater_modules_v3`. The 2026-05-05 audit doc also refers to these files as `src/clearwater_modules_v2/processes/...`. The v2/v3 labeling is internally inconsistent and will confuse provenance during future audits. | Standardize the module-title line to "v3 NSM1 FloatingAlgae/BenthicAlgae Process" (or document the v2/v3 naming convention once, authoritatively, and apply it consistently). |
| A9 | MINOR | `benthic_algae.py:503-556` | maintainability | `_cache_benthic_mortality_rates` is dead code: `_change_with_components` (line 444) calls `_compute_balgae_mortality_components_from_death`, and the Phase 5 dedup docstring (lines 382-389) states `_cache_benthic_mortality_rates` was the pre-Phase-5 path. The two bodies are now duplicated math (lines 489-502 vs 543-556) with no caller for the latter. Drift risk: a future formula fix applied to one and not the other. | Delete `_cache_benthic_mortality_rates`, or if retained for a test, add a comment naming the test and asserting equivalence with `_compute_balgae_mortality_components_from_death`. |
| A10 | OBSERVATION | `floating_algae.py:328-334` | algorithm-parity (needs verification) | `light_attenuation_coefficient` defaults to the scalar `1.0` (kwarg, audit O4). v1/Fortran compute lambda from the POM/Chla sum in `modGlobalParam`. The header comment correctly flags this as a wiring-only fallback per audit O4. Not an algae-module defect; cross-module (global-vars) scope. | No action in the algae module. Confirm the global-vars review tracks the lambda computation. |
| A11 | OBSERVATION | `floating_algae.py:414`, `benthic_algae.py:314` | algorithm-parity (needs verification) | `solar_radiation` is read from the registry and passed directly as `surface_light_intensity` into `limit_light`/`limit_light` with no `Fr_PAR=0.47` shortwave-to-PAR scaling. Audit F9 flagged this as needing verification against the registry-side `solar_radiation` variable contract: if the registry value is already PAR, no fix is needed; if it is total shortwave, effective irradiance is overstated by ~2x. | Confirm the registry contract for `solar_radiation` (PAR vs total shortwave). If shortwave, apply `Fr_PAR` before `limit_light`. This is the same item as audit F9 and should be resolved in the registry/global-vars review, not silently in the algae module. |
| A12 | OBSERVATION | `algae.py:56-58` | algorithm-parity (intentional, documented) | v3 `mu_max_20=2.0`, `kdp_20=0.05`, `krp_20=0.10` deliberately differ from v1 (`1.0`, `0.15`, `0.2`). The `algae.py` docstring (lines 31-58) documents this with three independent literature citations (Bowie et al. 1985; Chapra/Pelletier/Tao 2008 QUAL2K; Cole & Wells CE-QUAL-W2). Correctly-deferred / intentional-improvement; not a defect. The wiring test anchors on the v3 DEFAULTS values, not v1 literals, so it verifies wiring not v1 numeric parity (correct test design). | None. Listed for completeness so a future reviewer does not re-flag the v1 disagreement as a regression. |
| A13 | OBSERVATION | `balgae.py:36`, corrections doc 1.13 | algorithm-parity (intentional, documented) | v3 `BWa=1000` differs from v1 `BWa=3500` and Fortran `5000`. Documented in `balgae.py:9-22` and `parameter_defaults_corrections.md` Section 1.13 with the WASP7 canonical Chla:DW = 10 mg-Chla/g-DW derivation. Intentional correction; not a defect. Affects the registry-level `Chlb` derived variable (out of algae-Process scope). | None. Listed for completeness. |
| A14 | OBSERVATION | `floating_algae.py:735`, `partitioning.py:24-50` | maintainability | Two flag names gate the same physical decision: `FloatingAlgae.limit_phosphorus` gates on `self.use_phosphate`; the `fdp` partitioning utility gates on `self.use_TIP`. Both default `True`, so masked today (audit F11). v1's `FP` gates on `use_TIP` only. Two switch names for one decision is a maintainability hazard if a future config sets one but not the other. | Document the relationship, or collapse to a single flag. Not a correctness defect under current defaults. |

## 3. Bug #4 / #13 / #14 / #15 / #16 resolution-status table

| Bug | Claim | Actually fixed? | Evidence (current v3 code) |
|---|---|---|---|
| #4 broken multiplicative integrator `algae * rate * dt * 86400` | header says "replace ... with additive Forward Euler `algae + rate * dt_days`" | YES | `floating_algae.py:438-440`: `dt_days = self.time_step.total_seconds() / 86400.0`; `algae_new = algae + rate * dt_days`. No `* 86400` multiplicative term anywhere. Matches v1 `Ap = Ap + dApdt * dt` (`processes.py:685-697`) and `dApdt = ApGrowth - ApRespiration - ApDeath - ApSettling` (`processes.py:666-681`); v3 `rate()` at `floating_algae.py:628-639` is `growth - death - respiration - settling`. Benthic parallel at `benthic_algae.py:342-344` (no settling term, matching v1 `dAbdt`, `processes.py:1069-1082`). |
| #13 `ammonium_respiration()` returns 0 | header says "implement `ammonium_respiration` (was returning 0)" | YES | `floating_algae.py:917-925`: `rna = self.AWn / self.AWa; return rna * self.algal_respiration_rate`. Matches v1 `NH4_ApRespiration = rna * ApRespiration` (`processes.py:1472-1486`). Benthic at `benthic_algae.py:769-782`: `rnb * balgae_respiration_rate * Fb / depth`, matching v1 `NH4_AbRespiration = rnb * AbRespiration * Fb / depth` (`processes.py:1506-1525`). Not a 0-returning stub. |
| #14 `ammonium_growth()` returns 0 | header says "implement `ammonium_growth` (was returning 0)" | YES | `floating_algae.py:927-935`: `return self.algal_nh4_uptake_fraction * rna * self.algal_growth_rate`. Matches v1 `NH4_ApGrowth = ApUptakeFr_NH4 * rna * ApGrowth` (`processes.py:1488-1504`). Benthic at `benthic_algae.py:784-796`: `balgae_nh4_uptake_fraction * rnb * Fb * balgae_growth_rate / depth`, matching v1 `NH4_AbGrowth = AbUptakeFr_NH4 * rnb * Fb * AbGrowth / depth` (`processes.py:1527-1547`). Not a 0-returning stub. |
| #15 hard-coded `phosphate_fraction_dissolved=0.5` | header says "replace ... with the v3 `fdp` partitioning utility" | YES | `floating_algae.py:489-496`: imports `fdp as fdp_partition` and calls `fdp_partition(use_TIP=self.use_TIP, Solid=self.Solid, kdpo4=self.kdpo4)`. No literal `0.5` for the dissolved fraction. Benthic parallel at `benthic_algae.py:397-403`. The v3 `fdp` (`utils/partitioning.py:50`) also corrects the v1 unit-factor inversion (documented Phase 9.B). |
| #16 missing `set_at_time` persistence | header says "re-add `set_at_time` persistence after the integrator step" | YES | `floating_algae.py:450-451`: `registry.set_at_time("algae_floating", time, algae_new)` after the Forward Euler update and `clip_negative_state`. Benthic at `benthic_algae.py:353-354`: `registry.set_at_time("benthic_algae", time, algae_new)`. |
| NaN-guard `rate == np.nan` (always False) | header says "replace ... with `rate.isnull()` / `np.isnan`" | YES | `floating_algae.py:748-751` (`limit_phosphorus`), `774-777` (`limit_nitrogen`): `if isinstance(rate_raw, xr.DataArray): xr.where(rate_raw.isnull(), 0.0, rate_raw) else: xr.where(np.isnan(rate_raw), 0.0, rate_raw)`. `benthic_algae.py:757-761` (`limit_density`) same pattern. No `== np.nan` remains in scope. |

Supporting (audit-flagged formula bugs, independent of the Bug #N list):

| Audit ID | Defect described in 2026-05-05 audit | Actually fixed? | Evidence |
|---|---|---|---|
| F5 | `limit_light` option 1 misplaced parenthesis (np.log argument was only `(KL+PAR)`) | YES | `floating_algae.py:797-808`: `np.log((KL + PAR) / (KL + PAR * np.exp(-(L*depth))))` -- the ratio is now fully inside `np.log`. Matches v1 `processes.py:447`. Inline comment at lines 792-796 documents the fix. |
| F14 | harmonic-mean zero-guard fired on `FP == 1` instead of `FP == 0` | YES | `floating_algae.py:681-685`: `xr.where((limit_nitrogen == 0.0) | (limit_phosphorus == 0.0), 0, rate_raw)`. Matches v1 `mu` `np.select` condition `(FN == 0.0) | (FP == 0.0)` (`processes.py:587`). Inline comment lines 668-673 documents the fix. |
| B6 | benthic Steele `limit_light` option 3 used `x / exp(1-x)` (= `x*exp(x-1)`) instead of `x*exp(1-x)` | YES | `benthic_algae.py:730-735`: `surface * coef / KLb * np.exp(1.0 - surface*coef/KLb)` -- this is `x * exp(1 - x)`. Matches v1 `processes.py:841` `PAR*KEXT/KLb * np.exp(1.0 - PAR*KEXT/KLb)`. Inline comment lines 718-721 documents the fix. |
| F1-F4, B1-B3, B7-B9 | DEFAULTS keys silently shadowed by legacy v2 kwargs (default-instantiated rates = 0, theta = 1.0, wrong KsNb/KsPb/Ksb) | YES | `floating_algae.py:271-326` Phase 9.A.1 bridge: legacy kwargs default to `None`; when `None`, `setattr(self, legacy_name, getattr(self, defaults_name))` copies the DEFAULTS-merged value onto the legacy attribute the rate methods read. `benthic_algae.py:131-142` overrides `_LEGACY_TO_DEFAULTS` to `mub_max_20`/`kdb_20`/`krb_20`/`KsNb`/`KsPb`/`Ksb`, and `__init__` runs the same bridge via `FloatingAlgae.__init__`. Verified by `tests/v3/nsm1/test_phase9a1_algae_wiring.py` (e.g., `test_default_respiration_matches_v1`, `test_default_KsNb_KsPb_Ksb_match_v3_defaults` asserting `ba.KsNb == 0.25`, `ba.Ksb == 10.0`). |

Conclusion: all five named bugs and all four audit-flagged formula/wiring defects are resolved in current `streaming` code. The audit doc (`clearwater_modules_v3_nsm1_audit_algae.md`) is a snapshot of the pre-fix state and should be read as historical.

## 4. Algorithm parity matrix (v3 vs v1)

| Kinetic term | v1 reference | v3 implementation | Verdict |
|---|---|---|---|
| Temperature correction | `arrhenius_correction(T, k20, theta) = k20 * theta**(T-20)` (`processes.py:365-412`) | `utils/conversions.py:18-42`, identical form | match |
| `rna`/`rpa`/`rca`/`rda` | `AWn/AWa`, `AWp/AWa`, `AWc/AWa`, `AWd/AWa` (`processes.py:308-361`) | `floating_algae.py:588-590`, `601` (`self.AWn/self.AWa` etc.) | match |
| Floating growth rate (Arrhenius) | `mu_max_tc = arrhenius(T, mu_max_20, mu_max_theta)` (`processes.py:365-378`) | `rate_growth` `arrhenius_correction(T, self.growth_rate_max, self.growth_rate_correction)` with growth_rate_max bridged to `mu_max_20` (`floating_algae.py:651-655`, `291-302`) | match (intentional default-value change per A12) |
| Light limitation FL opt 1 (half-sat) | `(1/(L*d)) * log((KL+PAR)/(KL+PAR*exp(-Ld)))` (`processes.py:447`) | `floating_algae.py:797-808` | match (F5 fixed) |
| Light limitation FL opt 2 (Smith) | `processes.py:448-449`; abs(KL)<1e-10 -> 1 | `floating_algae.py:810-852`; abs(KL)<1e-10 -> 1 | match |
| Light limitation FL opt 3 (Steele) | `(2.718/(L*d)) * (exp(-PAR/KL*exp(-Ld)) - exp(-PAR/KL))` (`processes.py:451`); abs(KL)<1e-10 -> 0 | `floating_algae.py:854-871` | match |
| FL guards (Ap<=0, L*d<=0) | unified `np.select` guard incl. `PAR<=0` (`processes.py:436-437`) | `floating_algae.py:877-880`; no explicit `PAR<=0` guard | match in value (with F5 fixed, log(KL/KL)=0 at PAR=0); audit F8 minor, not re-raised because the missing guard is now benign |
| Nitrogen limitation FN | four-branch `(N_active)/(KsN+N_active)`; NaN->0; clamp>1->1 (`processes.py:474-527`) | `floating_algae.py:757-781`: `n/(KsN+n)` with use-flag gating; isnull->0; clamp>1->1 | match (algebra equivalent in all four sub-cases) |
| Phosphorus limitation FP | `fdp*TIP/(KsP+fdp*TIP)`; NaN->0; clamp>1->1 (`processes.py:530-561`) | `floating_algae.py:727-755` | match (gating flag name differs, see A14) |
| Growth combination opt 1/2/3 | `np.select` over multiplicative / min / harmonic (`processes.py:582-602`) | `floating_algae.py:657-687` | match (F14 harmonic guard fixed) |
| `ApGrowth`/`ApRespiration`/`ApDeath`/`ApSettling` | `mu*Ap`, `krp_tc*Ap`, `kdp_tc*Ap`, `vsap/depth*Ap` (`processes.py:606-662`) | `rate_growth` returns `rate*algae`; `rate_respiration` `algae*krp`; `rate_death` `algae*kdp`; `rate_settling` `algae/depth*vsap` (`floating_algae.py:689,704,715,723`) | match |
| `dApdt` and Ap update | `growth - resp - death - settling`; `Ap + dApdt*dt` (`processes.py:666-697`) | `floating_algae.py:628-639`, `438-440` | match (Bug #4/#16) |
| `ApUptakeFr_NH4` | NH4-only->1, NO3-only->0, neither->0.5, both->`PN*NH4/(PN*NH4+(1-PN)*NO3)`, NaN->PN (`processes.py:1206-1247`) | `floating_algae.py:887-915` | match |
| `NH4_ApRespiration` | `rna*ApRespiration` (`processes.py:1472-1486`) | `floating_algae.py:917-925` | match (Bug #13) |
| `NH4_ApGrowth` | `ApUptakeFr_NH4*rna*ApGrowth` (`processes.py:1488-1504`) | `floating_algae.py:927-935` | match (Bug #14) |
| Mortality routing OrgN/OrgP | `rna*ApDeath`, `rpa*ApDeath` (`processes.py:1347-1360`, `1897-1912`) | `floating_algae.py:593-594` | match |
| Mortality routing POC/DOC | `f_pocp*rca*ApDeath`, `(1-f_pocp)*rca*ApDeath` (`processes.py:2484-2502`, `2565-2583`) | `floating_algae.py:595-596` | match in form; `f_pocp` default differs (A2) |
| POM from algal settling | `vsap*Ap*rda/h2` (`processes.py:2200-2218`) | `floating_algae.py:600-602` `self.vsap*algae*(self.AWd/self.AWa)/self.h2` | match |
| Benthic growth Arrhenius | `mub_max_tc = arrhenius(T, mub_max_20, mub_max_theta)` (`processes.py:701-713`) | `benthic_algae.py:627-631` via bridged `mub_max_20`/`mub_max_theta` (`_LEGACY_TO_DEFAULTS` override) | match (wiring fixed) |
| FLb opt 1/2/3 | `processes.py:825-861` | `benthic_algae.py:684-748` | match (B6 Steele sign fixed) |
| FNb / FPb | inherited form with KsNb/KsPb | `benthic_algae.py` inherits `limit_nitrogen`/`limit_phosphorus`; `nitrogen_michaelis_menton_constant` bridged to `KsNb`, `phosphorus_michaelis_menton_constant` to `KsPb` (`benthic_algae.py:131-142`) | match (B7/B8 wiring fixed; verified by wiring test asserting `ba.KsNb==0.25`, `ba.KsPb==0.125`) |
| FSb density limitation | `1 - Ab/(Ab+Ksb)`; NaN->0; >1->1 (`processes.py:953-982`) | `benthic_algae.py:750-763` with `density_michaelis_menton_constant` bridged to `Ksb` | match (B9 wiring fixed; verified test asserts `ba.Ksb==10.0`) |
| mub combination opt 1/2 | `np.select` multiplicative / min, default 0 (`processes.py:1006-1022`) | `benthic_algae.py:633-650`; opt 3 raises ValueError | match (v1 also lacks opt 3) |
| `dAbdt` and Ab update | `AbGrowth - AbRespiration - AbDeath` (no settling); `Ab + dAbdt*dt` (`processes.py:1069-1101`) | `benthic_algae.py:603-614`, `342-344` | match |
| `AbUptakeFr_NH4` | same shape as `ApUptakeFr_NH4` with PNb (`processes.py:1263-1302`) | `benthic_algae.py:558-576` | match |
| `NH4_AbRespiration`/`NH4_AbGrowth` | `rnb*AbRespiration*Fb/depth`, `AbUptakeFr_NH4*rnb*Fb*AbGrowth/depth` (`processes.py:1506-1547`) | `benthic_algae.py:769-796` | match |
| Benthic mortality routing OrgN/OrgP/POC/DOC | `rnb*Fw*Fb*AbDeath/depth`, `rpb*Fw*Fb*AbDeath/depth`, `(1/depth)*f_pocb*Fb*Fw*rcb*AbDeath`, `(1/depth)*(1-f_pocb)*Fb*Fw*rcb*AbDeath` (`processes.py:1362-1381`, `1914-1935`, `2505-2529`, `2586-2610`) | `benthic_algae.py:490-497` | match in form; `f_pocb` default differs (A2) |
| Benthic mortality routing POM | `Ab*kdb_tc*Fb*(1-Fw)/h2 = AbDeath*Fb*(1-Fw)/h2` (`processes.py:2257-2277`) | `benthic_algae.py:500-502` | match |

No genuine algorithm-parity discrepancies were found. The only v1/v3 numerical differences are the intentional, documented default-value changes (`mu_max_20`/`kdp_20`/`krp_20`, `BWa`) and the undocumented `f_pocp`/`f_pocb` default (A2).

## 5. Stale-comment list

1. `floating_algae.py:1-46` (module header) -- MAJOR (A1). Entire docstring is a forward-looking plan ("Phase 2.A ... apply the FloatingAlgae bug fixes", "Bug #4: replace the broken multiplicative integrator", "Bug #13: implement `ammonium_respiration` (was returning 0)", "Bug #16: re-add `set_at_time` persistence"). Every listed item is implemented in the same file. The header reads as a TODO/plan; it must read as an "as-implemented / changelog" record.

2. `benthic_algae.py:1-20` (module header) -- MAJOR (same class as A1). "apply the parallel set of bug fixes", "Integrator fix: replace the broken multiplicative integrator", "Implement `ammonium_respiration` and `ammonium_growth`", "Persistence: persist updated state via `registry.set_at_time`". All implemented; header is stale.

3. `floating_algae.py:376` -- MINOR (A3). `# TODO: implement` directly above `self.use_phosphate = True`, which is implemented and consumed. The benthic sibling implements the same concept with `model.has_process` (`benthic_algae.py:253-276`); the floating header even contradicts itself by describing the diagnostics surface as complete.

4. `floating_algae.py:3` and `benthic_algae.py:1` -- MINOR (A8). Title line says "v2 NSM1 ... Process" for files in `clearwater_modules_v3`.

5. `benthic_algae.py:503-556` -- MINOR (A9). `_cache_benthic_mortality_rates` docstring describes itself as the active mortality-routing path ("Compute and cache benthic-algae mortality routing rates"), but it is dead code superseded by `_compute_balgae_mortality_components_from_death`; the only honest signal that it is dead is the Phase 5 dedup note in a *different* method's docstring.

Note on inline fix comments: the inline comments at `floating_algae.py:668-673` (F14), `floating_algae.py:792-796` (F5), and `benthic_algae.py:718-721` (B6) are correctly written in past/explanatory tense ("the previous form ... ", "v1/Fortran fire on ...") and accurately describe a completed fix. These are good and should be preserved (see Section 7).

## 6. Correctly-deferred list

The following are explicitly deferred or are intentional, documented deviations, and are NOT findings against the algae code:

1. Full Registry-side rate-variable plumbing (`Registry.set_rate_variable` / `clear_rate_variables`). The step-scoped `self.<name>` cache pattern is an explicit Phase 2.A workaround pending Phase 2.A.1 (`floating_algae.py:41-45`, `336-339`). Downstream Nitrogen.run reading `floating_algae_process.<rate_name>` directly is the documented interim contract.

2. `mu_max_20=2.0`, `kdp_20=0.05`, `krp_20=0.10` (A12). Intentional move toward the published mesotrophic-river consensus, documented in `algae.py:31-58` with three independent literature citations and an explicit statement that site-specific values belong in YAML overrides.

3. `BWa=1000` vs v1 `3500`/Fortran `5000` (A13). Intentional WASP7-canonical harmonization, documented in `balgae.py:7-22` and `parameter_defaults_corrections.md` Section 1.13.

4. `kdpo4=0.0` (TIP partitioning effectively disabled). `parameter_defaults_corrections.md` Section 2.2 explicitly defers proper TIP/DIP partitioning to the NSM2 path. The algae module correctly consumes the v3 `fdp` utility, which returns `fdp = 1/(1 + kdpo4*Solid*1e-6)` -> 1.0 at `kdpo4=0.0`, matching the intent.

5. `light_attenuation_coefficient` scalar fallback (A10) and `solar_radiation` / PAR coupling (A11). Both are cross-module (global-vars / registry-contract) items per audit O4 and F9; the algae module's header correctly labels the lambda scalar as a "wiring-only fallback per audit O4". These are out of algae-Process scope and should be tracked in the global-vars review, not patched in the algae module.

6. Benthic growth-rate option 3 (harmonic mean) raising `ValueError`. v1 and Fortran also implement only options 1 and 2 for benthic algae; the explicit reject is correct parity, not a missing feature.

## 7. Positive notes (preserve during future refactors)

1. The Phase 9.A.1 legacy-kwarg-to-DEFAULTS bridge (`floating_algae.py:271-326`, `benthic_algae.py:131-142`) cleanly resolved the entire audit "wiring" defect class without breaking v2 back-compat, and is covered by a purpose-built default-instantiation regression test. Do not remove the bridge or the test when v2 kwargs are eventually retired.

2. The inline fix-provenance comments at `floating_algae.py:668-673`, `792-796` and `benthic_algae.py:718-721` are exemplary: they state the prior wrong form, the correct form, and the v1/Fortran anchor. Preserve this convention.

3. The `_sanitize_cache` rationale docstring (`floating_algae.py:62-85`) correctly explains the mass-conservation argument for NaN -> 0 at the cache source (rate * 0 algae = 0). The reasoning is sound; only the inf gap (A6) needs closing.

4. Algorithm parity with v1 is faithful across all kinetic terms (Section 4). The Forward Euler integrator in days, the operator-split additive rate composition, and the clip-with-log contract are correctly and consistently implemented in both Floating and Benthic.

## 8. Recommended follow-up

1. Resolve A1/A8: rewrite both module headers as as-implemented records before any external (LimnoTech) review packet, since a reviewer reading the current headers would conclude the module is still broken.
2. Resolve A2: either set `f_pocp`/`f_pocb` inline fallback to 0.9 or add a corrections-doc entry; this is the only undocumented numeric deviation from v1 in the algae path.
3. Close A6 by routing the settling cache through `sanitize_rate` (inf-safe) rather than the NaN-only `_sanitize_cache`.
4. Add a multi-cell (>=2 cell, coord-bearing DataArray) regression that exercises `limit_phosphorus`/`limit_nitrogen` with a use-flag off and a sibling Process consuming the result, to convert A4/A5 from latent to verified-safe or to surface a real broadcasting defect.
5. Delete or test-anchor the dead `_cache_benthic_mortality_rates` (A9).

## 9. Open questions

1. Is registry `solar_radiation` PAR or total shortwave (A11 / audit F9)? This is the one item that could be a real ~2x irradiance error; it cannot be resolved from the algae source alone.
2. Is the `f_pocp`/`f_pocb` = 0.5 inline default an intentional v3 choice or an oversight (A2)? Needs author input; if intentional it needs a corrections-doc entry, if not it needs to be 0.9.
