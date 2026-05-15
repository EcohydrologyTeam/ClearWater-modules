# v3 NSM1 DOX Source-Code and Science-Correctness Review

Review date: 2026-05-15
Reviewer: water-quality-model-source-reviewer agent
Repository: ClearWater-modules-streaming, branch `streaming`
Scope reviewed (read in full):

- `src/clearwater_modules_v3/processes/dox.py` (846 LOC)
- `src/clearwater_modules_v3/parameters/dox.py`
- `src/clearwater_modules_v3/utils/reaeration.py`

Reference material consulted (targeted ranges only):

- v1 `src/clearwater_modules/nsm1/processes.py` DOX block lines 2876-3135 and reaeration/SOD helpers lines 50-236
- v1 `src/clearwater_modules/nsm1/constants.py` (DOX/SOD/reaeration/algal defaults)
- `design/clearwater_modules_v3_nsm1_audit_c_dox.md` (DOX and Carbon audit sections)
- `src/clearwater_modules_v3/parameter_defaults_corrections.md` (items 1.1-1.6)
- `src/clearwater_modules_v3/processes/nitrogen.py` (nitrification flux provenance) and `src/clearwater_modules_v3/utils/sediment.py` (SOD primitive)

Language: Python (NumPy / xarray). No source was modified. No benchmarks were run.

---

## 1. Summary verdict

The v3 DOX Process is a faithful, well-structured v3-native re-implementation of the v1 NSM1 dissolved-oxygen kinetics. The xarray refactor is essentially complete: there are no cell loops, no array-truthiness `if` over DataArrays, no `== np.nan` comparisons, and per-sub-flux NaN sanitization is applied before summation. The Phase 9.B atomic-weight resolution fix did reach dox.py: all four algal/benthic stoichiometric ratios are computed as resolved ratios (`rca = AWc/AWa`, `rcb = BWc/BWd`) rather than raw weights, and `roc` / `ron` are resolved-ratio constants. This is the opposite of the alkalinity.py defect the sibling reviewer found.

The single behaviorally significant issue is the documented-but-physically-incorrect default: with the v3 corrected defaults (`kah_20_user = 0`, `kaw_20_user = 0`) plus `wind_reaeration_option = 1`, the wind branch contributes zero and, although `hydraulic_reaeration_option` now defaults to 5 (Cover 1976), the explicit zero-reaeration short-circuit in `_change_with_components` only fires when *both* options are 1; with the shipped defaults (hydraulic option 5, wind option 1) reaeration is computed normally from hydraulics, so the C5 audit concern is now mitigated by the Phase 9.E option-5 default. The residual concern is narrower than C5 originally stated and is recorded below as a Major finding scoped to the user-override configuration. A documentation drift in the module docstring (it still describes the SOD-Monod factor and reaeration short-circuit in ways that no longer fully match Phase 9.E) and one stale marker in `utils/reaeration.py` are Minor.

Findings by severity: Critical 0, Major 1, Minor 4, Observation 4.

Confidence: high for algorithm parity (every term traced line-for-line against v1), high for the Phase 9.B verdict (all four ratio sites read), medium for the salinity-omission and reaeration-default findings (these are documented design choices whose downstream behavior was not executed here).

---

## 2. Findings table

| ID | Severity | file:line | Category | Description | Fix |
|----|----------|-----------|----------|-------------|-----|
| F1 | MAJOR | `processes/dox.py:754-761` | Default behavior / science | The zero-reaeration short-circuit fires only when *both* `hydraulic_reaeration_option == 1` and `wind_reaeration_option == 1` with both user values 0. With the shipped defaults (hydraulic option 5, wind option 1, both user values 0) this branch is *not* taken, so hydraulic reaeration is computed correctly. However, if a user sets `hydraulic_reaeration_option = 1` (the v1/Fortran convention) without supplying `kah_20_user`, the corrected default `kah_20_user = 0.0` yields identically zero atmospheric reaeration with no warning. This is the residual of audit finding C5. | Emit a one-time `logger.warning` when the short-circuit path is taken, or when `hydraulic_reaeration_option == 1 and kah_20_user == 0.0`, stating that atmospheric reaeration is disabled. Document the interaction in the DOX docstring. |
| F2 | MINOR | `processes/dox.py:36-40, 161-201` | Documentation-to-code fidelity / science | Neither `dox_sat_apha` nor the module docstring notes that the APHA saturation formula omits the salinity correction (`exp(-S*(0.017674 - 10.754/Tk + 2140.7/Tk^2))`) present in Fortran `modDOX.f90:97-99`. This is audit finding C6. Zero impact for fresh water (S = 0 gives factor 1.0) but silently overstates O2sat for brackish coupling. | Add a docstring note in `dox_sat_apha` stating the freshwater assumption and that salinity correction is deliberately omitted (matching v1), referencing the audit C6 deferral. |
| F3 | MINOR | `processes/dox.py:63-69, 607-633` | Documentation-to-code fidelity | The module docstring (lines 63-69) states `utils.sediment.SOD_tc` is "pure Arrhenius (no DOX-Monod attenuation)" and that DOX re-applies the Monod factor. The code does this correctly (`_sod_flux` lines 630-633), but the docstring's Fortran citation `modGlobalParam.f90:254` and v1 `shared.processes.SOD_tc:200` should be reconciled with the actual v1 location read here (`shared` reaeration/SOD helper at `processes.py:216-236`, `SOD_tc` Monod at line 234). The algebra is correct; only the citation precision is at issue. | Verify and correct the line citations in the docstring to the exact v1/Fortran lines, or generalize the reference to avoid stale line numbers. |
| F4 | MINOR | `utils/reaeration.py:168-169` | Stale comment | Comment at lines 168-169 ("See `kah_20` for the same `np.select` dim-stripping fix") references a fix that is implemented and working; the comment phrasing reads as a pending cross-reference rather than a completed rationale. Not a defect, but the only marker in this file per the task brief. | Rephrase to a settled rationale ("The `np.select` result loses xarray dims; reattach `wind_speed` metadata, as in `kah_20`."). No code change. |
| F5 | MINOR | `parameters/dox.py:26-27` | Documentation / parameter provenance | The inline comment for `hydraulic_reaeration_option = 5` labels option 5 as "Cover 1976 / Internal". The corrections doc and `utils/reaeration.py:64` Notes attribute option 5 to "Covar 1976" / "Cover (1976): depth-piecewise blend of options 2-4". The author name is inconsistently spelled "Cover" vs "Covar" across files. | Standardize the citation spelling (the QUAL2K manual uses "Covar 1976") across `parameters/dox.py`, `utils/reaeration.py`, and the corrections doc. |
| F6 | OBSERVATION | `processes/dox.py:672` | Robustness (needs verification) | Zero-NH4 fallback builds `xr.zeros_like(dox)` only when `dox` has `.dims`; otherwise scalar `0.0`. Correct for the nitrification path (which is gated on Nitrogen being wired) but the branch was not exercised by a test here. | Confirm via a Tier 1 unit test that DOX with no `ammonium` in the registry and no Nitrogen process integrates without shape errors. |
| F7 | OBSERVATION | `processes/dox.py:494-542` | Numerical robustness (needs verification) | Benthic algae photosynthesis/respiration divide by `depth` with no thin/dry-cell guard inside `_benthic_algae_growth_flux` / `_benthic_algae_respiration_flux`; `_sod_flux` similarly divides by `depth`. v1 has the same structure (no guard). Per-sub-flux `sanitize_rate` in `_change_with_components` catches the resulting inf/NaN at thin cells, so this is defended downstream, not at the division. | Confirm `sanitize_rate` maps `inf` (from `x/0`) to a finite value (e.g. 0) and that this is the intended thin-cell behavior; document it. This matches v1 and is not a regression. |
| F8 | OBSERVATION | `processes/dox.py:595-605` | Parity note | `_cbod_oxidation_flux` returns `cbod_oxidation_rate` directly (1 mg-O2 per mg-CBOD). v1 `DOX_CBOD_oxidation` (lines 3019-3029) likewise returns `CBOD_oxidation` directly and ignores its `roc` argument. v3 matches v1 exactly. Preserve this; do not "fix" the unused-`roc` asymmetry to match the DIC side. | None. Positive parity note. |
| F9 | OBSERVATION | `processes/dox.py:807-818` | Conservation | Net-rate assembly term order and signs are identical to v1 `dDOXdt` (line 3119): reaeration + ApGrowth - ApResp - Nitrif - DOC - CBOD + AbGrowth - AbResp - SOD. No conservation defect; DOX is a terminal integrator and does not publish a rate cache consumed elsewhere. | None. |

---

## 3. Algorithm parity matrix (v3 vs v1)

| Term | v1 reference | v3 location | Verdict |
|------|--------------|-------------|---------|
| O2 saturation (Benson-Krause 4-coeff log polynomial + pressure/pwv/alpha correction) | `processes.py:2901-2923` (`DOX_sat`), `pwv` 2878-2886, `DOs_atm_alpha` 2890-2898 | `dox_sat_apha` 161-201, `_pwv_atm` 146-148, `_do_atm_alpha` 151-158 | MATCH. Coefficients, `pressure_mb * 0.000986923`, and the full correction ratio are byte-identical. Salinity term omitted in both (see F2; documented audit C6 deferral). |
| Atmospheric reaeration `ka_tc * (O2sat - DOX)` | `processes.py:2927-2939` (`Atm_O2_reaeration`) | `_atm_reaeration_flux` 433-440 | MATCH. |
| Reaeration coefficient `ka_tc = kaw_tc/depth + kah_tc` | `processes.py:202-214` (`ka_tc`) | `utils/reaeration.py:217-246` | MATCH (Arrhenius-corrected kah/kaw, wind velocity divided by depth, hydraulic added). |
| `kah_20` (9 hydraulic options, depth/flow piecewise) | `processes.py:50-106` | `utils/reaeration.py:26-128` | MATCH. All 14 condlist branches and choicelist expressions identical, including option-5 depth piecewise and option-9 Froude shear form. v3 adds dim/coord reattachment after `np.select` (correct fix, not a behavior change). |
| `kaw_20` (13 wind options, Uw10 = wind*(10/2)^0.143) | `processes.py:125-182` | `utils/reaeration.py:131-214` | MATCH. All 17 condlist branches and choicelist coefficients identical. |
| Floating-algae photosynthesis `ApGrowth * rca * roc * (138/106 - 32/106 * ApUptakeFr_NH4)` | `processes.py:2942-2959` (`DOX_ApGrowth`) | `_floating_algae_growth_flux` 442-474 | MATCH. v3 uses resolved `rca = AWc/AWa` (Phase 9.B). v1 form `138/106 - 32*ApUptakeFr_NH4/106` is algebraically identical to v3 `138/106 - 32/106*ApUptakeFr_NH4`. |
| Floating-algae respiration `ApRespiration * rca * roc` | `processes.py:2962-2977` (`DOX_ApRespiration`) | `_floating_algae_respiration_flux` 476-492 | MATCH. Resolved `rca`. |
| Benthic-algae photosynthesis `(138/106 - 32/106*AbUptakeFr_NH4) * roc * rcb * AbGrowth * Fb / depth` | `processes.py:3032-3054` (`DOX_AbGrowth`) | `_benthic_algae_growth_flux` 494-525 | MATCH. Resolved `rcb = BWc/BWd`. Depth-normalized. |
| Benthic-algae respiration `roc * rcb * AbRespiration * Fb / depth` | `processes.py:3057-3078` (`DOX_AbRespiration`) | `_benthic_algae_respiration_flux` 527-542 | MATCH. Resolved `rcb`. Depth-normalized. |
| Nitrification O2 sink `ron * (1 - exp(-KNR*DOX)) * knit_tc * NH4` | `processes.py:2980-2999` (`DOX_Nitrification`) | `_nitrification_flux` 544-575 (delegated to `Nitrogen.nitrification_flux_rate`) | MATCH (architectural improvement). v3 reads Nitrogen's cached `nitrification_flux_rate` and multiplies by `ron`. Nitrogen's `ammonium_nitrification` (`nitrogen.py:724`) computes `ammonium * rate_corrected * nitrification_inhibition(DOX)` with `nitrification_inhibition = 1 - exp(-KNR*DOX)` (`nitrogen.py:885`). Composition reproduces the v1 closed form exactly with single ownership of the inhibition term. |
| DOC oxidation O2 sink `roc * DOC_DIC_oxidation` | `processes.py:3002-3015` (`DOX_DOC_oxidation`) | `_doc_oxidation_flux` 577-593 | MATCH. Gated on `use_carbon and use_DOC`. |
| CBOD oxidation O2 sink `CBOD_oxidation` (1:1) | `processes.py:3019-3029` (`DOX_CBOD_oxidation`) | `_cbod_oxidation_flux` 595-605 | MATCH. v1 also returns `CBOD_oxidation` directly and ignores `roc`. |
| SOD O2 sink `SOD_tc * DOX/(DOX+KsSOD) / depth` | `processes.py:3081-3092` (`DOX_SOD`) + Monod in `SOD_tc` 216-236 | `_sod_flux` 607-633 | MATCH (Phase 1.1 + 9.B re-composition). v3 pulls pure-Arrhenius `SOD_tc` from `utils.sediment` and re-applies the DOX-Monod factor in the consumer, gated on `use_DOX`. Net result identical to v1's `xr.where(use_DOX, SOD_tc*DOX/(DOX+KsSOD), SOD_tc)/depth`. |
| dDOX/dt assembly (term order and signs) | `processes.py:3095-3119` (`dDOXdt`) | `_change_with_components` 807-818 | MATCH. |
| Forward Euler `DOX + dDOXdt*dt` (dt in days) | `processes.py:3123-3135` (`DOX`) | `run` 700-701 + `_change_with_components` 825-826 | MATCH. `dt_days = time_step.total_seconds()/86400.0`; rate is per-day; delta is per-substep. Unit handling correct. |

No DISCREPANCY findings. No undocumented algorithmic deviation from v1.

---

## 4. Reaeration and DO-saturation method-option note

The reaeration menu (`kah_20` options 1-9, `kaw_20` options 1-13) is reproduced exactly from v1, including all piecewise depth/flow/Uw10 thresholds and all empirical coefficients. The v3 port adds a necessary correctness fix: `np.select` strips xarray dim/coord metadata, so `kah_20` and `kaw_20` reattach the template DataArray's dims/coords via `_first_dataarray`. Without this, downstream broadcasting against `oxygen_dissolved` would produce a spurious `cell x dim_0` outer product. This is an xarray-refactor improvement, not a behavior change, and is correctly implemented.

The combined coefficient `ka_tc = kaw_tc/depth + kah_tc` (units: m/d / m + 1/d = 1/d) is dimensionally consistent. `kaw_20` returns a velocity (m/d) and the depth division to convert to a rate occurs only once, in `ka_tc`, matching v1.

The DO-saturation formula is the Benson-Krause four-coefficient log polynomial with the QUAL2E/APHA pressure-vapor-alpha correction, identical to v1 to the digit. The empirical `DOs_atm_alpha` uses Celsius (the v3 docstring at lines 151-158 correctly notes that the v1 docstring mislabels the input as Kelvin while the formula and call site use Celsius; this is a documentation correction, not a behavior change). The salinity correction present in Fortran `modDOX.f90:97-99` is omitted in both v1 and v3 (F2; documented audit C6 deferral, correct for freshwater NSM1).

Default reaeration configuration: shipped defaults are `hydraulic_reaeration_option = 5` (Phase 9.E, Covar 1976 depth-piecewise blend), `wind_reaeration_option = 1`, `kah_20_user = 0.0`, `kaw_20_user = 0.0`. With these defaults the dox.py short-circuit at lines 754-761 does NOT fire (it requires hydraulic option 1), so atmospheric reaeration IS computed from cell hydraulics. The Phase 9.E option-5 default substantially mitigates the original C5 audit concern. The residual risk (F1) is confined to users who explicitly select `hydraulic_reaeration_option = 1` without supplying `kah_20_user`, which now silently yields zero reaeration.

---

## 5. Phase 9.B-in-dox.py verdict

VERDICT: the Phase 9.B atomic-weight resolution fix DID reach dox.py. Unlike alkalinity.py (which the sibling reviewer found unfixed), dox.py uses resolved ratios at every coupling site:

- `_floating_algae_growth_flux` (line 468): `rca = self.AWc / self.AWa` — resolved C:Chla ratio, not raw `AWc`.
- `_floating_algae_respiration_flux` (line 491): `rca = self.AWc / self.AWa` — resolved.
- `_benthic_algae_growth_flux` (line 517): `rcb = self.BWc / self.BWd` — resolved C:dry-weight ratio, not raw `BWc`.
- `_benthic_algae_respiration_flux` (line 541): `rcb = self.BWc / self.BWd` — resolved.
- Nitrification O2 stoichiometry uses `self.ron` (line 575), composed from `DOX_DEFAULTS` as `2.0*32.0/14.0 = 4.5714...` — a resolved O2:N ratio, not a raw atomic weight.
- Carbon oxidation O2 stoichiometry uses `self.roc` (lines 472, 492, 521, 542, 593), composed from `CARBON_DEFAULTS` as `32.0/12.0` — a resolved O2:C ratio.

The `__init__` composition (lines 336-349) correctly stores the raw weights `AWc=40`, `AWa=1000`, `BWc=40`, `BWd=100` and the helpers derive `rca = 40/1000 = 0.04` mg-C/ug-Chla and `rcb = 40/100 = 0.4` mg-C/mg-D at run time. The inline comments at lines 462-468, 489-491, 511-517, and 539-541 explicitly document the Phase 9.B audit fix and cite the Fortran `modDOX.f90` and v1 lines that derive the same way. This matches the v1 `DOX_ApGrowth`/`DOX_AbGrowth` contracts where the caller passes `rca = AWc/AWa`. No raw-weight-as-ratio defect exists in dox.py.

---

## 6. Parameter-default verification

| Param | v3 value (`parameters/dox.py`) | v1 (`constants.py`) | Corrections doc | Verdict |
|-------|-------------------------------|---------------------|-----------------|---------|
| `ron` | `2.0*32.0/14.0` = 4.5714 | `2.0*32.0/14.0` (line 189) | n/a | MATCH. Resolved O2:N ratio. |
| `KsSOD` | `1.0` | `1` (line 190) | n/a | MATCH. |
| `SOD_20` | `1.0` | `999` (sentinel, line 322) | item 1.3: v3 = 1.0 (Fortran 0.2; defensible midpoint, flagged for LimnoTech) | CORRECTED as documented. v1 sentinel 999 drives DOX negative on any wet timestep; v3 fix is correct and documented. |
| `SOD_theta` | `1.060` | `999` (sentinel, line 323) | item 1.4: v3 = 1.060 (matches Fortran/Chapra 1997) | CORRECTED as documented. v1 sentinel 999 is catastrophic (orders-of-magnitude blow-up per degree). |
| `kaw_20_user` | `0.0` | `999` (sentinel) | item 1.5: v3 = 0.0 (matches Fortran) | CORRECTED as documented. |
| `kah_20_user` | `0.0` | `999` (sentinel) | item 1.6: v3 = 0.0; Fortran 1.0 | CORRECTED as documented. Behavioral interaction with option default is F1 (Major). |
| `kaw_theta` | `1.024` | `1.024` (line 329) | n/a | MATCH. |
| `kah_theta` | `1.024` | `1.024` (line 330) | n/a | MATCH. |
| `hydraulic_reaeration_option` | `5` | `1` | item 1.6 / Phase 9.E: v3 = 5 (Covar 1976 / QUAL2K Internal default) | CORRECTED as documented. This is the key mitigation for audit C5: shipped default no longer routes through the zero-yielding user branch. |
| `wind_reaeration_option` | `1` | `1` | item 1.6: appropriate for stream/river NSM1 (QUAL2K wind option 1 "omitted") | MATCH and documented. |

All ten DOX parameter defaults are either an exact match to v1 or a documented correction with rationale in `parameter_defaults_corrections.md` Section 1. No undocumented default deviation. The single behavioral concern (F1) is the residual of audit C5 narrowed by the Phase 9.E option-5 default.

---

## 7. Stale-comment list

- `utils/reaeration.py:168-169` — comment phrased as a pending cross-reference for a completed `np.select` dim-stripping fix (F4, Minor; the only marker in this file per the task brief). Not a fixed-but-still-labeled-broken case; the fix is present and the comment is rationale, merely awkwardly worded.
- `processes/dox.py:63-69` — module docstring SOD-Monod / Phase 1.1 / Phase 9.B narrative is accurate to the code but carries Fortran/v1 line citations (`modGlobalParam.f90:254`, `shared.processes.SOD_tc:200`) that should be reconciled with the actual v1 helper location verified here (`processes.py:216-236`) (F3, Minor).
- No `TODO`, `FIXME`, `BUG`, or `XXX` markers were found in any of the three reviewed files. The v1 `processes.py:2876` `#TODO: make sure np.exp will work here...` is in v1, not v3, and is not in scope.

No fixed-but-still-labeled-broken comment exists in the reviewed v3 files.

---

## 8. Correctly-deferred list

- Salinity correction on O2sat (Fortran `modDOX.f90:97-99`): omitted in both v1 and v3; correctly deferred for freshwater NSM1. Recorded as audit C6 and as F2 (documentation note only). Not a finding against parity since v3 matches v1.
- SedFlux-coupled sediment diagenesis (NSM2 territory): the DOX SOD term uses the Arrhenius + DOX-Monod water-column formulation, not a coupled sediment-diagenesis flux. Correctly deferred to NSM2 / v3.x per project scope; not a finding.
- DIC-side CBOD source term and non-SedFlux DIC sediment release (audit C3, C11): these are Carbon-process scope, not DOX, and are out of scope for this review. Noted only to confirm the DOX side of CBOD coupling (F8) is correct and complete.
- Phase 5.5 semi-implicit source/sink split: `self.dox_rate` is cached for a future opt-in split (lines 372-375, 697). Correctly deferred infrastructure, not a defect.
