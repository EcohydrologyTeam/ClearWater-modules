# v3 NSM1 Carbon & Alkalinity — Line-Level Source & Science Review

Reviewer: water-quality model source-code reviewer (Claude)
Date: 2026-05-15
Branch: `streaming`
Repo: `/Users/todd/GitHub/ecohydrology/ClearWater-modules-streaming`

Scope (read line-by-line):

- `src/clearwater_modules_v3/processes/carbon.py` (836 lines)
- `src/clearwater_modules_v3/processes/alkalinity.py` (592 lines)
- `src/clearwater_modules_v3/parameters/carbon.py` (21 lines)
- `src/clearwater_modules_v3/parameters/alkalinity.py` (21 lines)

v1 reference (cross-read): `src/clearwater_modules/nsm1/processes.py`
(`kpoc_tc`, `POC_hydrolysis`, `POC_settling`, `POC_*_mortality`,
`dPOCdt`, `kdoc_tc`, `DOC_DIC_oxidation`, `dDOCdt`, `Henrys_k`,
`Atmospheric_CO2_reaeration`, `DIC_*`, `dDICdt`, `Alk_*`, `dAlkdt`,
`rca`, `rcb`), `dynamic_variables.py`, `static_variables.py`,
`constants.py`. Supporting v3 reads: `processes/floating_algae.py`,
`processes/benthic_algae.py`, `processes/nitrogen.py`,
`utils/numerics.py`, `utils/reaeration.py`, `utils/conversions.py`.
Audit docs cross-referenced: `clearwater_modules_v3_nsm1_audit_c_dox.md`,
`clearwater_modules_v3_nsm1_audit_simple_constituents.md`,
`clearwater_modules_v3_nsm1_audit_summary.md`,
`clearwater_modules_v3_nsm1_design_specification.md`,
`parameter_defaults_corrections.md`.

Method note: code/science/documentation correctness review. Benchmarks
were not re-run (per task instruction). All `file:line` citations were
read directly.

---

## 1. Summary verdict

The v3 Carbon Process is, in its current (post Phase 9.B / Phase 9.E)
state, a faithful and improved port of v1. The Phase 9.B `rca`/`rcb`
stoichiometric fix and the Phase 9.E DIC unit reconciliation are both
correctly applied in `carbon.py`, the missing CBOD→DIC source (audit
finding C3) and the spurious POC-hydrolysis DOX-Monod (C4) are both
resolved, and the carbonate system is correctly and explicitly deferred
to NSM2 as a simple-tracer (`FCO2`, `KH(T)`) placeholder consistent with
v1. Carbon is in good shape.

The v3 Alkalinity Process carries one Critical defect: it uses the raw
stoichiometric weights `self.AWc` (= 40) and `self.BWc` (= 40) as the
algal C:Chla and C:dry-weight ratios in all four algal/benthic-algae
alkalinity coupling terms, where v1 uses `rca = AWc/AWa = 0.04` and
`rcb = BWc/BWd = 0.4`. This is the **same root-cause defect** that the
Phase 9.B audit flagged as Critical finding C1 for Carbon and DOX and
that was subsequently fixed in `carbon.py` and `dox.py` — but the fix
sweep did not reach `alkalinity.py`, and the simple-constituents audit
(`clearwater_modules_v3_nsm1_audit_simple_constituents.md` §21–§23)
erroneously marked the v3 Alkalinity `rca = self.AWc` form as "Match"
with v1, masking the gap. Floating-algae alkalinity terms are 1000x too
large; benthic-algae terms are 100x too large, whenever algal coupling
is active.

Findings by severity: 1 Critical, 0 Major, 4 Minor, 4 Observations.

The reason the existing v1-vs-v3 benchmarks and the Santiam-Salem case
study still pass despite the Critical defect is the same masking pattern
the audit summary documents for Carbon/DOX: parity tests instantiate
processes with explicit kwargs and feed the same wrong stoichiometric
ratio into both sides, and the Santiam-Salem alkalinity signal is
dominated by the nitrification/denitrification flux terms (which are
correct in v3) rather than by the algal coupling terms (where the defect
lives). The defect is latent under low-algal-productivity conditions and
becomes dominant under bloom conditions.

---

## 2. Findings table

| ID | Severity | File:line | Category | Description | Recommended fix |
|----|----------|-----------|----------|-------------|-----------------|
| CA-1 | Critical | `alkalinity.py:362, 386, 411, 441` | Scientific correctness / algorithm parity | v3 Alkalinity uses raw weights `rca = self.AWc` (40) and `rcb = self.BWc` (40) in the four algal/benthic-algae alkalinity coupling terms. v1 binds the resolved dynamic variable `rca = AWc/AWa = 0.04` (`processes.py:337-347`, `dynamic_variables.py:121-128`) and `rcb = BWc/BWd = 0.4` (`processes.py:776-786`). Floating-algae alkalinity flux is 1000x too large; benthic-algae 100x too large, whenever `use_Algae`/`use_Balgae` and the algae processes are wired. Identical root cause to audit C1 (Carbon/DOX), which was fixed in `carbon.py:495-496` and `dox.py` but not here. | Derive `rca = self.AWc / self.AWa` and `rcb = self.BWc / self.BWd` once (mirroring `carbon.py:495-496`); compose `AWa` and `BWd` into `Alkalinity.DEFAULTS`. Replace `self.AWc`→`rca` at lines 362, 386 and `self.BWc`→`rcb` at lines 411, 441. Add a Fortran-anchored numerical regression test (not same-error parity). |
| CA-2 | Minor | `carbon.py:379-381, 466-473`; `alkalinity.py:540-542` | Documentation-to-code fidelity / stale comment | Docstrings describe a "companion shadow `_change_legacy_inline`" method and `tests/v3/nsm1/test_carbon_helper_vs_inline.py` / `test_alkalinity_helper_vs_inline.py` parity tests "used through Phase 10". Neither the `_change_legacy_inline` methods nor the test files exist anywhere in the tree. The documented helper-vs-inline parity guarantee has no enforcing code. | Remove the stale docstring paragraphs, or implement the shadow method + test if the parity invariant is still wanted. At minimum, delete the references to non-existent infrastructure. |
| CA-3 | Minor | `carbon.py:659-670` (`_ka_tc` docstring) | Documentation-to-code fidelity / stale comment | The `_ka_tc` docstring states it "re-wrap[s] the combined `ka_tc` result with `depth`'s dims if the value count matches; otherwise we fall through to scalar broadcasting." The method body (lines 671-693) does no such re-wrapping — it returns `ka_tc(...)` directly. The dim-stripping issue is now handled inside `utils/reaeration.py:127-128, 213-214` (the utility reattaches `template.coords/dims`). The docstring describes a workaround that was moved upstream. | Update the `_ka_tc` docstring to state that dim preservation is handled in `utils.reaeration`; drop the obsolete "re-wrap … fall through to scalar broadcasting" paragraph. |
| CA-4 | Minor | `carbon.py:33-34, 41, 75-79` (module docstring) | Documentation-to-code fidelity | The module-level dDIC/dt pseudo-equation block (lines 32-42) still shows `KH * pCO2 / 1e6 * 12000` and `+ JDIC / depth` while the prose at lines 75-79 says v1 derives `JDIC` from `SOD_tc / roc` and "a Phase 5.B sediment integration may rewire this." The current code (`carbon.py:594-597`) gates `dic_sed_release` on `use_SedFlux` and yields `0.0` otherwise; the Fortran/v1 non-SedFlux SOD-derived fallback is not implemented. The docstring's `+ JDIC / depth` line does not flag that the term is identically zero under defaults. | Add a one-line note in the docstring dDIC/dt block that `dic_sed_release` is `0.0` unless `use_SedFlux` is True and `JDIC > 0` (the SOD-derived fallback is not ported in 1.0.0; this is audit finding C11, Minor/scope). |
| CA-5 | Minor | `alkalinity.py:46-47, 50-51` (module docstring); `parameters/alkalinity.py:14-15, 18-19` | Documentation-to-code fidelity / units | The Alkalinity module docstring and `parameters/alkalinity.py` comment `r_alkaa`/`r_alkan`/`r_alkba`/`r_alkbn` as `eq/mg-C`. With CA-1 unfixed, the code multiplies these by `self.AWc` (mg-C per **stoichiometric unit**, not per ug-Chla), so the dimensional chain stated in the docstrings (`ug-Chla/L/d * mg-C/ug-Chla * eq/mg-C * mg-CaCO3/eq`) does not actually hold in the code. The unit comments are correct for `r_alk*` themselves; the mis-citation is the implicit claim (e.g., `alkalinity.py:359-361`, `_floating_algae_growth_alk_flux` comment) that `rca = self.AWc` makes the product `mg-C/L/d`. v1's `Alk_algal_growth` docstring (`processes.py:3334-3335`) compounds the confusion by labeling `r_alkaa` as `eq/ug-Chla`; v3's `eq/mg-C` label is the correct one. | Fix CA-1; then the docstring dimensional chains become true. Independently, correct the inline comment at `alkalinity.py:359-361` once `rca` is derived. Note the v1 docstring unit error (`eq/ug-Chla`) as a known upstream inaccuracy that v3's `eq/mg-C` correctly supersedes. |

---

## 3. Algorithm parity matrix

Legend: MATCH = v3 math equals v1 math; IMPROVED(doc) = intentional v3
improvement over v1, documented in an audit/corrections doc;
DISCREPANCY = unintended divergence (finding).

### Carbon — POC / DOC

| v3 term (`carbon.py`) | v1 reference (`processes.py`) | Verdict |
|---|---|---|
| `kpoc_tc = arrhenius_correction(T, kpoc_20, kpoc_theta)` (476-478) | `kpoc_tc` (2439-2451) | MATCH |
| `kdoc_tc = arrhenius_correction(T, kdoc_20, kdoc_theta)` (479-481) | `kdoc_tc` (2614-2626) | MATCH |
| `poc_hydrolysis = kpoc_tc * poc` (501) — no DOX-Monod | `POC_hydrolysis = kpoc_tc * POC` (2455-2465) | MATCH (audit C4 fix applied; earlier v3 had a spurious DOX factor, now removed) |
| `poc_settling = vsoc / depth * poc` (502) | `POC_settling = vsoc/depth*POC` (2469-2481) | MATCH |
| `poc_algal_mortality` ← FloatingAlgae cache `algal_poc_from_mortality_rate` (715-729) | `POC_algal_mortality = f_pocp*kdp_tc*rca*Ap` (2484-2502) | MATCH (cache bakes in correct `rca=AWc/AWa`, `f_pocp`; verified `floating_algae.py:590, 595`) |
| `poc_balgae_mortality` ← BenthicAlgae cache (768-785) | `POC_benthic_algae_mortality` (2505-2529) | MATCH (cache bakes in `rcb=BWc/BWd`, `Fb`, `Fw`, `/depth`) |
| `d_poc = mort_a + mort_b - hydrolysis - settling` (506-511) | `dPOCdt` (2532-2546) | MATCH |
| `doc_oxidation = kdoc_tc * doc * dox/(KsOxmc+dox)` (514) | `DOC_DIC_oxidation` (2629-2647) | MATCH (Monod form identical; v3 always applies the DOX factor — v1 only when `use_DOX`, but v3 default `use_DOX=True`) |
| `doc_algal_mortality` / `doc_balgae_mortality` ← caches (731-803) | `DOC_algal_mortality` / `DOC_benthic_algae_mortality` (2565-2610) | MATCH |
| `pom_doc_source` ← POM cache (526-531) | (v1 `dDOCdt` has no POM term) | IMPROVED(doc) — POM→DOC coupling; consumer-ready cache documented in `carbon.py:517-525`; degrades to 0 when POM absent |
| `d_doc = hydrolysis + mort_a + mort_b + pom - oxidation` (533-539) | `dDOCdt` (2651-2667) | MATCH (+ documented POM addition) |

### Carbon — DIC

| v3 term (`carbon.py`) | v1 reference (`processes.py`) | Verdict |
|---|---|---|
| `henrys_k_co2 = 10**(2385.73/Tk + 0.0152642*Tk - 14.0184)` (136-148) | `Henrys_k` (2687-2695) | MATCH (formula identical, Tk = T+273.15) |
| `co2_reaeration = 0.923*ka_tc*(KH*pCO2/1e6*MG_C_PER_MOL_C - FCO2*dic)` (566-569) | `Atmospheric_CO2_reaeration = 0.923*ka_tc*(K_H*pCO2/1e6 - FCO2*DIC)` (2698-2714) | IMPROVED(doc) — Phase 9.E mg-C/L unit reconciliation; v1's mol-C-vs-mg-C inconsistency corrected by `*MG_C_PER_MOL_C=12000` on the Henry term and removal of `/12000` elsewhere. Documented `parameter_defaults_corrections.md` §1.11 |
| `dic_algal_resp = algae_resp * (AWc/AWa)` (495, 579) | `DIC_algal_respiration = ApRespiration*rca/12000`, `rca=AWc/AWa` (2717-2731) | IMPROVED(doc) — `rca` correctly derived (Phase 9.B C1 fix) AND `/12000` removed (Phase 9.E). Net v3 = v1×12000, intentional unit fix |
| `dic_algal_photo = algae_growth * (AWc/AWa)` (495, 580) | `DIC_algal_photosynthesis` (2734-2748) | IMPROVED(doc) — same as above |
| `dic_balgae_resp = balgae_resp*(BWc/BWd)*Fb/depth` (496, 588) | `DIC_benthic_algae_respiration = AbResp*rcb*Fb/depth/12000` (2751-2769) | IMPROVED(doc) — `rcb=BWc/BWd` derived + `/12000` removed |
| `dic_balgae_photo` (496, 589) | `DIC_benthic_algae_photosynthesis` (2772-2790) | IMPROVED(doc) — same |
| `dic_cbod_oxidation = cbod_oxidation_rate / roc` (605-613) | `DIC_CBOD_oxidation = (1/roc)*(DOX/(KsOxbod+DOX))*kbod_tc*CBOD/12000` (2793-2814) | IMPROVED(doc) — audit C3 fix: v3 reads the pre-attenuated CBOD cache (Monod baked upstream) and drops `/12000` (Phase 9.E). v1 `dDICdt` *did* include this term; current v3 restores it |
| `dic_sed_release = JDIC/depth` if `use_SedFlux` else `0.0` (594-597) | `DIC_sed_release = SOD_tc/roc/depth/12000` (unconditional) (2817-2830) | DISCREPANCY (Minor/scope, = audit C11) — v3 omits the v1 SOD-derived non-SedFlux fallback; identically 0 under defaults. Documented as Phase 5.A scope; see CA-4 |
| `doc_oxidation` added to `d_dic` (616) | v1 `dDICdt` does **not** add `DOC_DIC_oxidation` (2834-2854) | IMPROVED(doc) — v3 restores Fortran `modCarbon.f90:268` DOC→DIC coupling that v1 dropped; audit finding C10. mg-C/L/d basis after Phase 9.E |
| `d_dic = doc_ox + co2 + alg_resp - alg_photo + bal_resp - bal_photo + cbod + sed` (615-624) | `dDICdt` (2834-2854) | MATCH in structure/sign; differs by the documented Phase 9.E unit basis and the +DOC-oxidation / +CBOD restorations |
| Forward Euler `state + rate*dt_days`, clip, set (424-437) | `DIC`/`POC`/`DOC` (`= X + dXdt*dt`, `dt` in days, `static_variables.py:402-408`) | MATCH |

### Alkalinity

| v3 term (`alkalinity.py`) | v1 reference (`processes.py`) | Verdict |
|---|---|---|
| `nitr_sink = r_alkn * nitrification_flux_rate * 50000` (277-300) | `Alk_nitrification = r_alkn*(1-exp(-KNR*DOX))*knit_tc*NH4*50000` (3284-3319) | IMPROVED(doc) — v3 reads pre-attenuated `nitrification_flux_rate` from Nitrogen (Monod baked in `ammonium_nitrification`, `nitrogen.py:706-725`); equivalent under matched params; matches Fortran single-source pattern. `parameter_defaults_corrections.md` §3.3 |
| `denit_source = r_alkden * denitrification_flux_rate * 50000` (302-325) | `Alk_denitrification = r_alkden*(1-DOX/(DOX+KsOxdn))*kdnit_tc*NO3*50000` (3246-3281) | IMPROVED(doc) — same architectural deviation; `nitrate_denitrification` (`nitrogen.py:740-767`) bakes in the O2-inhibition factor |
| `(r_alkaa*fNH4 - r_alkan*(1-fNH4)) * ApGrowth * rca * 50000`, `rca = self.AWc` (327-369) | `Alk_algal_growth` with `rca = AWc/AWa` injected (3322-3342; `dynamic_variables.py:121-128, 1319-1326`) | **DISCREPANCY (Critical, CA-1)** — formula/sign correct, but `rca = self.AWc` (40) vs v1 `AWc/AWa` (0.04): 1000x error |
| `ApRespiration * r_alkaa * self.AWc * 50000` (371-386) | `Alk_algal_respiration = ApResp*r_alkaa*50000*rca`, `rca=AWc/AWa` (3345-3361) | **DISCREPANCY (Critical, CA-1)** — 1000x error from `self.AWc` |
| `(1/depth)*(r_alkba*fbNH4 - r_alkbn*(1-fbNH4))*AbGrowth*Fb*rcb*50000`, `rcb=self.BWc` (388-420) | `Alk_benthic_algae_growth`, `rcb=BWc/BWd` injected (3364-3388) | **DISCREPANCY (Critical, CA-1)** — `rcb=self.BWc` (40) vs `BWc/BWd` (0.4): 100x error |
| `(1/depth)*r_alkba*AbResp*self.BWc*Fb*50000` (422-444) | `Alk_benthic_algae_respiration`, `rcb=BWc/BWd` (3391-3410) | **DISCREPANCY (Critical, CA-1)** — 100x error |
| `rate = denit - nitr - alg_growth + alg_resp - bal_growth + bal_resp` (558-565) | `dAlkdt` (3413-3431) | MATCH (sign convention and operand order verbatim; growth-as-sink with internal NH4/NO3 sign flip) |
| Forward Euler `alk + rate*dt_days`, clip, set (497-504) | `Alk = Alk + dAlkdt*dt` (3435-3447) | MATCH |
| Constant `EQ_TO_MG_CACO3 = 50000` (99) | literal `50000` in every `Alk_*` (3274 etc.) | MATCH |

### Parameter defaults

`parameters/carbon.py` vs v1 `DEFAULT_CARBON` (`constants.py:156-167`):
all ten entries (`f_pocp=0.9`, `kdoc_20=0.01`, `kdoc_theta=1.047`,
`f_pocb=0.9`, `kpoc_20=0.005`, `kpoc_theta=1.047`, `KsOxmc=1.0`,
`pCO2=383.0`, `FCO2=0.2`, `roc=32/12`) MATCH exactly.

`parameters/alkalinity.py` vs v1 `DEFAULT_ALKALINITY`
(`constants.py:54-60`): all six ratios (`r_alkaa=14/106/12/1000`,
`r_alkan=18/106/12/1000`, `r_alkn=2/14/1000`, `r_alkden=4/14/1000`,
`r_alkba=14/106/12/1000`, `r_alkbn=18/106/12/1000`) MATCH exactly.

The parameter files are correct. CA-1 is purely a usage defect in
`alkalinity.py` (treating the raw weight `AWc`/`BWc` as the ratio
`rca`/`rcb`), not a defaults defect.

---

## 4. Carbonate-equilibrium correctness note

v3 1.0.0 deliberately implements **no** carbonate equilibrium solver.
There is no K1, K2, KW, ionic-strength, or [H+]/pH computation in either
scoped file, and none is expected. The design specification
(`clearwater_modules_v3_nsm1_design_specification.md:356, 645`),
the Section 14 resolved-Q "Alkalinity simple-tracer" decision, and the
`alkalinity.py` module docstring (lines 7-12) all consistently defer the
full carbonate-pH solver (carbonate speciation, NH3/NH4+ partitioning,
free-CO2 fraction, alpha fractions) to NSM2 in v3 1.1+.

What the carbon code does retain is dimensionally and scientifically
correct for a simple-tracer model:

1. `henrys_k_co2(T)` (`carbon.py:136-148`) reproduces the v1 empirical
   Henry's-law constant `10^(2385.73/Tk + 0.0152642*Tk - 14.0184)`
   exactly (v1 `Henrys_k`, `processes.py:2687-2695`). Units mol/L/atm.
   No temperature-range guard, but this matches v1 (Observation OBS-3).
2. The atmospheric CO2 exchange is `0.923 * ka_tc * ([CO2*]_eq - FCO2*DIC)`
   where `[CO2*]_eq = KH * pCO2 / 1e6` (mol-C/L) is converted to mg-C/L
   via `MG_C_PER_MOL_C = 12000`. This is the Phase 9.E unit
   reconciliation: the v1/Fortran formula mixed a mol-C/L/d rate with a
   mg-C/L state (a 12000x scaling error that froze DIC). v3's correction
   is internally consistent (every dDIC/dt term is mg-C/L/d) and is the
   correct fix; it is documented in `parameter_defaults_corrections.md`
   §1.11 and flagged here as IMPROVED(doc), not a discrepancy.
3. `FCO2` is a constant from `parameters/carbon.py` (default 0.2),
   correctly documented at `carbon.py:63-66` as a simple-tracer
   placeholder for the pH-dependent free-CO2 fraction that NSM2 will
   compute. This matches v1 (v1 also uses a constant `FCO2`, no alpha
   solve). MATCH and correctly-deferred.
4. The CO2 gas-transfer coefficient reuses the O2 reaeration menu
   (`_ka_tc`, `carbon.py:659-693`), per the customary assumption that
   CO2 and O2 share the transfer velocity. This matches v1
   (`Atmospheric_CO2_reaeration` takes `ka_tc`, the same O2 coefficient)
   and is documented at `carbon.py:54-57`. MATCH.

No iterative carbonate-system solver exists, so there is no
vectorization or per-cell-loop hazard in the carbonate path: there is no
carbonate path. The `henrys_k_co2`, `co2_reaeration`, and all DIC terms
are pure elementwise array arithmetic over the cell dimension and
broadcast correctly. When the NSM2 carbonate/pH solver lands (v3 1.1+),
the reviewer's earlier guidance about vectorizing the [H+] iteration
over cells will become applicable; it is not applicable to 1.0.0.

Conclusion: the carbonate-equilibrium handling is correct for the
declared simple-tracer scope, faithfully matches v1's level of fidelity,
and the one substantive deviation (the Henry-term unit conversion) is a
deliberate, documented correction of a v1/Fortran inconsistency rather
than a defect.

---

## 5. xarray-refactoring assessment

The xarray refactor of both Processes is complete and free of the
classic hazards:

- No Python-level per-cell loops in either `_change_with_components`;
  all term arithmetic is elementwise over `ArrayLike` operands.
- No array-truthiness `if` on state/rate arrays. The only `if`
  statements branch on scalar Python coupling flags
  (`self.use_pom`, `self.use_cbod`, `self.use_floating_algae`, etc.)
  and on `is None` sentinels, never on array contents.
- No `== np.nan`. NaN/inf handling uses `sanitize_rate`
  (`utils/numerics.py:210-252`), which uses `.isnull()` /
  `np.isinf` / `np.where` and preserves container type. `clip_negative_state`
  is likewise container-type-aware and vectorized.
- Coupling reads via `getattr(..., default)` degrade to scalar 0/0.0
  when a sibling Process is absent; `_zeros_like` (`carbon.py:829-833`)
  and `_dox_from_registry` (`carbon.py:811-826`) correctly broadcast
  to `xr.full_like` / `xr.zeros_like` when given a DataArray template,
  scalar otherwise. Multi-cell-safe.
- The historical `np.select` `dim_0` anonymous-axis hazard in the
  reaeration menu is resolved inside `utils/reaeration.py` (it reattaches
  `template.coords/dims`, lines 127-128, 213-214); `_ka_tc` no longer
  needs the workaround its docstring still describes (see CA-3).

No multi-cell breakage was found in the carbon/alkalinity arithmetic
itself. The CA-1 defect is a scalar stoichiometric error that
mis-scales every cell equally; it is not an xarray/broadcasting bug.

---

## 6. Stale-comment list

| Location | Comment claims | Reality | Action |
|---|---|---|---|
| `carbon.py:379-381` | "See `_change_legacy_inline` for the pre-refactor inline composition (retained through Phase 10 for the helper-vs-inline parity test under §11.3)." | No `_change_legacy_inline` method exists in `carbon.py`. | CA-2: remove or implement. |
| `carbon.py:466-473` | "The companion shadow `_change_legacy_inline` returns just the deltas and is used by `tests/v3/nsm1/test_carbon_helper_vs_inline.py` to verify this helper produces bit-identical deltas through Phase 10." | Neither method nor test file exists. | CA-2. |
| `alkalinity.py:540-542` | "The companion shadow `_change_legacy_inline` … used by `tests/v3/nsm1/test_alkalinity_helper_vs_inline.py`." | Neither method nor test file exists. | CA-2. |
| `carbon.py:659-670` (`_ka_tc` docstring) | The method "re-wrap[s] the combined `ka_tc` result with `depth`'s dims if the value count matches; otherwise … fall[s] through to scalar broadcasting." | Method body returns `ka_tc(...)` unmodified; dim preservation moved into `utils/reaeration.py`. | CA-3: update docstring. |
| `alkalinity.py:359-361` (inline comment) | "rca: algal C:Chla. ApGrowth (ug-Chla/L/d) * rca (mg-C/ug-Chla) = mg-C/L/d." | Code sets `rca = self.AWc` (mg-C per stoichiometric unit, not per ug-Chla); the stated dimensional identity does not hold until CA-1 is fixed. | Fix as part of CA-1 / CA-5. |
| `carbon.py:75-79` (docstring `JDIC` note) | "v1 derives this from SOD (`SOD_tc / roc`); a Phase 5.B sediment integration may rewire this." | Accurate as scope intent, but the dDIC/dt equation block at lines 41 shows `+ JDIC/depth` without noting it is identically 0 under defaults (v1's SOD fallback not ported). | CA-4: minor doc clarification (not strictly stale; incomplete). |

No `TODO` / `FIXME` / `XXX` / `HACK` / `BUG#N` markers were found in
either scoped process file or either parameter file. The
`#TODO: make sure np.exp will work here...` at v1 `processes.py:2876` is
in the v1 reference, not in scope, and is not carried into v3.

---

## 7. Correctly-deferred list

These are genuinely open relative to v1/Fortran but are explicitly and
correctly deferred; they are **not** findings:

1. Full carbonate-pH solver (K1, K2, KW, ionic strength, [H+]/pH, alpha
   fractions, pH-dependent `FCO2`, NH3/NH4+ partitioning). Deferred to
   NSM2 in v3 1.1+. Documented at
   `clearwater_modules_v3_nsm1_design_specification.md:356, 645`,
   `alkalinity.py:7-12`, `carbon.py:63-66`,
   `parameters/carbon.py:18`. v1 also has no pH solver, so the
   simple-tracer scope is parity-preserving. Correctly deferred.
2. DIC sediment-release SOD-derived non-SedFlux fallback
   (`SOD_tc / roc / depth`). v3 only supports the `use_SedFlux` +
   user-`JDIC` branch and yields 0 otherwise; v1/Fortran have an
   unconditional SOD-derived release. Documented as Phase 5.A scope at
   `carbon.py:75-79` and audit finding C11 (Minor). Acceptable as a
   documented scope deferral, with the doc-clarity caveat in CA-4.
3. Alkalinity nitrification/denitrification DOX-Monod attenuation is
   sourced from the upstream Nitrogen flux cache rather than recomputed
   locally. Numerically equivalent to v1 under matched parameters and
   matches the Fortran single-source-of-truth pattern. Documented
   `parameter_defaults_corrections.md` §3.3 and `alkalinity.py:28-38,
   284-289, 307-311`. Correctly deferred/accepted as an architectural
   improvement.
4. POM→DOC dissolution source in dDOC/dt is a v3 completion (not in v1
   `dDOCdt`); reads a consumer-ready cache and degrades to 0 when POM
   absent. Documented at `carbon.py:28-29, 517-525`. Acceptable v3
   improvement.

---

## 8. Observations (not defects)

- OBS-1 (positive). The Carbon Process correctly carries the Phase 9.B
  C1 fix (`carbon.py:495-496` derive `rca = self.AWc / self.AWa` and
  `rcb = self.BWc / self.BWd`) and composes `AWa`/`BWd` into
  `Carbon.DEFAULTS` (`carbon.py:286, 288`). This is exactly the pattern
  Alkalinity needs for CA-1 and should be copied verbatim. Worth
  preserving through future refactors.
- OBS-2 (positive). `sanitize_rate` is applied both to the integrated
  deltas and to every cached component (`carbon.py:627-651`,
  `alkalinity.py:568`), with a clear rationale comment
  (`carbon.py:637-641`) about NaN propagation through the DOX rate sum.
  This is sound defense-in-depth for the producer→consumer cache
  contract.
- OBS-3 (needs verification — low priority). Neither `henrys_k_co2`
  (`carbon.py:136-148`) nor the algal/benthic coupling applies any
  temperature-range or domain guard. At extreme `t_water_c`, `KH(T)`
  and `arrhenius_correction` can produce very large values; `sanitize_rate`
  only catches NaN/inf, not large-but-finite excursions. This matches
  v1 behavior exactly (v1 `Henrys_k` is unguarded), so it is parity-
  preserving and not a v3 regression. Flagged only as a latent
  numerical-robustness item to consider when the NSM2 carbonate solver
  lands.
- OBS-4 (parity quirk, v3 is more correct). v1 `Alk_nitrification`'s
  second `np.select` branch (`use_NH4` true, `use_DOX` false) is
  `knit_tc * NH4 * 50000` — it drops the `r_alkn` factor (v1
  `processes.py:3312-3313`), an apparent v1 bug. v3 always applies
  `r_alkn` (`alkalinity.py:300`), so v3 is more correct in the
  `use_DOX=False` configuration. Under the default `use_DOX=True` the
  two agree. Not a finding; noted so it is not mistaken for a v3 defect
  during future parity work.

---

## 9. Recommended follow-up tests

1. Fortran-anchored numerical regression for Alkalinity algal coupling
   (closes CA-1 and the same-error masking the audit summary describes).
   Example anchor: `AWc=40`, `AWa=1000`, `ApGrowth=0.5 ug-Chla/L/d`,
   `fNH4=1.0`, `r_alkaa=14/106/12/1000` should yield
   `(r_alkaa)*0.5*(40/1000)*50000 ≈ 1.10e-2 mg-CaCO3/L/d`, not the
   `*40` value `≈ 11.0 mg-CaCO3/L/d` the current code produces. Assert
   the explicit reference value, not v1 called with the same wrong `rca`.
2. Closed-system Tier 1 alkalinity + carbon conservation test with
   nonzero floating- and benthic-algae growth/respiration, asserting
   `model.diagnostics.clip_events == 0` and that algal-coupling
   alkalinity flux magnitude is within the physical band (this would
   have caught CA-1: 1000x inflation forces clip events / mass blowup).
3. Helper-vs-inline parity test (or removal of the docstring claim) to
   resolve CA-2 — either implement `_change_legacy_inline` + the two
   `test_*_helper_vs_inline.py` files, or delete the references.
4. DIC unit-reconciliation regression that asserts the mg-C/L/d basis
   of every dDIC/dt term and the `*MG_C_PER_MOL_C` Henry conversion,
   pinning the Phase 9.E correction against accidental reversion.

---

## 10. Open questions

1. Was Alkalinity intentionally excluded from the Phase 9.B rca/rcb fix
   sweep, or was it an oversight propagated by the simple-constituents
   audit's incorrect "Match" verdict
   (`clearwater_modules_v3_nsm1_audit_simple_constituents.md` §21–§23)?
   The audit summary line items enumerate the fix for Carbon (#19) and
   DOX (#22) but contain no Alkalinity entry. Author input needed to
   confirm CA-1 is an oversight (recommended treatment) rather than a
   deliberate but undocumented convention.
2. Should the simple-constituents audit document be corrected? Its
   §21/§22/§23 "Match" verdicts for v3 Alkalinity `rca = self.AWc`
   are factually wrong (v1 `rca = AWc/AWa`) and actively mask CA-1.
   Recommend annotating it once CA-1 is fixed so the audit trail is
   consistent.
3. Is the DIC SOD-derived non-SedFlux fallback (audit C11) intended for
   a Phase 5.5 follow-up, or is the SedFlux-only behavior the final
   1.0.0 contract? This affects whether CA-4 is a doc-only fix or a
   tracked scope item.
