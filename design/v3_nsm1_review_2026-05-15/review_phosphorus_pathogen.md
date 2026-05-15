# v3 NSM1 Phosphorus & Pathogen — Line-Level Source & Science Review

Reviewer: water-quality model source-code reviewer (automated)
Date: 2026-05-15
Branch: `streaming`  Commit: `54f2b1209f07a1eb7ac5d26238b4ca18d5ba7777`

Scope (read line-by-line):

- `src/clearwater_modules_v3/processes/phosphorus.py`
- `src/clearwater_modules_v3/processes/pathogen.py`
- `src/clearwater_modules_v3/parameters/phosphorus.py`
- `src/clearwater_modules_v3/parameters/pathogen.py`
- `src/clearwater_modules_v3/utils/partitioning.py`

v1 cross-read: `src/clearwater_modules/nsm1/processes.py`,
`dynamic_variables.py`, `static_variables.py`, `constants.py`,
`src/clearwater_modules/shared/processes.py`.

Authoritative docs consulted: `design/clearwater_modules_v3_nsm1_audit_n_p.md`,
`design/clearwater_modules_v3_nsm1_audit_simple_constituents.md`,
`design/clearwater_modules_v3_nsm1_design_specification.md` (Sections 6, 7,
11), `src/clearwater_modules_v3/parameter_defaults_corrections.md`
(Sections 1.1, 1.2, 1.15, 1.16, 2.1, 2.2).

---

## 1. Summary verdict

The v3 phosphorus and pathogen kinetics are faithful, science-correct
ports of the v1 NSM1 formulas. Every kinetic term — OrgP→TIP hydrolysis,
OrgP settling, OrgP from algal/benthic mortality, TIP particulate
settling, TIP sediment release, TIP algal uptake/release, the dOrgP/dt
and dTIP/dt budgets, pathogen natural death, light-induced decay, and
settling — matches v1 algebraically (parity matrix in Section 3). Units,
depth normalization, and per-day rate conventions are correct throughout.
Forward Euler with `dt_days = total_seconds / 86400` is consistent with
v1's per-day `dt`.

The xarray refactoring is essentially complete. There are no cell loops,
no `== np.nan` comparisons, no array-truthiness `if` on per-cell data,
and the NaN/inf guard uses `.isnull()` / `xr.where`. Several scalar
fallbacks exist (`0.0` returns from algal-coupling helpers and the
`use_TIP`/`use_OrgP` False branches); these are broadcast-safe because
they enter additive rate expressions, with one narrow exception flagged
as a minor observation (F6).

The one substantive science topic — TIP sorption partitioning (`fdp`) —
needs a precise statement because the source, the v1 reference, and the
audit doc disagree with each other:

1. The **current v3 source** (`utils/partitioning.py:50`) implements the
   dimensionally correct sorption isotherm
   `fdp = 1 / (1 + kdpo4 * Solid * 1e-6)`. This is the *fixed* form.
2. The **v1 NSM1 model as actually run** wires `nsm1/processes.py:290`
   `fdp`, which returns `xr.where(use_TIP, 1, 0)` — a constant 1.0 with
   no isotherm at all. (The isotherm copy in `shared/processes.py:271`
   exists but is *not* wired into v1 NSM1's `dynamic_variables.py`.)
3. The **n_p audit P6** and its summary still quote the *old buggy* v3
   form `/ 0.000001` and label v3 critical-but-gated. That audit text is
   now stale; the source has been corrected.

At the shipped default `kdpo4 = 0`, all three converge to `fdp = 1.0`,
so `(1 - fdp)·TIP = 0`, TIP settling is zero, and v3 reproduces v1 NSM1
and the benchmarks exactly. The divergence is entirely latent and only
manifests if a user sets `kdpo4 > 0`. Findings F1 and F2 capture the
documentation/parity consequences. Section 4 gives the full analysis.

Findings: 0 Critical, 0 Major, 5 Minor, 4 Observation. No defect blocks
v3 1.0.0 at the validated default-parameter regime. The Minor items are
documentation/parity-traceability and code-hygiene issues.

---

## 2. Findings table

| ID | Severity | Location | Category | Description | Recommended fix |
|----|----------|----------|----------|-------------|-----------------|
| F1 | Minor | `design/clearwater_modules_v3_nsm1_audit_n_p.md:20,45-48,401` | Stale audit doc vs source | Audit P6 + top-concerns item 5 + summary line 20 state the v3 `fdp` divides by `0.000001` (factor-1E12 bug) and is "critical, gated." The current source `utils/partitioning.py:50` already uses the corrected `* 1.0e-6` form. The audit is stale relative to the fixed code. | Update P6 / summary to record the fix landed: v3 `fdp` is now dimensionally correct; reclassify from "critical (latent)" to "resolved; v3 deliberately diverges from v1-NSM1-as-run and from Fortran scaling — see F2." |
| F2 | Minor | `src/clearwater_modules_v3/utils/partitioning.py:11-19,38-40`; `src/clearwater_modules_v3/processes/phosphorus.py:36-41,383-390` | Algorithm parity / documentation | The v3 `fdp` utility implements the sorption isotherm; the v1 NSM1 model as actually run (`nsm1/processes.py:290`, wired by `dynamic_variables.py:99`) returns a constant `xr.where(use_TIP,1,0)` and never computes the isotherm. The partitioning docstring compares only against v1 `shared/processes.py` (which v1 NSM1 does not use) and against Fortran; it does not state that v3 *adds* sorption behavior absent from v1 NSM1's runtime. At `kdpo4>0` v3 and v1-NSM1 diverge (v1 stays `fdp=1`; v3 partitions). | This is a deliberate, defensible improvement (the v1 NSM1 wired `fdp` is the documented bug). Add an explicit note in `partitioning.py` and the phosphorus module docstring: "v1 NSM1 as run wires the degenerate `nsm1/processes.py:fdp` (constant 1.0); v3 restores the intended linear-equilibrium isotherm. The two agree only at `kdpo4=0`, the validated regime." Cite as an intentional deviation, not silent. |
| F3 | Minor | `src/clearwater_modules_v3/processes/phosphorus.py:12-13` | Documentation | Module docstring claims kinetics "mirrors v1 `processes.py:1833-2168`." The relevant v1 phosphorus block ends at the `OrgP` integrator (line 2125); 2168 overshoots into unrelated functions. The per-function line citations in the helper docstrings (e.g. `processes.py:308-348`, `1937-1956`) are accurate, but several inline citations point at v2/Appendix-A line numbers not present in scope (e.g. `floating_algae.py:321`, `benthic_algae.py:236`) and cannot be verified from the v1 reference. | Correct the `1833-2168` range to `1833-2125`. Tag the v2/Appendix-A line citations as "v2 overlay reference" so a reader does not expect them in v1 `nsm1/processes.py`. |
| F4 | Minor | `src/clearwater_modules_v3/processes/phosphorus.py:83` | Code hygiene | `import numpy as np` is unused in `phosphorus.py` (only `xr` is used; the `xr.zeros_like` guards at lines 445/457 do not need numpy). | Remove the unused `import numpy as np`. |
| F5 | Minor | `src/clearwater_modules_v3/processes/phosphorus.py:445,457` | xarray robustness | The `use_TIP`/`use_OrgP` False branches return `xr.zeros_like(tip) if hasattr(tip,"dims") else 0.0`. The `hasattr(tip,"dims")` duck-type test is weaker than v1's `xr.where(use_TIP, ..., 0)` pattern: a numpy `ndarray` lacks `.dims`, so a multi-cell numpy input with `use_TIP=False` collapses the rate to a Python scalar `0.0`. Forward Euler `tip + 0.0*dt` still broadcasts correctly, so there is no incorrect result today, but the type contract is inconsistent with the rest of the module (which preserves container type via `sanitize_rate`). | Use `xr.zeros_like` only after a robust array check, or return `0.0` uniformly and rely on broadcasting (simpler and equivalent here since the branch result is purely additive). Document the chosen contract. |
| F6 | Observation (needs verification) | `src/clearwater_modules_v3/processes/phosphorus.py:417-420` | Broadcasting | `dip_from_bed` returns `rpo4_tc / depth` (array) when `use_TIP`, else scalar `0.0`. `rpo4_tc = arrhenius_correction(water_temperature, rpo4_20, rpo4_theta)` — at default `rpo4_20=0` this is the array `0.0 * theta**(T-20)`, i.e. an all-zero DataArray, which is fine. No defect at default; flagged only because the term becomes active and depth-normalized the instant a calibrator sets `rpo4_20>0`, and depth==0 cells would yield inf there (mitigated downstream by `sanitize_rate` at lines 460-461). Confirm `sanitize_rate` runs before the value is persisted (it does — lines 460-461 precede the components dict). | No change required; recorded for traceability. The `sanitize_rate` net-rate guard covers the depth==0 inf path. |
| F7 | Observation | `src/clearwater_modules_v3/parameters/phosphorus.py:17`; `parameter_defaults_corrections.md:45-72` | Default-value deviation (documented) | v3 `vsop = 0.1` m/d deliberately deviates from Fortran/v1-intended `0.01` m/d (v1's literal default is the `999` sentinel). The deviation is documented with a physical-consistency rationale (OrgP detritus should settle near the `vsap=0.15` algal rate) and has regression coverage. This is a defensible, documented choice, not a defect; recorded so it is not lost in a future refactor. | None. Preserve the corrections-doc Section 1.1 rationale and the `test_phase9e_vsop_consistent_with_vsap` regression test. |
| F8 | Observation | `design/clearwater_modules_v3_nsm1_audit_simple_constituents.md:160-169` | Stale audit doc vs source | Simple-constituents audit item 10 states v3 pathogen "replaces `q_solar` with `i0 = PAR(q_solar, Fr_PAR) = q_solar * Fr_PAR`." The current source `pathogen.py:398` uses `i0 = q_solar` directly; the Phase 9.F.B revert is documented in the source docstring (`pathogen.py:385-396`). The audit text is stale; the source now matches v1 exactly (broadband `q_solar`, no `Fr_PAR` scaling). | Update audit item 10 to record that Phase 9.F.B reverted the PAR substitution; v3 pathogen light decay now matches v1 with no calibration-target offset. |
| F9 | Observation | `src/clearwater_modules_v3/processes/pathogen.py:235-237` | Registry key casing | Optional light-extinction inputs are read with mixed casing: `"Solid"` (capitalized) vs `"poc"` / `"ap"` (lowercase). `L()` and the `_LIGHT_DEFAULTS` use `Solid`. If the v3 registry standardizes constituent keys to lowercase (`solid`), the `"Solid"` lookup would silently miss and default to 0 with a one-time warning, dropping the ISS extinction term. Not a defect in the Tier 1 standalone harness (term defaults to 0 by design there); flagged for Phase 5 integration verification. | During Phase 5 integration, confirm the registry key for suspended solids and align the `_get_optional` lookup name with the registered constituent name. |

No Critical or Major findings. At the shipped default parameter set the
phosphorus and pathogen Processes reproduce v1 NSM1 exactly, consistent
with the passing v1-vs-v3 benchmark and the Willamette Santiam-Salem
case study.

---

## 3. Algorithm parity matrix

Notation: MATCH = v3 algebra equals v1 NSM1 as run. IMPROVED(doc) =
deliberate, documented improvement over v1. v1 references are
`src/clearwater_modules/nsm1/processes.py` unless noted.

### Phosphorus

| v3 term (phosphorus.py) | v1 reference | Verdict |
|---|---|---|
| `kop_tc = arrhenius_correction(T, kop_20, kop_theta)` (376-378) | `kop_tc` (1847) | MATCH (arg order `(T, k20, theta)` matches `conversions.arrhenius_correction`) |
| `rpo4_tc = arrhenius_correction(T, rpo4_20, rpo4_theta)` (379-381) | `rpo4_tc` (1864) | MATCH |
| OrgP→TIP hydrolysis `kop_tc * orgp` (394) | `OrgP_DIP_decay` (1879) | MATCH (use_OrgP gate preserved) |
| TIP settling `vs/depth*(1-fdp)*tip` (401) | `TIP_Settling` (1988) | MATCH (formula); fdp source differs — see Section 4 |
| OrgP settling `vsop/depth*orgp` (408) | `OrgP_Settling` (1895) | MATCH (raw `vsop`, no Arrhenius — correct per v1/Fortran); default `vsop` differs by design (F7) |
| Sediment P release `rpo4_tc/depth` (418) | `DIPfromBed` (1969) | MATCH (default `rpo4_20=0` ⇒ silently zero on both) |
| TIP←floating-algae uptake `rpa*growth` (515) | `DIP_ApGrowth` (2018) | MATCH (`rpa = AWp/AWa`) |
| TIP←floating-algae respiration `rpa*resp` (527) | `DIP_ApRespiration` (2003) | MATCH |
| TIP←benthic-algae uptake `rpb*Fb*growth/depth` (543) | `DIP_AbGrowth` (2056) | MATCH |
| TIP←benthic-algae respiration `rpb*Fb*resp/depth` (557) | `DIP_AbRespiration` (2037) | MATCH |
| OrgP←floating-algae mortality (cached `algal_orgp_from_mortality_rate`) (567-571) | `ApDeath_OrgP = rpa*ApDeath` (1912) | MATCH (routing delegated to FloatingAlgae cache; algebra equivalent) |
| OrgP←benthic-algae mortality (cached, incl. `Fw*Fb/depth`) (583-587) | `AbDeath_OrgP = rpb*Fw*Fb*AbDeath/depth` (1935) | MATCH (cache documented to carry `Fw*Fb/depth`) |
| `dtip_dt = orgp_to_tip - tip_settling + dip_from_bed - ap_uptake + ap_release - ab_uptake + ab_release` (435-443) | `dTIPdt` (2091) | MATCH (sign-for-sign identical to v1's `-TIP_Settling + DIPfromBed + OrgP_DIP_decay + DIP_ApRespiration - DIP_ApGrowth + DIP_AbRespiration - DIP_AbGrowth`) |
| `dorgp_dt = ap_mort + ab_mort - orgp_to_tip - orgp_settling` (450-455) | `dOrgPdt` (1956) | MATCH |
| Forward Euler `tip + dtip_dt*dt_days` (329-330) | `TIP`/`OrgP` integrators (2108, 2124) | MATCH (`dt_days = total_seconds/86400` ≡ v1 per-day `dt`) |

### Pathogen

| v3 term (pathogen.py) | v1 reference | Verdict |
|---|---|---|
| `kdx_tc = arrhenius_correction(T, kdx_20, kdx_theta)` (356-358) | `kdx_tc` (3155) | MATCH |
| Natural decay `kdx_tc * px` (359) | `PathogenDeath` (3170) | MATCH |
| Light decay `apx*i0*(1-e^{-kd})/kd*px`, `i0=q_solar` (398-424) | `PathogenDecay = apx*q_solar/(L*depth)*(1-exp(-L*depth))*PX` (3190) | MATCH (Phase 9.F.B reverted the PAR substitution; broadband `q_solar` now used, matching v1; audit doc stale — F8). NaN guard at `kd→0` is a documented robustness improvement, IMPROVED(doc) |
| Settling `vx/depth*px` (433) | `PathogenSettling` (3206) | MATCH |
| `rate = -(natural + light + settling)` (308) | `dPXdt = -Death - Decay - Settling` (3224) | MATCH |
| Forward Euler `px + rate*dt_days` (256-257) | `PX = PX + dt*dPXdt` (3241) | MATCH |
| `apx = 0.017`, `vx = 1.38` defaults | v1 `apx=1`, `vx=1` placeholders | IMPROVED(doc) — Auer & Niehaus (1993) / Chapra (1997) canonical values; documented in `parameters/pathogen.py:7-26` and corrections doc 1.15/1.16 |

Every kinetic term matches v1 algebraically. The only behavioral
divergence from v1 NSM1 is the `fdp` partitioning utility, which is a
deliberate restoration of the intended isotherm and is fully latent at
the validated default `kdpo4=0` (Section 4).

---

## 4. Phosphorus partitioning (fdp) correctness note

**Convention.** In both v1 and v3, `fdp` is the *dissolved* fraction of
total inorganic phosphorus (TIP). The particulate (sorbed) fraction is
`(1 - fdp)`. TIP settling acts only on the particulate fraction:
`TIP_Settling = vs/depth · (1 - fdp) · TIP`
(`phosphorus.py:401`, v1 `processes.py:1988`). This convention is
internally consistent and correctly applied in v3 — the dissolved
fraction does not settle, the sorbed fraction does. The docstring at
`partitioning.py:29-48` correctly labels `fdp` as the dissolved
fraction in `[0,1]`.

**Three-way state of the formula.**

1. v3 source, `utils/partitioning.py:50`:
   `xr.where(use_TIP, 1.0/(1.0 + kdpo4*Solid*1.0e-6), 0.0)`.
   Dimensional check: `kdpo4 [L/kg] · Solid [mg/L] = mg/kg`; multiplying
   by `1e-6 [kg/mg]` yields a dimensionless mass ratio. This is the
   dimensionally correct linear-equilibrium isotherm and matches the
   Fortran scaling (`modGlobalParam.f90` `/ 1.0E6`).
2. v1 NSM1 *as actually run*: `dynamic_variables.py:99` wires
   `processes.fdp`, defined at `nsm1/processes.py:290-304` as
   `return xr.where(use_TIP, 1, 0)` — a constant 1.0, no isotherm,
   `kdpo4` and `Solid` ignored. The isotherm copy in
   `shared/processes.py:271` (with the inverted `/0.000001` unit factor)
   exists but is **not** referenced by v1 NSM1's variable graph.
3. n_p audit P6 / summary: quotes the *old* v3 form `/0.000001` and
   marks it "critical, gated." That text predates the source fix and is
   now stale (F1).

**Consequence.**

- At the shipped default `kdpo4 = 0`: all three give `fdp = 1.0` (for
  `use_TIP=True`), so `(1 - fdp)·TIP = 0`. TIP settling is exactly zero.
  v3 reproduces v1 NSM1 and the benchmarks bit-for-bit. No defect.
- At `kdpo4 > 0`: v3 computes a physically sensible dissolved fraction
  decreasing toward 0 as `kdpo4·Solid` grows (more sorption ⇒ more
  particulate ⇒ more settling). v1 NSM1 *as run* would still return
  `fdp = 1.0` (no settling) because it never evaluates the isotherm.
  v3 therefore *diverges from v1 NSM1's runtime behavior at `kdpo4>0`*,
  but in the correct direction — v1 NSM1's wired `fdp` is the documented
  degenerate stub, and v3 restores the intended physics with the
  dimensionally correct scaling.

**Assessment.** The v3 `fdp` is scientifically correct and an
improvement over both the v1-NSM1-as-run stub and the inverted-unit
`shared/processes.py` copy. The remaining issues are documentation
traceability: (a) the audit doc still describes the pre-fix buggy form
(F1), and (b) neither the partitioning docstring nor the phosphorus
module docstring states that v3 deliberately departs from v1 NSM1's
runtime `fdp ≡ 1.0` for `kdpo4 > 0` (F2). The `kdpo4 = 0` default keeps
the divergence fully latent, so v3 1.0.0 is safe at the validated
regime. The corrections-doc Section 2.2 disposition ("kept at zero by
design; full DIP-solid partitioning is NSM2 territory") remains the
correct framing, but Section 2.2 also mislabels the helper location as
`utils/phosphorus.py` (it is `utils/partitioning.py`) — minor doc fix
folded into F2's recommendation.

---

## 5. Stale-comment / marker list

Markers in scoped files and their fixed-vs-open disposition:

| Marker | Location | Disposition |
|---|---|---|
| `FIXME(phase1-audit): TIP partitioning feature disabled (NSM2 territory)` | `parameters/phosphorus.py:16` (the "1 marker") | NOT stale — correctly open. Per corrections-doc Section 2.2 this tag is an intentional cross-reference to the NSM2-scope deferral, explicitly "not a defect to be fixed in 1.0.0." The `kdpo4=0` default and NSM2 deferral are the right disposition. Keep, but it would read more clearly as `NSM2-DEFERRED(phase1-audit):` rather than `FIXME:` to avoid implying an unaddressed bug. |
| `FIXME(phase1-audit)` reference in module docstring | `processes/phosphorus.py:35` (marker 1 of 2) | NOT stale — accurate cross-reference to the `kdpo4=0` deferral. The surrounding prose correctly states the term is silently zero at default. Keep. |
| `FIXME(phase1-audit)` reference in module docstring | `processes/phosphorus.py:38` (marker 2 of 2) | NOT stale, but the adjacent claim "With `kdpo4=0` … `fdp = 1.0`" is correct only because the source was fixed to the `*1e-6` form. The text is accurate for the current source. Keep; no change needed. |
| n_p audit P6 + summary line 20 + top-concern 5 | `design/clearwater_modules_v3_nsm1_audit_n_p.md:20,45-48,401` | STALE — describes the pre-fix `/0.000001` form and marks v3 critical. Source is fixed. Captured as F1. |
| Simple-constituents audit item 10 | `design/clearwater_modules_v3_nsm1_audit_simple_constituents.md:160-169` | STALE — describes a `q_solar*Fr_PAR` substitution that Phase 9.F.B reverted. Source uses raw `q_solar`. Captured as F8. |
| Corrections-doc Section 2.2 helper path | `parameter_defaults_corrections.md:583` | MINOR STALE — refers to `utils/phosphorus.py`; the helper is `utils/partitioning.py`. Folded into F2. |

No in-scope *source* comment falsely labels a fixed issue as broken. The
stale labels are all in the audit/design docs (F1, F8) and are corrected
relative to the shipped source.

---

## 6. Correctly-deferred list

The following are genuinely deferred and correctly so; they are NOT
findings:

1. **TIP solid-dissolved partitioning at `kdpo4>0`** — deferred to NSM2
   (corrections-doc Section 2.2). Requires a multi-class suspended-solids
   model and coupling to NSM2 sediment diagenesis. v3 1.0.0 keeps
   `kdpo4=0`, matching v1 NSM1 and Fortran behavior at default. Correctly
   deferred. The v3 `fdp` utility is nonetheless implemented with the
   correct isotherm so the NSM2 path inherits a sound primitive.
2. **Full SedFlux sediment-P budget / `use_SedFlux=True`** — deferred to
   NSM2. The `Phosphorus.__init__` guard (lines 217-225) raises
   `NotImplementedError` on `use_SedFlux=True`, preventing the historical
   silent-partial behavior. The `rpo4_20=0` default is the documented
   de-facto gate (corrections-doc Section 2.1). Correctly deferred and
   defensively guarded.
3. **DIP derived post-step variable (`DIP = TIP·fdp`)** — out of scope
   for the Phosphorus Process (n_p audit P14, a minor Fortran-vs-v1
   reciprocal disagreement). Not exercised by the in-scope kinetics.
   Correctly out of scope.
4. **Pathogen light-extinction inputs (`Solid`/`poc`/`ap`)** — optional
   registry reads that default to 0 with a one-time warning so Pathogen
   runs standalone in the Tier 1 harness. This is the intended Phase 3.1
   standalone behavior; full coupling is a Phase 5 integration concern
   (see F9 for the registry-key-casing verification item). Correctly
   deferred to Phase 5 integration.

---

## Appendix: method and evidence

- Read every line of the five scoped v3 files.
- Cross-read v1 phosphorus functions `fdp` (290-304), `kop_tc` (1847),
  `rpo4_tc` (1864), `OrgP_DIP_decay` (1879), `OrgP_Settling` (1895),
  `ApDeath_OrgP` (1912), `AbDeath_OrgP` (1935), `dOrgPdt` (1956),
  `DIPfromBed` (1969), `TIP_Settling` (1988), `DIP_ApRespiration` (2003),
  `DIP_ApGrowth` (2018), `DIP_AbRespiration` (2037), `DIP_AbGrowth`
  (2056), `dTIPdt` (2091), `TIP`/`OrgP` integrators (2108, 2124); and
  v1 pathogen functions `kdx_tc` (3155), `PathogenDeath` (3170),
  `PathogenDecay` (3190), `PathogenSettling` (3206), `dPXdt` (3224),
  `PX` (3241).
- Verified the v1 NSM1 fdp wiring via `dynamic_variables.py:93-100`
  (`process=processes.fdp`) and the `shared/processes.py:257-271`
  isotherm copy is not referenced by v1 NSM1's variable graph.
- Verified v3 utility signatures: `arrhenius_correction(water_temperature,
  reaction_kinetics, theta)` (conversions.py:18), `L(...)` (light.py:13),
  `sanitize_rate` (numerics.py:210, uses `.isnull()|np.isinf`),
  `clip_negative_state` (numerics.py:71).
- Did NOT re-run benchmarks (per scope); relied on the stated passing
  v1-vs-v3 and Willamette Santiam-Salem results plus algebraic parity.
