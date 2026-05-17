# ClearWater Modules v3 — Gold-Standard Release Specification (NSM1 + TSM)

**Date:** 2026-05-16
**Status:** design specification — plan of record; no code changed by this document.
**Scope:** every open TSM and NSM1 issue required to declare v3 the gold standard, with v1 deprecated.

## Sources of record (this spec consolidates, it does not re-derive)

- `design/clearwater_modules_v3_nsm1_review_findings_2026-05-15.md` (NSM1 findings, determinations)
- `design/v3_nsm1_review_2026-05-15/` (per-domain `review_*.md`, `scival_*.md`, `review_SUMMARY.md`)
- `design/clearwater_modules_v3_tsm_line_review_2026-05-15.md` (TSM line review)
- `src/clearwater_modules_v3/parameter_defaults_corrections.md` (16 default corrections + audit dispositions)
- `docs/clearwater_modules_v3_nsm1_research_2_3_ksbod.md` (CBOD-settling literature research)
- Independent QUAL2Kw v6 cross-check: `qual2kw-ver6/design/qual2kw_vs_nsm1_verification_report.md` (corroborates SCI-CB1 and the salinity item from a different baseline)

## 1. Definition of "gold standard"

v3 is the gold standard when **all** hold:

1. **Correctness:** no CRITICAL or science-MAJOR defect open; every kinetic term either matches the literature/CE-QUAL-W2 ground truth or is a *documented, reference-anchored* intentional divergence.
2. **v3 ⊇ v1:** v3 reproduces every v1 capability and is verifiably ≥ v1 in correctness (the parity-≠-correctness payoff items are retained as evidence). v1 is deprecated and archived only after this gate passes.
3. **Documentation:** no stale/forward-looking comment, header, or README that would mislead an external (LimnoTech/sponsor) reviewer; intentional divergences from v1 are explicitly flagged.
4. **Tests:** every gate fix has a regression test that does not share the constant/path under test between the v1-mirror and v3 sides; the highest-value deferred conservation tests are landed or explicitly risk-accepted.
5. **Deferrals are clean:** NSM2-scope items are guarded, documented as deferrals, and not mislabeled as defects.

## 2. Status snapshot

| Component | State | Gate work |
|---|---|---|
| **TSM** | Near-clean. Review verdict: nothing blocks current v3 TSM; algorithms match v1 or improve deliberately; no live TODO in `temperature.py`; 88 tests pass. | Doc hygiene + the deferred MMS conservation test + wind_c divergence callout. **No code defect.** |
| **NSM1** | 1 CRITICAL, a MAJOR science set, a documentation-staleness cluster. v1→v3 port faithful; v3 independently fixes ~9 genuine v1 science defects. | Workstreams A–D below. |

v1-only defects (the v1 `fdp` stub, v1 `SOD_20=999` sentinels, v1 θ-transposition, v1 DIC unit error, v1 `vb` 1460× bug, v1 pathogen placeholders) are **Closed by v1 deprecation** — v3 is already correct; retained as v3 ⊇ v1 evidence (findings doc §5).

---

## 3. Workstream A — Correctness blockers (hard gate)

### A1 · NSM1-CA-1 (CRITICAL) — alkalinity algal/benthic coupling uses raw weights

- **Source:** findings §3; `review_SUMMARY.md` §2; scival_carbon_alkalinity.
- **Current state:** `processes/alkalinity.py:362,386,411,441` bind `rca = self.AWc` (=40), `rcb = self.BWc` (=40); composition block `:202-206` never pulls `AWa`/`BWd`; docstring `:61-62` and comments `:202-203` rationalize the wrong convention as intentional.
- **Root cause:** the Phase 9.B fix `rca = AWc/AWa` reached `carbon.py:495-496` and `dox.py` but not `alkalinity.py`. Magnitude is pure arithmetic: floating-algae alkalinity flux ×1000 (40/0.04), benthic ×100 (40/0.4). External ground truth: CE-QUAL-W2 `water-quality.f90:3164-3166` routes via the intensive carbon fraction; v1 is correct, v3 is wrong.
- **Design / fix:**
  1. In `alkalinity.py`, derive `rca = self.AWc / self.AWa` and `rcb = self.BWc / self.BWd` once (mirror `carbon.py:495-496`).
  2. Add `AWa` (from `ALGAE_DEFAULTS`) and `BWd` (from `BALGAE_DEFAULTS`) to the composition block `:202-206`.
  3. Rewrite the root-cause docstring `:61-62` and comments `:202-203` to state the correct intensive-ratio convention.
  4. Correct `clearwater_modules_v3_nsm1_audit_simple_constituents.md` §21–§23 (currently marks the wrong v3 form "Match").
- **Affected files:** `processes/alkalinity.py`; `parameters/alkalinity.py` (no value change — `r_alkaa`/`r_alkba` already correct); audit doc.
- **Tests:** new regression test that feeds the v1-mirror and v3 paths **independently computed** constants (not a shared symbol); assert floating-algae alkalinity flux ≈ 9.172e-4 mg-CaCO₃/L/d (not 9.172e-1) at default stoichiometry; an alkalinity mass-conservation assertion under a synthetic bloom (no silent negative-clip).
- **Acceptance:** v3 alkalinity flux matches the resolved-ratio value to machine precision; no negative-clip events in the bloom test; audit doc corrected.
- **Depends on:** none. **Do first.**

### A2 · NSM1-SCI-N1 (MAJOR, science) — denitrification alkalinity coefficient 4× too high

- **Source:** findings §4 SCI-N1; scival_nitrogen.
- **Current state:** `r_alkden = 4/14/1000` eq/mg-N. Shared by Fortran (`modAlkalinity.f90:54`), v1, and v3 — the canonical "wrong at all stages" case the parity pass structurally could not catch.
- **Correct value:** `r_alkden = 1/14/1000` (1 eq alkalinity per mol N denitrified) — CE-QUAL-W2 `water-quality.f90:3157`; Stumm & Morgan.
- **Design / fix:**
  1. **Pre-req investigation:** confirm what NSM1 `denitrification_flux_rate` physically represents (mol-N basis, sign, whether it already nets sediment denitrification) so the 1:1 eq:mol mapping is applied to the right quantity.
  2. Set `r_alkden = 1/14/1000` in `parameters/alkalinity.py`, with an inline citation block (CE-QUAL-W2 line ref + Stumm & Morgan) marking it a **deliberate, reference-anchored divergence from upstream NSM1 Fortran**.
  3. File/report the upstream Fortran error (`modAlkalinity.f90:54`) to the NSM1 Fortran maintainers.
- **Affected files:** `parameters/alkalinity.py`; `processes/alkalinity.py` (citation/comment); upstream-defect note.
- **Tests:** closed-system benchmark — a box with only denitrification active; assert Δalkalinity = 1 eq per mol N denitrified (not 4); regression guard on `r_alkden`.
- **Acceptance:** closed-system alkalinity balance correct to 1 eq/mol-N; divergence-from-Fortran documented; upstream report filed.
- **Depends on:** the `denitrification_flux_rate`-semantics investigation (cheap; do alongside A1).

---

## 4. Workstream B — Investigation gate (must resolve before sign-off)

### B1 · NSM1-SCI-A3 — light-limitation irradiance basis (PAR vs total shortwave)

- **Source:** findings §4 SCI-A3 + §6 (Santiam–Salem determination).
- **Issue:** `limit_light` compares `KL` (default 10 W/m², a **PAR**-scale value) directly against `registry.get_at_time("solar_radiation")` in the depth-averaged Beer–Lambert half-saturation form. If the coupling driver supplies **total shortwave** while `KL` is PAR-based (PAR ≈ 0.43–0.47·SW), light under-limits and algal growth is over-predicted **~30–60% wherever light limitation is active**. Shared by v1 and v3 → invisible in the v1↔v3 benchmark. **Leading hypothesis for the Santiam–Salem algae mismatch** (CA-1 is *ruled out* as the cause — NSM1 algae depends on N/P/light/T, not alkalinity).
- **Decisive check (required, not yet performed):** trace the Santiam–Salem case-study input that is written into `solar_radiation` (variable + units), and the `KL`, `light_limitation_option`, `growth_rate_option` settings on that same basis. The **direction** of the observed miss discriminates: over-prediction ⇒ light/PAR path; under-prediction ⇒ multiplicative-nutrient form.
- **Design / fix (conditional on the trace):** if the driver supplies total shortwave, either (a) scale to PAR at the `limit_light` boundary (`PAR = Fr_PAR · SW`, `Fr_PAR ≈ 0.45`, single documented constant) or (b) re-express `KL`'s default/units on the shortwave basis — pick one convention, document it, and make the irradiance basis explicit in the docstring and a registry-units assertion. Re-run the Santiam–Salem case study.
- **Affected files:** `processes/floating_algae.py` / `benthic_algae.py` `limit_light`; the coupling driver / case-study input mapping; docs.
- **Acceptance:** irradiance convention proven and made explicit (docstring + runtime units guard); Santiam–Salem algae re-run and the determination (confirm/refute SCI-A3) recorded in the findings doc §6.
- **Depends on:** access to the Santiam–Salem run configuration. **Gate item because the model's largest potential growth bias is unresolved until the trace is done.**

---

## 5. Workstream C — MAJOR science / quality (gate)

### C1 · NSM1-SCI-A2 — `f_pocp` / `f_pocb = 0.5` too low

- **Source:** findings §4 SCI-A2; corrections-doc precedent (`mu_max_20`).
- **Issue:** fraction of algal/benthic mortality carbon routed to POC defaults to 0.5; CE-QUAL-W2 `APOM` ≈ 0.6–0.9 (~0.8), QUAL2K/Bowie/Chapra treat dead-algal C as predominantly particulate. v1's 0.9 is defensibly correct; v3's 0.5 mis-routes ~40% of mortality C to DOC, perturbing DOC→DIC and DO demand (not algal biomass).
- **Design / fix — author decision required (D-track):** restore `f_pocp = f_pocb ≈ 0.8–0.9` (v1 was 0.9; CE-QUAL-W2 ~0.8), **or** keep a value with a cited, documented deliberate deviation entered in `parameter_defaults_corrections.md` (the `mu_max_20` precedent).
- **Tests:** carbon-partition unit test asserting the chosen split; DOC/DO sensitivity regression.
- **Acceptance:** value set with literature citation in the corrections doc; no undocumented inline default.

### C2 · NSM1-SCI-CB1 — `ksbod_theta = 1.047` + 1/depth form mismatch  *(= QUAL2Kw cross-check F3)*

- **Source:** findings §4 SCI-CB1; `parameter_defaults_corrections.md` §2.3 (RESOLVED for the *default value* only); `docs/clearwater_modules_v3_nsm1_research_2_3_ksbod.md` (fix already specified); independently corroborated by the QUAL2Kw v6 cross-check (report F3).
- **Issue:** CBOD settling applies θ = 1.047 (the *oxidation* Arrhenius coefficient); Bowie/QUAL2E specify **θ = 1.024** for settling (CE-QUAL-W2 applies no Arrhenius to `CBODS`). Separately, `processes/cbod.py:240` divides by depth (`ksbod_tc/depth·cbod`), treating the parameter as a velocity (m/d); Fortran `modCBOD.f90:114` and QUAL2E use a **1/d first-order rate** (no depth division). Dormant at the shipped default `ksbod_20 = 0`; silently wrong the moment a user calibrates nonzero CBOD settling (a usage the corrections doc documents as legitimate).
- **Design / fix (per the research doc's recommendation):**
  1. `ksbod_theta = 1.024` in `parameters/cbod.py` (or a documented deliberate deviation).
  2. Change `processes/cbod.py` to `settling_rate = ksbod_tc * cbod` (no depth division); docstring units → "1/d at 20 °C".
  3. Replace any residual unit-label "m/d" with "1/d"; cross-reference `parameter_defaults_corrections.md` §2.3 and the research doc.
- **Tests:** nonzero-`ksbod_20` regression vs the Fortran 1/d form (must agree, no `1/depth` factor); θ regression guard.
- **Acceptance:** nonzero CBOD-settling result equals the Fortran/QUAL2E 1/d convention; θ = 1.024; docs consistent.

### C3 · NSM1-DOX-F1 — DO-saturation salinity omission  *(= QUAL2Kw cross-check F7 / audit C6)*

- **Source:** findings §4 DOX-F1; audit_summary C6; review_dox F2 (MINOR).
- **Issue:** the DO-sat salinity correction is silently omitted — correct for fresh water (factor 1.0) but with no salinity input, no guard, no docstring statement of the freshwater assumption. CE-QUAL-W2 includes it (APHA-exact); ~18% DOsat overstatement at 35 ppt.
- **Design / fix:** add an explicit freshwater-assumption statement to `dox_sat_apha` and the module docstring; add a runtime guard/warning if a nonzero salinity is present in the registry (deliberate-deferral note referencing audit C6). No numeric change for freshwater.
- **Acceptance:** freshwater assumption documented at the function; brackish input cannot pass silently.

### C4 · NSM1-DOX-F2 — silent zero atmospheric reaeration

- **Source:** findings §4 DOX-F2; `review_SUMMARY.md` §3.
- **Issue:** with `hydraulic_reaeration_option = 1` and `kah_20_user = 0.0`, the user gets **silent zero atmospheric reaeration**. v3-divergent from CE-QUAL-W2, which enforces a `MINKL` minimum-reaeration floor on every branch (`gas-transfer.f90`). Not triggered by the shipped default (option 5).
- **Design / fix:** emit a one-time warning when option 1 is selected with no coefficient; add an opt-in minimum-reaeration floor following the W2 `MINKL` precedent (off by default to preserve parity; documented).
- **Acceptance:** the silent-zero path warns; opt-in floor available and documented; default behavior unchanged.

### C5 · NSM1-SCI-A1 — alkalinity↔algae coupling is carbon-routed, process is N-driven  *(stage with NSM2)*

- **Source:** findings §4 SCI-A1.
- **Issue:** even after CA-1's magnitude is fixed, the alkalinity-algae term is a carbon-routed approximation of an intrinsically nitrogen-driven process (CE-QUAL-W2 / Stumm & Morgan compute it from the N flux); self-consistent only at exact Redfield N:C, and NSM1's `AWn/AWc` is ~2% off Redfield.
- **Design / disposition:** a true formulation limitation beyond CA-1. The deeper fix is to reimplement on the N-flux basis as W2 does. **Reasonable to stage with NSM2** — record as a documented known limitation now (docstring + corrections doc), schedule the N-flux reimplementation in the NSM2 work plan. Not a v3-1.0 gate blocker provided it is explicitly documented as a known limitation.

---

## 6. Workstream D — Documentation & engineering gate

### D1 · NSM1 documentation-staleness cluster (explicit review ask)

- **A1 module headers:** rewrite `floating_algae.py:1-46` and `benthic_algae.py:1-20` from forward-looking work plans to **as-implemented records** (state each bug *corrected*, cite the fixing phase).
- **README:** rewrite `src/clearwater_modules_v3/README.md` — remove the v2-re-export claims (Riverine/BenthicAlgae/FloatingAlgae/Nitrogen/Process/ProcessFactory are all v3-native in-tree), fix "2 of 18 MAJOR resolved" → the authoritative 17 of 18.
- **Minor doc-fidelity:** clear the stale shadow-method/test refs in `carbon.py`; the N2 docstring naming non-existent `denitrification_rate`; `nitrate_uptake_*` docstrings leading with the old broken implementation; `utils/reaeration.py:168-169` stale marker; stale audit-doc entries (`fdp` pre-fix; simple-constituents §21–§23). Full inventory: `review_SUMMARY.md` §6.
- **Acceptance:** **before any external (LimnoTech/sponsor) review packet ships**, no header/README/comment states resolved work as outstanding; minor cluster cleared.

### D2 · TSM documentation hygiene

- **Source:** TSM review §3, §4.
- Annotate `clearwater_modules_v3_tsm_audit_2026-05-05.md` (Q6, lines 31, 222–224) and the TSM gap analysis with a one-line "superseded by the wind-function specification (`wind_c = 2.0`)" where they state the old `0.3/1.5/3.0` defaults.
- **wind_c divergence callout:** add an explicit, prominent statement in sponsor/LimnoTech-facing TSM materials that **v3 `wind_c = 2.0` ≠ v1 `3.0` by design** (CE-QUAL-W2 `CFW=2.0`; QUAL2K Brady-Graves-Geyer 2.0; all seven W2 example cases; Santiam–Salem validated). Not a defect — a reviewer-awareness requirement.
- **Acceptance:** no sponsor-facing TSM document cites the superseded `3.0` default without the supersession note.

### D3 · Test gate

- **TSM MMS energy-conservation test (highest-value deferred test):** land the method-of-manufactured-solutions end-to-end energy-conservation test (TSM review §5.2 — it would have caught F2 before it shipped), **or** explicitly risk-accept with sign-off and a tracked v3.x ticket.
- **NSM1:** the A1 (CA-1) and A2 (SCI-N1) regression tests must not share the constant/path under test between the v1-mirror and v3 sides (the structural reason CA-1/SCI-N1 escaped). Add a closed-system alkalinity benchmark covering CA-1 + SCI-N1 jointly.
- **Inherited-zero-defaults user note:** add one consolidated note (params docs) that `kdpo4=0`, `ksbod_20=0`, `rpo4_20=0`, `rnh4_20=0`, `vno3_20=0` are intentional inherited-lineage defaults that silently disable P sorption, CBOD settling, and sediment-flux release — so users are not surprised.
- **Acceptance:** gate fixes have non-shared-path regression tests; MMS test landed or risk-accepted with a ticket.

---

## 7. Workstream E — Decisions required (author/PI)

| ID | Decision | Options | Recommendation |
|---|---|---|---|
| E1 | SCI-A2 `f_pocp/f_pocb` | restore 0.9 (v1) · set ~0.8 (CE-QUAL-W2) · keep 0.5 with cited rationale | Set ~0.8–0.9 with a corrections-doc citation |
| E2 | SCI-A1 timing | fix N-flux basis now · stage with NSM2 | Stage with NSM2; document as a known limitation now |
| E3 | SCI-CB1 form | match Fortran 1/d (recommended by research doc) · keep v3 m/d with rationale | Match Fortran 1/d; θ = 1.024 |
| E4 | TSM 273.15 vs 273.16 | hold 273.15 · reconcile with LimnoTech v2 (273.16) | Decide only if v2 numerical parity is renegotiated (TSM review §5.3) |
| E5 | MMS test | land now · risk-accept to v3.x | Land before "gold standard" claim if feasible |

---

## 8. Explicitly out of gate — deferred to NSM2 (document only, NOT defects)

Guarded and documented as deferrals (findings §5; corrections §2.1/§2.2): `Nitrogen use_SedFlux=True` → `NotImplementedError`; full carbonate/pH solver; SOD-derived DIC sediment release; alkalinity DOX-Monod single-source routing; POM→DOC source; multi-group algae; `kdpo4` phosphorus partitioning (TIP sorption — needs multi-class solids + NSM2 sediment-flux coupling). **No action beyond confirming each remains guarded and labeled a deferral, not a defect.** SCI-A1's deeper N-flux reimplementation (C5) also stages here.

## 9. Future feature backlog (NOT a gold-standard gate)

From the QUAL2Kw v6 cross-check (`qual2kw-ver6/design/qual2kw_vs_nsm1_verification_report.md` §4) — capabilities present in QUAL2K/literature, absent in NSM1, offered as future candidates with benefit:cost. **Lacking these is not a deficiency.** Near-term high-B:C: selectable per-process O₂-attenuation form; inorganic-carbon (CO₂) growth limitation. Larger/scope-dependent: slow/fast CBOD split; Droop variable-quota benthic algae; mechanistic Di Toro SOD (NSM2). Trivial optional polish: `dox_sat_apha` leading constant `−139.34410 → −139.34411` for APHA-citation fidelity (numerically ~1e-7, immaterial; v1↔v3↔Fortran are byte-identical — not a regression, do only as cosmetic literature-fidelity).

## 10. Sequencing & exit criteria

**Order** (extends findings-doc §7): A1 (CA-1) → A2 (SCI-N1) + the `denitrification_flux_rate` investigation → B1 (Santiam–Salem irradiance trace; apply SCI-A3 fix if confirmed) → C1 (SCI-A2) → C2/C3/C4 (SCI-CB1, DOX-F1, DOX-F2) → C5 documented (SCI-A1, stage NSM2) → D1/D2/D3 documentation + tests (D1 before any external packet) → E decisions resolved.

**Exit criteria (gold standard declared):** §1 (1)–(5) all satisfied; CA-1 and SCI-N1 fixed with non-shared-path tests; SCI-A3 resolved (fixed or refuted with evidence); SCI-A2/CB1/DOX-F1/DOX-F2 fixed; SCI-A1 documented as a scheduled NSM2 limitation; NSM1 doc cluster + TSM doc hygiene cleared; MMS test landed or risk-accepted; v1 archived.

---

*This specification is the plan of record for the v3 gold-standard release. It consolidates and does not supersede the per-domain analysis in the sources of record; where it cites a fix already researched (SCI-CB1), follow that research doc's detail.*
