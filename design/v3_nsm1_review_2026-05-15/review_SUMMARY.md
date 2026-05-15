# v3 NSM1 Line-Level Code & Science Review — Consolidated Summary

- **Date:** 2026-05-15
- **Branch:** `streaming` @ `54f2b12`
- **Scope:** Full v3 NSM1 process kinetics (~9,300 LOC), parameters, utilities, model
  framework, config. Line-level read of every in-scope file, cross-referenced
  against the v1 reference (`src/clearwater_modules/nsm1`) and the authoritative
  design/audit doc trail.
- **Method:** Six parallel domain reviewers (water-quality-model-source-code-reviewer,
  Opus) + lead synthesis & independent verification of the CRITICAL finding.
- **Not in scope:** numerical re-validation. The v1-vs-v3 benchmark and the
  Willamette Santiam–Salem case study already pass; this review is code,
  science, and documentation **correctness**.

Per-domain reports (read these for full findings tables and parity matrices):

| Domain | Report |
|---|---|
| Nitrogen + N2 | `review_nitrogen_n2.md` |
| Floating + Benthic Algae | `review_algae.md` |
| Carbon + Alkalinity | `review_carbon_alkalinity.md` |
| Dissolved Oxygen | `review_dox.md` |
| CBOD + POM | `review_cbod_pom.md` |
| Phosphorus + Pathogen | `review_phosphorus_pathogen.md` |
| Framework / utils / config / temperature | `review_framework_utils.md` |

---

## 1. Severity rollup

| Domain | CRITICAL | MAJOR | MINOR | OBS | Confidence |
|---|---:|---:|---:|---:|---|
| Nitrogen + N2 | 0 | 0 | 7 | 8 | High |
| Algae | 0 | 2 | 7 | 6 | High |
| Carbon + Alkalinity | **1** | 0 | 4 | 4 | High |
| DOX | 0 | 1 | 4 | 4 | High |
| CBOD + POM | 0 | 1 | 2 | 3 | High |
| Phosphorus + Pathogen | 0 | 0 | 5 | 4 | High |
| Framework / utils | 0 | 0 | 5 | 9 | High |
| **Total** | **1** | **4** | **34** | **38** | — |

**Bottom line:** The xarray refactor is complete and the v1→v3 algorithm port
is faithful across the entire model. Every traced kinetic term is `MATCH` or a
*documented* intentional improvement; the historically buggy modules
(FloatingAlgae Bug #4/#13/#14/#15/#16, Nitrogen Phase 9.A/B/E, Carbon/DOX
Phase 9.B) are verifiably fixed. A second pass validated the kinetics against
**CE-QUAL-W2 v2026.02 (local source), QUAL2K/QUAL2E, WASP, and the rate
literature (Bowie et al. 1985, Chapra 1997, Stumm & Morgan)** — because
v1-parity is not proof of correctness if v1 is wrong. That pass **confirmed
CA-1 against external ground truth** (v1 is correct here, v3 is wrong),
surfaced **one new defect shared by v1 *and* v3 that the parity pass
structurally could not catch** (SCI-N1, the 4× denitrification-alkalinity
coefficient), upgraded three latent/calibration-regime issues, and — the
parity-≠-correctness payoff — **independently verified that v3 fixed ~9 genuine
v1 science defects** (see §9). Net: v3 is materially *more* correct than v1.
**One CRITICAL (CA-1) and a set of MAJOR latent issues are open**, plus a
documentation-staleness cluster that would mislead an external reviewer.

---

## 2. CRITICAL — must fix before any algal-active production run

### CA-1 — `alkalinity.py`: algal/benthic alkalinity coupling uses raw stoichiometric weights, not resolved ratios (1000× / 100× flux error)

- **Files:** `src/clearwater_modules_v3/processes/alkalinity.py:362` (`rca = self.AWc`),
  `:411` (`rcb = self.BWc`); also `:386`, `:441` and the composition block `:202-206`.
- **Defect:** v1 binds the *resolved concentration ratios* `rca = AWc/AWa = 0.04`
  mg-C/µg-Chla and `rcb = BWc/BWd = 0.4` mg-C/mg-D (v1 `processes.py:337`/`:776`,
  registered as dynamic variables; used identically by v1 alkalinity
  `processes.py:3340/3359/3386/3408`). v3 `alkalinity.py` instead uses the **raw
  weights** `AWc = 40` and `BWc = 40` directly, and never composes `AWa`/`BWd`.
- **Independent verification (lead):** the alkalinity stoichiometric constants
  `r_alkaa`/`r_alkba` in `parameters/alkalinity.py:14-19` are **byte-identical**
  to v1 `constants.py:55-60`, so the *only* divergence is the rca/rcb factor.
  With v3 defaults `AWc=40, AWa=1000` → error factor = 40 / 0.04 = **1000×**
  (floating-algae growth + respiration alkalinity terms); `BWc=40, BWd=100`
  → 40 / 0.4 = **100×** (benthic-algae terms).
- **Why it escaped:** (a) The Phase 9.B audit fix that introduced
  `rca = self.AWc/self.AWa` reached `carbon.py:495-496` and `dox.py` but **not**
  `alkalinity.py`. (b) The `alkalinity.py` inline comments (`:202-203`,
  `:61-62`, docstring `:18/:24`) *document the wrong convention as intentional*
  ("rca = AWc, mg-C/µg-Chla"), so it reads as deliberate. (c) The parity test
  feeds the same wrong constant into both the v1-mirror and v3 paths, and the
  simple-constituents audit doc (§21–§23) marked the v3 form "Match" — masking it.
- **Impact:** Alkalinity (hence pH / carbonate system when coupled) is correct
  in the Santiam–Salem validation because that signal is dominated by the
  nitrification/denitrification terms (which use the correct `r_alkn`/`r_alkden`).
  Under productive / algal-bloom conditions the algal coupling term dominates and
  alkalinity is over-driven by 2–3 orders of magnitude. **Latent, not benign.**
- **Fix:** Mirror the `carbon.py:495-496` pattern — derive
  `rca = self.AWc / self.AWa` and `rcb = self.BWc / self.BWd` once; add `AWa`
  (from `ALGAE_DEFAULTS`) and `BWd` (from `BALGAE_DEFAULTS`) to the
  `alkalinity.py` composition block `:202-206`; correct the misleading comments
  and docstring; add a regression test that does **not** share the constant
  between sides; correct the simple-constituents audit doc §21–§23.

---

## 3. MAJOR

| ID | File | Issue | Action |
|---|---|---|---|
| A1 | `floating_algae.py:1-46`, `benthic_algae.py:1-20` | **Stale module headers.** Written as forward-looking work plans ("apply the fixes", "Bug #13: implement (was returning 0)"). The work is *done and verified*. A reader concludes the module is still broken — directly the failure mode the user asked to catch. | Rewrite headers as as-implemented records (state bugs are *corrected*, cite the fixing phase). |
| A2 | `floating_algae.py` / `benthic_algae.py` (mortality C-split) | **Undocumented `f_pocp`/`f_pocb` = 0.5 vs v1's 0.9.** Silently routes ~40% of algal/benthic mortality carbon from POC to DOC. No entry in `parameter_defaults_corrections.md`. | Author decision: restore 0.9, or document it like the `mu_max_20` precedent. |
| DOX-F1 | `dox.py:754-761` | Zero-reaeration short-circuit + corrected `kah_20_user=0.0` default: a user who sets `hydraulic_reaeration_option=1` (v1/Fortran convention) without a value gets **silent zero atmospheric reaeration**. | Emit a one-time warning when the option is selected with no coefficient. |
| CP-F1 | `pom.py:400`, `REGISTRY_DIAGNOSTICS` | `"pom_settling_rate"` registry slot is populated with `rate_burial = vb*pom/h2` (sediment burial), not water-column settling. Budget/diagnostic consumers misattribute the loss pathway. | Rename slot to `"pom_burial_rate"` (grep consumers first). |

---

## 4. The user's five questions — answered

### 4.1 Is the xarray refactoring complete?
**Yes.** Across all seven domains: no per-cell Python loops, no `== np.nan`
(all NaN/inf logic uses `.isnull()` / `xr.where` / `np.isnan`), no
array-truthiness `if` on per-cell data, container-aware clip/sanitize guards.
Minor latent broadcasting hazards only (scalar early-returns that drop
DataArray coords — Algae A4/A5, Phosphorus F5; not triggered in current runs).
One inf-safety gap: Algae A6 (`_sanitize_cache` is NaN-only; `algae/depth` →
`+inf` at `depth==0` passes through).

### 4.2 Do the algorithms match v1 (or improve on it intentionally)?
**Yes, with one exception (CA-1).** Every traced kinetic term across N, P,
algae, C, DOX, CBOD, POM, pathogen is algebraically `MATCH` to v1 sign-for-sign,
or a *documented* improvement (cite the audit / corrections doc). Verified
intentional improvements include: Phase 9.E nitrogen θ transposition (v3 matches
Fortran `modNitrogen.f90`; v1 had the pairs swapped — v3 is *more* correct),
Phase 9.B carbon/DOX rca/rcb stoichiometry, Phase 9.E DIC mg-C/L unit
reconciliation, restored CBOD→DIC source, removed spurious POC-hydrolysis
DOX-Monod, canonical pathogen `apx`/`vx` defaults, the 16 parameter-default
corrections. The Forward-Euler-in-days integrator equals v1's per-day additive
step everywhere. **Only CA-1 is an unintended discrepancy.**

### 4.3 Is the source well, appropriately, correctly documented?
**Mostly yes, with a documentation-staleness cluster** (the user's specific
concern). Resolved bugs are *not* falsely labeled broken in the live kinetic
**source comments** — except where the comment rationalizes CA-1. The
staleness is concentrated in: module headers (A1), the package `README.md`
(Sec. 5 below), and several stale cross-references to deleted shadow
methods/test files (Carbon CA-2/CA-3/CA-4; Nitrogen F1/F4/F6 — incl. an N2
docstring naming a non-existent attribute `denitrification_rate`).

### 4.4 Are there fixed TODO/error comments not updated to current status?
**Yes — these are findings, list below (Sec. 6).** The largest are A1 (algae
headers) and the README. Most numbered `Bug #N` / `Phase X` inline markers in
the live code are accurate as-implemented records; the stale ones are
enumerated per-domain.

### 4.5 Any unaddressed TODOs / known issues? (NSM2 deferral is OK)
**Yes, and they are correctly classified.** Genuinely-deferred-to-NSM2 items —
`Nitrogen use_SedFlux=True` NotImplementedError, full carbonate/pH solver,
SOD-derived DIC sediment release, sediment diagenesis generally — are properly
guarded and documented as deferrals, **not defects** (matches the user's "OK if
saved for NSM2"). Genuine open framework items needing a decision: `riverine.py:101`
registers `wetted_surface_area` under the name `depth` (m² used as m) with an
honest TODO — only reached inside `if model.has_process("FloatingAlgae")`;
potential blocker for the NSM1+Riverine coupled path (open question). CBOD-F3:
`use_DOX=True` hardcoded post-merge, not overridable via `parameters`.

---

## 5. README staleness (documentation finding)

`src/clearwater_modules_v3/README.md` predates NSM1 pattern-alignment Phases 6–10
and is now internally contradictory:

- States Riverine / BenthicAlgae / FloatingAlgae / Nitrogen and Process /
  ProcessFactory are **v2 re-exports** — they are all v3-native in-tree; the
  top-level `__init__.py` says v2 was removed entirely.
- Claims "2 of 18 MAJOR resolved, Phase 5 in progress" — the authoritative
  triage records 17 of 18 resolved (M4 deferred) and the code confirms it.
- Stale `review_findings.md` M4 "STILL DEFERRED" label and a utilities-audit
  "re-exports 273.16" claim (code uses SI 273.15, runtime-verified).

---

## 6. Stale-comment / stale-doc inventory (the user's explicit ask)

| Where | Stale claim | Reality |
|---|---|---|
| `floating_algae.py:1-46` | "apply the fixes", "Bug #13/#14: implement (was returning 0)" | All fixed & regression-tested |
| `benthic_algae.py:1-20` | forward-looking bug plan | All fixed |
| `floating_algae.py:376` | `# TODO: implement` | Annotates an already-implemented line |
| `alkalinity.py:202-203,61-62` | "rca = AWc, mg-C/µg-Chla" | **Wrong convention rationalized** (CA-1) |
| `nitrogen.py` N2 docstring (`n2.py:133`) | reads `nitrogen_process.denitrification_rate` | Code reads `denitrification_flux_rate`; the named attr doesn't exist |
| `nitrogen.py` class docstring | "v2 … Phase 2.B fixes applied" | Phase 9.A.2 critical fixes not reflected at class level |
| `nitrate_uptake_*` docstrings | lead with the *old broken* implementation | Current behavior is correct; doc buries it |
| `carbon.py` (CA-2/3/4) | refs `_change_legacy_inline`, `test_*_helper_vs_inline.py`, `_ka_tc` workaround | Methods/files deleted; workaround moved upstream |
| `README.md` | v2 re-exports; "2 of 18 MAJOR" | All v3-native; 17 of 18 |
| `utils/reaeration.py:168-169` | pending cross-reference | The `np.select` dim fix is done |
| audit docs (P F1/F8, simple-constituents §21-23) | pre-fix `fdp` form; "Match" for CA-1 | Stale relative to fixed source / wrong |

Dead-but-callable back-compat code (drift risk, not wrong today): Nitrogen F3
(4 dead methods incl. phantom `ammonium_decay_nitrate`), Algae A9
(`_cache_benthic_mortality_rates` duplicate), CBOD/POM F2 (legacy `POM.rate()`
side-effect duplicate).

---

## 7. Correctly deferred to NSM2 / future phases (NOT defects)

`Nitrogen use_SedFlux=True` → `NotImplementedError` (NSM2 diagenesis); full
carbonate/pH solver (NSM2); SOD-derived DIC sediment-release fallback;
alkalinity DOX-Monod single-source routing; POM→DOC source; multi-group algae
(separate design spec). All explicitly guarded and documented.

---

## 8. Science validation vs CE-QUAL-W2 / QUAL2K / QUAL2E / WASP / literature

Second pass. Question: *is the formulation/coefficient scientifically correct
against validated references, regardless of v1?* CE-QUAL-W2 v2026.02 source was
read line-by-line locally (authoritative & checkable); QUAL2K (Chapra, Pelletier
& Tao 2008), QUAL2E (Brown & Barnwell 1987), WASP, Bowie et al. 1985, Chapra
1997, Stumm & Morgan cited from domain knowledge as corroborating references.
Per-domain reports: `scival_algae.md`, `scival_nitrogen.md`, `scival_dox.md`,
`scival_carbon_alkalinity.md`, `scival_om_cbod_pom_p_path.md`,
`scival_temperature_defaults.md`.

### 8.1 CA-1 — authoritatively confirmed (v1 correct, v3 wrong)

CE-QUAL-W2 `ENTRY ALKALINITY` (`water-quality.f90:3164-3166`) routes algal
alkalinity coupling through nitrogen fluxes derived from biomass via the
*intensive* carbon fraction, never a raw weight; Stumm & Morgan Table 4.5
gives the same stoichiometry W2 quotes verbatim (`:3152-3159`). Correct
floating-algae flux ≈ **9.172e-4 mg-CaCO₃/L/d**; v3 produces **9.172e-1**
(**1000×**; benthic **100×**). v1's resolved `rca=AWc/AWa=0.04` recovers the
right order; v3's raw `AWc=40` does not. **Ruling: CA-1 is a true v3 science
error that v1 gets right.** It also violates alkalinity mass conservation under
bloom conditions (forces silent negative-clips).

### 8.2 NEW — defects the v1-parity pass structurally could not catch

| ID | Sev | Finding | Evidence | Class |
|---|---|---|---|---|
| **SCI-N1** | **MAJOR** | `r_alkden = 4/14/1000` eq/mg-N overstates **denitrification** alkalinity production by **4×**. Correct = `1/14/1000` (1 eq alk per mol N). | CE-QUAL-W2 `water-quality.f90:3157` (`+1.*NO3D...`); Stumm & Morgan. Upstream NSM1 Fortran `modAlkalinity.f90:54` error → propagated to v1 **and** v3 identically. | **v1 & v3 BOTH wrong** (parity passed it because v1==v3). Carbon reviewer wants author confirmation of what `denitrification_flux_rate` represents before patching; nitrogen reviewer rates it confirmed on W2 line-level evidence. |
| **SCI-A1** | **MAJOR** | NSM1's alkalinity↔algae term is **carbon-routed**; the process is intrinsically **nitrogen-driven**. Self-consistent only at exact Redfield N:C — and NSM1's own `AWn/AWc` is ~2 % off Redfield, so any non-Redfield calibration silently desyncs alkalinity from the N balance. | W2 + Stumm & Morgan compute it from the N flux. | Conceptual; *additional* to CA-1's magnitude bug. Deeper fix = reimplement on N-flux basis like W2. |
| **SCI-CB1** | **MAJOR (latent)** | `ksbod_theta = 1.047` applies the *oxidation* Arrhenius θ to CBOD *settling*; Bowie/QUAL2E use **1.024** for settling, W2 applies no Arrhenius to `CBODS`. Also a `1/depth` form mismatch. | Bowie et al. 1985; QUAL2E; W2. | **v1 & v3 BOTH wrong.** Dormant at shipped `ksbod_20=0`; real silent error the instant a user calibrates nonzero CBOD settling (a documented legitimate usage). |

### 8.3 Upgraded by the external reference (parity pass under-rated these)

- **A2 → SCI-A2 (MAJOR, decided):** `f_pocp/f_pocb = 0.5` is not merely
  "undocumented" — it is scientifically **low**. W2 `APOM`≈0.6–0.9 (~0.8),
  QUAL2K, Bowie/Chapra: dead-algal C is predominantly particulate. **v1's 0.9
  is defensibly correct; v3's 0.5 is wrong.** Restore ≈0.8–0.9.
- **A11 → SCI-A3 (MAJOR, conditional):** `KL=10 W/m²` is a PAR-scale constant;
  `solar_radiation` reaches `limit_light` with no PAR scaling. If the registry
  feeds total shortwave, light limitation is overstated and growth
  over-predicted **~30–60 %** under light limitation. Shared v1+v3. Resolve by
  confirming the registry irradiance convention at the coupling driver.
- **DOX-F1 (MAJOR, conditional):** DO-saturation salinity correction silently
  omitted (correct for freshwater; ~18 % DOsat overstatement at 35 ppt). W2
  includes it with APHA-exact constants. Add a guard + freshwater docstring.
- **DOX-F2 (MAJOR):** silent zero reaeration under
  `hydraulic_reaeration_option=1` + `kah_20_user=0.0` is **V3-DIVERGENT** from
  W2, which enforces a `MINKL` minimum-reaeration floor on every branch. Adopt
  W2's precedent (warning + opt-in floor).

### 8.4 Parity-≠-correctness payoff: v3 fixed ~9 real v1 science defects

Independently verified against the external references — v3 is materially
*more* correct than v1:

- **Nitrogen θ transposition** (Phase 9.E): v1 `constants.py:134-137` has the
  within-pair θ values swapped (kon↔rnh4, kdnit↔vno3); v3 matches Fortran *and*
  the physical convention (OM-hydrolysis θ≈1.047, sediment-exchange ≈1.074–1.08,
  denitrification ≈1.045). **Confirmed real v1 bug, v3 correct.**
- **v1 `SOD_20=999`, `SOD_theta=999` sentinels** → first wet step yields
  `999·999^(T−20)` SOD, driving DOX catastrophically negative. v3 `1.0/1.060`
  (within Bowie). **Critical v1 defect, fixed in v3.**
- **v1 `fdp` degenerate stub** (returns 1.0 for all params, zeroing TIP
  sorbed-settling) vs v3's correct `1/(1+kdpo4·Solid·1e-6)` isotherm,
  structurally identical to W2 `PARTP`. **Genuine v1 defect, fixed in v3.**
- **v1 DIC mol-C/mg-C unit error** + omitted DOC→DIC and CBOD→DIC terms — v3
  is correct on all three (Phase 9.E `*12000`).
- **v1 `vb` POM settling unit bug (~1460×)** — corrected in v3.
- **v1 pathogen `apx`/`vx` = 1.0 placeholders** → v3 canonical Auer & Niehaus
  1993 / Chapra 1997 values.
- **The 16 parameter-default corrections**: all 16 independently checked move
  toward the literature/W2-consistent value; **none moves the wrong way**.

### 8.5 Cleared as defensible method differences (NOT errors)

θ^(T−20) Arrhenius vs W2's 4-point Thornton-Lessem rising/falling curve
(acceptable; shared with QUAL2E/QUAL2K; **one documented limitation:** monotone,
cannot reproduce W2 falling-limb thermal inhibition, so >28–30 °C algal growth
over-predicts vs W2); half-saturation light model (matches QUAL2E; W2 uses
Steele); multiplicative nutrient limitation (matches QUAL2K; W2 uses Liebig
minimum); exponential low-DO nitrification inhibition `1−e^(−KNR·DO)` (matches
QUAL2K, not W2's Monod); QUAL2E/EUTRO simple-ratio ammonia preference;
carbonate K1/K2 deferred to NSM2 (parity-preserving); `pCO2=383 ppm` lags
current ~420 (minor, parity-preserving).

---

## 9. Recommended action order

1. **Fix CA-1** (`alkalinity.py` `rca=AWc/AWa`, `rcb=BWc/BWd`; compose
   `AWa`/`BWd`; correct the root-cause docstring `:61-62` / comments
   `:202-206`) + a regression test that does **not** share the constant between
   sides + correct the simple-constituents audit doc §21–§23. **Only CRITICAL;
   externally confirmed; alkalinity-mass-conservation risk under blooms.**
2. **SCI-N1** — confirm what `denitrification_flux_rate` physically represents,
   then set `r_alkden = 1/14/1000` as a documented, W2/Stumm-Morgan-anchored
   divergence from upstream NSM1; report the Fortran `modAlkalinity.f90:54`
   error upstream; add a closed-system denitrification-alkalinity benchmark.
3. **SCI-A2** — restore `f_pocp/f_pocb ≈ 0.8–0.9` (or document a defensible
   deliberate deviation with a citation). **SCI-A3** — confirm the registry
   irradiance convention (PAR vs shortwave) at the coupling driver; this is the
   largest *potential* growth-bias and needs a cross-module answer.
4. **SCI-CB1** — set `ksbod_theta = 1.024` (or document the deviation) before
   any nonzero-CBOD-settling application. **DOX-F1 / DOX-F2** — salinity guard
   + reaeration floor/warning, following W2 precedent.
5. **SCI-A1** — schedule the N-flux reimplementation of alkalinity↔algae
   coupling (deeper than CA-1; aligns with W2). Reasonable to stage with NSM2.
6. **A1 / README** — rewrite stale forward-looking module headers and the
   package README as as-implemented **before any external (LimnoTech/sponsor)
   review packet ships**. Then the MINOR doc-fidelity cluster (§6).

**Shipped-default safety:** every MAJOR latent item (SCI-N1 magnitude in
alkalinity only, SCI-CB1, fdp, rpo4) is gated to zero or dominated by correct
terms at the validated Santiam–Salem defaults — which is why the benchmarks
pass. The risk surface is **productive/bloom, brackish, and user-calibrated**
regimes, plus alkalinity/pH whenever algae or denitrification are active.
