# v3 NSM1 Review — Findings, Determinations, and Recommendations

- **Date:** 2026-05-15
- **Subject:** Code-correctness, scientific-correctness, and documentation review of v3 NSM1
- **Branch / commit reviewed:** `streaming` @ `54f2b12` (review artifacts committed at `b9230fe`)
- **Status of code:** **Review only. No source code was modified.** All items below are open unless explicitly marked resolved.
- **Detailed backing reports:** `design/v3_nsm1_review_2026-05-15/` — `review_SUMMARY.md` (consolidated), `review_*.md` (6 parity reports), `scival_*.md` (6 science-validation reports).

---

## 1. Purpose and scope

A two-pass review of the full v3 NSM1 implementation (~9,300 LOC of process
kinetics plus parameters, utilities, model framework, and configuration):

1. **Pass 1 — parity, refactor, documentation.** Line-level read of every
   in-scope v3 file against the v1 reference (`clearwater_modules/nsm1`) and the
   design/audit-doc trail. Questions: is the xarray refactor complete; do the
   algorithms match v1 (or improve on it intentionally); is the code correctly
   documented; are fixed TODO/bug comments updated; are any issues unaddressed.

2. **Pass 2 — scientific correctness.** Because v1-parity is not proof of
   correctness if v1 itself is wrong, the kinetics were re-validated against
   **CE-QUAL-W2 v2026.02** (read line-by-line from local ERDC source —
   authoritative and checkable), **QUAL2K** (Chapra, Pelletier & Tao 2008),
   **QUAL2E** (Brown & Barnwell 1987), **WASP/EUTRO**, and the rate literature
   (**Bowie et al. 1985**, **Chapra 1997**, **Stumm & Morgan**). v1 was treated
   as a suspect, not an oracle.

Twelve specialist reviewer agents (six per pass), partitioned by scientific
domain, with lead synthesis and independent verification of the CRITICAL
finding.

---

## 2. Headline determinations

1. **The xarray refactor is complete.** No per-cell Python loops, no `== np.nan`
   antipatterns (NaN/inf logic uses `.isnull()` / `xr.where` / `np.isnan`), no
   array-truthiness branching on per-cell data, container-aware guards
   throughout. Residual items are latent broadcasting hazards (scalar
   early-returns that drop DataArray coordinates) not triggered in current runs.

2. **The v1→v3 port is faithful**, and the historically buggy modules are
   verifiably fixed: FloatingAlgae Bug #4/#13/#14/#15/#16, Nitrogen Phase
   9.A/9.B/9.E, Carbon/DOX Phase 9.B. Every traced kinetic term is `MATCH` to v1
   or a *documented* intentional improvement — with one exception (CA-1).

3. **v3 is materially *more* scientifically correct than v1.** The
   science-validation pass independently confirmed that v3 fixes ~9 genuine v1
   science defects (Section 5). Blind v1-parity would have propagated all of
   them.

4. **One CRITICAL defect (CA-1) is open and externally confirmed**, plus one
   new MAJOR defect shared by v1 *and* v3 that the parity pass structurally
   could not catch (SCI-N1), plus a set of MAJOR latent/calibration-regime
   issues and a documentation-staleness cluster.

5. **Shipped-default safety.** Every MAJOR latent item is gated to zero or
   dominated by correct terms at the validated Santiam–Salem default
   configuration — which is why the v1-vs-v3 benchmark and the case study
   passed. The risk surface is **productive/bloom, brackish, and
   user-calibrated** regimes, plus alkalinity/pH whenever algae or
   denitrification are active.

### Severity rollup (parity pass)

| Domain | CRITICAL | MAJOR | MINOR | OBS |
|---|---:|---:|---:|---:|
| Nitrogen + N2 | 0 | 0 | 7 | 8 |
| Algae | 0 | 2 | 7 | 6 |
| Carbon + Alkalinity | **1** | 0 | 4 | 4 |
| DOX | 0 | 1 | 4 | 4 |
| CBOD + POM | 0 | 1 | 2 | 3 |
| Phosphorus + Pathogen | 0 | 0 | 5 | 4 |
| Framework / utils | 0 | 0 | 5 | 9 |
| **Total** | **1** | **4** | **34** | **38** |

The science-validation pass then confirmed CA-1 against external ground truth,
added SCI-N1, and upgraded several parity-pass MINOR/observation items to MAJOR.

---

## 3. CRITICAL finding

### CA-1 — `alkalinity.py`: algal/benthic alkalinity coupling uses raw stoichiometric weights (1000× / 100× flux error)

- **Location:** `src/clearwater_modules_v3/processes/alkalinity.py:362, 386, 411,
  441`; composition block `:202-206`; root-cause docstring `:61-62`.
- **Defect:** v3 binds `rca = self.AWc` (= 40) and `rcb = self.BWc` (= 40) — the
  *raw* stoichiometric weights — where the scientifically correct quantities are
  the *resolved intensive ratios* `rca = AWc/AWa = 0.04` mg-C/µg-Chla and
  `rcb = BWc/BWd = 0.4` mg-C/mg-D. The alkalinity stoichiometric constants
  `r_alkaa`/`r_alkba` (eq/mg-C) are byte-identical to v1, so the ratio factor is
  the sole divergence.
- **Magnitude:** floating-algae alkalinity flux inflated **1000×**
  (40 / 0.04); benthic-algae flux inflated **100×** (40 / 0.4). Correct
  floating-algae flux ≈ 9.172e-4 mg-CaCO₃/L/d; v3 produces 9.172e-1.
- **External confirmation:** CE-QUAL-W2 `ENTRY ALKALINITY`
  (`water-quality.f90:3164-3166`) routes algal alkalinity coupling through the
  *intensive* carbon fraction, never a raw weight; Stumm & Morgan Table 4.5
  gives the same stoichiometry W2 quotes verbatim. **Here v1 is correct and v3
  is wrong.**
- **Why it escaped:** the Phase 9.B fix that introduced `rca = self.AWc/self.AWa`
  reached `carbon.py:495-496` and `dox.py` but **not** `alkalinity.py`; the
  `alkalinity.py` docstring and comments rationalize the wrong convention as
  intentional; the parity test feeds the same wrong constant to both sides; the
  simple-constituents audit doc (§21–§23) marked it "Match."
- **Impact:** alkalinity (and any downstream pH) is correct in Santiam–Salem
  because that signal is dominated by the nitrification/denitrification terms
  (which use the correct `r_alkn`/`r_alkden` stoichiometry). Under productive /
  bloom conditions the algal coupling term dominates and over-drives alkalinity
  by 2–3 orders of magnitude, including an alkalinity mass-conservation
  violation (silent negative-clip events). **Latent, not benign.**
- **Recommended fix:** mirror the `carbon.py:495-496` pattern — derive
  `rca = self.AWc / self.AWa`, `rcb = self.BWc / self.BWd`; add `AWa`
  (`ALGAE_DEFAULTS`) and `BWd` (`BALGAE_DEFAULTS`) to the `alkalinity.py`
  composition block; correct the root-cause docstring `:61-62` and comments
  `:202-206`; add a regression test that does **not** share the constant
  between the v1-mirror and v3 paths; correct the simple-constituents audit doc.

---

## 4. MAJOR findings

### SCI-N1 — denitrification alkalinity coefficient 4× too high (v1 *and* v3)

`r_alkden = 4/14/1000` eq/mg-N overstates denitrification alkalinity production
by **4×**; the correct value is `1/14/1000` (1 eq alkalinity per mol N
denitrified — CE-QUAL-W2 `water-quality.f90:3157`; Stumm & Morgan). This is an
**upstream NSM1 Fortran error** (`modAlkalinity.f90:54`) propagated unchanged to
v1 and v3 — the canonical case the parity pass *structurally could not catch*
because v1 == v3. It does not affect nitrogen mass balance (it lives in the
Alkalinity process) but biases simulated alkalinity and derived pH wherever
denitrification is significant. **Determination:** confirmed wrong on
W2 line-level evidence; before patching, confirm what NSM1
`denitrification_flux_rate` physically represents, then set `r_alkden =
1/14/1000` as a documented, reference-anchored divergence from upstream NSM1,
and report the Fortran error upstream.

### SCI-A1 — alkalinity↔algae coupling is carbon-routed; the process is nitrogen-driven

Even after CA-1's magnitude is fixed, NSM1's alkalinity-algae term is a
carbon-routed approximation of an intrinsically *nitrogen*-driven process
(CE-QUAL-W2 and Stumm & Morgan compute it from the N flux). It is self-consistent
only at exact Redfield N:C; NSM1's own `AWn/AWc` is ~2% off Redfield, so any
non-Redfield calibration silently desyncs alkalinity from the N balance.
**Determination:** a true formulation limitation, additional to CA-1. The deeper
fix is to reimplement on the N-flux basis as W2 does; reasonable to stage with
NSM2.

### SCI-A2 — `f_pocp` / `f_pocb` = 0.5 is scientifically low

The inline default 0.5 (fraction of algal/benthic mortality carbon routed to
POC) is not merely undocumented — it is low against every validated reference:
CE-QUAL-W2 `APOM` ≈ 0.6–0.9 (typically ~0.8), QUAL2K, Bowie et al., and Chapra
all treat dead-algal carbon as predominantly particulate. **v1's 0.9 is
defensibly correct; v3's 0.5 is wrong.** It biases ~40% of mortality carbon into
DOC, perturbing the DOC→DIC pathway and DO demand (not algal biomass directly).
**Determination:** restore ≈0.8–0.9, or document a defensible deliberate
deviation with a citation.

### SCI-A3 — light limitation: PAR vs total-shortwave irradiance basis

`limit_light` compares the light half-saturation constant `KL` (default
10 W/m², a PAR-scale value) directly against `registry.get_at_time(
"solar_radiation")` inside the standard depth-averaged Beer–Lambert
half-saturation form. If the coupling driver supplies **total shortwave** while
`KL` is on a PAR basis (PAR ≈ 0.43–0.47 of shortwave), light limitation is
understated and algal growth is over-predicted by **~30–60% wherever light
limitation is active**. Shared by v1 and v3 (so it does not appear in the
v1↔v3 benchmark — only against observed data). **Determination:** highest
*potential* growth bias in the model; **conditional** on the registry
irradiance convention at the coupling driver, which must be confirmed
(decisive check — see Section 6).

### SCI-CB1 — `ksbod_theta = 1.047` misapplies the oxidation θ to CBOD settling (v1 *and* v3)

CBOD settling uses θ = 1.047 (the *oxidation* Arrhenius coefficient); Bowie et
al. and QUAL2E specify θ = 1.024 for settling, and CE-QUAL-W2 applies no
Arrhenius to `CBODS`. There is also a `1/depth` form mismatch versus
QUAL2E/Fortran. Shared by v1 and v3. **Dormant at the shipped default
`ksbod_20 = 0`**; a real, silent science error the instant a user calibrates
nonzero CBOD settling (a usage the corrections doc itself documents as
legitimate). **Determination:** set `ksbod_theta = 1.024` (or document a
deliberate deviation) before any nonzero-CBOD-settling application.

### DOX-F1 / DOX-F2 — DO-saturation salinity omission; silent zero reaeration

- **DOX-F1:** the DO-saturation salinity correction is silently omitted —
  correct for fresh water (factor = 1.0) but with no salinity input, no guard,
  and no docstring statement of the freshwater assumption. CE-QUAL-W2 includes
  the term with APHA-exact constants; ~18% DOsat overstatement at 35 ppt.
  Relevant to any brackish/estuarine coupling.
- **DOX-F2:** with `hydraulic_reaeration_option = 1` and `kah_20_user = 0.0`, a
  user gets **silent zero atmospheric reaeration**. This is **V3-DIVERGENT from
  CE-QUAL-W2**, which enforces a `MINKL` minimum-reaeration floor on every
  branch (`gas-transfer.f90`). Not triggered by the shipped default (option 5);
  conditional on a user override to the v1/Fortran convention.

**Determination:** add a freshwater-assumption docstring + salinity guard
(DOX-F1); follow W2's precedent with a one-time warning plus an opt-in
minimum-reaeration floor (DOX-F2).

### Documentation-staleness cluster (the explicit review ask)

- **Algae A1:** `floating_algae.py:1-46` and `benthic_algae.py:1-20` module
  headers are written as forward-looking work plans ("apply the fixes",
  "Bug #13: implement (was returning 0)"). The work is done and regression-tested
  — a reader concludes the module is still broken. This is precisely the failure
  mode the review was asked to catch.
- **README:** `src/clearwater_modules_v3/README.md` predates NSM1 pattern
  alignment Phases 6–10 and is internally contradictory — it states Riverine /
  BenthicAlgae / FloatingAlgae / Nitrogen and Process / ProcessFactory are v2
  re-exports (they are all v3-native in-tree; the top-level `__init__.py` says
  v2 was removed), and claims "2 of 18 MAJOR resolved, Phase 5 in progress"
  (the authoritative triage records 17 of 18).
- **Minor doc-fidelity:** stale references to deleted shadow methods/test files
  in `carbon.py`; an N2 docstring naming a non-existent attribute
  `denitrification_rate`; `nitrate_uptake_*` docstrings that lead with the old
  broken implementation; a stale `utils/reaeration.py:168-169` marker; stale
  audit-doc entries (`fdp` pre-fix form; simple-constituents §21–§23 "Match" for
  CA-1). Full inventory in `review_SUMMARY.md` §6.

**Determination:** rewrite the stale forward-looking module headers and the
package README as as-implemented records **before any external
(LimnoTech/sponsor) review packet ships**; then clear the minor doc-fidelity
cluster. None of these is a logic defect; all would mislead a reviewer.

---

## 5. Parity-≠-correctness payoff — v1 defects that v3 fixes

Independently verified against the external references:

- **Nitrogen θ transposition (Phase 9.E):** v1 `constants.py:134-137` has the
  within-pair θ values swapped (kon↔rnh4, kdnit↔vno3). v3 matches the Fortran
  source *and* the physical convention (OM-hydrolysis θ ≈ 1.047,
  sediment-exchange ≈ 1.074–1.08, denitrification ≈ 1.045). Confirmed real v1
  bug; v3 correct.
- **v1 `SOD_20 = 999`, `SOD_theta = 999` sentinels** → first wet step yields
  `999 · 999^(T−20)` SOD, driving DOX catastrophically negative. v3's
  `1.0 / 1.060` is within Bowie ranges. Critical v1 defect, fixed in v3.
- **v1 `fdp` degenerate stub** (returns 1.0 for all parameters, zeroing TIP
  sorbed settling) vs v3's correct `1/(1 + kdpo4·Solid·1e-6)` isotherm,
  structurally identical to CE-QUAL-W2 `PARTP`. Genuine v1 defect, fixed in v3.
- **v1 DIC mol-C/mg-C unit error** plus omitted DOC→DIC and CBOD→DIC terms — v3
  correct on all three (Phase 9.E `*12000`).
- **v1 `vb` POM settling unit bug (~1460×)** — corrected in v3.
- **v1 pathogen `apx`/`vx` = 1.0 placeholders** → v3 canonical Auer & Niehaus
  1993 / Chapra 1997 values.
- **The 16 parameter-default corrections** — all 16 independently checked move
  toward the literature/W2-consistent value; none moves the wrong way.

### Cleared as defensible method differences (NOT errors)

θ^(T−20) Arrhenius vs CE-QUAL-W2's 4-point Thornton–Lessem rising/falling curve
(acceptable; shared with QUAL2E/QUAL2K; one documented limitation: monotone,
cannot reproduce W2 falling-limb thermal inhibition, so >28–30 °C algal growth
over-predicts vs W2); half-saturation light model (matches QUAL2E; W2 uses
Steele); multiplicative nutrient limitation (matches QUAL2K; W2 uses Liebig
minimum); exponential low-DO nitrification inhibition (matches QUAL2K, not W2's
Monod); QUAL2E/EUTRO simple-ratio ammonia preference; carbonate K1/K2 deferred
to NSM2 (parity-preserving); `pCO2 = 383 ppm` lags current ~420 (minor,
parity-preserving).

### Correctly deferred to NSM2 (NOT defects)

`Nitrogen use_SedFlux=True` → `NotImplementedError`; full carbonate/pH solver;
SOD-derived DIC sediment-release fallback; alkalinity DOX-Monod single-source
routing; POM→DOC source; multi-group algae. All explicitly guarded and
documented.

---

## 6. Determination on the Santiam–Salem algae mismatch

Question posed: could the issues found explain algae not matching observed data
well in the Willamette River Santiam–Salem case study?

**Rule out the CRITICAL and the alkalinity/nitrogen-alkalinity findings.** CA-1,
SCI-N1, and SCI-A1 all live in the **Alkalinity** process. In NSM1, algal growth
depends on N, P, light, and temperature — **not** alkalinity — and v3 defers the
carbonate/pH solver to NSM2, so there is no pH→algae feedback path. **CA-1 is
not the cause of the algae mismatch.** This is stated explicitly so the fix
effort is not misdirected: CA-1 corrupts simulated alkalinity/pH, not
chlorophyll.

**Leading hypothesis — SCI-A3 (light: PAR vs total shortwave).** The code path
compares `KL` (PAR-scale, default 10 W/m²) directly against `solar_radiation`
inside the depth-averaged Beer–Lambert half-saturation light limitation. If the
case study feeds **total shortwave** while `KL` is PAR-based, the
light-limitation factor is biased high, light *under*-limits, and algal growth
is **over-predicted by ~30–60% wherever light limitation is active** — which, on
the Willamette mainstem with depth, turbidity, and self-shading, is most of the
reach most of the time. The bias is systematic and sign-consistent, and because
it is **shared by v1 and v3** it would not surface in the v1↔v3 benchmark (which
looked good) — only against observed data. This matches the reported symptom
(other constituents acceptable, algae poor) closely. **This is the most
probable contributor**, pending the irradiance-convention trace.

**Secondary contributors (smaller, possibly offsetting):**

- **Nutrient-limitation form.** If the run used `growth_rate_option = 1`
  (multiplicative `FN·FP`), that is *more* restrictive than the Liebig minimum
  CE-QUAL-W2 uses; with both N and P moderately limiting it can *under*-predict
  growth, partially opposing the light bias. Net effect is reach- and
  season-specific.
- **θ^(T−20) monotone (no falling limb)** → over-prediction at sustained warm
  summer temperatures (>~28–30 °C); plausible as a late-summer secondary effect
  on the lower Willamette.
- **SCI-A2 (`f_pocp/f_pocb = 0.5`)** does *not* change algal biomass directly
  (it mis-routes dead-algal carbon between POC and DOC); at most a weak indirect
  feedback via DO/DOC. Not a primary algae driver.
- **Ordinary calibration** of `mu_max_20` / `krp_20` / `kdp_20` — all within
  literature ranges, but reach-specific tuning still matters.

**Decisive check (not yet performed).** Confirm the irradiance convention at the
coupling driver: which variable/units the Santiam–Salem run writes into
`solar_radiation`, and whether `KL`, `light_limitation_option`, and
`growth_rate_option` were set on that same basis. The **direction** of the
observed miss further discriminates: over-prediction strongly implicates the
light/PAR path; under-prediction points to the multiplicative-nutrient form.

**Determination:** the algae mismatch is **plausibly explained by SCI-A3**, not
by the CRITICAL alkalinity defect. Confidence is *moderate and conditional*
pending the irradiance-convention trace and the miss-direction.

---

## 7. Recommendations (action order)

1. **Fix CA-1.** Only CRITICAL; externally confirmed; alkalinity
   mass-conservation risk under blooms. Code + non-shared-constant regression
   test + docstring/comment correction + audit-doc correction.
2. **SCI-N1.** Confirm the `denitrification_flux_rate` definition, set
   `r_alkden = 1/14/1000` as a documented W2/Stumm-Morgan-anchored divergence,
   report the upstream Fortran error, add a closed-system benchmark.
3. **Resolve the Santiam–Salem algae question.** Trace the case-study
   `solar_radiation` input and the `KL` / `light_limitation_option` /
   `growth_rate_option` settings to confirm or refute SCI-A3. If confirmed,
   apply the PAR scaling (or correct `KL`'s basis) and re-run the case study.
4. **SCI-A2.** Restore `f_pocp/f_pocb ≈ 0.8–0.9` (or document a cited deliberate
   deviation).
5. **SCI-CB1 / DOX-F1 / DOX-F2.** `ksbod_theta = 1.024`; DO-sat freshwater
   docstring + salinity guard; reaeration warning + opt-in floor (W2 precedent).
   Required before nonzero-CBOD-settling or brackish applications.
6. **SCI-A1.** Schedule the N-flux reimplementation of alkalinity↔algae coupling
   (deeper than CA-1; aligns with CE-QUAL-W2). Reasonable to stage with NSM2.
7. **Documentation.** Rewrite the stale forward-looking module headers and the
   package README as as-implemented records **before any external review packet
   ships**; then clear the minor doc-fidelity cluster (`review_SUMMARY.md` §6).

---

## 8. Confidence and caveats

- High confidence on all CE-QUAL-W2-anchored conclusions: the W2 v2026.02 source
  was read line-by-line locally; CA-1's 1000×/100× magnitude is pure arithmetic
  and needs no model run.
- Medium confidence on QUAL2K/QUAL2E/WASP/Bowie/Chapra numeric *range
  endpoints*: cited from domain knowledge, not local copies of the source
  documents. The directional verdicts are robust to this because the canonical
  conventions (θ ≈ 1.047 OM hydrolysis, ≈ 1.024 reaeration/settling, ≈ 1.06–1.08
  SOD/nitrification) are stable across editions, and every severity rests on the
  locally verified CE-QUAL-W2 source plus arithmetic, never solely on a
  knowledge-cited reference.
- The Santiam–Salem determination (Section 6) is conditional pending the
  irradiance-convention trace and the miss-direction; no benchmarks were re-run
  as part of this review.
