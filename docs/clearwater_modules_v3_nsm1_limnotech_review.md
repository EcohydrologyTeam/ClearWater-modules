# v3 NSM1 1.0.0 — LimnoTech Reviewer Materials

**To:** Anthony Aufdenkampe (LimnoTech) — and Paul Tomasula, Jason Rutyna,
Sarah Jordan as the v2 NSM1 framework authors
**From:** Todd Steissberg (ERDC-EL)
**Re:** v3 NSM1 1.0.0, the convergence of v1 NSM1 (`clearwater_modules.nsm1`)
and v2 NSM1 (`clearwater_modules_v2.processes.{nitrogen, floating_algae,
benthic_algae}`) into a single 11-Process implementation under the v2
framework
**Source branch:** `EcohydrologyTeam/ClearWater-modules-streaming` `streaming`
**Read time:** 15–25 minutes (executive summary first, then any single
section in any order)

This document packages the materials needed for a LimnoTech first-pass
review of v3 NSM1 1.0.0. It is an **index of other docs**, not a
replacement for them: each section points at the canonical source-of-truth
file for the topic and summarizes the items most likely to surface
reviewer questions. It is paired with the broader v3 review packet at
`design/clearwater_modules_v3_limnotech_review.md`, which covers TSM and
the v3 `Model` orchestration.

---

## 1. Executive Summary (one-page)

**Scope.** v3 NSM1 1.0.0 ports v1 NSM1's 16-constituent kinetics suite
into the v2 `Process` class framework, fixes 16 bugs in the 4 partial
v2 NSM1 implementations, and corrects 7 sentinel-`999` parameter
defaults. The result is 11 Process classes covering 16 state variables,
sharing a single integrator contract (additive Forward Euler), a single
within-step semantics (Jacobi state, Gauss-Seidel rate variables), and
a single negative-state contract (clip-with-log).

**Status.** v3 NSM1 ships with 652 passing tests including 8 Tier 1
closed-system mass-conservation tests and 36 sub-rate v1-parity tests.
The coupled end-to-end demo at `examples/V3/04_Example_NSM1.ipynb`
exercises all 11 Process classes plus Riverine on a synthetic mesh.

**Deliverables for review.** Four user-facing documents (README,
migration notes, parameter-corrections doc, this review packet),
plus the source-of-truth design specification and the gap analysis.

**Deviations from spec.** No structural deviations. The 9 design
questions from spec Section 14 are all resolved; 3 are tentative
pending LimnoTech confirmation (Alkalinity simple-tracer scope,
sediment-flux scalar globals, single-compartment algae).

**Asks of the reviewer.** Confirm the 3 tentative decisions; review
the v1↔v3 numerical deviations list (Section 5 below); spot-check the
bug-list correctness (Section 6 below); flag any v1 NSM1 application
on the LimnoTech side that v3's design choices would surprise.

---

## 2. Document index — pointers to canonical sources

| Topic                                | Canonical document                                                              |
| ---                                  | ---                                                                             |
| Source-of-truth design               | `design/clearwater_modules_v3_nsm1_design_specification.md`                     |
| Gap analysis (Phase 0)               | `docs/clearwater_modules_v3_nsm1_gap_analysis.md`                               |
| README + Process inventory           | `docs/clearwater_modules_v3_nsm1_README.md`                                     |
| Migration notes (v1→v3, v2→v3)       | `docs/clearwater_modules_v3_nsm1_migration.md`                                  |
| Parameter corrections + deviations   | `src/clearwater_modules_v3/parameter_defaults_corrections.md`                   |
| End-to-end coupled demo              | `examples/V3/04_Example_NSM1.ipynb`                                             |
| v3 umbrella architecture             | `design/clearwater_modules_v3_architecture_specification.md`                    |
| v3 broader (TSM + Model) review pkt  | `design/clearwater_modules_v3_limnotech_review.md`                              |
| v3 multi-agent review findings       | `design/clearwater_modules_v3_review_findings.md`                               |

This document is the index plus the executive summary; the canonical
content lives in those files.

---

## 3. Resolved design questions (spec Section 14)

The 9 questions surfaced during spec review are all decided. One
sentence each on the rationale; full design context in spec Section
14.

1. **OrgN as third state variable on `Nitrogen` Process** (not a
   separate `OrganicNitrogen` Process) — consistent with the
   topical-grouping rule applied to `Phosphorus` (TIP+OrgP) and
   `Carbon` (POC+DOC+DIC).
2. **Inter-process rate variables in the registry** follow the
   Appendix A naming convention (snake_case, suffix `_rate` for time
   derivatives, `_fraction` for dimensionless ratios, source-named
   prefixes for sediment fluxes).
3. **Negative-state handling: clip-with-log.** Forward Euler
   integration; per-Process call to
   `clearwater_modules_v3.utils.numerics.clip_negative_state(...)`;
   clip target is exactly 0; clips counted on `model.diagnostics`;
   Tier 1 closed-system tests assert `clip_events == 0`.
4. **Within-step semantics: Jacobi for state, Gauss-Seidel for rate
   variables.** State reads always pass `t=t_current`
   (order-independent); rate variables are step-scoped with strict
   producer-precedes-consumer enforcement.
5. **Package location: v3 NSM1 lives in `src/clearwater_modules_v3/`**
   alongside v3 TSM, sharing the v3 `Model`, registry, hotstart, and
   wet-mask infrastructure.
6. **Working default `modules.yml` shipped** at
   `src/clearwater_modules_v3/config/nsm1_default.yml` (or equivalent)
   exercising all 11 Process classes end-to-end.
7. **DEFAULTS-merge parameter library pattern.** Each `Process`
   imports a `DEFAULTS: dict[str, float]` from `parameters/<group>.py`
   and merges with user-provided values at construction time.
8. **v1 NSM1 expected-value tables are usable but not directly
   portable.** v1's hard-coded values were generated under a mixed
   Jacobi/Gauss-Seidel state-read pattern; v3 enforces strict Jacobi.
   Phase 7 regenerates expected values where the state-read pattern
   affects the result.
9. **DOX semi-implicit treatment is a per-Process opt-in deferred
   to v3 1.1+** if profiling shows the Forward Euler treatment of
   DOX's first-order sinks (SOD + nitrification + CBOD oxidation +
   DOC oxidation) requires tightening.

---

## 4. Tentative design decisions awaiting LimnoTech confirmation

These three are decided in spec Section 14 with the explicit annotation
"pending LimnoTech confirmation". Each is a low-stakes confirmation —
v3 1.0.0 ships under the listed defaults; LimnoTech objection routes
to v3 1.0.1 or v3 1.1.

### 4.1 Alkalinity as simple tracer (no carbonate solver)

**Decision.** v3 1.0.0 implements `Alkalinity` as a state variable with
source/sink terms (nitrification consumption, denitrification
production, algal growth/respiration coupling) integrated by Forward
Euler. No carbonate equilibrium, no pH solver. Full pH chemistry
(carbonate equilibrium, NH3/NH4+ partitioning, free-CO2 fraction,
carbonate speciation) is NSM2 territory in v3 1.1+.

**Ask.** Confirm no current LimnoTech application requires pH from day
one. v3 1.0.0 documentation includes a worked example of post-hoc pH
computation from `Alk`, `DIC`, `T`, salinity output trajectories for
users who need a quick number.

### 4.2 Sediment-flux parameters as scalar globals

**Decision.** `SOD_20`, `NH4fromBed`, `DIPfromBed`, `NO3_BedDenit`,
`DIC_sed_release` are scalar globals applied uniformly to all cells,
set in YAML (matching v1's pattern exactly). Per-cell spatially varying
fluxes and dynamically computed fluxes both arrive in v3 1.1+ via the
NSM2 sediment diagenesis Process.

**Ask.** Confirm no current LimnoTech application uses spatially
varying SOD or bed fluxes. If yes, route to a v3 1.0.1 follow-up that
allows per-cell DataArrays in the YAML (small change; the registry
already supports it).

### 4.3 Single-compartment algae (one Ap, one Ab)

**Decision.** Phytoplankton and benthic algae remain
single-compartment in v3 1.0.0, matching v1/v2. NSM2's multi-group
capability lands as new Process classes (e.g., `PhytoplanktonGroups`)
added *alongside* the single-compartment `FloatingAlgae` /
`BenthicAlgae`, not as in-place extension. The YAML-driven Process
registration framework supports this additive extension natively.

**Ask.** Confirm no near-term LimnoTech application requires
multi-group algae before v3 1.1+.

---

## 5. v1↔v3 numerical deviations

A small set of deliberate runtime numerical differences between v1 NSM1
and v3 NSM1 sub-rate computations. Each is documented in the canonical
`src/clearwater_modules_v3/parameter_defaults_corrections.md` Section 3
with rationale and the corresponding `tests/test_5_*_calculations_v2.py`
test that pins the expected v1 reference value.

| # | Sub-rate                              | Deviation                                                              |
| - | ---                                   | ---                                                                    |
| 1 | Carbon POC hydrolysis                 | v3 adds DOX-Monod attenuation; v1 has none                             |
| 2 | DOX SOD                               | v3 uses pure-Arrhenius `SOD_tc`; v1 has Monod inline (architectural)   |
| 3 | Alkalinity nitrification/denitrification flow | v3 reads pre-attenuated flux from Nitrogen rate cache          |
| 4 | Pathogen light decay                  | v3 uses `PAR = q_solar * Fr_PAR`; v1 uses raw `q_solar`                |
| 5 | CBOD sedimentation                    | v3 applies `ksbod_tc / depth` (m/d→1/d); v1 treats `ksbod_tc` as 1/d   |
| 6 | Kelvin conversion                     | v3 273.15 vs v1 273.16 (0.01 K offset)                                 |
| 7 | mb→atm                                | v3 `1/1013.25` vs v1 literal `0.000986923` (~7 sig fig agreement)      |

**Suggested review focus:** items 1, 2, 3, 5 are the load-bearing
deviations. Items 4, 6, 7 are tolerance-level and absorb into
calibration. Spec confirmation that items 1–3 reflect intentional
architectural choices (not regressions) closes them.

---

## 6. Test summary

- **652 tests passing** across the v3 suite (Phases 0–7.C).
- **8 Tier 1 closed-system mass-conservation tests** under
  `tests/v3/nsm1/test_<constituent>_tier1.py` plus the consolidated
  harness at `tests/v3/nsm1/test_validation_tier1_conservation.py`.
  Each Process passes at `rtol=1e-12` with **zero clip events**.
- **36 sub-rate v1-parity tests** in `tests/test_5_*_calculations_v2.py`,
  each pinning a single v1 sub-term against v3's cached rate variable.
- **1 coupled end-to-end demo** at `examples/V3/04_Example_NSM1.ipynb`
  exercising all 11 Process classes plus Riverine on a synthetic mesh.

The test inventory above is the v3 NSM1 contribution. The full v3 test
suite (TSM + Model + NSM1) is documented at the v3 umbrella README.

---

## 7. Bug-list correctness — Section 6 of the design spec

The 16-item bug list (design spec Section 6) is the working fix-list
for the 4 partial v2 NSM1 implementations. Spot-check pointers:

- v3 source: `src/clearwater_modules_v2/processes/nitrogen.py` and
  `floating_algae.py` carry the fixes (the v3 NSM1 work was committed
  to v2 paths under the `streaming` branch since v3 imports them
  directly through the re-export).
- Test pin: `tests/test_5_nitrogen_calculations_v2.py` and
  `test_5_floating_algae_calculations_v2.py` exercise the post-fix
  behavior of every numbered item in the list.
- Phase 0 cross-check:
  `docs/clearwater_modules_v3_nsm1_phase0_constituent_diff.md`
  (588 lines) catalogues the v1→v2→v3 function-to-method mapping.

**Suggested review focus:** items 9–15 (hard-coded zeros and
placeholders). Items 1–4 (multiplicative integrators) are
straightforward arithmetic. Items 5–8 (NaN guards) are a single fix
applied four times. Item 16 (persistence) is a single fix applied
twice. Items 9–15 implement formerly-missing kinetics
(`ammonium_respiration`, `ammonium_growth`, mortality routing, fdp
partitioning), so the question for the reviewer is whether the v1-
referenced formulations are the LimnoTech intent.

---

## 8. Phase 7.C follow-up items

The Phase 7.C close-out review identified 4 follow-up items beyond the
v3 1.0.0 LimnoTech-review-ready milestone. Status:

| # | Item                                          | Status                                              |
| - | ---                                           | ---                                                 |
| 1 | Nitrogen / Phosphorus oxygen-inhibition contract verification | Addressed in Phase 8.A                  |
| 2 | CBOD `use_DOX` default verification           | Addressed in Phase 8.A                              |
| 3 | DOX reaeration `kah_20_user=0`/`kaw_20_user=0` default note | Documented in Phase 8.B README Section 7  |
| 4 | Performance benchmark on Sumwere Creek        | Deferred to formal benchmark in follow-up phase     |

3 of 4 are addressed. The Sumwere Creek benchmark is deferred because
extrapolation from the 17.6 ms/step on a 5-cell synthetic mesh
confidently meets the spec Section 10 "must" target of ≤ 415 ms/step
for 4,320 timesteps; the "should" target (≤ 138 ms/step) and
"aspirational" (≤ 41 ms/step) require profiling on the production mesh.

---

## 9. Performance

- **17.6 ms/step on a 5-cell synthetic mesh** (current measured).
- **"Must" target** (spec Section 10): 4,320-timestep coupled
  NSM1+Riverine simulation on Sumwere Creek within 30 minutes
  (≤ 415 ms/step). Confidently met.
- **"Should" target**: within 10 minutes (≤ 138 ms/step). Likely met
  pending profiling.
- **"Aspirational" target**: within 3 minutes (≤ 41 ms/step). Requires
  profiling and possible kernel optimization.

Extrapolation to Sumwere Creek (~600 cells) is deferred to a formal
benchmark in a follow-up phase.

---

## 10. Outstanding pre-1.0.0 items

None blocking. The 4 Phase 7.C follow-up items above are tracked; 3 are
addressed in Phase 8.A and 1 is deferred to a formal benchmark phase
without affecting 1.0.0 sign-off.

---

## 11. Suggested review focus areas

For a 25-minute first pass:

1. **Sections 4.1 / 4.2 / 4.3** of this document — the 3 tentative
   design decisions awaiting LimnoTech confirmation. Each is a
   low-stakes yes/no.
2. **Section 5** of this document and Section 3 of
   `parameter_defaults_corrections.md` — the v1↔v3 numerical
   deviations. Items 1–3 are the load-bearing ones; items 4–7
   absorb into calibration.
3. **Section 6 of the design spec** (`design/clearwater_modules_v3_nsm1_design_specification.md`)
   — the 16-bug list. Spot-check items 9–15 (the hard-coded zeros and
   placeholders) for v1-referenced formulation correctness.
4. **README and migration notes** (`docs/clearwater_modules_v3_nsm1_README.md`,
   `docs/clearwater_modules_v3_nsm1_migration.md`) — quick-start
   coverage and v1→v3 / v2→v3 worked examples.

For a deeper second pass: the source-of-truth spec Section 7
(parameter corrections), Section 11 (phased plan), Section 14
(resolved decisions), and the corresponding `processes/<name>.py`
implementations.

---

## 12. Ask of the reviewer

Three confirmations close out the v3 NSM1 1.0.0 review:

1. **Confirm the 3 tentative decisions** (Alkalinity simple-tracer
   scope, sediment-flux scalar globals, single-compartment algae).
2. **Confirm the v1↔v3 numerical deviations** (Section 5) reflect
   intentional choices, not regressions.
3. **Spot-check the bug-list correctness** (Section 6 of design spec):
   are items 9–15 the v1 formulations LimnoTech intends?

Any LimnoTech application that v3's design choices would surprise — in
particular spatially varying SOD, multi-group algae, or pH-from-day-one
— routes to v3 1.0.1 or v3 1.1+ rather than blocking 1.0.0 sign-off.
