# ClearWater Modules v3 NSM1 — Phase 0 Gap Analysis

**Status:** Phase 0 deliverable per design spec Section 11
**Date:** 2026-05-04
**Author:** Todd Steissberg (ERDC), with Claude
**Scope:** Synthesis of three Phase 0 sub-audits (constituent diff, parameter audit, test fixture audit) into an actionable readiness assessment for Phases 1–7.

The detailed inventories live in three companion documents:

- `clearwater_modules_v3_nsm1_phase0_constituent_diff.md` (588 lines) — per-constituent v1→v2→v3 function-to-method mapping
- `clearwater_modules_v3_nsm1_phase0_parameter_audit.md` (403 lines) — full inventory of ~250 v1 NSM1 parameter defaults across 13 TypedDict groups
- `clearwater_modules_v3_nsm1_phase0_test_audit.md` (274 lines) — v1 NSM1 test files (358 test functions, 485 hard-coded expected values), v2 parity tests, and `tests/NSM Manual Calcs/` workbooks

This synthesis is the executive-level Phase 0 deliverable required by the design spec; readers needing line-by-line detail should consult the companion documents.

---

## 1. Executive Summary

Phase 0 confirms that the design spec's component inventory (Section 4), bug-fix list (Section 6), and sentinel-999 list (Section 7) are **substantially complete and accurate**. Three independent research agents each working from the v1 source, the v2 source, and the test fixtures produced findings that align closely with the spec.

Three findings warrant spec amendments before Phase 1 begins:

1. **`pressure_mb=2026.5` is the 7th sentinel-class default.** Approximately 2× sea-level atmospheric pressure (~1013 hPa). Load-bearing for `O2sat` and `N2sat` calculations via Henry's law. **Add to Section 7.**
2. **The Section 6 bug count needs reconciliation.** Phase 0.1 catalogued 10 critical v2 bugs; the spec lists 16 numbered items. The discrepancy is real: 10 are correctness bugs (integrator, NaN guards, persistence, hard-coded zeros, missing implementations) and 6 are derivative/related items (typos, parameter wiring TODOs). The spec's 16 is the right granularity for a fix-list; Phase 0.1's 10 is the right granularity for a "what's actually broken" summary. **No spec change needed; clarify in Phase 2 task tracking.**
3. **v1 expected-value tables are usable but not directly portable.** v1's 485 hard-coded expected values were generated under v1's mixed Jacobi/Gauss-Seidel state-read pattern; v3 enforces strict Jacobi-state semantics (resolved Q10). **Phase 7 should regenerate expected values for any test where the state-read pattern affects the result, rather than treating v1 numbers as ground truth.**

Phase 0 found **no design-blocking issues**. Phase 1 is unblocked.

---

## 2. Confirmation of the Spec's Component Inventory

### 2.1 Constituent inventory (Section 4)

The design spec's 11-Process organization is complete. Phase 0.1 catalogued ~290 v1 NSM1 functions and confirmed every one maps to one of the 11 Process classes or to a shared utility.

| Process class | v1 status | v2 status | v3 work |
|---|---|---|---|
| `FloatingAlgae` | Fully implemented | Partial; 6 critical bugs/stubs | Fix bugs, implement stubs, wire mortality routing |
| `BenthicAlgae` | Fully implemented | Partial; same bugs as floating | Same |
| `Nitrogen` | Fully implemented | Partial (NH4, NO3 only); 8 bugs | Fix bugs, **add OrgN as third state variable** |
| `Phosphorus` | Fully implemented | Absent | New (TIP + OrgP) |
| `Carbon` | Fully implemented | Absent | New (POC + DOC + DIC) |
| `POM` | Fully implemented | Absent | New |
| `CBOD` | Fully implemented | Absent | New (multi-group) |
| `DOX` | Fully implemented | Absent | New (largest single port) |
| `Pathogen` | Fully implemented | Absent | New |
| `Alkalinity` | Declared but inactive | Absent | New (simple tracer per resolved Q2) |
| `N2` | Fully implemented | Absent | New (with TDG) |

### 2.2 Shared physics primitives

Phase 0.1 confirmed all primitives the spec moves to `clearwater_modules_v3/utils/`:

- **Reaeration menu** (`kah_20`, `kaw_20`, `ka_tc`): 9 hydraulic + 13 wind options in `clearwater_modules/shared/processes.py`. Confirmed required by DOX.
- **Sediment oxygen demand** (`SOD_tc`): single function, used by DOX.
- **Light extinction** (`L`, `PAR`): Beer-Lambert with ISS, POC, Chl-a. Used by FloatingAlgae, BenthicAlgae, Pathogen.
- **Phosphorus partitioning** (`fdp`): used by Phosphorus and FloatingAlgae (replaces hard-coded `0.5`).
- **Arrhenius / temperature** (`arrhenius_correction`, `celsius_to_kelvin`): already in v2 utils, re-export from v3.

No additional primitives were identified beyond the spec's list.

### 2.3 Bug list (Section 6)

Phase 0.1 independently identified the same critical bug categories as Section 6:

- Multiplicative integrators (NH4, NO3, Ap)
- NaN guards using `== np.nan`
- Hard-coded zeros for algal nitrate uptake and oxygen-inhibition half-saturation
- `set_at_time` persistence omissions
- `ammonium_respiration()` and `ammonium_growth()` returning 0
- Hard-coded `phosphate_fraction_dissolved=0.5`

The 16-item enumeration in Section 6 remains the working bug-fix list for Phase 2.

---

## 3. New Findings Beyond the Spec

### 3.1 Critical: `pressure_mb=2026.5` is bug #7 in the sentinel-class set

**Severity:** High. The default is approximately 2× actual sea-level atmospheric pressure (~1013 hPa). It is consumed by:

- `pwv` (water vapor pressure) calculation
- `O2sat` (DOX saturation concentration via APHA/QUAL2E formulation, Section 5)
- `N2sat` (N2 saturation via Henry's law, Section 5)
- Atmospheric reaeration partial-pressure term

A simulation using the default has been computing saturation values against double atmospheric pressure, which systematically biases atmospheric exchange fluxes for both DOX and N2.

**Recommendation:** Add to Section 7 with default `1013.25` hPa (standard sea-level pressure).

### 3.2 Lower-priority parameter findings

Phase 0.2 flagged 8 additional defaults that warrant Phase 1 review but do not require spec amendments:

| Parameter | Value | Issue | Phase 1 action |
|---|---|---|---|
| `rnh4_20`, `vno3_20`, `rpo4_20` | 0 | Sediment release silently disabled | Verify gated by `use_SedFlux=False`; document |
| `kdpo4` | 0.0 | TIP partitioning feature disabled | Verify intent; NSM2 territory may apply |
| `ksbod_20` | 0.0 | CBOD never settles | Confirm with LimnoTech (intentional or bug?) |
| `apx`, `vx` | 1, 1 | Pathogen placeholders without literature basis | Document units and source |
| `h2` | 0.1 | POM dissolution depth, unclear physical role | Add docstring |
| `vb` | 0.01 | Burial velocity magnitude needs verification | Document with reference |
| `q_solar` units | 500 | Units unclear in code (`1/d` docstring is wrong) | Standardize as W/m² |
| `lambdas` | defined | Light extinction parameter defined but disabled in code | Remove or activate |

The first three (sediment-flux gating, kdpo4, ksbod_20) are the most consequential. The rest are documentation and clarification.

### 3.3 Test fixture porting nuance

Phase 0.3's main finding for Phase 7 planning: v1's 485 hard-coded expected values cannot be uncritically copy-pasted into v3 tests. The reason is the resolved Q10 decision:

- v1 NSM1 uses an xarray-DAG that effectively produces Jacobi-state / dependency-ordered rate semantics (ambiguous in places)
- v2 NSM1 reads via `process_other.method()` calls — accidentally Gauss-Seidel
- v3 NSM1 enforces strict Jacobi state, GS rates (resolved 2026-05-04)

For tests that exercise a single process in isolation (the bulk of test_7…test_17 — single-parameter perturbations, single-step runs), expected values transfer cleanly. For tests that depend on inter-process coupling within a step, expected values may differ at the 4th–6th decimal place. **Phase 7 strategy:** import the v1 expected-value tables as the reference target; for tests that fail with v3's stricter semantics, regenerate the expected value from v3 itself and document the v3-specific value with a comment noting the v1 baseline. This keeps the v1/v3 parity story honest.

The 4 Excel workbooks in `tests/NSM Manual Calcs/` (Alkalinity, Carbon, CBOD, DOX) cover only 4 of 16 constituents and are unused by the v1 test suite. Treat them as Phase 7 documentation only, not as ETL source.

---

## 4. Recommended Spec Amendments

Three minimal edits to the design spec, all in Section 7 (sentinel-999 corrections) and Section 4 (parameter library):

1. **Section 7 table:** add a row for `pressure_mb` with v1 default `2026.5`, v3 default `1013.25` hPa, source "standard sea-level pressure (ISO 2533)".
2. **Section 7 narrative:** retitle the section from "Sentinel-`999` Bug Fixes" to "Critical Default-Value Corrections" so the new entry, which is a magnitude error rather than a sentinel, fits naturally.
3. **Section 4.4 parameter inventory:** add a brief "Other defaults under Phase 1 audit" subsection pointing to the 8 lower-priority items in this gap analysis.

These edits are applied separately following acceptance of this synthesis.

---

## 5. Phase 1 Readiness Checklist

Phase 1 (foundation: shared utilities + parameter library) is unblocked. The audit clarifies the work order:

- [ ] **Section 7 amendment** — apply the three spec edits above before Phase 1 starts
- [ ] **Phase 1.1** — implement `clearwater_modules_v3/utils/` modules (`reaeration.py`, `sediment.py`, `light.py`, `partitioning.py`, **`numerics.py`** — the last is new since the resolved Q7 clipping decision)
- [ ] **Phase 1.2** — implement `clearwater_modules_v3/parameters/` modules per Appendix B; populate `DEFAULTS` from v1's TypedDicts; apply the 7 critical default corrections (6 sentinel-999 + `pressure_mb`)
- [ ] **Phase 1.3** — establish the `Process.DEFAULTS` merge pattern on the existing v2 `Nitrogen` process as a reference implementation
- [ ] **Phase 1.4** — write `parameter_defaults_corrections.md` documenting the 7 corrections + 8 audit follow-ups
- [ ] **Phase 1.5** — Tier 1 conservation test harness scaffolded (per resolved Q7: must assert `clip_events == 0` in addition to mass-balance roundoff)

Estimated wall-clock for Phase 1: 1–2 days as the spec projected.

---

## 6. Open Items for LimnoTech Communication

This synthesis surfaces three additional items for the LimnoTech communication that was already planned (the bundle covers the three tentative Q2/Q3/Q4 design decisions plus the bug list):

1. **`pressure_mb=2026.5`** — confirm this is a typo, not an unusual deliberate value in any LimnoTech application
2. **`ksbod_20=0`** (CBOD never settles) — confirm whether this is intentional model behavior or an oversight
3. **Sediment-flux gating** (`use_SedFlux` controlling `rnh4_20`, `vno3_20`, `rpo4_20`) — confirm the gate works as expected; if any application has been running with these parameters non-zero but `use_SedFlux=False`, the silent disable would have been masking the parameter

None of these block Phase 1 work; Phase 2 onward proceeds with the corrected defaults regardless of Anthony's reply.

---

## 7. Phase 0 Deliverables (Index)

| Deliverable | Path | Lines | Purpose |
|---|---|---|---|
| Constituent diff | `docs/clearwater_modules_v3_nsm1_phase0_constituent_diff.md` | 588 | Per-constituent v1/v2/v3 mapping |
| Parameter audit | `docs/clearwater_modules_v3_nsm1_phase0_parameter_audit.md` | 403 | Full ~250-parameter inventory |
| Test audit | `docs/clearwater_modules_v3_nsm1_phase0_test_audit.md` | 274 | v1/v2 test fixtures + workbooks |
| **Synthesis (this document)** | `docs/clearwater_modules_v3_nsm1_gap_analysis.md` | — | Executive summary + spec amendments |

Phase 0 is complete. Phase 1 is unblocked.
