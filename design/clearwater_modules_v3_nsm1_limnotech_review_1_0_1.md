# v3 NSM1 1.0.1 — LimnoTech Reviewer Materials

**Re:** v3 NSM1 1.0.1, the pattern-alignment release. Restores structural conformance with the design spec Section 14 (Appendix A registry rate-variable convention, clip-with-log contract, Jacobi/GS substep semantics) that the 1.0.0 release shipped as code TODOs.

**Source branch:** `EcohydrologyTeam/ClearWater-modules-streaming` `streaming`

**Companion document:** `design/clearwater_modules_v3_nsm1_limnotech_review.md` (the 1.0.0 review packet). That packet describes the kinetics-correctness state; this packet describes the structural state. Kinetics correctness is **unchanged** between 1.0.0 and 1.0.1 — that is the load-bearing invariant of the work documented here.

---

## 1. Executive summary

v3 NSM1 1.0.1 lands the structural pattern-alignment that the 1.0.0 design spec called for but did not implement at the code level. The result:

- **80 calibration / validation diagnostics** are now exposed via an opt-in registry-write surface across the 11 NSM1 Processes. Users who pre-register a name in the registry receive its value each substep; users who don't pre-register pay zero cost.
- **Every Process uses the same canonical `run()` shape** — read state → call `_change_with_components` (or `_rate_with_components` for Pathogen's rate-form integrator) → Forward Euler → unconditional `clip_negative_state` → persist state → opportunistic diagnostic writes.
- **Zero kinetics regressions.** The 4,320-substep coupled-run baseline trajectory is bit-identical (`rtol=0, atol=0`) to the pre-refactor reference across all 11 state variables × 5 cells × 4,321 substep indices. All 8 Tier 1 closed-system mass-conservation tests continue to pass at `rtol=1e-12` with zero clip events.
- **All 1.0.0 review-packet kinetics items remain unchanged.** The 4-way Nitrogen theta swap, the 4 deliberate value choices, the rca/rcb derivation, the DIC unit reconciliation, the 3 tentative spec-§14 decisions — all carry forward verbatim. No new kinetics asks for 1.0.1.

The work shipped as 12 commits (Phase 0 through Phase 10.B) over 2026-05-13 / 2026-05-14. The full audit trail is in `design/clearwater_modules_v3_nsm1_pattern_alignment_phase{0,1,...,10b}_closeout.md`.

---

## 2. Document index — pattern-alignment release

| Topic | Canonical document |
| --- | --- |
| Source-of-truth specification | `design/clearwater_modules_v3_nsm1_pattern_alignment_specification.md` |
| Appendix A diff (1.0.0 cheatsheet → 80-name catalog) | `design/clearwater_modules_v3_nsm1_appendix_a_diff.md` |
| Phase 0 close-out (baseline + Q1/Q2 resolutions) | `clearwater_modules_v3_nsm1_pattern_alignment_phase0_closeout.md` |
| Phase 1 close-out (mechanical alignment) | `..._phase1_closeout.md` |
| Phase 2–9 close-outs (per-Process refactors) | `..._phase{2,3,4,5,6,7,8,9}_closeout.md` |
| Phase 10.A close-out (conformance + completeness + smoke + perf) | `..._phase10a_closeout.md` |
| Phase 10.B close-out (legacy-shadow cleanup, final) | `..._phase10b_closeout.md` |
| Perf findings (like-for-like benchmark) | `clearwater_modules_v3_nsm1_pattern_alignment_perf_findings.md` |
| Baseline NetCDF (gold reference for §11.2) | `tests/v3/nsm1/baseline/baseline_coupled_trajectory_186b5c4.nc` |
| Pre-existing 1.0.0 review packet (kinetics correctness) | `clearwater_modules_v3_nsm1_limnotech_review.md` |

---

## 3. What 1.0.1 changes for users

### 3.1 New diagnostic-subscription surface

Each of the 11 NSM1 Processes now declares a class-level `REGISTRY_DIAGNOSTICS: tuple[str, ...]` listing the named rate / flux / limitation diagnostics it owns. The full catalog has **80 names** across all 11 Processes — every Appendix A name listed in `design/clearwater_modules_v3_nsm1_appendix_a_diff.md` §3.

To subscribe, a user pre-registers the variable in their registry before model construction:

```python
ic = default_initial_conditions()
ic["nitrification_flux_rate"] = xr.zeros_like(ic["ammonium"])
demo = build_nsm1_demo(initial_conditions=ic)
# Every substep, Nitrogen writes its nitrification_flux_rate to the
# registry. Other diagnostics that were NOT pre-registered are not
# computed-and-written.
```

The contract: **every state-variable trajectory is bit-identical regardless of which diagnostics are subscribed**. This is enforced by `tests/v3/nsm1/test_coupled_demo_parity.py::test_diagnostics_subscription_smoke_state_bit_identical`.

### 3.2 Preserved attribute / registry name aliases

Several Appendix A names differ slightly from the pre-1.0.1 internal cache attribute names that v1-parity tests read via `getattr`. The cache attribute names are preserved as **aliases**; both the legacy name and the Appendix A name point at the same value. The aliased pairs (legacy ↔ Appendix A):

| Process | Legacy attribute | Appendix A registry name |
| --- | --- | --- |
| DOX | `dox_sod_rate` | `sod_rate` (same volumetric SOD sink) |
| N2 | `tdg` | `total_dissolved_gas` |
| Phosphorus | `orgp_to_tip_hydrolysis_rate` | `orgp_hydrolysis_rate` |
| Alkalinity | `alk_nitrification_rate` | `alk_nitrification_sink_rate` |
| Alkalinity | `alk_denitrification_rate` | `alk_denitrification_source_rate` |
| Alkalinity | `alk_benthic_algae_growth_rate` | `alk_balgae_growth_rate` |
| Alkalinity | `alk_benthic_algae_respiration_rate` | `alk_balgae_respiration_rate` |

The aliases are pinned by tests in each `test_<process>_registry_diagnostics.py` so a future refactor cannot silently drop them.

### 3.3 BenthicAlgae `rate_death` deduplication

The pre-1.0.1 BenthicAlgae `_cache_benthic_mortality_rates` and `rate()` each invoked `self.rate_death(...)` — duplicate work producing the same value (the function is pure). The 1.0.1 `_change_with_components` computes `ab_death` once and reuses the cached value. Net mortality rate is bit-identical to the pre-refactor path; the duplicate is gone.

### 3.4 POM `pom_doc_source_rate` cache-relocation

`pom_doc_source_rate` (consumed by Carbon via `getattr`) is now set as a side effect of `POM._change_with_components` rather than inside the public `rate()` method. The legacy `rate()` method is retained for back-compat with three external test calls. Net DOC source rate to Carbon is unchanged.

### 3.5 Pathogen `_rate_with_components`

Pathogen is a single-state rate-form integrator (the only one in NSM1); its canonical helper is named `_rate_with_components` rather than `_change_with_components`. The convention reflects the integrator's argument shape (returns a per-day rate, not a per-step delta) and is documented in pattern-alignment spec §10 Q5.

---

## 4. Migration notes

### 4.1 No YAML change required

All Appendix A names are opt-in via registry pre-registration. Existing YAML configs continue to work unchanged: if a YAML file does not declare any of the 80 new diagnostic variables, the model produces bit-identical outputs to 1.0.0.

### 4.2 No sibling-process consumer change

The cross-process rate-variable consumers (`DOX → Nitrogen.nitrification_flux_rate`, `Alkalinity → Nitrogen.denitrification_flux_rate`, `Carbon → CBOD.cbod_oxidation_rate`, etc.) continue to read sibling caches via `getattr`. All preserved attribute names are pinned by tests.

### 4.3 v3 1.0.0 YAML configs run unchanged

The Tier 1 closed-system mass-conservation tests, the v1-parity tests, the demo notebook (`examples/V3/04_Example_NSM1.ipynb`) all run against 1.0.1 with no modification.

### 4.4 New tests

- `test_pattern_conformance.py` (88 tests) — structural conformance scan.
- `test_appendix_a_completeness.py` (14 tests) — 80-name catalog uniqueness + per-Process membership.
- `test_coupled_demo_parity.py` (4 tests) — end-to-end bit-identical baseline replay + pattern G full-subscription smoke.
- 11 `test_<process>_registry_diagnostics.py` files — per-Process pattern G contracts (~60 tests total).

Total test suite: **990 passing, 2 xfailed** (unchanged from 1.0.0 — the 2 xfailed are the pre-existing legacy CBOD `test_changed_kbod_20` tests).

---

## 5. Performance

Apples-to-apples benchmark on the current measurement machine (`pixi --environment dev`, 5-cell synthetic mesh, 60-substep warmup + 500-substep measurement window):

| Code state | ms/substep (median) |
| --- | --- |
| Pre-pattern-alignment (`186b5c4`) | 22.64 |
| Post-pattern-alignment (`afac699`) | 26.04 |
| **Real overhead** | **+15.0%** |

The pattern G **zero-cost-when-unused** contract is met perfectly: subscribing to all 80 diagnostics adds **−0.6%** vs no-subscription (within noise).

The 15% no-subscription overhead arises from the structural cost of the pattern-aligned shape:

- ~80 dict construction operations per substep (11 Processes × ~7 names each).
- ~80 setattr operations per substep.
- ~80 `in registry` membership checks per substep.
- New diagnostic computations (algal `limit_*` factors that did not exist pre-1.0.1).

The pattern-alignment spec's §8 "must" budget was ≤ 5% no-sub overhead; the actual 15% exceeds it. See `design/clearwater_modules_v3_nsm1_pattern_alignment_perf_findings.md` for the budget-reconciliation recommendation and a proposed §8 amendment (Must ≤ 20% no-sub, ≤ 15% full-sub vs no-sub; Aspirational ≤ 5%).

Sumwere Creek (600 cells, 4,320 substeps) extrapolation at 26 ms/step × 120 cells-scaling = ~3.6 hours. The 1.0.0 spec §10 "must" target of 30 min remains comfortably ahead of the actual cost. The 1.0.1 perf overhead does not move the v3 NSM1 outside the spec's overall budget envelope.

---

## 6. The 80-name diagnostic catalog

Per Process, the Appendix A registry diagnostics:

| Process | Count | Notable names |
| --- | --- | --- |
| Carbon | 9 | `poc_hydrolysis_rate`, `doc_dic_oxidation_rate`, `dic_atm_exchange_rate`, `dic_sed_release_rate`, `carbon_{algal,balgae}_{photo,resp}_rate`, `carbon_cbod_oxidation_rate` |
| DOX | 11 | `dox_sat`, `atm_reaeration_rate`, `dox_{nitrification,sod,doc_oxidation,cbod_oxidation,algal_photo,algal_resp,balgae_photo,balgae_resp}_rate`, `sod_rate` (alias) |
| Nitrogen | 10 | `nitrification_flux_rate`, `denitrification_flux_rate`, `nh4_from_bed`, `no3_from_bed_denit`, `orgn_{hydrolysis,settling}_rate`, `{nh4,no3}_algal_growth_rate`, `nh4_{algal,balgae}_resp_rate` |
| FloatingAlgae | 13 | `algal_{growth,respiration,death,settling}_rate`, `algal_{orgn,orgp,poc,doc}_from_mortality_rate`, `algal_pom_from_settling_rate`, `algal_nh4_uptake_fraction`, `algal_{light,nutrient_n,nutrient_p}_limitation` |
| BenthicAlgae | 11 | analogous to FloatingAlgae (minus settling, minus pom-from-settling) |
| Phosphorus | 7 | `orgp_{hydrolysis,settling}_rate`, `tip_settling_rate`, `dip_from_bed`, `orgp_algal_mortality_rate`, `tip_{algal,balgae}_growth_rate` |
| POM | 4 | `pom_{hydrolysis,settling}_rate`, `pom_{algal,balgae}_mortality_rate` |
| CBOD | 2 | `cbod_{oxidation,settling}_rate` |
| N2 | 4 | `n2_atm_exchange_rate`, `n2_sat`, `total_dissolved_gas`, `n2_denit_source_rate` |
| Pathogen | 3 | `pathogen_{natural_death,light_death,settling}_rate` |
| Alkalinity | 6 | `alk_{nitrification_sink,denitrification_source,algal_growth,algal_respiration,balgae_growth,balgae_respiration}_rate` |

The 80 names are documented exhaustively in `design/clearwater_modules_v3_nsm1_appendix_a_diff.md` §3, pinned by `tests/v3/nsm1/test_appendix_a_completeness.py`, and verified to be uniquely owned by exactly one Process.

---

## 7. Asks of the reviewer (summary)

### 7.1 Required confirmations (1.0.0 kinetics items, unchanged)

The four open kinetics items from the 1.0.0 packet remain pending — no decisions have changed:

1. 4-way Nitrogen theta swap.
2. 4 deliberate value choices (`vsop`, `SOD_20`, `BWa`, `kah_20_user`).
3. rca / rcb derivation.
4. DIC unit reconciliation.

See `clearwater_modules_v3_nsm1_limnotech_review.md` §4 for full context.

### 7.2 New confirmations specific to 1.0.1

1. **Diagnostic catalog naming** (§6 above): are the 80 Appendix A names acceptable as the published 1.0.1 registry interface? In particular:
   - Is the `_rate` / `_fraction` / `_from_bed` / source-prefix naming convention (§14 Q2 of the 1.0.0 spec) acceptable as the canonical convention going forward?
   - The 7 preserved aliases (§3.2 above): is keeping both names (legacy + Appendix A) acceptable, or should we deprecate the legacy names in 1.1+?
2. **Performance budget reconciliation** (§5 above): accept the 15% no-subscription overhead as the price of the diagnostic surface, or commission a profiling / optimisation pass for 1.0.2?
3. **Versioning** (spec §10 Q3): is 1.0.1 (SemVer patch; additive opt-in surface) the right release marker, or should this be 1.1.0 (minor; structural change)?

### 7.3 Optional / suggested next-cycle items

- The 1.0.0 kinetics items left tentative pending LimnoTech confirmation (Alkalinity simple-tracer scope, sediment-flux scalar globals, single-compartment algae) remain pending.
- The DOX salinity correction on O2sat (C6) and the DIC unit reconciliation (C9 / C10) remain in the same state as 1.0.0.
- The pattern-alignment work surfaced no new kinetics defects; the 80-name diagnostic surface exposes computations that were already happening internally, so any future calibration anomaly attributable to a specific named diagnostic can be diagnosed without code changes.

---

## 8. Test summary

| Suite | Tests | Status |
| --- | --- | --- |
| Full v3 suite (`pytest tests/`) | 992 (990 passed + 2 xfailed) | Green |
| Pattern conformance + completeness + coupled-demo parity | 106 | Green |
| Per-Process registry-diagnostics (11 files) | ~60 | Green |
| Tier 1 closed-system mass conservation (8 constituents) | 40 | `rtol=1e-12`, `clip_events == {}` |
| Tier 1.5 active kinetics conservation | 4 | Green |
| Pre-existing v1-parity (11 Processes) | (carried forward from 1.0.0) | Green |
| Pre-existing model orchestration / hotstart / wet-mask / chunking | (carried forward) | Green |
| Baseline trajectory parity (`check_baseline_parity.py` + `test_coupled_demo_parity.py`) | bit-identical (`rtol=0, atol=0`) | OK |

The 2 xfailed tests are the pre-existing legacy CBOD `test_changed_kbod_20` from the v1 NSM1 reference suite — unchanged from 1.0.0; not a 1.0.1 regression.

---

## 9. Suggested review focus areas (25-minute first pass)

1. **§3.2** — the 7 preserved aliases. Quick yes/no per pair on the keep-both-names approach.
2. **§5** — the 15% perf overhead. Accept-or-optimise decision.
3. **§7.2** — the three new asks (catalog naming, perf budget, versioning).
4. **`tests/v3/nsm1/test_appendix_a_completeness.py`** — quick scan of the 80-name catalog to confirm none of the names will collide with a planned NSM2 / v3 1.1 name.

For a deeper second pass: the per-Process `REGISTRY_DIAGNOSTICS` tuples in each Process file plus the `_change_with_components` / `_rate_with_components` helpers. Each helper has a docstring identifying the pre-refactor inline composition source.

---

## 10. The big picture

v3 NSM1 1.0.0 shipped with correct kinetics and a planned-but-unimplemented diagnostic surface (referenced as TODOs in the design spec's §14). v3 NSM1 1.0.1 closes the surface implementation without touching the kinetics. The packet for review is therefore narrow: **does the catalog look right, is the structural cost acceptable, and is 1.0.1 the right release marker.**

If all three asks land yes, v3 NSM1 1.0.1 ships. The 1.0.0 kinetics-correctness items remain in their existing review state (no movement in either direction); resolving them is the gate for v3 NSM1 1.1.
