# Phase 10.B Close-Out — Legacy-Inline-Shadow Cleanup

**Date:** 2026-05-14
**Branch:** `streaming`
**Baseline:** Phase 0 commit `d862d68`
**Phase 10.A base:** commit `7ead050`
**Spec reference:** `design/clearwater_modules_v3_nsm1_pattern_alignment_specification.md` §6 Phase 10

Phase 10.B is the **final commit** of the pattern-alignment work. It deletes the `_change_legacy_inline` / `_rate_legacy_inline` shadow methods from every Process and removes the 11 `test_<process>_helper_vs_inline.py` parity-test files. The pattern-aligned canon is now the single source of truth.

---

## What changed

### Shadow methods deleted

10 `_change_legacy_inline` methods (Alkalinity, BenthicAlgae, Carbon, CBOD, DOX, FloatingAlgae, N2, Nitrogen, Phosphorus, POM) and 1 `_rate_legacy_inline` method (Pathogen). Total: ~720 lines removed.

| File | Lines removed |
|---|---|
| `alkalinity.py` | 46 |
| `benthic_algae.py` | 59 |
| `carbon.py` | 125 |
| `cbod.py` | 45 |
| `dox.py` | 93 |
| `floating_algae.py` | 54 |
| `n2.py` | 74 |
| `nitrogen.py` | 56 |
| `pathogen.py` | 34 |
| `phosphorus.py` | 93 |
| `pom.py` | 71 |

### Test files deleted

11 `test_<process>_helper_vs_inline.py` files (~270 parity tests). The helper-vs-inline contract is no longer relevant because the shadow methods are gone; the canonical `_change_with_components` / `_rate_with_components` helpers are the sole rate-composition path.

| File | Tests |
|---|---|
| `test_alkalinity_helper_vs_inline.py` | 10 |
| `test_benthic_algae_helper_vs_inline.py` | 18 |
| `test_carbon_helper_vs_inline.py` | 15 |
| `test_cbod_helper_vs_inline.py` | 19 |
| `test_dox_helper_vs_inline.py` | 21 |
| `test_floating_algae_helper_vs_inline.py` | 17 |
| `test_n2_helper_vs_inline.py` | 14 |
| `test_nitrogen_helper_vs_inline.py` | 27 |
| `test_pathogen_helper_vs_inline.py` | 14 |
| `test_phosphorus_helper_vs_inline.py` | 24 |
| `test_pom_helper_vs_inline.py` | 14 |
| **Total** | **193** |

### Docstring stragglers cleaned

The `_change_with_components` / `run` docstrings in every Process previously had two cross-references to the shadow methods:

- "See `_change_legacy_inline` for the pre-Phase-N inline composition…"
- "The companion shadow `_change_legacy_inline` returns just…"

Both references are stripped. The canonical helpers' docstrings now stand on their own.

---

## §11 Zero-Regression Contract Verification

| Clause | Required | Achieved |
|---|---|---|
| §11.2 bit-identical state trajectory | `rtol=0, atol=0` against `baseline_coupled_trajectory_186b5c4.nc` | **OK** |
| §11.4 full test suite | 1,187 − 197 deleted helper-vs-inline tests = 990 passed, 2 xfailed | **OK** — 990 passed, 2 xfailed, 239.92 s wall |
| §11.4 Tier 1 conservation | `rtol=1e-12`, `clip_events == {}` (40 tests) | **OK** — 40 passed |

The bit-identical baseline trajectory and Tier 1 conservation confirm that the canonical `_change_with_components` helpers carry the full kinetics-correctness contract. The helper-vs-inline parity contract that the shadow methods supported was the temporary scaffold that proved the refactor was numerically perfect at each step; it has served its purpose and is now retired.

---

## Test suite after Phase 10.B

- **Pattern conformance** (`test_pattern_conformance.py`, 88 tests, retained): structural scan over every Process.
- **Appendix A completeness** (`test_appendix_a_completeness.py`, 14 tests, retained): catalog uniqueness + 80-name pin.
- **Coupled demo parity + diagnostics smoke** (`test_coupled_demo_parity.py`, 4 tests, retained): end-to-end §11.2 baseline replay + pattern G full-subscription smoke.
- **11 registry-diagnostics tests** (`test_<process>_registry_diagnostics.py`, retained): per-Process pattern G + sibling-consumer contracts.
- **8 Tier 1 conservation tests** (retained): closed-system mass-balance closure at `rtol=1e-12`.
- **Pre-existing test files** (retained): `test_<process>_v1_parity_v3.py`, Tier 1.5 active kinetics, model orchestration, etc.

Total: 990 passing tests, 2 xfailed (both the pre-existing legacy CBOD `test_changed_kbod_20` from the v1 NSM1 reference tests).

---

## Files changed

- `src/clearwater_modules_v3/processes/alkalinity.py` — shadow + docstring refs deleted.
- `src/clearwater_modules_v3/processes/benthic_algae.py` — same.
- `src/clearwater_modules_v3/processes/carbon.py` — same.
- `src/clearwater_modules_v3/processes/cbod.py` — same.
- `src/clearwater_modules_v3/processes/dox.py` — same.
- `src/clearwater_modules_v3/processes/floating_algae.py` — same.
- `src/clearwater_modules_v3/processes/n2.py` — same.
- `src/clearwater_modules_v3/processes/nitrogen.py` — same.
- `src/clearwater_modules_v3/processes/pathogen.py` — same (`_rate_legacy_inline` deleted).
- `src/clearwater_modules_v3/processes/phosphorus.py` — same.
- `src/clearwater_modules_v3/processes/pom.py` — same.
- 11 `tests/v3/nsm1/test_<process>_helper_vs_inline.py` — deleted.

---

## Pattern-alignment release: 10 phases complete

| Phase | Date | Commit | Summary |
|---|---|---|---|
| 0 | 2026-05-13 | `d862d68` | Baseline capture, Q1/Q2 resolutions, `Diagnostics.current_step` |
| 1 | 2026-05-13 | `481532f` | Mechanical alignment pass (8 Process files) |
| 2 | 2026-05-13 | `92cfa58` | Carbon `_change_with_components` (canonical template) |
| 3 | 2026-05-13 | `b0f2af6` | DOX |
| 4 | 2026-05-13 | `78134a1` | Nitrogen (preserved flux-rate attribute names) |
| 5 | 2026-05-13 | `7725f1a` | FloatingAlgae + BenthicAlgae (with `rate_death` dedup) |
| 6 | 2026-05-13 | `89ffe38` | Phosphorus (with `orgp_hydrolysis_rate` alias) |
| 7 | 2026-05-13 | `a15bc75` | POM + CBOD (Phase 1 `pom_doc_source_rate` relocation resolved) |
| 8 | 2026-05-13 | `b19fdfd` | N2 + Pathogen (Pathogen `_rate_with_components`) |
| 9 | 2026-05-13 | `bfec53e` | Alkalinity (with 4 legacy aliases) |
| 10.A | 2026-05-14 | `7ead050` | Conformance + completeness + smoke + perf |
| 10.B | 2026-05-14 | this commit | Legacy-shadow cleanup |

**80 Appendix A registry-diagnostic names** exposed across 11 pattern-aligned Processes. **Pattern G zero-cost-when-unused verified** (full-sub vs no-sub: −0.6%). **Bit-identical baseline trajectory** verified at every commit since Phase 0.

The pattern-alignment work as scoped by `clearwater_modules_v3_nsm1_pattern_alignment_specification.md` is complete.

---

## Recommended Phase 11 / v3 1.0.1 follow-ups

1. **Like-for-like perf benchmark** on the current machine, comparing the pre-pattern-alignment commit (e.g. `186b5c4`) to the current branch tip. Resolves whether the +50% no-subscription overhead vs the Phase 0 documented 17.6 ms/step is a hardware artifact, a real overhead, or both.
2. **If real overhead**: profile the hot path. Likely candidates per Phase 10.A closeout: setattr-loop overhead, duplicate sub-flux computations in algal-coupling components, per-cell membership-check overhead.
3. **LimnoTech review packet refresh** — `clearwater_modules_v3_nsm1_limnotech_review.md` updates to document the 80-name diagnostic surface and the pattern-alignment as part of the 1.0.1 release content.
4. **Versioning decision** (spec §10 Q3) — 1.0.1 (patch; additive opt-in surface) vs 1.1.0 (minor; structural change). Recommended in the spec: 1.0.1 with 1.1.0 reserved for first kinetics-additive work.
