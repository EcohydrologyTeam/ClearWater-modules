# Phase 7 Close-Out — POM & CBOD `_change_with_components`

**Date:** 2026-05-13
**Branch:** `streaming`
**Baseline:** Phase 0 commit `d862d68`
**Phase 6 base:** commit `89ffe38`
**Spec reference:** `design/clearwater_modules_v3_nsm1_pattern_alignment_specification.md` §6 Phase 7

Phase 7 applies the canonical refactor template to POM and CBOD together. **POM finally relocates the `pom_doc_source_rate` cache from `rate()` into `_change_with_components`** (Phase 1's deferred item).

---

## What changed

### POM

**`REGISTRY_DIAGNOSTICS`** — 4 names per Appendix A diff §3:

- `pom_hydrolysis_rate` — POM dissolution rate (mg-D/L_sed/d), the raw `kpom_tc * pom`. Distinct from `pom_doc_source_rate` (the consumer-ready unit-converted DOC source).
- `pom_settling_rate` — POM burial sink (mg/L/d). Spec uses "settling" but the term is `vb * pom / h2`; documented in code.
- `pom_algal_mortality_rate` — floating-algae settling source.
- `pom_balgae_mortality_rate` — benthic-algae mortality source.

**`pom_doc_source_rate`** is **NOT** in `REGISTRY_DIAGNOSTICS` — it remains a sibling-cache attribute that Carbon reads via `getattr` (NOT exposed to the registry). Pinned by `test_g_pom_doc_source_rate_remains_on_self_after_run`.

**Phase 7 cache-relocation (deferred from Phase 1):** the `self.pom_doc_source_rate = ...` assignment moves from inside `rate()` (line 304-306 pre-Phase-7) into `_change_with_components`. The new helper takes `depth` as an explicit argument (no longer reading it from the registry inside the helper, which `rate()` did). The `rate()` method is no longer called by `run`; `_change_with_components` inlines what `rate()` did and adds the components dict.

The legacy `rate(...)` method retains its prior behaviour for back-compat with external test calls (`tests/v3/nsm1/test_pom_v1_parity_v3.py`); it still reads `depth` from the registry and sets `self.pom_doc_source_rate` as a side effect. Its arithmetic is unchanged.

**`_change_with_components(...)`** returns `(rate, components)`. Operand order, sub-flux helpers, and the use-flag gating preserved verbatim.

**`_change_legacy_inline(...)`** shadow returns just the net rate. Verbatim copy of pre-Phase-7 `rate()` body. Deleted in Phase 10.

### CBOD

**`REGISTRY_DIAGNOSTICS`** — 2 names per Appendix A diff §3:

- `cbod_oxidation_rate` (preserved attribute name; consumed by DOX and Carbon via `getattr`).
- `cbod_settling_rate` (preserved attribute name; was already cached on self).

**`_change_with_components(...)`** returns `(rate, components)`. The net rate is `-oxidation_rate - settling_rate` (the existing convention; CBOD's state decreases via both terms). Per-sub-flux `sanitize_rate` calls preserved verbatim from Phase 1.E.

**`_change_legacy_inline(...)`** shadow returns just the net rate. Verbatim copy of pre-Phase-7 `run` body. Deleted in Phase 10.

---

## New tests

| File | Tests | Retained / deleted |
|---|---|---|
| `tests/v3/nsm1/test_pom_helper_vs_inline.py` | 14 (6 scenarios × parity + 6 components-set + Phase 7 cache-relocation parity tests + zero-state finiteness) | Deleted in Phase 10 |
| `tests/v3/nsm1/test_pom_registry_diagnostics.py` | 6 (pattern G full contract + pom_doc_source_rate-stays-off-registry) | Retained |
| `tests/v3/nsm1/test_cbod_helper_vs_inline.py` | 19 (5 scenarios × 2 use_DOX modes parity + components-set + positive-magnitudes + zero-state finiteness + hypoxic monotonicity) | Deleted in Phase 10 |
| `tests/v3/nsm1/test_cbod_registry_diagnostics.py` | 6 (pattern G full contract + cbod_oxidation_rate sibling-consumer visibility) | Retained |

---

## §11 Zero-Regression Contract Verification

| Clause | Required | Achieved |
|---|---|---|
| §11.2 bit-identical state trajectory | `rtol=0, atol=0` against `baseline_coupled_trajectory_186b5c4.nc` | **OK** |
| §11.3 helper-vs-inline parity | bit-identical helper-vs-shadow across the input matrices | **OK** — 33/33 passed (14 POM + 19 CBOD) |
| §11.4 full test suite | 978 + 45 new = 1,023 passed, 2 xfailed | **OK** — 1,023 passed, 2 xfailed, 123.42 s wall |
| §11.4 Tier 1 conservation | `rtol=1e-12`, `clip_events == {}` (40 tests) | **OK** — 40 passed |

The bit-identical baseline parity confirms the **Phase 7 cache-relocation** (`pom_doc_source_rate` from `rate()` to the helper) is numerically perfect — Carbon continues to read identical values via `getattr(pom_process, 'pom_doc_source_rate')`.

No regressions. No silent tolerance relaxation.

---

## Files changed

- `src/clearwater_modules_v3/processes/pom.py` — POM refactor + Phase 1 cache-relocation.
- `src/clearwater_modules_v3/processes/cbod.py` — CBOD refactor.
- `tests/v3/nsm1/test_pom_helper_vs_inline.py` — new (14 tests).
- `tests/v3/nsm1/test_pom_registry_diagnostics.py` — new (6 tests).
- `tests/v3/nsm1/test_cbod_helper_vs_inline.py` — new (19 tests).
- `tests/v3/nsm1/test_cbod_registry_diagnostics.py` — new (6 tests).

---

## Phase 1 deferred items resolved

Phase 7 closes out **2 of the 3 Phase 1 deferred items**:

- ✅ POM cache-relocation (moved `pom_doc_source_rate` from `rate()` to `_change_with_components` per Phase 1.C / Phase 7 spec).
- ✅ BenthicAlgae `rate_death` dedup (resolved in Phase 5).
- N/A CBOD `init_process` sibling discovery (closed in Phase 1 as not-applicable; CBOD reads DOX from registry, not via sibling-process cache).

The full Phase 1 deferral list from `clearwater_modules_v3_nsm1_pattern_alignment_phase1_closeout.md` is now fully resolved.

---

## What changes for Phase 8

Phase 8 covers N2 and Pathogen together (~1 day). N2 already has the opportunistic-write loop for `total_dissolved_gas` (the only pre-existing example of pattern G); Phase 8 extends it to the full Appendix A set. Pathogen is a single-state Process; Phase 8 introduces `_rate_with_components` (the per-Process naming-by-integrator-form convention from spec §10 Q5).
