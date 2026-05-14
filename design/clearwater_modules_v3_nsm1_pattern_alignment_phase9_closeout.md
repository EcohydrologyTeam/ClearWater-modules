# Phase 9 Close-Out — Alkalinity

**Date:** 2026-05-14
**Branch:** `streaming`
**Baseline:** Phase 0 commit `d862d68`
**Phase 8 base:** commit `b19fdfd`
**Spec reference:** `design/clearwater_modules_v3_nsm1_pattern_alignment_specification.md` §6 Phase 9

Phase 9 applies the canonical refactor template to Alkalinity — the **last per-Process phase**. After this commit, all 11 v3 NSM1 Processes carry the canonical `_change_with_components` / `_rate_with_components` shape with `REGISTRY_DIAGNOSTICS` tuples and pattern G opportunistic-write loops.

---

## What changed

### `Alkalinity.REGISTRY_DIAGNOSTICS`

Six names per Appendix A diff §3:

- `alk_nitrification_sink_rate` — **renamed** from legacy `alk_nitrification_rate`; legacy name preserved as alias.
- `alk_denitrification_source_rate` — **renamed** from legacy `alk_denitrification_rate`; legacy name preserved as alias.
- `alk_algal_growth_rate` (already matched Appendix A naming).
- `alk_algal_respiration_rate` (already matched Appendix A naming).
- `alk_balgae_growth_rate` — **renamed** from legacy `alk_benthic_algae_growth_rate`; legacy name preserved as alias.
- `alk_balgae_respiration_rate` — **renamed** from legacy `alk_benthic_algae_respiration_rate`; legacy name preserved as alias.

The four legacy attribute names are read by `test_alkalinity_v1_parity_v3.py` and `test_alkalinity_tier1.py` via `getattr`; `_change_with_components` sets all four as side-effect aliases so the existing tests continue to pass without modification. Pinned by `test_g_legacy_attribute_aliases_remain_set` and `test_legacy_attribute_aliases_match_components`.

### `_change_with_components(...)`

Takes only `depth` as a kwarg — the per-source / per-sink sub-fluxes are read internally via sub-flux helpers that consult sibling-Process caches (`Nitrogen.nitrification_flux_rate`, `Nitrogen.denitrification_flux_rate`, `FloatingAlgae.algal_growth_rate / algal_respiration_rate`, `BenthicAlgae.balgae_growth_rate / balgae_respiration_rate`). Returns `(rate, components)`. Operand order and the v1 `dAlkdt` sign convention preserved verbatim.

Side effect: sets the four legacy attribute aliases plus `self.alk_rate` (used by some test fixtures) so the shadow and helper leave identical self state.

### `_change_legacy_inline(...)`

Returns just the net rate. Verbatim copy of pre-Phase-9 `run` body. Sets all six pre-Phase-9 cache attribute names as side effects (matches the exact pre-Phase-9 self state). Deleted in Phase 10.

### `Alkalinity.run` rewired

Reads alkalinity / water_temperature / depth → `_change_with_components(depth=depth)` → setattr loop (writes the 6 Appendix A names + the legacy `alk_rate`) → Forward Euler → unconditional `clip_negative_state` → persist → opportunistic-write loop.

---

## New tests

### `tests/v3/nsm1/test_alkalinity_helper_vs_inline.py` (10 tests; deleted in Phase 10)

- 4-scenario depth matrix (`uniform`, `thin`, `deep`, `randomised`) × bit-identical parity (no siblings wired).
- 4 components-set pinning tests across the same depth matrix.
- Zero-state finiteness.
- Legacy-attribute-alias parity (the four `alk_*_rate` legacy names equal their Appendix A counterparts after the helper runs).
- Helper-vs-shadow parity with mock-injected Nitrogen sibling caches (exercises the sub-flux helper path).

### `tests/v3/nsm1/test_alkalinity_registry_diagnostics.py` (7 tests; retained)

Pattern G full contract: pre-registration writes, no-subscription reproducibility, zero-cost-when-unused state invariant, attribute / registry value parity, partial subscription, and the legacy-attribute-aliases-stay-set-on-self contract.

---

## §11 Zero-Regression Contract Verification

| Clause | Required | Achieved |
|---|---|---|
| §11.2 bit-identical state trajectory | `rtol=0, atol=0` against `baseline_coupled_trajectory_186b5c4.nc` | **OK** |
| §11.3 helper-vs-inline parity | bit-identical across the depth matrix and with mock siblings | **OK** — 10/10 passed |
| §11.4 full test suite | 1,064 + 17 new = 1,081 passed, 2 xfailed | **OK** — 1,081 passed, 2 xfailed, 141.47 s wall |
| §11.4 Tier 1 conservation | `rtol=1e-12`, `clip_events == {}` (40 tests) | **OK** — 40 passed |

No regressions. The four legacy attribute aliases are preserved exactly; existing parity / Tier 1 tests pass without modification.

---

## Files changed

- `src/clearwater_modules_v3/processes/alkalinity.py` — Alkalinity refactor.
- `tests/v3/nsm1/test_alkalinity_helper_vs_inline.py` — new (10 tests).
- `tests/v3/nsm1/test_alkalinity_registry_diagnostics.py` — new (7 tests).

---

## All 11 Processes pattern-aligned

After Phase 9 the full v3 NSM1 1.0.x pattern-alignment is structurally complete:

| Process | Phase | Helper | Diagnostics | Legacy-name aliases |
|---|---|---|---|---|
| Temperature (TSM) | — | (canonical exemplar) | n/a (single state) | — |
| Carbon | 2 | `_change_with_components` | 9 | — |
| DOX | 3 | `_change_with_components` | 11 | `sod_rate` ↔ `dox_sod_rate` |
| Nitrogen | 4 | `_change_with_components` | 10 | `nitrification_flux_rate` / `denitrification_flux_rate` (preserved) |
| FloatingAlgae | 5 | `_change_with_components` | 13 | 10 preserved |
| BenthicAlgae | 5 | `_change_with_components` (with `rate_death` dedup) | 11 | preserved sibling-consumer names |
| Phosphorus | 6 | `_change_with_components` | 7 | `orgp_hydrolysis_rate` ↔ `orgp_to_tip_hydrolysis_rate` |
| POM | 7 | `_change_with_components` (with `pom_doc_source_rate` relocation) | 4 | `pom_doc_source_rate` off-registry |
| CBOD | 7 | `_change_with_components` | 2 | preserved |
| N2 | 8 | `_change_with_components` | 4 | `self.tdg` ↔ `self.total_dissolved_gas` |
| Pathogen | 8 | `_rate_with_components` (rate-form per §10 Q5) | 3 | — |
| Alkalinity | 9 | `_change_with_components` | 6 | 4 preserved (`alk_*_rate` legacy names) |

**Total**: 11 Processes × ~7 diagnostics each = ~80 registry-diagnostic names exposed via pattern G opportunistic writes; ~270 helper-vs-inline parity tests across Phases 2–9.

---

## What changes for Phase 10

Phase 10 is the final phase (~1.5 days per spec). Deliverables:

1. **Conformance test** (`tests/v3/nsm1/test_pattern_conformance.py`) — single test that iterates every Process class and asserts canonical shape: `_change_with_components` (or `_rate_with_components`) exists; `REGISTRY_DIAGNOSTICS` is a non-empty tuple; `init_process` captures `self.diagnostics`; no `isinstance/diagnostics-not-None` guards remain in clip-with-log call sites.
2. **Appendix A completeness test** — every Appendix A name maps to exactly one Process's `REGISTRY_DIAGNOSTICS`.
3. **Final baseline parity run** — replays Phase 0 baseline; verifies bit-identical across all 11 state variables × 4,320 substeps × 5 cells.
4. **Diagnostics-subscription smoke test** — all Appendix A names pre-registered; state-variable subset still bit-identical to the no-subscription baseline.
5. **Perf check** — within 5% (no-subscription) / 15% (full-subscription) of the pre-refactor baseline.
6. **Legacy-inline-cleanup commit** — delete the `_change_legacy_inline` / `_rate_legacy_inline` shadow methods from every Process and delete the corresponding `test_*_helper_vs_inline.py` files. One final §11 contract run after the cleanup.
