# Phase 6 Close-Out — Phosphorus `_change_with_components`

**Date:** 2026-05-13
**Branch:** `streaming`
**Baseline:** Phase 0 commit `d862d68`
**Phase 5 base:** commit `7725f1a`
**Spec reference:** `design/clearwater_modules_v3_nsm1_pattern_alignment_specification.md` §6 Phase 6

Phase 6 applies the canonical refactor template to Phosphorus — two state variables (TIP, OrgP), 7 Appendix A diagnostics, Phase 1 `step=0` placeholders cleaned up.

---

## What changed

### `Phosphorus.REGISTRY_DIAGNOSTICS`

Seven names per Appendix A diff §3:

- `orgp_hydrolysis_rate` — **alias** for the existing `orgp_to_tip_hydrolysis_rate` attribute (which `test_phosphorus_v1_parity_v3.py` reads via `getattr`); both names are set to the same value as side effects in the helper.
- `orgp_settling_rate` (preserved attribute name).
- `tip_settling_rate` (preserved attribute name).
- `dip_from_bed` — new diagnostic; sediment P release.
- `orgp_algal_mortality_rate` — new; sum of floating + benthic algal mortality contributions to OrgP.
- `tip_algal_growth_rate` — new; floating-algae TIP uptake.
- `tip_balgae_growth_rate` — new; benthic-algae TIP uptake.

### `_change_with_components(...)`

Returns `(dtip_dt, dorgp_dt, components)`. Per-day rates (mg-P/L/d). Code-motion-only refactor of the pre-Phase-6 `run` body — operand order, intermediate names, kinetic-helper calls, and the `use_TIP` / `use_OrgP` gating preserved verbatim. Side effect: sets `self.orgp_to_tip_hydrolysis_rate`, `self.tip_settling_rate`, `self.orgp_settling_rate` (preserved attribute names that the Tier 1 tests, the v1 parity tests, and the v3 parity tests all depend on). The pattern F `setattr` loop in `run` is idempotent on these names and additionally writes the four new diagnostic names.

### `_change_legacy_inline(...)`

Shadow returning just `(dtip_dt, dorgp_dt)`. Verbatim copy of the pre-Phase-6 `run` body. Deleted in Phase 10.

### `Phosphorus.run` rewired

Reads → `_change_with_components` → setattr loop → Forward Euler → unconditional `clip_negative_state` (no `step=0` placeholders) → persist → opportunistic-write loop.

### Phase 1 cleanup

The two `clip_negative_state(..., step=0)` placeholder calls (`phosphorus.py:382` and `:385`) are removed. Step attribution is automatic via `diagnostics.current_step` (Phase 0.6 Q1). Phosphorus was not in Phase 1's harmonisation list; Phase 6 picks it up as part of the rewire.

---

## New tests

### `tests/v3/nsm1/test_phosphorus_helper_vs_inline.py` (24 tests; deleted in Phase 10)

Six scenarios — `zero_state`, `uniform_state`, `randomised_state`, `cold_water` (0.5 °C), `hot_water` (35 °C), `thin_depth` (0.05 m) — bit-identical between helper and shadow on `(dtip_dt, dorgp_dt)`, separately under `(use_TIP=True, use_OrgP=True)`, `(use_TIP=True, use_OrgP=False)`, and `(use_TIP=False, use_OrgP=True)` (18 parity tests). Plus 6 components-set pinning, 1 alias-equality test (`orgp_hydrolysis_rate == self.orgp_to_tip_hydrolysis_rate`), 1 zero-state finiteness.

### `tests/v3/nsm1/test_phosphorus_registry_diagnostics.py` (8 tests; retained)

Pattern G full contract: pre-registration writes, no-subscription reproducibility, zero-cost-when-unused state invariant, attribute / registry value parity, partial subscription, and the alias contract enforced through the registry path (`orgp_hydrolysis_rate` registry value matches `self.orgp_to_tip_hydrolysis_rate`).

---

## §11 Zero-Regression Contract Verification

| Clause | Required | Achieved |
|---|---|---|
| §11.2 bit-identical state trajectory | `rtol=0, atol=0` against `baseline_coupled_trajectory_186b5c4.nc` | **OK** |
| §11.3 helper-vs-inline parity | bit-identical between helper and shadow across the input matrix | **OK** — 24/24 passed |
| §11.4 full test suite | 946 + 32 new = 978 passed, 2 xfailed | **OK** — 978 passed, 2 xfailed, 110.25 s wall |
| §11.4 Tier 1 conservation | `rtol=1e-12`, `clip_events == {}` (40 tests) | **OK** — 40 passed |

No regressions. The alias contract for `orgp_hydrolysis_rate` ↔ `orgp_to_tip_hydrolysis_rate` keeps `test_phosphorus_v1_parity_v3.py` passing without modification.

---

## Files changed

- `src/clearwater_modules_v3/processes/phosphorus.py` — Phosphorus refactor.
- `tests/v3/nsm1/test_phosphorus_helper_vs_inline.py` — new (24 tests).
- `tests/v3/nsm1/test_phosphorus_registry_diagnostics.py` — new (8 tests).

---

## What changes for Phase 7

Phase 7 covers POM and CBOD together (~1.5 days). POM finally relocates the `pom_doc_source_rate` cache from `rate()` into `run` via `_change_with_components` (Phase 1's deferred item), per spec §6 Phase 7. CBOD multi-group diagnostics aggregate per-group oxidation/settling.
