# Phase 8 Close-Out — N2 & Pathogen

**Date:** 2026-05-14
**Branch:** `streaming`
**Baseline:** Phase 0 commit `d862d68`
**Phase 7 base:** commit `a15bc75`
**Spec reference:** `design/clearwater_modules_v3_nsm1_pattern_alignment_specification.md` §6 Phase 8

Phase 8 applies the canonical refactor template to N2 and Pathogen. Pathogen introduces the **`_rate_with_components` naming convention** per spec §10 Q5 (rate-form integrator).

---

## What changed

### N2

**`REGISTRY_DIAGNOSTICS`** — 4 names per Appendix A diff §3:

- `n2_atm_exchange_rate` (preserved attribute name)
- `n2_sat` (preserved attribute name)
- `total_dissolved_gas` (the only pre-Phase-8 example of pattern G in v3 NSM1 1.0.0; **extended** to cover the full Appendix A set)
- `n2_denit_source_rate` — new

**`self.tdg` ↔ `self.total_dissolved_gas` aliases**: the back-compat attribute name `self.tdg` (read by `test_n2_v1_parity_v3.py` and `test_n2_tier1.py`) is preserved alongside the Appendix A name `self.total_dissolved_gas`. Both point at the same value; pinned by `test_g_tdg_alias_to_total_dissolved_gas`.

**`_change_with_components(...)`** returns `(rate, components)`. The TDG diagnostic is computed *post-integrator-step* inside `run` because it depends on `n2_new`, not `n2_state`. The helper covers 3 of 4 REGISTRY_DIAGNOSTICS names; `run` adds `total_dissolved_gas` to the components dict after computing TDG, then runs the pattern F setattr loop and the pattern G opportunistic-write loop over the complete set.

**`_change_legacy_inline(...)`** shadow returns just the net rate. Verbatim copy of pre-Phase-8 `run` body (rate composition portion). Deleted in Phase 10.

### Pathogen — `_rate_with_components` (rate-form per §10 Q5)

**`REGISTRY_DIAGNOSTICS`** — 3 names per Appendix A diff §3:

- `pathogen_natural_death_rate` (new)
- `pathogen_light_death_rate` (new)
- `pathogen_settling_rate` (new)

All three are positive-magnitude sinks; the integrator applies them via `rate = -(natural + light + settling)`. Pinned by `test_components_are_positive_magnitudes`.

**`_rate_with_components(...)`** — first use of the rate-form naming convention. Returns `(rate, components)` where `rate` is the net per-day rate (cfu/100mL/d). The name distinguishes from `_change_with_components(...)` used by multi-state delta-form Processes; the convention is documented in spec §10 Q5.

**`_rate_legacy_inline(...)`** shadow returns just the net rate. Verbatim copy of pre-Phase-8 `run` body invoking `self.rate(...)` plus the `sanitize_rate` post-step. Deleted in Phase 10.

**Legacy `rate()` method** is retained for back-compat with external callers; `_rate_with_components` delegates to its three sub-helpers (`_rate_natural_decay`, `_rate_light_decay`, `_rate_settling`) directly without going through `rate()`.

### Phase 1 cleanup on Pathogen

The `clip_negative_state(..., step=0)` placeholder call at `pathogen.py:247` is removed (Pathogen was not in Phase 1's harmonisation list; Phase 8 picks it up).

---

## New tests

| File | Tests | Retained / deleted |
|---|---|---|
| `tests/v3/nsm1/test_n2_helper_vs_inline.py` | 14 (6 scenarios × bit-identical parity + 6 components-set + zero-state finiteness + supersaturated-sign monotonicity) | Deleted in Phase 10 |
| `tests/v3/nsm1/test_n2_registry_diagnostics.py` | 7 (pattern G full contract + pre-existing `total_dissolved_gas` path preserved + `self.tdg` / `self.total_dissolved_gas` alias contract) | Retained |
| `tests/v3/nsm1/test_pathogen_helper_vs_inline.py` | 14 (6 scenarios × bit-identical parity + 6 components-set + positive-magnitudes + zero-state finiteness + high-light monotonicity) | Deleted in Phase 10 |
| `tests/v3/nsm1/test_pathogen_registry_diagnostics.py` | 5 (pattern G full contract) | Retained |

---

## §11 Zero-Regression Contract Verification

| Clause | Required | Achieved |
|---|---|---|
| §11.2 bit-identical state trajectory | `rtol=0, atol=0` against `baseline_coupled_trajectory_186b5c4.nc` | **OK** |
| §11.3 helper-vs-inline parity | bit-identical across the input matrices | **OK** — 28/28 passed (14 N2 + 14 Pathogen) |
| §11.4 full test suite | 1,023 + 41 new = 1,064 passed, 2 xfailed | **OK** — 1,064 passed, 2 xfailed, 135.70 s wall |
| §11.4 Tier 1 conservation | `rtol=1e-12`, `clip_events == {}` (40 tests) | **OK** — 40 passed |

The bit-identical baseline parity is particularly meaningful for N2 because the pre-existing `total_dissolved_gas` opportunistic-write path had to be preserved exactly while the loop was extended to the full Appendix A set.

No regressions. No silent tolerance relaxation. No tests changed status.

---

## Files changed

- `src/clearwater_modules_v3/processes/n2.py` — N2 refactor.
- `src/clearwater_modules_v3/processes/pathogen.py` — Pathogen refactor (with `_rate_with_components` naming convention).
- `tests/v3/nsm1/test_n2_helper_vs_inline.py` — new (14 tests).
- `tests/v3/nsm1/test_n2_registry_diagnostics.py` — new (7 tests).
- `tests/v3/nsm1/test_pathogen_helper_vs_inline.py` — new (14 tests).
- `tests/v3/nsm1/test_pathogen_registry_diagnostics.py` — new (5 tests).

---

## What changes for Phase 9

Phase 9 is the last per-Process phase: **Alkalinity** (one state variable, 6 diagnostics). After Phase 9, all 9 pattern-aligned Processes (Carbon, DOX, Nitrogen, FloatingAlgae, BenthicAlgae, Phosphorus, POM, CBOD, N2, Pathogen, Alkalinity) carry the canonical structure. Phase 10 then runs the final end-to-end verification, the conformance test that scans every Process for the canonical shape, and the legacy-inline-shadow cleanup commit.
