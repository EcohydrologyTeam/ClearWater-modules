# Phase 1 Close-Out — Mechanical Alignment Pass

**Date:** 2026-05-13
**Branch:** `streaming`
**Baseline:** Phase 0 commit `d862d68` (NetCDF `baseline_coupled_trajectory_186b5c4.nc`)
**Spec reference:** `design/clearwater_modules_v3_nsm1_pattern_alignment_specification.md` §6 Phase 1

Phase 1 lands the auxiliary-convention harmonisation across all 11 v3 NSM1 Processes, without touching kinetics. Every change is code-motion-only per §11.6 and is verified against the Phase 0 baseline NetCDF bit-identically.

---

## Changes by Process

| Process | Phase 1 change(s) |
|---|---|
| **Nitrogen** (`nitrogen.py`) | Lifted `clip_negative_state` import to module level; removed `_clip` helper (`isinstance/diagnostics-not-None` wrapper); replaced all three call sites with direct `clip_negative_state(state, name, self.diagnostics)`. |
| **Alkalinity** (`alkalinity.py`) | Removed the `isinstance(... DataArray) and self.diagnostics is not None` guard branch around `clip_negative_state`; replaced with unconditional call. |
| **DOX** (`dox.py`) | Same as Alkalinity. |
| **N2** (`n2.py`) | Same as Alkalinity. |
| **FloatingAlgae** (`floating_algae.py`) | Lifted `clip_negative_state` import to module level; replaced `step=0` placeholder with no-arg form (step comes from `diagnostics.current_step` per Phase 0.6 Q1). |
| **POM** (`pom.py`) | Lifted `clip_negative_state` import to module level; removed `step=0` placeholder. Cache-relocation deferred to Phase 7 — see "Deferred" below. |
| **BenthicAlgae** (`benthic_algae.py`) | Replaced hardcoded `self.use_nitrate / use_ammonium / use_phosphate = True` in `init_process` with `model.has_process("<Class>")` sibling discovery (mirrors `FloatingAlgae.init_process`). Falls back to the previous defaults when `model` does not expose `has_process` (legacy fixture path). Lifted `clip_negative_state` import; removed `step=0` placeholder. `rate_death` de-duplication deferred to Phase 5 — see "Deferred" below. |
| **CBOD** (`cbod.py`) | Replaced `try/except KeyError` DOX-registry read with `if "oxygen_dissolved" in registry:` (matches FloatingAlgae / DOX / Pathogen idiom). Replaced two inline `xr.where(isnull, 0, ...)` / `np.where(np.isnan, 0, ...)` NaN scrubs with the canonical `sanitize_rate(...)` helper. Removed `step=0` placeholder. |
| **Carbon** | No Phase 1 change (already conformant for Phase 1 items; refactor lands in Phase 2). |
| **Phosphorus** | No Phase 1 change (already conformant; refactor lands in Phase 6). |
| **Pathogen** | No Phase 1 change (already conformant). |

**Aggregate change footprint:** 8 process files modified, 87 insertions, 67 deletions (`git diff --stat`).

---

## Deferred items

Two items in the spec §6 Phase 1 list were re-scoped after analysis:

### 1. POM cache relocation (deferred to Phase 7)

**Spec text:** "Move POM cache write from `rate()` to `run()`."

**What's actually there:** `POM.rate()` sets `self.pom_doc_source_rate` as a side effect at `pom.py:297-299`. The cache attribute is read by `Carbon.run` for DOC sourcing.

**Why deferred:** The cache assignment uses an intermediate value (`fcom * rate_dissolution * h2 / depth`) that is computed inside `rate()` from other locals. Moving the assignment to `run()` requires either (a) re-computing the same intermediate after `rate()` returns — which violates §11.6 rule 1 (code motion not code rewrite) — or (b) changing `rate()` to return a tuple `(rate, pom_doc_source_rate)` — which changes the public-ish `rate()` signature used by 3 external tests (`tests/v3/nsm1/test_pom_v1_parity_v3.py`).

**Resolution:** Phase 7 introduces `POM._change_with_components` returning `(d_pom, components)` and the cache assignment moves to `run()` as a natural part of that refactor. Spec §6 Phase 7 already anticipates this: "Phase 1 already moved it, but Phase 7 wraps the consolidated rate-and-cache flow into `_change_with_components`." The spec text "Phase 1 already moved it" is updated here to reflect that the relocation properly belongs to Phase 7.

### 2. BenthicAlgae redundant `rate_death` de-duplication (deferred to Phase 5)

**Spec text (§6 Phase 1):** "Fix BenthicAlgae redundant `rate_death` re-invocation."

**What's actually there:** `BenthicAlgae.rate(...)` calls `self.rate_death(algae, water_temperature)` inline (`benthic_algae.py:433`), and `_cache_benthic_mortality_rates(...)` calls it again (`benthic_algae.py:355`). `FloatingAlgae` has the analogous duplicate (lines 459 and 509).

**Why deferred:** De-duplicating without changing `rate()`'s public signature requires introducing an `_ab_death` cache attribute that both helpers read with a `getattr` fallback. Either approach is a code rewrite, not code motion, and the two calls *do* produce identical values today (deterministic, side-effect-free arithmetic) — so the duplicate work is a correctness no-op and the rewrite carries risk without functional benefit.

**Resolution:** Phase 5 consolidates `_cache_benthic_mortality_rates` into `_change_with_components` as a single-call rate-and-cache flow; the duplicate disappears there. Spec §6 Phase 5 already names this: "BenthicAlgae: consolidate the `_cache_benthic_mortality_rates` helper into the components-dict path (one call, no duplicate `rate_death`)." Same for FloatingAlgae in the same phase.

### 3. CBOD `init_process` sibling discovery (no-op; not applicable)

**Spec text (§5 table):** "J — **add** sibling discovery for DOX (currently `try/except`)"

**What's actually there:** CBOD reads DOX from the registry, not from a sibling-process cache. The `try/except` was a registry-membership check expressed in an off-pattern idiom; replacing it with the `if "oxygen_dissolved" in registry:` form (already done above) is the actual harmonisation. There is no `BenthicAlgae`/`DOX`/`Nitrogen`-style "consume from sibling cache" pattern for CBOD because no sibling Process produces a CBOD-relevant cache.

**Resolution:** No `init_process` change needed. The §5 table conflated `try/except` removal with sibling-process discovery; the closeout disambiguates.

---

## §11 Zero-Regression Contract Verification

| Clause | Required | Achieved |
|---|---|---|
| §11.2 bit-identical state trajectory | `xr.testing.assert_identical`, `rtol=0, atol=0` against `baseline_coupled_trajectory_186b5c4.nc` | **OK** — `check_baseline_parity.py` reports `OK: bit-identical parity` |
| §11.4 full test suite | 805 passed + 14 Phase 0.6 = 819 passed, 2 xfailed, no new skips | **OK** — 819 passed, 2 xfailed, 71.55 s wall (`pytest tests/`) |
| §11.4 Tier 1 conservation | `rtol=1e-12`, `clip_events == {}` (40 tests) | **OK** — 40 passed |

No regressions. No silent tolerance relaxation. The only behavior change vs the Phase 0 baseline is structural (guard branches removed, NaN scrub style harmonised, `init_process` sibling discovery on BenthicAlgae); state trajectories are byte-identical.

---

## Notes on the CBOD `sanitize_rate` adoption

The CBOD inline scrub previously caught only `NaN`. The canonical `sanitize_rate(...)` catches both `NaN` and `inf`. The change therefore widens the scrub to also catch `inf` — a behavior change in principle.

**Why bit-identicality still holds:** the 4,320-step baseline run does not produce any `inf` in `cbod_oxidation_rate` or `cbod_settling_rate`. The widening is a no-op on the baseline path. If a future scenario does generate `inf` in these rates, the new code zeroes them at CBOD instead of letting them propagate to DOX's `sanitize_rate(rate)` call site, which is the safer (and intended) behavior. This was verified empirically: the parity check passed unchanged after the swap.

---

## What changes for Phase 2

Carbon is the first Process to receive the `_change_with_components` refactor. Phase 2 introduces:

- `Carbon._change_with_components(...)` returning `(d_poc, d_doc, d_dic, components: dict)`.
- `Carbon.REGISTRY_DIAGNOSTICS: tuple[str, ...]` listing the Carbon-owned Appendix A names.
- `Carbon._change_legacy_inline(...)` — a shadow method holding the pre-refactor inline composition verbatim, retained through Phase 10 for the helper-vs-inline parity test.
- `tests/v3/nsm1/test_carbon_helper_vs_inline.py` (new) — bit-identical parity matrix.
- `tests/v3/nsm1/test_carbon_registry_diagnostics.py` (new) — opportunistic-write semantics.

The §11 test gate becomes the canonical template for Phases 3–9.
