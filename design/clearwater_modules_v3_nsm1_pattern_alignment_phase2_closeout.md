# Phase 2 Close-Out — Carbon `_change_with_components`

**Date:** 2026-05-13
**Branch:** `streaming`
**Baseline:** Phase 0 commit `d862d68` (NetCDF `baseline_coupled_trajectory_186b5c4.nc`)
**Phase 1 base:** commit `481532f`
**Spec reference:** `design/clearwater_modules_v3_nsm1_pattern_alignment_specification.md` §6 Phase 2

Phase 2 lands the canonical pattern-alignment refactor on Carbon — the first Process to receive a fused `_change_with_components` helper, a class-level `REGISTRY_DIAGNOSTICS` tuple, and the opportunistic-write loop. Sets the template for Phases 3–9.

---

## What changed

### 1. `Carbon.REGISTRY_DIAGNOSTICS` class attribute

Nine names — the full Carbon-owned subset of the Appendix A diff catalog (§3 of `clearwater_modules_v3_nsm1_appendix_a_diff.md`):

- `poc_hydrolysis_rate`
- `doc_dic_oxidation_rate`
- `dic_atm_exchange_rate` (new — `co2_reaeration`)
- `dic_sed_release_rate` (new — was the `dic_from_bed` sediment-global; now owned by Carbon)
- `carbon_algal_resp_rate` (new)
- `carbon_balgae_resp_rate` (new)
- `carbon_algal_photo_rate` (new)
- `carbon_balgae_photo_rate` (new)
- `carbon_cbod_oxidation_rate` (new)

### 2. `Carbon._change_with_components(...)`

New canonical helper. Returns `(d_poc, d_doc, d_dic, components)` where `components` is a `dict[str, ArrayLike]` containing every Appendix A diagnostic. Computes every sub-rate exactly once and is the single source of truth for both the integrator step and the diagnostic exposure.

Operand order, intermediate variable names, and arithmetic are preserved verbatim from the pre-refactor `run` body per §11.6 rule 1 (code motion, not code rewrite).

### 3. `Carbon._change_legacy_inline(...)`

Shadow method holding the pre-Phase-2 inline composition. Returns only `(d_poc, d_doc, d_dic)`. Used exclusively by `test_carbon_helper_vs_inline.py` to prove that `_change_with_components` produces bit-identical deltas. Deleted in Phase 10.

### 4. `Carbon.run` rewired

`run` now:
1. Reads state and forcings from the registry (unchanged).
2. Calls `self._change_with_components(...)` once.
3. Caches each `REGISTRY_DIAGNOSTICS` name on `self.<name>` via `setattr` from the components dict (pattern F).
4. Applies Forward Euler + unconditional `clip_negative_state` (pattern C+D — replaces the `_clip` wrapper, which is deleted).
5. Persists POC/DOC/DIC to the registry (pattern E).
6. Iterates `REGISTRY_DIAGNOSTICS` and writes each pre-registered name to the registry (pattern G).

### 5. `Carbon._clip` removed

The Phase-0.6-Q2-superseded wrapper (`isinstance/diagnostics-not-None` guard with `xr.where` fallback) is deleted. The three former call sites use `clip_negative_state(state, name, self.diagnostics)` directly. Step attribution is automatic via `diagnostics.current_step` (Phase 0.6 Q1).

---

## New tests

### `tests/v3/nsm1/test_carbon_helper_vs_inline.py` (15 tests; deleted in Phase 10)

Parametrised matrix of state/forcing scenarios:

- `zero_state` — all state at 0.
- `uniform_state` — typical mid-range values.
- `randomised_state` — fixed-seed (20260513) sample of realistic ranges.
- `cold_water` (0.5 °C) and `hot_water` (35 °C).
- `thin_depth` (0.05 m) — edge case for the Arrhenius / depth divisions.
- `low_dox` (0.1 mg/L) — DOX-Monod attenuation regime.

Each scenario is asserted bit-identical (`np.testing.assert_array_equal`, `rtol=0, atol=0`) between `_change_with_components` and `_change_legacy_inline`. The components dict structure (key set == `REGISTRY_DIAGNOSTICS`) is also pinned.

### `tests/v3/nsm1/test_carbon_registry_diagnostics.py` (5 tests; retained through Phase 10)

Pattern G contract:

- `test_g_diagnostics_written_when_pre_registered` — pre-register every name, run for 60 substeps, assert finite per-cell values are written each step.
- `test_g_diagnostics_skipped_when_not_pre_registered` — un-subscribed run is bit-reproducible.
- `test_g_state_bit_identical_with_and_without_diagnostics` — the zero-cost-when-unused invariant: subscribing to diagnostic outputs does NOT change the POC/DOC/DIC trajectory.
- `test_g_diagnostics_attribute_caches_match_registry_writes` — the value at `self.<name>` matches `registry.get(name)` after each substep (single source of truth between F and G).
- `test_g_partial_subscription_writes_only_requested_names` — subscribing to a subset of names does not silently add the rest to the registry.

---

## §11 Zero-Regression Contract Verification

| Clause | Required | Achieved |
|---|---|---|
| §11.2 bit-identical state trajectory | `xr.testing.assert_identical`, `rtol=0, atol=0` against `baseline_coupled_trajectory_186b5c4.nc` | **OK** — `check_baseline_parity.py` reports `OK: bit-identical parity` |
| §11.3 helper-vs-inline parity | `_change_with_components` vs `_change_legacy_inline`, `rtol=0, atol=0`, parametrised matrix | **OK** — 15/15 passed (7 bit-identical scenarios + 7 component-set scenarios + 1 zero-state finiteness) |
| §11.4 full test suite | 819 + 20 new = 839 passed, 2 xfailed | **OK** — 839 passed, 2 xfailed, 78.20 s wall (`pytest tests/`) |
| §11.4 Tier 1 conservation | `rtol=1e-12`, `clip_events == {}` (40 tests) | **OK** — 40 passed |

No regressions. No silent tolerance relaxation. No tests changed status.

---

## Files changed

- `src/clearwater_modules_v3/processes/carbon.py` — Carbon refactor (one file, +~330 / -~220 lines net; same arithmetic, restructured into helpers).
- `tests/v3/nsm1/test_carbon_helper_vs_inline.py` — new (15 tests).
- `tests/v3/nsm1/test_carbon_registry_diagnostics.py` — new (5 tests).

---

## Template for Phases 3–9

Phase 2 establishes the canonical refactor template. Each subsequent Process follows the same structure:

1. Add `<Process>.REGISTRY_DIAGNOSTICS: tuple[str, ...]` to the class body, sourced from `clearwater_modules_v3_nsm1_appendix_a_diff.md` §3.
2. Move the existing inline rate composition from `run` into a new `_change_with_components(...)` helper that returns `(delta_state(s)..., components: dict)`. Populate the components dict from the same intermediates the integrator consumes. No algebraic rearrangement.
3. Create the shadow `_change_legacy_inline(...)` returning just the deltas — verbatim copy of the pre-refactor body, used only by `test_<process>_helper_vs_inline.py`.
4. Rewire `run` to call `_change_with_components`, cache via `setattr(self, name, components[name])`, integrate, clip, persist, then run the opportunistic-write loop.
5. Add `tests/v3/nsm1/test_<process>_helper_vs_inline.py` (deleted in Phase 10).
6. Add `tests/v3/nsm1/test_<process>_registry_diagnostics.py` (retained).
7. Verify all §11 gates: full suite, bit-identical baseline parity, helper-vs-inline parity, Tier 1.

The naming conventions, the test-file templates, the docstring shape, and the integration with `Phase 0.6` (no `_clip` wrappers, no `step=` placeholders, no `isinstance/diagnostics-not-None` guards) all carry forward verbatim.

---

## What changes for Phase 3

DOX (`dox.py`) is next. DOX has eight sub-fluxes (atmospheric reaeration, nitrification sink, SOD sink, DOC oxidation sink, CBOD oxidation sink, four algal photo/resp terms) and is the second-heaviest integrator after Carbon. Its `REGISTRY_DIAGNOSTICS` tuple covers 11 names per the Appendix A diff.

The `_change_with_components` helper-vs-inline parity matrix for DOX must additionally cover the hypoxic regime (`DOX → 0`) where SOD attenuation kicks in, and the DOX-saturated regime where atmospheric reaeration reverses sign. The §11.6 code-motion-only rules apply identically.
