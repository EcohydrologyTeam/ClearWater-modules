# Phase 3 Close-Out — DOX `_change_with_components`

**Date:** 2026-05-13
**Branch:** `streaming`
**Baseline:** Phase 0 commit `d862d68` (NetCDF `baseline_coupled_trajectory_186b5c4.nc`)
**Phase 2 base:** commit `92cfa58`
**Spec reference:** `design/clearwater_modules_v3_nsm1_pattern_alignment_specification.md` §6 Phase 3

Phase 3 applies the Phase 2 canonical refactor template to DOX — the second-heaviest integrator (eight sub-fluxes) and the second of nine pattern-aligned Processes.

---

## What changed

### 1. `DOX.REGISTRY_DIAGNOSTICS` class attribute

Eleven names per the Appendix A diff §3:

- `dox_sat` (already cached pre-Phase 3)
- `atm_reaeration_rate` (already cached)
- `dox_nitrification_rate` (already cached)
- `dox_sod_rate` (already cached)
- `sod_rate` — **new**; an alias for `dox_sod_rate` exposing the same volumetric SOD sink under the 1.0.0 "sediment-globals" registry name
- `dox_doc_oxidation_rate` — **new**
- `dox_cbod_oxidation_rate` — **new**
- `dox_algal_photo_rate` — **new**
- `dox_algal_resp_rate` — **new**
- `dox_balgae_photo_rate` — **new**
- `dox_balgae_resp_rate` — **new**

### 2. `DOX._change_with_components(...)`

New canonical helper returning `(delta_dox, rate, components)`:

- `delta_dox` is the Forward Euler increment (mg-O2/L per substep).
- `rate` is the net rate (mg-O2/L/d) — caller writes it to `self.dox_rate` for downstream debugging.
- `components` is the `dict[str, ArrayLike]` indexed by `REGISTRY_DIAGNOSTICS`.

Operand order, per-sub-flux sanitisation, intermediate names, and the `is_user_*_zero` fast-path short-circuit are all preserved verbatim from the pre-refactor body per §11.6 rule 1 (code motion, not code rewrite).

### 3. `DOX._change_legacy_inline(...)`

Shadow returning just `(delta_dox, rate)`. Verbatim pre-Phase-3 body. Used only by `test_dox_helper_vs_inline.py`. Deleted in Phase 10.

### 4. `DOX.run` rewired

`run` now:
1. Reads state and optional forcings (unchanged).
2. Calls `self._change_with_components(...)` once.
3. Caches every `REGISTRY_DIAGNOSTICS` name on `self.<name>` via `setattr` from the components dict (pattern F). The four pre-existing cache attributes (`dox_sat`, `atm_reaeration_rate`, `dox_nitrification_rate`, `dox_sod_rate`) and the new seven all go through the same loop. `self.dox_rate` retains its non-Appendix-A meaning (the net rate).
4. Applies Forward Euler + unconditional `clip_negative_state` (no `_clip` wrapper in DOX — already harmonised in Phase 1; just absorbs the integrator step as `dox_new = dox + delta_dox`).
5. Persists `oxygen_dissolved` (pattern E).
6. Runs the opportunistic-write loop over `REGISTRY_DIAGNOSTICS` (pattern G).

### 5. `sod_rate` / `dox_sod_rate` alias contract

Both names map to the same sanitised volumetric SOD sink (mg-O2/L/d). The components dict has two keys pointing at one value. `tests/v3/nsm1/test_dox_helper_vs_inline.py::test_sod_and_dox_sod_rate_are_aliases` and `tests/v3/nsm1/test_dox_registry_diagnostics.py::test_g_sod_and_dox_sod_rate_match_in_registry` both pin this contract so a future refactor can't silently drop the alias.

---

## New tests

### `tests/v3/nsm1/test_dox_helper_vs_inline.py` (21 tests; deleted in Phase 10)

Parametrised matrix across nine scenarios — the Phase 2 seven plus two DOX-specific edges:

- `zero_state` — all state at 0 (DOX clamped to 1e-12 to avoid 0/0 in SOD Monod).
- `uniform_state` — typical mid-range values.
- `randomised_state` — fixed-seed sample.
- `cold_water` (0.5 °C) and `hot_water` (35 °C).
- `thin_depth` (0.05 m).
- **`hypoxic`** (DOX = 0.05 mg/L) — exercises the SOD-Monod attenuation regime.
- **`supersaturated`** (DOX = 15 mg/L > O2sat) — exercises atmospheric-reaeration sign reversal.
- `high_nh4` (NH4 = 5 mg/L) — exercises the nitrification sink coupling.

Per scenario, bit-identical assertions (`np.testing.assert_array_equal`, `rtol=0, atol=0`) on both `delta_dox` and `rate` between `_change_with_components` and `_change_legacy_inline`. Components-dict structure pinned to `REGISTRY_DIAGNOSTICS`. Plus three contract tests: SOD-attenuation monotonicity (hypoxic SOD ≤ normoxic SOD), zero-state finiteness, sod_rate / dox_sod_rate alias equality.

### `tests/v3/nsm1/test_dox_registry_diagnostics.py` (6 tests; retained)

Pattern G contract:

- `test_g_diagnostics_written_when_pre_registered` — all 11 names finite after 60 substeps.
- `test_g_diagnostics_skipped_when_not_pre_registered` — un-subscribed DOX trajectory reproducible.
- `test_g_state_bit_identical_with_and_without_diagnostics` — subscribing to diagnostics does not change `oxygen_dissolved`.
- `test_g_diagnostics_attribute_caches_match_registry_writes` — cached `self.<name>` matches registry value.
- `test_g_partial_subscription_writes_only_requested_names`.
- `test_g_sod_and_dox_sod_rate_match_in_registry` — alias contract enforced through the registry path too.

---

## §11 Zero-Regression Contract Verification

| Clause | Required | Achieved |
|---|---|---|
| §11.2 bit-identical state trajectory | `rtol=0, atol=0` against `baseline_coupled_trajectory_186b5c4.nc` | **OK** — `check_baseline_parity.py` reports `OK` |
| §11.3 helper-vs-inline parity | `_change_with_components` vs `_change_legacy_inline`, `rtol=0, atol=0`, 9-scenario matrix | **OK** — 21/21 passed |
| §11.4 full test suite | 839 + 27 new = 866 passed, 2 xfailed | **OK** — 866 passed, 2 xfailed, 83.25 s wall |
| §11.4 Tier 1 conservation | `rtol=1e-12`, `clip_events == {}` (40 tests) | **OK** — 40 passed |

No regressions. No silent tolerance relaxation. No tests changed status.

---

## Files changed

- `src/clearwater_modules_v3/processes/dox.py` — DOX refactor.
- `tests/v3/nsm1/test_dox_helper_vs_inline.py` — new (21 tests).
- `tests/v3/nsm1/test_dox_registry_diagnostics.py` — new (6 tests).

---

## What changes for Phase 4

Phase 4 is Nitrogen (`nitrogen.py`). Three state variables (NH4, NO3, OrgN), and Nitrogen owns the `nitrification_flux_rate` / `denitrification_flux_rate` caches that DOX and Alkalinity consume — the helper must preserve those instance-attribute names exactly so sibling consumers continue to read the same values. The parity matrix must cover `use_OrgN={True,False}` and the NH4+NO3 mass-balance closure regime. Phase 4 commit will trigger a re-run of the Phase 3 DOX baseline parity (DOX is the downstream consumer of Nitrogen's flux caches).
