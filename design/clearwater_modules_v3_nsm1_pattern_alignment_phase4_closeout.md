# Phase 4 Close-Out — Nitrogen `_change_with_components`

**Date:** 2026-05-13
**Branch:** `streaming`
**Baseline:** Phase 0 commit `d862d68`
**Phase 3 base:** commit `b0f2af6`
**Spec reference:** `design/clearwater_modules_v3_nsm1_pattern_alignment_specification.md` §6 Phase 4

Phase 4 applies the Phase 2 canonical refactor template to Nitrogen — three state variables (NH4, NO3, OrgN) and the upstream producer of two flux caches that DOX, Alkalinity, and N2 already consume.

---

## What changed

### 1. `Nitrogen.REGISTRY_DIAGNOSTICS` class attribute

Ten names per Appendix A diff §3:

- `nitrification_flux_rate` — **preserved** attribute name; DOX/Alkalinity/N2 already read this via `getattr`.
- `denitrification_flux_rate` — **preserved**; same consumers.
- `nh4_from_bed` — sediment NH4 release.
- `no3_from_bed_denit` — sediment NO3 denitrification sink.
- `orgn_hydrolysis_rate` — OrgN → NH4 hydrolysis.
- `orgn_settling_rate` — OrgN settling sink.
- `nh4_algal_growth_rate` — total algal NH4 uptake (floating + benthic sum).
- `no3_algal_growth_rate` — total algal NO3 uptake (floating + benthic sum).
- `nh4_algal_resp_rate` — floating-algae respiration NH4 source.
- `nh4_balgae_resp_rate` — benthic-algae respiration NH4 source.

The two `_flux_rate`-suffixed names are the renamed-from-1.0.0 entries that Appendix A documented as part of the migration. They were already published as `self.nitrification_flux_rate` / `self.denitrification_flux_rate` instance attributes pre-Phase-4, so the rename is documentation-only — no code consumer change.

### 2. `Nitrogen._change_with_components(...)`

New canonical helper returning `(ammonium_rate, nitrate_rate, orgn_rate, components)`. All rates are per-day (mg-N/L/d); `run` applies `dt_days` (matches the Carbon convention).

The helper:
- Computes `ammonium_nitrification` and `nitrate_denitrification` once each for the cache (= `nitrification_flux_rate`, `denitrification_flux_rate`).
- Calls `change_ammonium`, `change_nitrate`, `change_organic_nitrogen` exactly as `run` did pre-Phase-4 (these methods internally re-invoke the same helpers; the recomputation is bit-identical and matches the pre-refactor behaviour exactly).
- Computes the per-source diagnostic sub-fluxes (bed releases, algal couplings, OrgN hydrolysis/settling) for the components dict.

Code-motion-only per §11.6: operand order, intermediate names, and arithmetic preserved verbatim.

### 3. `Nitrogen._change_legacy_inline(...)`

Shadow returning just `(ammonium_rate, nitrate_rate, orgn_rate)`. Verbatim copy of the pre-Phase-4 `run` body's rate computation. Intentionally also calls `ammonium_nitrification` and `nitrate_denitrification` once each (discarding their values via `_ = ...`) so the call counts match the helper exactly. Used only by `test_nitrogen_helper_vs_inline.py`. Deleted in Phase 10.

### 4. `Nitrogen.run` rewired

`run` now reads → delegates to `_change_with_components` → caches every `REGISTRY_DIAGNOSTICS` name via `setattr` loop → applies Forward Euler + unconditional `clip_negative_state` for all three states → persists NH4 / NO3 / OrgN (with the existing `if "organic_nitrogen" in registry` guard for legacy v2 registries) → runs the opportunistic-write loop.

---

## New tests

### `tests/v3/nsm1/test_nitrogen_helper_vs_inline.py` (27 tests; deleted in Phase 10)

Eight scenarios — Phase 2/3 standards plus Nitrogen-specific edges:

- `zero_state`, `uniform_state`, `randomised_state` (fixed seed).
- `cold_water` (0.5 °C), `hot_water` (35 °C), `thin_depth` (0.05 m).
- **`high_nh4_low_no3`** — NH4+NO3 mass-balance closure regime (Phase 4 deliverable).
- **`hypoxic`** (DOX = 0.05 mg/L) — nitrification inhibition + denitrification enhancement regime.

For each scenario, bit-identical assertions on `(ammonium_rate, nitrate_rate, orgn_rate)` between `_change_with_components` and `_change_legacy_inline`, **separately under `use_OrgN ∈ {True, False}`** (16 parity tests). Plus 8 components-dict structure tests, 1 preserved-attribute-name finiteness test, 1 high-NH4 nitrification-monotonicity test, 1 zero-state finiteness test.

### `tests/v3/nsm1/test_nitrogen_registry_diagnostics.py` (6 tests; retained)

Pattern G contract:

- `test_g_diagnostics_written_when_pre_registered` — all 10 names finite after 60 substeps.
- `test_g_diagnostics_skipped_when_not_pre_registered` — un-subscribed NH4 / NO3 / OrgN reproducible.
- `test_g_state_bit_identical_with_and_without_diagnostics` — subscription does NOT change NH4 / NO3 / OrgN trajectories.
- `test_g_diagnostics_attribute_caches_match_registry_writes` — cached `self.<name>` matches `registry.get(name)`.
- `test_g_partial_subscription_writes_only_requested_names`.
- `test_g_preserved_names_are_consumer_visible` — `nitrification_flux_rate` / `denitrification_flux_rate` populated on `self` after `run`, finite, non-negative (the F contract for sibling consumers).

---

## §11 Zero-Regression Contract Verification

| Clause | Required | Achieved |
|---|---|---|
| §11.2 bit-identical state trajectory | `rtol=0, atol=0` against `baseline_coupled_trajectory_186b5c4.nc` | **OK** |
| §11.3 helper-vs-inline parity | `_change_with_components` vs `_change_legacy_inline`, `rtol=0, atol=0`, 8-scenario × 2-OrgN-mode matrix | **OK** — 27/27 passed |
| §11.4 full test suite | 866 + 33 new = 899 passed, 2 xfailed | **OK** — 899 passed, 2 xfailed, 99.17 s wall |
| §11.4 Tier 1 conservation | `rtol=1e-12`, `clip_events == {}` (40 tests) | **OK** — 40 passed |

DOX (Phase 3) and Alkalinity / N2 (still pre-refactor) all read Nitrogen's preserved-name caches via `getattr` — the bit-identical baseline parity confirms those sibling reads still resolve to the same values.

---

## Files changed

- `src/clearwater_modules_v3/processes/nitrogen.py` — Nitrogen refactor.
- `tests/v3/nsm1/test_nitrogen_helper_vs_inline.py` — new (27 tests).
- `tests/v3/nsm1/test_nitrogen_registry_diagnostics.py` — new (6 tests).

---

## What changes for Phase 5

Phase 5 covers FloatingAlgae and BenthicAlgae — the upstream producers of the eight `algal_*_from_mortality_rate` caches that Nitrogen / Phosphorus / Carbon / POM consume. Phase 5 also consolidates the BenthicAlgae `_cache_benthic_mortality_rates` into the components-dict path, eliminating the redundant `rate_death` re-invocation that Phase 1 deferred to this phase. Phase 5 commit will trigger re-runs of the Phase 2 / 3 / 4 baseline parity (Carbon, DOX, Nitrogen all read algal caches).
