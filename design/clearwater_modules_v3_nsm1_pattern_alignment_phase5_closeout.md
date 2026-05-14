# Phase 5 Close-Out — FloatingAlgae & BenthicAlgae `_change_with_components`

**Date:** 2026-05-13
**Branch:** `streaming`
**Baseline:** Phase 0 commit `d862d68`
**Phase 4 base:** commit `78134a1`
**Spec reference:** `design/clearwater_modules_v3_nsm1_pattern_alignment_specification.md` §6 Phase 5

Phase 5 applies the canonical refactor template to both algae Processes — the upstream producers of every `algal_*` and `balgae_*` cache that Carbon, DOX, Nitrogen, Phosphorus, and POM consume. **Includes the BenthicAlgae `rate_death` deduplication that Phase 1 deferred** (spec §6 Phase 5 explicit deliverable).

---

## What changed

### FloatingAlgae

**`REGISTRY_DIAGNOSTICS`** — 13 names per Appendix A diff §3:

- `algal_growth_rate` (preserved; consumed by Carbon, DOX, Nitrogen, Phosphorus)
- `algal_respiration_rate` (preserved; consumed by Carbon, DOX, Nitrogen)
- `algal_death_rate` (preserved)
- `algal_settling_rate` (preserved)
- `algal_orgn_from_mortality_rate` (preserved; consumed by Nitrogen)
- `algal_orgp_from_mortality_rate` (preserved; consumed by Phosphorus)
- `algal_poc_from_mortality_rate` (preserved; consumed by Carbon)
- `algal_doc_from_mortality_rate` (preserved; consumed by Carbon)
- `algal_pom_from_settling_rate` (preserved; consumed by POM)
- `algal_nh4_uptake_fraction` (preserved; consumed by Nitrogen, DOX)
- `algal_light_limitation` — **new** diagnostic
- `algal_nutrient_limitation_n` — **new** diagnostic
- `algal_nutrient_limitation_p` — **new** diagnostic

**Helpers:**
- `_change_with_components(...)` returns `(rate, components)`. Calls `rate()` and `_cache_mortality_rates()` to populate the existing side-effect caches, plus pure-function recomputes of `limit_light` / `limit_nitrogen` / `limit_phosphorus` for the three new diagnostic entries. Sets `self.algal_nh4_uptake_fraction` as a side effect (matches pre-Phase-5 `run`).
- `_change_legacy_inline(...)` shadow returns just `rate`; verbatim copy of the pre-Phase-5 `run` body. Deleted in Phase 10.

**`run` rewired** to read forcings → call helper → setattr loop → Forward Euler → unconditional clip → persist → opportunistic-write loop.

### BenthicAlgae

**`REGISTRY_DIAGNOSTICS`** — 11 names per Appendix A diff §3 (analogous to FloatingAlgae, minus `algal_settling_rate` and `algal_pom_from_settling_rate` which BenthicAlgae doesn't have):

- `balgae_growth_rate`, `balgae_respiration_rate`, `balgae_death_rate` (preserved sibling-consumer caches)
- `balgae_orgn_from_mortality_rate`, `balgae_orgp_from_mortality_rate`, `balgae_poc_from_mortality_rate`, `balgae_doc_from_mortality_rate` (preserved)
- `balgae_nh4_uptake_fraction` (preserved; consumed by Nitrogen)
- `balgae_light_limitation`, `balgae_nutrient_limitation_n`, `balgae_nutrient_limitation_p` (new)

**Note:** `balgae_pom_from_mortality_rate` is consumed by POM but is **not** in REGISTRY_DIAGNOSTICS (not in Appendix A's BenthicAlgae list). It is still set on `self` as a side effect for sibling reads but is not exposed via the opportunistic-write loop; `test_g_pom_routing_cache_remains_set` pins this contract.

**Helpers (with `rate_death` dedup):**
- `_change_with_components(...)` computes `ab_death = self.rate_death(algae, water_temperature)` **once** and reuses the cached value for both the rate composition (inlines `rate()`'s body using `ab_death` instead of re-invoking `rate_death`) AND the mortality routing (via the new private `_compute_balgae_mortality_components_from_death(ab_death, depth)` helper). Net rate is bit-identical to the pre-Phase-5 path because `rate_death` is pure.
- `_change_legacy_inline(...)` shadow invokes `rate()` (which calls `rate_death` once) and `_cache_benthic_mortality_rates(...)` (which calls `rate_death` again) — verbatim pre-Phase-5 behaviour. The two-call vs one-call dedup contract is pinned by:
  - `test_dedup_calls_rate_death_exactly_once` — helper invokes `rate_death` exactly once.
  - `test_shadow_calls_rate_death_twice` — shadow invokes `rate_death` exactly twice.

Both tests use `unittest.mock.patch.object(... wraps=...)` to count invocations without changing return values.

**`_compute_balgae_mortality_components_from_death(ab_death, depth)`** — new private helper. Verbatim copy of `_cache_benthic_mortality_rates`'s body MINUS the `self.rate_death(...)` call. Takes the pre-computed `ab_death` as an argument. The original `_cache_benthic_mortality_rates` is retained for back-compat with any external caller.

**`run` rewired** identically to FloatingAlgae's pattern.

---

## New tests

| File | Tests | Retained or deleted? |
|---|---|---|
| `tests/v3/nsm1/test_floating_algae_helper_vs_inline.py` | 17 (7 scenarios × parity + 7 components-set + zero-state finiteness + limit-bounds + preserved-cache parity) | Deleted in Phase 10 |
| `tests/v3/nsm1/test_floating_algae_registry_diagnostics.py` | 6 (pattern G full contract + preserved-names finiteness on sibling consumers) | Retained |
| `tests/v3/nsm1/test_benthic_algae_helper_vs_inline.py` | 18 (6 scenarios × parity + 6 components-set + dedup-call-count contracts (×2) + zero-state finiteness + limit-bounds + preserved-cache parity + pom-routing finiteness) | Deleted in Phase 10 |
| `tests/v3/nsm1/test_benthic_algae_registry_diagnostics.py` | 6 (pattern G full contract + pom-routing-not-in-registry-but-on-self) | Retained |

---

## §11 Zero-Regression Contract Verification

| Clause | Required | Achieved |
|---|---|---|
| §11.2 bit-identical state trajectory | `rtol=0, atol=0` against `baseline_coupled_trajectory_186b5c4.nc` | **OK** |
| §11.3 helper-vs-inline parity | bit-identical between helper and shadow across full input matrix | **OK** — 35/35 passed (17 FA + 18 BA) |
| §11.4 full test suite | 899 + 47 new = 946 passed, 2 xfailed | **OK** — 946 passed, 2 xfailed, 103.92 s wall |
| §11.4 Tier 1 conservation | `rtol=1e-12`, `clip_events == {}` (40 tests) | **OK** — 40 passed |

The bit-identical baseline parity confirms the **Phase 5 dedup is numerically perfect**: BenthicAlgae now invokes `rate_death` once per substep (down from twice) without changing any downstream state. Carbon (Phase 2), DOX (Phase 3), and Nitrogen (Phase 4) all read FloatingAlgae / BenthicAlgae sibling caches via `getattr` and continue to see identical values.

No regressions. No silent tolerance relaxation. No tests changed status.

---

## Files changed

- `src/clearwater_modules_v3/processes/floating_algae.py` — FloatingAlgae refactor.
- `src/clearwater_modules_v3/processes/benthic_algae.py` — BenthicAlgae refactor + dedup.
- `tests/v3/nsm1/test_floating_algae_helper_vs_inline.py` — new (17 tests).
- `tests/v3/nsm1/test_floating_algae_registry_diagnostics.py` — new (6 tests).
- `tests/v3/nsm1/test_benthic_algae_helper_vs_inline.py` — new (18 tests).
- `tests/v3/nsm1/test_benthic_algae_registry_diagnostics.py` — new (6 tests).

---

## What changes for Phase 6

Phase 6 is Phosphorus (`phosphorus.py`) — two state variables (TIP, OrgP), 7 Appendix A diagnostics. Reads `algal_growth_rate` / `algal_orgp_from_mortality_rate` from FloatingAlgae and the analogous BenthicAlgae caches, both of which Phase 5 has now formalised under their `REGISTRY_DIAGNOSTICS` tuples. Smaller refactor than Phase 5 (1 day).
