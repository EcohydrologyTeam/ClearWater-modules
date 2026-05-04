# ClearWater Modules v3 — TSM Design Specification

**Status:** Approved for implementation
**Author:** Todd Steissberg (ERDC), with Claude
**Date:** 2026-05-04 (revised 2026-05-04 with all open questions resolved)
**Scope:** TSM (Temperature Simulation Module) within `clearwater_modules_v3`.

**Read this with the umbrella spec.** Quick start, env setup, branch conventions, package architecture, integrator-pattern contract, retirement plan, and umbrella risks are documented in `clearwater_modules_v3_architecture_specification.md`. This document covers TSM-specific design only.

---

## 1. TSM Background

### v1 TSM (`clearwater_modules.tsm`)

- Function-style framework with ~150 dynamic variables registered through `dynamic_variables.py` and process functions in `processes.py` (~520 lines)
- Hardened by the author in April–May 2026 with:
  - Latent-heat unit fix (Celsius polynomial vs Kelvin input) — commit `d9505c6`
  - Thin-water stability guard (depth ramp + `dTdt_max_per_hour` rate cap) — commit `d9505c6`
  - Hotstart from xr.Dataset checkpoint — commit `1a226dd`
  - Multi-cell-safe debug print removal — commit `4a95da4`
  - 418× kernel optimization (orchestration layer, benefits TSM) — commit `6daa65e`
  - Wet-mask gating (orchestration layer) — commit `3d18965`
- Comprehensive test coverage: 15 calculation tests + 4 latent-heat tests + 6 stability-ramp tests + hotstart roundtrip tests

### v2 TSM (`clearwater_modules_v2.processes.temperature.Temperature`)

- Class-based framework with `Process` subclass, `ProcessFactory.register`, per-process `time_step`, `registry.get_at_time`/`set_at_time`
- `Temperature` class is 645 lines as a single class
- LimnoTech contributions:
  - `mixing_ratio_air` divide-by-zero correction (Anthony, April 2026)
  - "Skip first time step" coupling logic (Paul, March 2026)
  - Sediment-temperature optional flag (Anthony, February 2026)
- Verified working end-to-end: coupled TSM+Riverine on Sumwere Creek runs in 89s for 4,320 timesteps as of 2026-05-04
- Missing relative to v1: latent-heat unit fix, thin-water stability guard, hotstart, wet-mask, debug print removal (some commented out, some still live)

### What v3 TSM accomplishes

A single coherent TSM implementation that combines v2's class-based framework with all of v1's recent corrections and capabilities. v3 TSM uses the v3 `Model` (developed in this work) which inherits the kernel optimization, wet-mask, and hotstart from v1.

---

## 2. TSM Goals and Non-Goals

### Goals

1. v3 TSM matches or exceeds v1 TSM's numerical correctness on the regression test suite (latent-heat fix applied; thin-water stability guard applied).
2. v3 TSM uses the v2 class-based framework (Process subclass, YAML config, per-process substepping).
3. v3 TSM inherits hotstart support from the v3 `Model`, including v3's per-process opt-in mechanism for substep-internal state preservation.
4. End-to-end coupled TSM+Riverine demo runs cleanly on v3 with no notebook code changes beyond import statements.
5. Test coverage at least matching v1's TSM regression suite (25 tests + hotstart roundtrip).
6. Performance equal to or better than v2's verified Sumwere Creek baseline (89 s for 4,320 timesteps).
7. Numerical outputs match v1's corrected TSM within floating-point tolerance for the regression test suite.

### Non-Goals

1. v3 TSM will not introduce new heat-flux components (ice cover, riparian shading, atmospheric coupling beyond met forcing remain out of scope; v1 doesn't have them either).
2. v3 TSM will not change the meteorological forcing schema. The v2 YAML data-source pattern is preserved.

---

## 3. Component Inventory

### 3.1 v3 TSM source: `src/clearwater_modules_v3/processes/temperature.py`

A full file, not an overlay (because TSM is the first merged module). It contains the v2 `Temperature` class structure with the following modifications applied:

| Modification | Source | Action |
|---|---|---|
| Latent-heat unit fix | v1 commit `d9505c6` (`mf_latent_heat_vaporization`) | Apply to `latent_heat_vaporization`: convert Kelvin to Celsius before applying the polynomial coefficients |
| Thin-water stability guard | v1 commit `d9505c6` (`dTdt_water_c`) | Apply to `temperature_change`: depth ramp on net flux + `dTdt_max_per_hour` rate cap, with the same parameter defaults that disable both effects when set to 0 / +inf |
| Debug print removal | v1 commit `4a95da4` and v2 commented-out prints | Delete the commented-out print blocks throughout the file (do not just leave them commented) |
| `mixing_ratio_air` divide-by-zero guard | v2 `memory-refactor-pytestUpdate` | Carry forward as-is (already correct on v2 latest) |
| Sediment-temperature optional flag | v2 `memory-refactor-pytestUpdate` | Carry forward as-is |
| "Skip first time step" logic | v2 `memory-refactor-pytestUpdate` | Carry forward as-is, with one-line documentation comment explaining the v1-coupling-compat reason |
| Richardson `-1` factor — **Resolved 2026-05-04: REMOVE** | v2 `memory-refactor-pytestUpdate` | Delete the commented-out `-1` line and the TODO comment at `temperature.py:597`. Rationale: Jason Rutyna investigated the v1↔v2 Richardson_number diff in commits `8218962` (2026-01-23) and `7f4166a` (2026-01-27) and concluded `-1` should not be there; v1 `tsm/processes.py:150-179` does not have it. v3 formula matches v1: `gravity * (density_air - density_air_sat) * 2.0 / (density_air * wind_speed**2.0)`. |
| `flux_sediment / 0.5` | v2 `memory-refactor-pytestUpdate` | Resolve TODO: document as the sediment active-layer half-thickness convention. Matches v1 legacy code at `processes.py:411`. Replace TODO comment with one-line docstring explanation; do not change the formula. |

### 3.2 v3 Model: `src/clearwater_modules_v3/model.py`

The v3 `Model` class extends v2's `Model` with three orchestration-level capabilities ported from v1. **These improvements are part of the v3 framework, not TSM-specific**, but TSM v3 is the first consumer and drives their initial design.

| Capability | v1 source | v3 design |
|---|---|---|
| 418× kernel optimization | v1 commit `6daa65e` (`_iter_computations_fast`, cached compute plan) | Adapt to v2's `Process.run`-based dispatch: cache the (process, time_step, variables) tuples once per Model instance; in the per-step loop, dispatch directly without re-checking `current_time_seconds % process.time_step_seconds` (precompute next-fire time per process). **Design note (to be flagged in v3 code at the kernel-optimization site):** this adapts v1's variable-level cached-plan optimization to v2's process-level dispatch. v2's process-level loop is coarser-grained than v1's variable-level loop, so the speedup may not match v1's 418× directly. **Follow-up optimization recommended after Phase 3:** profile against v2 baseline; if v3 is materially slower than v1, add a second optimization layer at the variable level inside `Process.run` (e.g., cache the registry-name lookups, batch `get_at_time`/`set_at_time` calls per process, or build a per-Process compute plan that bypasses the registry for hot variables). The recommendation is documented at the relevant code site so future maintainers see the design intent. |
| Wet-mask gating | v1 commit `3d18965` (`wet_mask` kwarg) | Add a `wet_mask` registry variable that processes can opt into; `Model.run` honors the mask by skipping `set_at_time` for masked cells. Existing per-process mask logic (e.g., the `xr.where(volume > 0, ...)` in v2 `Temperature.run`) becomes redundant once the registry-level mask is in place; remove from process bodies |
| Hotstart from xr.Dataset | v1 commit `1a226dd` (`hotstart_dataset`, `hotstart_timestep`); v2 has no hotstart | Add `hotstart_dataset` and `hotstart_timestep` kwargs to `init_from_file`. `Model.__init__` honors them by seeding the registry from the dataset at `hotstart_timestep`. **Resolved 2026-05-04: substep state semantics** — hotstart preserves dataset-level (registry) state at the saved time. Per-process substep-internal state (e.g., v2 `Temperature.__skip_first_time_step`) defaults to "fresh start" semantics after hotstart. Each `Process` may optionally implement `to_hotstart() -> dict` and `from_hotstart(state: dict)` methods to preserve richer internal state; default is no-op. v3 `Temperature.from_hotstart` sets `__skip_first_time_step = False` (don't skip the first post-hotstart step — you're not starting from scratch). |
| Chunking execution — **Resolved 2026-05-04** | v2 `Model.__process_loop_chunked` (Paul's commit `d712c59`, 2026-03-11) | Inherit v2's `__process_loop_chunked` skeleton in v3. **Resolve the four TODOs Paul left in the method:** (1) "this need actual chunking logic", (2) "look at riverine's code and mirror where applicable", (3) "align with riverine", (4) "confirm if this is necessary to write out the last chunk". Resolution approach: study `clearwater_riverine.transport.py` and `clearwater_riverine.constituents.py` chunking patterns (since the author's streaming/chunking work lives in the riverine repo, not in `clearwater_modules` — searching the legacy modules package for "chunk" returns only one comment in `nsm1/model.py:51`). Mirror riverine's chunk-boundary handling, last-chunk write semantics, and dataset alignment. |

### 3.3 v3 Configuration: `src/clearwater_modules_v3/config/init.py`

`init_from_file` accepts the same YAML schema v2 accepts, plus three additional optional top-level keys:

```yaml
hotstart:                     # optional
  dataset_path: hotstart.nc
  timestep: '2022-05-13 12:00:00'

wet_mask:                     # optional
  variable: wetted_surface_area_threshold
  threshold: 1.0              # m^2; cells with surface_area below this are considered dry

# everything else identical to v2
```

If neither `hotstart` nor `wet_mask` keys are present, v3 behavior matches v2 exactly. This preserves backward compatibility with all existing v2 configs.

### 3.4 Migration table (v2 → v3 for TSM)

| v2 import | v3 equivalent |
|---|---|
| `import clearwater_modules_v2 as cwm` | `import clearwater_modules_v3 as cwm` |
| `from clearwater_modules_v2.config import init_from_file` | `from clearwater_modules_v3.config import init_from_file` |
| `from clearwater_modules_v2.processes.temperature import Temperature` | `from clearwater_modules_v3.processes.temperature import Temperature` |

A migration table will be included in v3's README.

---

## 4. Testing and Validation

### Test infrastructure

v3 inherits the test directory structure from v2 (`tests/`) and adds:

- `tests/v3/test_5_tsm_calculations_v3.py` — v3 equivalent of v1's `test_5_tsm_calculations.py`, using the same expected-value tables in `tests/NSM Manual Calcs/` but exercising the v3 `Temperature.run` flow
- `tests/v3/test_tsm_latent_heat_v3.py` — port of v1's `test_tsm_latent_heat.py`
- `tests/v3/test_tsm_stability_ramp_v3.py` — port of v1's `test_tsm_stability_ramp.py`
- `tests/v3/test_hotstart_roundtrip_v3.py` — port of v1's `test_hotstart_roundtrip.py`, adapted for v3's hotstart API
- `tests/v3/test_v2_v3_parity.py` — A/B comparison: load Sumwere Creek modules.yml, run v2 and v3 in the same Python session, assert that water temperature timeseries match within tolerance after applying the v3 corrections to a v2-equivalent control case (i.e., when v3's latent-heat fix and thin-water guard are disabled, outputs should be identical to v2's outputs)
- `tests/v3/test_coupled_tsm_riverine_v3.py` — programmatic version of `examples/03_Example_Coupled_TSM_and _Riverine.ipynb` with v3, asserting completion within a wall-time budget

### Validation tiers (per `Validation_strategy_no_reference.md`)

- **Tier 1 (conservation):** mass and energy conservation tests for TSM (closed-system test with no fluxes; total enthalpy must be conserved)
- **Tier 2 (analytical limits):** verification of individual heat-flux components against analytical Brutsaert and Brunt formulations
- **Tier 3 (steady-state):** constant-forcing equilibrium temperature should match the algebraic balance of net flux against zero
- **Tier 5 (sensitivity):** parameter sweeps confirm physically reasonable response signs (e.g., warmer air → warmer water, stronger wind → faster equilibration)

These validation tests may be skipped for the initial v3 release (which is primarily a port + merge) but are essential foundations for the NSM1 work that follows.

---

## 5. Performance Targets

The current v2 baseline on Sumwere Creek (4,320 timesteps, coupled TSM+Riverine) is:

| Metric | v2 baseline |
|---|---|
| Wall time | 89 seconds |
| CPU time | 87 seconds |
| Per-step | ~20 ms (includes Riverine transport) |

v3 targets:

- **Must:** wall time within 1.5× of v2 baseline (≤ 135 seconds)
- **Should:** wall time equal to or better than v2 baseline (≤ 90 seconds)
- **Aspirational:** 2-5× faster than v2 baseline by porting the 418× kernel optimization concept to v2's Process dispatch (≤ 18-45 seconds)

If the initial port without optimization meets the "must" target, optimization work is deferred to a later v3 release.

---

## 6. Phased Implementation Plan

### Phase 0 — Gap-analysis diff table (½ day)

Produce a Markdown table comparing every meaningful difference between v1 TSM (with author's recent fixes applied, on the streaming repo) and v2 TSM (on `memory-refactor-pytestUpdate`). For each row, classify: `Port v1→v3`, `Keep v2 in v3`, `Reconcile in v3`, `Resolve TODO in v3`. Saved at `design/clearwater_modules_v3_tsm_gap_analysis.md`.

**Deliverable:** the gap analysis table itself, used as the work-tracking artifact for Phases 1-5.

### Phase 1 — v3 scaffold (½ day)

Create the v3 directory structure (architecture spec Section 4), with overlay imports pointing at v2 for everything except TSM. Add v3 to the pixi dev env. Verify `import clearwater_modules_v3` works and resolves symbols correctly.

**Deliverable:** v3 package scaffold, importable but functionally equivalent to v2 (no v3-native code yet).

### Phase 2 — v3 Temperature class (1-2 days)

Implement `clearwater_modules_v3/processes/temperature.py` as the merged TSM, applying every modification in Section 3.1.

**Deliverable:** v3 TSM importable and runnable; coupled notebook executes against v3 successfully.

### Phase 3 — v3 Model orchestration (2-3 days)

Implement `clearwater_modules_v3/model.py` with the three orchestration-level capabilities in Section 3.2 (kernel optimization, wet-mask, hotstart). Update `clearwater_modules_v3/config/init.py` to honor the new YAML keys.

**Deliverable:** v3 supports hotstart and wet-mask end-to-end; coupled notebook with hotstart kwargs runs.

### Phase 4 — Tests and validation (2-3 days)

Port the v1 TSM regression suite to v3 (Section 4). Add the v2/v3 parity test. Run the coupled notebook test programmatically.

**Deliverable:** v3 TSM at full feature parity with v1, with passing test suite.

### Phase 5 — Documentation and review prep (½ day)

Write v3 README. Write migration notes for downstream users. Prepare materials for LimnoTech review.

**Deliverable:** v3 ready for LimnoTech sign-off.

**Total estimated wall-clock with Claude doing the coding: 6-9 working days.**

---

## 7. TSM-Specific Risks

For umbrella risks (e.g., LimnoTech objection to v3 directory pattern, three-way split), see the architecture spec.

| Risk | Probability | Impact | Mitigation |
|---|---|---|---|
| Hotstart/wet-mask API designed differently than LimnoTech expects, requires rework | Medium | Medium (1-2 days lost) | Design in clearly-isolated modules; make rework cheap; communicate design early |
| v2 TSM produces materially different numerical outputs than v1 (beyond floating-point noise), root-cause unclear | Medium | Medium-High (1-3 days investigation) | Build the v2/v3 parity test early so divergences surface immediately; consult Jason's January 2026 work for prior diff investigation |
| Chunking reconciliation reveals incompatible designs that cannot be merged cleanly | Low-Medium | Medium | Discuss with LimnoTech early; if irreconcilable, choose one and document the rationale |
| Performance significantly worse than v2 (>1.5× slower) | Medium | Medium-High | Profile early; if regressed, port the kernel optimization into v3 in Phase 3 rather than deferring |

### Coordination with LimnoTech (TSM-specific items)

For umbrella coordination notes, see the architecture spec.

TSM-specific items that benefit from LimnoTech input:

- Hotstart API design for v2's `Process.run`-based dispatch (architecture spec Section 4 establishes the pattern; TSM is the first consumer)
- Wet-mask design for v2's per-substep filter (same)
- Chunking reconciliation: confirm that mirroring riverine's chunking patterns is the intended design

Items that can proceed without LimnoTech input:

- Phase 0 (gap analysis): pure analysis, no commitment
- Phase 1 (scaffold): structural, easily reversible
- Phase 2 items: pure ports of v1 fixes (latent-heat, thin-water, debug print removal), no design ambiguity
- Phase 4 test ports: framework-independent

---

## 8. Open Questions

All open questions resolved 2026-05-04 from evidence review (commit history, source diffs). Resolutions are summarized here and incorporated into the spec body.

1. ~~**Should v3's chunking implementation be Paul's `__process_loop_chunked` (v2) or a new design that incorporates the author's streaming work?**~~ **Resolved 2026-05-04:** Adopt v2's `__process_loop_chunked` skeleton in v3. There is no competing implementation in `clearwater_modules` to merge from — searching the legacy modules package for "chunk" returns exactly one match (a comment in `nsm1/model.py:51`). The author's streaming/chunking work lives in the `clearwater_riverine` repo, not in `clearwater_modules`. v3 resolves the four TODOs in Paul's `__process_loop_chunked` by mirroring riverine's chunking conventions per `clearwater_riverine.transport.py` and `clearwater_riverine.constituents.py`. See Section 3.2.
2. ~~**What is the semantics of `wet_mask` at the registry level?**~~ **Resolved 2026-05-04:** `wet_mask` is a registry-level concept; processes opt in by consulting the mask in `run`. See Section 3.2.
3. ~~**Should hotstart preserve the substep state of each process, or only the time-aligned global state?**~~ **Resolved 2026-05-04:** Hotstart preserves dataset-level (registry) state at the saved time. Per-process substep-internal state defaults to "fresh start" semantics after hotstart. Each `Process` may optionally implement `to_hotstart() -> dict` and `from_hotstart(state: dict)` methods to preserve richer internal state; default is no-op. v3 `Temperature.from_hotstart` sets `__skip_first_time_step = False`. Evidence: v1 hotstart (`base.py:64-100`, `nsm1/model.py:40-80`) is dataset-level; v2 has no hotstart (verified by grep); the only mutating substep state in v2 `Temperature` is `__skip_first_time_step`. See Section 3.2.
4. ~~**Does the Richardson `-1` factor belong in v3?**~~ **Resolved 2026-05-04: NO.** Delete the commented `-1` and TODO comment at `temperature.py:597`. v3 formula matches v1 (no leading `-1`). Evidence: Jason Rutyna's January 2026 investigation (commits `8218962` and `7f4166a`) concluded `-1` should not be there; v1 `tsm/processes.py:150-179` does not have it; the current v2 state has `-1` commented out with a TODO note that explicitly says "not in v1 of code." See Section 3.1.
5. ~~**What is the v2 retirement trigger?**~~ **Resolved 2026-05-04: Feature-based.** Trigger is "v3 1.0.0 ships." May 31 is the schedule target, not the trigger. If v3 1.0.0 slips, retirement also slips. See architecture spec Section 5.

---

## 9. Approval Criteria

**Status as of 2026-05-04: APPROVED for implementation.**

The author has reviewed and accepted:

1. ✓ The TSM-specific motivation (Section 1).
2. ✓ The TSM-specific goals and non-goals (Section 2).
3. ✓ The component inventory (Section 3) for TSM.
4. ✓ The test infrastructure and validation plan (Section 4).
5. ✓ The performance targets (Section 5) calibrated to the verified Sumwere Creek baseline.
6. ✓ The phased plan (Section 6) — Phases 0 through 5, ~6-9 working days with Claude executing.
7. ✓ The TSM-specific risks (Section 7) and mitigations.
8. ✓ All five open questions (Section 8) resolved and incorporated into the spec body.

Implementation can begin with Phase 0 (Section 6). See architecture spec Section 1 for prerequisites, branch setup, and environment activation.

---

## Appendix A: TSM Gap Analysis (preview, full table TBD in Phase 0)

A representative subset of the differences that Phase 0 will catalog:

| v1 TSM (with fixes) | v2 TSM (memory-refactor-pytestUpdate) | v3 disposition |
|---|---|---|
| `mf_latent_heat_vaporization`: converts K→C before polynomial | `latent_heat_vaporization`: converts C→K before polynomial (incorrect) | Port v1 logic |
| `dTdt_water_c`: depth ramp + rate cap | `temperature_change`: simple energy balance, no guard | Port v1 logic |
| `q_net`: returns `(...) * 86400 * dt` (W·s = J) | `flux_net`: returns sum of fluxes in W/m² (no time scaling) | Reconcile: v3 follows v2 framework convention; the 86400*dt factor moves to `temperature_change` |
| Debug prints commented out in some functions | Debug prints commented out in different functions | Remove all commented prints |
| No `mixing_ratio_air` divide-by-zero guard | `mixing_ratio_air`: `if atmospheric_pressure == atmospheric_vapor_pressure: return 0.0` | Keep v2 logic |
| No "skip first time step" coupling logic | `__skip_first_time_step` flag in `__init__` and `run` | Keep v2 logic |
| Richardson formula does not have leading `-1` | Richardson formula has `-1` commented with TODO | **Resolved 2026-05-04:** Remove the `-1` and TODO comment; v3 matches v1 (no leading `-1`). Per Jason Rutyna's January 2026 investigation. |
| `q_sediment` uses `/ 0.5` factor undocumented | `flux_sediment` uses `/ 0.5` factor with TODO | Document the `/ 0.5` as the sediment active-layer half-thickness convention (matches v1 legacy); replace TODO with one-line docstring. Do not change the formula. |

Phase 0 produces the full table.
