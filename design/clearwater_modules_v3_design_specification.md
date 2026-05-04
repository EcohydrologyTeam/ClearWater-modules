# ClearWater Modules v3 — Design Specification

**Status:** Approved for implementation
**Author:** Todd Steissberg (ERDC), with Claude
**Date:** 2026-05-04 (revised 2026-05-04 with all open questions resolved)
**Scope:** Architectural specification for `clearwater_modules_v3`. TSM is the first module to be implemented; NSM1 will follow under a separate but compatible specification (`clearwater_modules_v3_nsm1_design_specification.md`).

---

## 0. Implementation Quick Start

This section gives the implementing agent everything needed to start work without re-reading the conversational history that produced this spec.

### Prerequisites

- Pixi installed (already done; binary at `/opt/homebrew/Cellar/pixi/0.67.2`)
- Pixi dev environment **already built in the modules repo** at `/Users/todd/GitHub/ecohydrology/ClearWater-modules/.pixi/envs/dev/` (verified working with the coupled TSM+Riverine demo). v3 development happens in the streaming repo (see below), so a **second pixi env** must be built there before any v3 code can run.
- All three baseline repos cloned and on the right branches:
  - `ClearWater-modules` on `memory-refactor-pytestUpdate` (LimnoTech-visible baseline; this branch is where the canonical v3 spec was committed yesterday — **do not push v3 development to this repo until the work is ready to share with LimnoTech**)
  - `ClearWater-modules-streaming` on `streaming` (private fork; v3 development happens here on a new branch)
  - `ClearWater-riverine` on `refactor-demo`
  - `ClearWater-data` on `main`
- Coupled TSM+Riverine demo verified to run on the modules-repo baseline (89s for 4,320 timesteps on Sumwere Creek as of 2026-05-04)

### Branch and working directory

- All v3 development happens in `EcohydrologyTeam/ClearWater-modules-streaming` (the private fork) on a new feature branch off `memory-refactor-pytestUpdate`. The streaming repo has `memory-refactor-pytestUpdate` mirrored from the modules repo via the `upstream` remote.
- Suggested branch name: `v3-tsm-merge`.
- v3 source goes in `src/clearwater_modules_v3/`.
- The new branch will be pushed to `EcohydrologyTeam/ClearWater-modules-streaming` (private fork), **not** to `EcohydrologyTeam/ClearWater-modules` (the LimnoTech-visible repo).

```bash
cd /Users/todd/GitHub/ecohydrology/ClearWater-modules-streaming
git fetch --all
git checkout -b v3-tsm-merge origin/memory-refactor-pytestUpdate
# (or upstream/memory-refactor-pytestUpdate; both point at the same commit)
```

### Pixi env setup (one-time, before any v3 code can run)

The existing pixi `dev` env in the modules repo has editable installs pointing at modules-repo paths and won't pick up new `src/clearwater_modules_v3/` files added to the streaming repo. A second pixi env must be built in the streaming repo:

```bash
cd /Users/todd/GitHub/ecohydrology/ClearWater-modules-streaming
pixi install -e dev
```

This creates `.pixi/envs/dev/` in the streaming repo with editable installs from streaming-repo paths. ~5 minutes the first time. After that, all subsequent `pixi shell -e dev` invocations from the streaming-repo cwd activate the streaming-repo env.

### Environment activation

```bash
cd /Users/todd/GitHub/ecohydrology/ClearWater-modules-streaming
pixi shell -e dev
```

All commands subsequently assume this shell (Python, pytest, jupyter all resolve to the streaming-repo dev env, which sees `src/clearwater_modules_v3/` once it exists).

### First action

**Start with Phase 0 (Section 8 below, "Phase 0 — Gap-analysis diff table"): produce the gap-analysis diff table comparing v1 TSM (with the author's recent fixes on the `streaming` branch of `ClearWater-modules-streaming`, file `src/clearwater_modules/tsm/processes.py`) against v2 TSM (on `memory-refactor-pytestUpdate`, file `src/clearwater_modules_v2/processes/temperature.py`).** Save the table at `/Users/todd/GitHub/ecohydrology/ClearWater-modules-streaming/design/clearwater_modules_v3_tsm_gap_analysis.md`. The gap-analysis table becomes the work-tracking artifact for Phases 1-5.

### Implementation invariants

The following must hold true at every commit:

1. **The coupled TSM+Riverine demo notebook (`examples/03_Example_Coupled_TSM_and _Riverine.ipynb`) continues to execute without error.** This is the smoke test. It currently runs against v2 TSM. After Phase 2, it should also run against v3 TSM via a config edit. After all phases, it produces equivalent results.
2. **All v2 tests continue to pass.** v3 work does not modify v2 source files.
3. **No new dependencies introduced** beyond what's already in the pixi `dev` environment.
4. **No mention of Claude or AI assistance in commit messages or PR descriptions.**
5. **Commit at every meaningful checkpoint** (after each Phase, or after completing a discrete fix). Push frequently so progress is visible.

### Parallel NSM1 work

NSM1 v3 will be developed in a sibling session, also against `memory-refactor-pytestUpdate`, on a separate branch (suggested: `v3-nsm1-merge`). NSM1 v3 depends on TSM v3 having completed Phase 1 (the v3 package scaffold) and Phase 3 items 1-2 (kernel optimization + wet-mask in v3 Model). NSM1 v3 Phase 0 can start immediately; NSM1 v3 Phase 1 (shared utilities) can start as soon as TSM v3 Phase 1 completes. Coordinate via git: NSM1 session pulls TSM session's commits before starting its Phase 1.

### Where to find supporting docs

- This spec: `streaming/design/clearwater_modules_v3_design_specification.md`
- NSM1 v3 spec: `streaming/design/clearwater_modules_v3_nsm1_design_specification.md`
- v1 vs v2 inventory: `streaming/design/TSM_NSM1_v1_vs_v2_inventory.md`
- Validation strategy: `Validation_strategy_no_reference.md` (in `~/Downloads/NSM_comparison/` or wherever Todd has it; provides Tier 1 conservation + Tier 2 analytical-limit test patterns)
- DOX bug investigation (legacy modules constants): `streaming/design/NSM1_DOX_rate_bug_investigation.md`

---

## 1. Background and Motivation

The ClearWater Modules codebase currently contains two parallel implementations of TSM and NSM1 that have evolved on independent tracks:

- **`clearwater_modules` (v1)** — the function-style framework originally developed by the ECOMOD team and substantially extended by the author with Claude in April–May 2026. Lives in the `streaming` branch of `EcohydrologyTeam/ClearWater-modules-streaming`.
- **`clearwater_modules_v2`** — the class-based framework developed by LimnoTech (Paul Tomasula, Anthony Aufdenkampe, Jason Rutyna, Sarah Jordan) from July 2025 through April 2026, with TSM and partial NSM1 implementations. Lives in the `memory-refactor-pytestUpdate` branch of `EcohydrologyTeam/ClearWater-modules`.

Each track produced substantive contributions the other lacks. v1 received orchestration-level optimization, numerical correctness fixes, and new capabilities (hotstart, wet-mask) but stayed in the function-style framework. v2 introduced a modern class-based framework, YAML configuration, per-process substepping, and a chunking execution path but did not receive v1's later corrections and capabilities. The two tracks diverged from a common ancestor in late 2024 and have not been synchronized since.

**v3 is the convergence of these two tracks** into a single coherent codebase that combines v2's framework with v1's optimization and correctness work, plus the synthesis required to make them coexist.

### Why a new package rather than modifying v2

1. v2 remains unmodified during v3's development, so neither LimnoTech's framework work nor the author's optimization work is at risk while the merge is in progress.
2. v3 and v2 can be imported side by side in the same Python session, enabling direct A/B regression comparison against the verified Sumwere Creek baseline.
3. The diff between v2 and v3 is exactly the merge work — no noise from copied-but-unchanged files.
4. The v3 directory is a concrete artifact LimnoTech can review at their own pace without coordination friction.
5. The pattern follows the precedent LimnoTech established when they created `clearwater_modules_v2` alongside `clearwater_modules` rather than modifying v1 in place.

### What makes this v3, not v2.x

| Source | Contribution |
|---|---|
| LimnoTech v2 | `Process` class abstraction with explicit composition |
| LimnoTech v2 | YAML-driven configuration via `init_from_file` |
| LimnoTech v2 | Per-process substepping (`time_step` on each `Process`) |
| LimnoTech v2 | `clearwater_data.VariableRegistry` integration |
| LimnoTech v2 | Chunking execution path (`Model.__process_loop_chunked`, March 2026) |
| LimnoTech v2 | `mixing_ratio_air` divide-by-zero correction (Anthony, April 2026) |
| LimnoTech v2 | TSM "skip first step" coupling logic (Paul, March 2026) |
| LimnoTech v2 | Sediment-temperature optional flag (Anthony, February 2026) |
| Author / Claude (v1) | 418× compute-plan kernel optimization with cached plan and direct array writes |
| Author / Claude (v1) | Streaming/chunking framework |
| Author / Claude (v1) | Wet-mask gating at the orchestration layer |
| Author / Claude (v1) | Hotstart from `xr.Dataset` checkpoint |
| Author / Claude (v1) | TSM latent-heat unit correction (Celsius polynomial vs Kelvin input) |
| Author / Claude (v1) | TSM thin-water stability guard (depth ramp + `dTdt_max_per_hour` rate cap) |
| Author / Claude (v1) | Multi-cell-safe debug print removal |
| v3 synthesis | Reconciliation of two chunking implementations |
| v3 synthesis | Resolution of Richardson `-1` factor TODO |
| v3 synthesis | Resolution of `flux_sediment / 0.5` TODO |
| v3 synthesis | A coherent design where wet-mask, hotstart, and per-process substepping coexist |
| v3 synthesis | Test infrastructure validating v3 against both v1 and v2 outputs |

The combined contribution is a meaningful version increment, not a bug-fix release.

---

## 2. Goals and Non-Goals

### Goals

1. Single coherent codebase that supersedes both v1 and v2 once mature.
2. Backward-compatible YAML configuration: existing v2 configs run on v3 without modification.
3. Performance equal to or better than v1's optimized kernel on the Sumwere Creek coupled benchmark (current v2 baseline: 89 s for 4,320 timesteps).
4. Numerical outputs match v1's corrected TSM within floating-point tolerance for the regression test suite.
5. End-to-end coupled TSM+Riverine demo runs cleanly on v3 with no notebook code changes beyond import statements.
6. Test coverage at least matching v1's TSM regression suite (15 calculation tests + 4 latent-heat + 6 stability ramp + hotstart roundtrip).
7. Clear retirement path for v2.

### Non-Goals

1. v3 will not introduce a new architectural framework. The class-based v2 framework remains the architectural baseline.
2. v3 will not break v2's public API. Process names, configuration keys, and registry names remain compatible.
3. v3 will not address NSM2 features (multi-pool organic matter, alkalinity, methane/sulfide, silica). Those are scoped for separate work after NSM1 is complete in v3.
4. v3 will not retire v1 immediately. v1 stays in the codebase as a deprecated reference for one release cycle, then is removed.
5. v3 will not require a different Python environment or dependency set than v2.

---

## 3. Architectural Approach

### Thin-overlay strategy

v3 starts as a thin overlay over v2. Components that have not been modified relative to v2 are imported from v2 at the package level; components that have been modified are implemented in v3 directly.

```
src/clearwater_modules_v3/
├── __init__.py                    # re-exports from v2 + v3-specific
├── README.md                      # the framing in Section 1
├── model.py                       # v3 Model with kernel optimization + wet-mask + hotstart + chunking
├── config/
│   ├── __init__.py                # re-exports init_from_file from v2 initially
│   └── init.py                    # v3 init_from_file (added wet-mask, hotstart kwargs)
├── processes/
│   ├── __init__.py                # imports unmodified from v2; overrides v3-native
│   ├── base.py                    # re-export from v2 (Process, ProcessFactory)
│   ├── temperature.py             # v3-native (the merged TSM)
│   ├── riverine.py                # re-export from v2 initially
│   ├── nitrogen.py                # re-export from v2 initially; v3-native after NSM1 merge
│   ├── floating_algae.py          # re-export from v2 initially
│   └── benthic_algae.py           # re-export from v2 initially
└── utils/
    ├── __init__.py
    ├── constants.py               # re-export from v2
    └── conversions.py             # re-export from v2
```

As more components migrate from v2-overlay to v3-native, the imports in `processes/__init__.py` are progressively replaced. When every component is v3-native, the overlay imports are deleted in a final cleanup PR and v3 stands alone.

### Why thin overlay rather than full copy

- Day-one code duplication is ~200 lines instead of ~2,300.
- The diff between v2 and v3 in version control reflects only the merge work, not file moves.
- LimnoTech's continued v2 work (if any) before v3 is mature can be pulled forward into v3 by replacing the overlay imports incrementally.
- Reduces the merge-conflict surface if v2 and v3 evolve in parallel during the transition.

### Package registration

v3 is registered as a separate editable Python package via the modules repo's `pyproject.toml`. The pixi `dev` environment installs both v2 and v3 editable, so both are importable in the same shell.

```toml
[tool.pixi.feature.dev.pypi-dependencies]
ClearWater-modules = { path = ".", editable = true }
clearwater_riverine = { path = "../ClearWater-riverine", editable = true }
clearwater_data = { path = "../ClearWater-data", editable = true }
```

The existing entry already installs both `clearwater_modules` (v1) and `clearwater_modules_v2` from the same package (the `path = "."` entry covers both subpackages). Adding `clearwater_modules_v3` requires no new pixi entry — it is part of the same package.

---

## 4. Component Inventory: TSM

This section enumerates what v3 TSM contains and where each piece originates. NSM1 will receive an analogous inventory under its own specification.

### 4.1 v3 TSM source: `src/clearwater_modules_v3/processes/temperature.py`

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

### 4.2 v3 Model: `src/clearwater_modules_v3/model.py`

The v3 `Model` class extends v2's `Model` with three orchestration-level capabilities ported from v1:

| Capability | v1 source | v3 design |
|---|---|---|
| 418× kernel optimization | v1 commit `6daa65e` (`_iter_computations_fast`, cached compute plan) | Adapt to v2's `Process.run`-based dispatch: cache the (process, time_step, variables) tuples once per Model instance; in the per-step loop, dispatch directly without re-checking `current_time_seconds % process.time_step_seconds` (precompute next-fire time per process). **Design note (to be flagged in v3 code at the kernel-optimization site):** this adapts v1's variable-level cached-plan optimization to v2's process-level dispatch. v2's process-level loop is coarser-grained than v1's variable-level loop, so the speedup may not match v1's 418× directly. **Follow-up optimization recommended after Phase 3:** profile against v2 baseline; if v3 is materially slower than v1, add a second optimization layer at the variable level inside `Process.run` (e.g., cache the registry-name lookups, batch `get_at_time`/`set_at_time` calls per process, or build a per-Process compute plan that bypasses the registry for hot variables). The recommendation is documented at the relevant code site so future maintainers see the design intent. |
| Wet-mask gating | v1 commit `3d18965` (`wet_mask` kwarg) | Add a `wet_mask` registry variable that processes can opt into; `Model.run` honors the mask by skipping `set_at_time` for masked cells. Existing per-process mask logic (e.g., the `xr.where(volume > 0, ...)` in v2 `Temperature.run`) becomes redundant once the registry-level mask is in place; remove from process bodies |
| Hotstart from xr.Dataset | v1 commit `1a226dd` (`hotstart_dataset`, `hotstart_timestep`); v2 has no hotstart | Add `hotstart_dataset` and `hotstart_timestep` kwargs to `init_from_file`. `Model.__init__` honors them by seeding the registry from the dataset at `hotstart_timestep`. **Resolved 2026-05-04: substep state semantics** — hotstart preserves dataset-level (registry) state at the saved time. Per-process substep-internal state (e.g., v2 `Temperature.__skip_first_time_step`) defaults to "fresh start" semantics after hotstart. Each `Process` may optionally implement `to_hotstart() -> dict` and `from_hotstart(state: dict)` methods to preserve richer internal state; default is no-op. v3 `Temperature.from_hotstart` sets `__skip_first_time_step = False` (don't skip the first post-hotstart step — you're not starting from scratch). |
| Chunking execution — **Resolved 2026-05-04** | v2 `Model.__process_loop_chunked` (Paul's commit `d712c59`, 2026-03-11) | Inherit v2's `__process_loop_chunked` skeleton in v3. **Resolve the four TODOs Paul left in the method:** (1) "this need actual chunking logic", (2) "look at riverine's code and mirror where applicable", (3) "align with riverine", (4) "confirm if this is necessary to write out the last chunk". Resolution approach: study `clearwater_riverine.transport.py` and `clearwater_riverine.constituents.py` chunking patterns (since the author's streaming/chunking work lives in the riverine repo, not in `clearwater_modules` — searching the legacy modules package for "chunk" returns only one comment in `nsm1/model.py:51`). Mirror riverine's chunk-boundary handling, last-chunk write semantics, and dataset alignment. |

### 4.3 v3 Configuration: `src/clearwater_modules_v3/config/init.py`

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

---

## 5. Migration Strategy

### From v2 to v3

For users who have v2 code:

| v2 import | v3 equivalent |
|---|---|
| `import clearwater_modules_v2 as cwm` | `import clearwater_modules_v3 as cwm` |
| `from clearwater_modules_v2.config import init_from_file` | `from clearwater_modules_v3.config import init_from_file` |
| `from clearwater_modules_v2.processes.temperature import Temperature` | `from clearwater_modules_v3.processes.temperature import Temperature` |

A migration table will be included in v3's README.

### v2 retirement plan

**Trigger (resolved 2026-05-04):** v2 retirement is **feature-based, not time-based.** The trigger event is "v3 1.0.0 ships" (which means TSM and NSM1 complete and validated). The May 31 date is the **target schedule**, not the trigger itself; if v3 1.0.0 slips, retirement also slips. The schedule controls communications and external commitments; the trigger controls the action.

Target schedule: **2026-05-31** for v3 1.0.0 with both TSM and NSM1 complete and validated. This is an aggressive 27-day window from spec date (2026-05-04) and assumes Claude executes the implementation autonomously with the author reviewing at each phase boundary.

Scope clarification for v3 1.0.0: v3 1.0.0 includes the *legacy NSM1 capability set* ported into the v3 framework with the merged corrections and orchestration improvements — i.e., the 16 NSM1 constituents, reaeration menu, light extinction, and other v1 NSM1 features. NSM2 features (multi-pool organic matter, alkalinity/pH, methane/sulfide, silica) are explicitly out of scope for v3 1.0.0 and will be added incrementally in subsequent v3 releases (1.1, 1.2, ...) over the weeks and months following.

1. **2026-05-04 to ~2026-05-13:** v3 0.1.0 — TSM complete, tested, coupled demo runs. v2 enters frozen state at this point: bug fixes only, no new features. All new development goes to v3.
2. **2026-05-13 to 2026-05-31:** v3 1.0.0 — NSM1 complete and validated. v2 marked deprecated in its README.
3. **Shortly after 2026-05-31:** v3 1.1.0 (or 1.0.1, depending on what surfaces during NSM1 validation) removes v2 from the source tree in a single cleanup PR.

If the May 31 target slips, the slip is communicated as it happens rather than discovered after the fact. The aggressive timeline is a goal, not a constraint that will be honored at the cost of correctness or test coverage.

### v1 retirement plan

v1 (`clearwater_modules`) remains in the repo as a reference implementation through the v3 1.0.0 release. v1 is removed in the same cleanup PR that retires v2 (target: shortly after 2026-05-31).

---

## 6. Testing and Validation

### Test infrastructure

v3 inherits the test directory structure from v2 (`tests/`) and adds:

- `tests/v3/test_5_tsm_calculations_v3.py` — v3 equivalent of v1's `test_5_tsm_calculations.py`, using the same expected-value tables in `tests/NSM Manual Calcs/` but exercising the v3 `Temperature.run` flow
- `tests/v3/test_tsm_latent_heat_v3.py` — port of v1's `test_tsm_latent_heat.py`
- `tests/v3/test_tsm_stability_ramp_v3.py` — port of v1's `test_tsm_stability_ramp.py`
- `tests/v3/test_hotstart_roundtrip_v3.py` — port of v1's `test_hotstart_roundtrip.py`, adapted for v3's hotstart API
- `tests/v3/test_v2_v3_parity.py` — A/B comparison: load Sumwere Creek modules.yml, run v2 and v3 in the same Python session, assert that water temperature timeseries match within tolerance after applying the v3 corrections to a v2-equivalent control case (i.e., when v3's latent-heat fix and thin-water guard are disabled, outputs should be identical to v2's outputs)
- `tests/v3/test_coupled_tsm_riverine_v3.py` — programmatic version of `examples/03_Example_Coupled_TSM_and _Riverine.ipynb` with v3, asserting completion within a wall-time budget

### Validation tiers (per the validation strategy doc)

- **Tier 1 (conservation):** mass and energy conservation tests for TSM (closed-system test with no fluxes; total enthalpy must be conserved)
- **Tier 2 (analytical limits):** verification of individual heat-flux components against analytical Brutsaert and Brunt formulations
- **Tier 3 (steady-state):** constant-forcing equilibrium temperature should match the algebraic balance of net flux against zero
- **Tier 5 (sensitivity):** parameter sweeps confirm physically reasonable response signs (e.g., warmer air → warmer water, stronger wind → faster equilibration)

These validation tests may be skipped for the initial v3 release (which is primarily a port + merge) but are essential foundations for the NSM1 work that follows.

---

## 7. Performance Targets

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

## 8. Phased Implementation Plan (TSM)

### Phase 0 — Gap-analysis diff table (½ day)

Produce a Markdown table comparing every meaningful difference between v1 TSM (with author's recent fixes applied, on the streaming repo) and v2 TSM (on `memory-refactor-pytestUpdate`). For each row, classify: `Port v1→v3`, `Keep v2 in v3`, `Reconcile in v3`, `Resolve TODO in v3`. Saved at `docs/clearwater_modules_v3_tsm_gap_analysis.md`.

**Deliverable:** the gap analysis table itself, used as the work-tracking artifact for Phase 1.

### Phase 1 — v3 scaffold (½ day)

Create the v3 directory structure (Section 3), with overlay imports pointing at v2 for everything except TSM. Add v3 to the pixi dev env. Verify `import clearwater_modules_v3` works and resolves symbols correctly.

**Deliverable:** v3 package scaffold, importable but functionally equivalent to v2 (no v3-native code yet).

### Phase 2 — v3 Temperature class (1-2 days)

Implement `clearwater_modules_v3/processes/temperature.py` as the merged TSM, applying every modification in Section 4.1.

**Deliverable:** v3 TSM importable and runnable; coupled notebook executes against v3 successfully.

### Phase 3 — v3 Model orchestration (2-3 days)

Implement `clearwater_modules_v3/model.py` with the three orchestration-level capabilities in Section 4.2 (kernel optimization, wet-mask, hotstart). Update `clearwater_modules_v3/config/init.py` to honor the new YAML keys.

**Deliverable:** v3 supports hotstart and wet-mask end-to-end; coupled notebook with hotstart kwargs runs.

### Phase 4 — Tests and validation (2-3 days)

Port the v1 TSM regression suite to v3 (Section 6). Add the v2/v3 parity test. Run the coupled notebook test programmatically.

**Deliverable:** v3 TSM at full feature parity with v1, with passing test suite.

### Phase 5 — Documentation and review prep (½ day)

Write v3 README. Write migration notes for downstream users. Prepare materials for LimnoTech review.

**Deliverable:** v3 ready for LimnoTech sign-off.

**Total estimated wall-clock with Claude doing the coding: 6-9 working days.**

---

## 9. Coordination with LimnoTech

### What to discuss with Anthony

1. The v3 directory exists as a merge proposal, not a fait accompli. LimnoTech is invited to review and propose changes.
2. v3 explicitly intends to converge with or replace v2 over time. v3 is not a permanent third version.
3. The orchestration-level additions (hotstart, wet-mask, kernel optimization) are designed to be backward-compatible with v2's process composition pattern. If LimnoTech wants different API surfaces for these, that is in scope for review.
4. The chunking reconciliation (v2's `__process_loop_chunked` versus the author's v1 streaming work) needs LimnoTech's input on the intended design.
5. The Richardson `-1` factor decision should reference Jason Rutyna's January 2026 investigation; if Jason has a documented conclusion, v3 should adopt it.

### What can proceed without LimnoTech input

- Phase 0 (gap analysis): pure analysis, no commitment
- Phase 1 (scaffold): structural, easily reversible
- Phase 2 items 1, 2, 3 (latent-heat, thin-water, debug print removal): pure ports of v1 fixes, no design ambiguity
- Phase 4 test ports: framework-independent

Items that benefit from LimnoTech input:

- Phase 2 Richardson `-1` resolution
- Phase 3 hotstart and wet-mask API design
- Chunking reconciliation
- v2 retirement timeline

---

## 10. Risks and Mitigations

| Risk | Probability | Impact | Mitigation |
|---|---|---|---|
| Hotstart/wet-mask API designed differently than LimnoTech expects, requires rework | Medium | Medium (1-2 days lost) | Design in clearly-isolated modules; make rework cheap; communicate design early |
| v2 TSM produces materially different numerical outputs than v1 (beyond floating-point noise), root-cause unclear | Medium | Medium-High (1-3 days investigation) | Build the v2/v3 parity test early so divergences surface immediately; consult Jason's January 2026 work for prior diff investigation |
| Permanent three-way split: v1, v2, v3 all maintained indefinitely | Low (with explicit retirement plan) | High (compounds maintenance burden) | Section 5 retirement plan; commit to v2 freeze when v3 0.1.0 ships |
| Chunking reconciliation reveals incompatible designs that cannot be merged cleanly | Low-Medium | Medium | Discuss with LimnoTech early; if irreconcilable, choose one and document the rationale |
| Performance significantly worse than v2 (>1.5× slower) | Medium | Medium-High | Profile early; if regressed, port the kernel optimization into v3 in Phase 3 rather than deferring |
| LimnoTech objects to v3 as a directory | Low (consistent with their own v2 precedent) | High (requires re-planning) | Frame as a working merge proposal in the conversation with Anthony; v3 is reversible |

---

## 11. Open Questions

All five open questions resolved 2026-05-04 from evidence review (commit history, source diffs). Resolutions are summarized here and incorporated into the spec body.

1. ~~**Should v3's chunking implementation be Paul's `__process_loop_chunked` (v2) or a new design that incorporates the author's streaming work?**~~ **Resolved 2026-05-04:** Adopt v2's `__process_loop_chunked` skeleton in v3. There is no competing implementation in `clearwater_modules` to merge from — searching the legacy modules package for "chunk" returns exactly one match (a comment in `nsm1/model.py:51`). The author's streaming/chunking work lives in the `clearwater_riverine` repo, not in `clearwater_modules`. v3 resolves the four TODOs in Paul's `__process_loop_chunked` by mirroring riverine's chunking conventions per `clearwater_riverine.transport.py` and `clearwater_riverine.constituents.py`. See Section 4.2.
2. ~~**What is the semantics of `wet_mask` at the registry level?**~~ **Resolved 2026-05-04:** `wet_mask` is a registry-level concept; processes opt in by consulting the mask in `run`. See Section 4.2.
3. ~~**Should hotstart preserve the substep state of each process, or only the time-aligned global state?**~~ **Resolved 2026-05-04:** Hotstart preserves dataset-level (registry) state at the saved time. Per-process substep-internal state defaults to "fresh start" semantics after hotstart. Each `Process` may optionally implement `to_hotstart() -> dict` and `from_hotstart(state: dict)` methods to preserve richer internal state; default is no-op. v3 `Temperature.from_hotstart` sets `__skip_first_time_step = False`. Evidence: v1 hotstart (`base.py:64-100`, `nsm1/model.py:40-80`) is dataset-level; v2 has no hotstart (verified by grep); the only mutating substep state in v2 `Temperature` is `__skip_first_time_step`. See Section 4.2.
4. ~~**Does the Richardson `-1` factor belong in v3?**~~ **Resolved 2026-05-04: NO.** Delete the commented `-1` and TODO comment at `temperature.py:597`. v3 formula matches v1 (no leading `-1`). Evidence: Jason Rutyna's January 2026 investigation (commits `8218962` and `7f4166a`) concluded `-1` should not be there; v1 `tsm/processes.py:150-179` does not have it; the current v2 state has `-1` commented out with a TODO note that explicitly says "not in v1 of code." See Section 4.1.
5. ~~**What is the v2 retirement trigger?**~~ **Resolved 2026-05-04: Feature-based.** Trigger is "v3 1.0.0 ships." May 31 is the schedule target, not the trigger. If v3 1.0.0 slips, retirement also slips. See Section 5.

---

## 12. Approval Criteria

**Status as of 2026-05-04: APPROVED for implementation.**

The author has reviewed and accepted:

1. ✓ The motivation for v3 (Section 1), including LimnoTech's contributions and the author's contributions.
2. ✓ The thin-overlay architecture (Section 3) as the day-one structure.
3. ✓ The component inventory (Section 4) for TSM.
4. ✓ The migration strategy (Section 5) for v2 and v1 retirement, with feature-based trigger and 2026-05-31 schedule target.
5. ✓ The phased plan (Section 8) — Phases 0 through 5, ~6-9 working days with Claude executing.
6. ✓ The risk register (Section 10) and mitigations.
7. ✓ All five open questions (Section 11) resolved and incorporated into the spec body.

Implementation can begin with Phase 0 (Section 8.0). See Section 0 (Implementation Quick Start) for prerequisites, branch setup, and environment activation.

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

---

## Appendix B: Pixi / install changes

v3 requires no new pixi entries. Verify after scaffolding:

```bash
cd ClearWater-modules
pixi run -e dev python -c "import clearwater_modules_v3; print(clearwater_modules_v3.__file__)"
```

Expected: prints the path to `src/clearwater_modules_v3/__init__.py`.
