# ClearWater Modules v3 — Architecture Specification

**Status:** Approved for implementation
**Author:** Todd Steissberg (ERDC), with Claude
**Date:** 2026-05-04
**Scope:** Umbrella architectural specification for the `clearwater_modules_v3` package. Defines the framework, package structure, environment setup, migration plan, and v2/v1 retirement criteria that apply to **every** module in v3. Module-specific design lives in companion specifications.

**Companion documents:**

- `clearwater_modules_v3_tsm_design_specification.md` — TSM (Temperature Simulation Module)
- `clearwater_modules_v3_nsm1_design_specification.md` — NSM1 (Nutrient Simulation Module v1)

This document and the module-specific specifications are intended to be read together. Anything that applies to all of v3 lives here; anything specific to one module lives in that module's spec.

---

## 1. Implementation Quick Start

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
- Suggested branch names: `v3-tsm-merge` for the TSM session, `v3-nsm1-merge` for the parallel NSM1 session.
- v3 source goes in `src/clearwater_modules_v3/`.
- Branches are pushed to `EcohydrologyTeam/ClearWater-modules-streaming` (private fork), **not** to `EcohydrologyTeam/ClearWater-modules` (the LimnoTech-visible repo).

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

### Implementation invariants (apply to every v3 commit)

1. **The coupled TSM+Riverine demo notebook (`examples/03_Example_Coupled_TSM_and _Riverine.ipynb`) continues to execute without error.** This is the smoke test. After modules migrate to v3, the demo runs against v3 via a config edit.
2. **All v2 tests continue to pass.** v3 work does not modify v2 source files.
3. **No new dependencies introduced** beyond what's already in the pixi `dev` environment.
4. **No mention of Claude or AI assistance in commit messages or PR descriptions.**
5. **Commit at every meaningful checkpoint** (after each Phase, or after completing a discrete fix). Push frequently so progress is visible.

### Parallel TSM and NSM1 sessions

TSM v3 and NSM1 v3 are developed in parallel sessions, each on its own branch off `memory-refactor-pytestUpdate`. Coordination via git:

- TSM session creates the v3 package scaffold (TSM Phase 1) and the v3 `Model` orchestration improvements (TSM Phase 3 items 1-2: kernel optimization, wet-mask).
- NSM1 session can start its Phase 0 (gap analysis) immediately, since gap analysis doesn't depend on v3 code yet.
- NSM1 session waits for TSM Phase 1 (v3 package scaffold) before starting NSM1 Phase 1 (shared utilities, which live in the v3 package).
- NSM1 session pulls TSM session's commits before starting any phase that depends on TSM-session output.
- Final integration: both branches merge into a v3 integration branch (e.g., `v3-integration`) for end-to-end testing.

### Where to find supporting docs

All under `/Users/todd/GitHub/ecohydrology/ClearWater-modules-streaming/design/`:

- This spec: `clearwater_modules_v3_architecture_specification.md` (umbrella)
- TSM v3 spec: `clearwater_modules_v3_tsm_design_specification.md`
- NSM1 v3 spec: `clearwater_modules_v3_nsm1_design_specification.md`
- v1 vs v2 inventory: `TSM_NSM1_v1_vs_v2_inventory.md`
- DOX bug investigation: `NSM1_DOX_rate_bug_investigation.md`
- Validation strategy: `Validation_strategy_no_reference.md` (in `~/Downloads/NSM_comparison/`; provides Tier 1 conservation + Tier 2 analytical-limit test patterns)

---

## 2. Background and Motivation

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
| Author / Claude (v1) | Streaming/chunking framework (in the riverine repo) |
| Author / Claude (v1) | Wet-mask gating at the orchestration layer |
| Author / Claude (v1) | Hotstart from `xr.Dataset` checkpoint |
| Author / Claude (v1) | TSM latent-heat unit correction (Celsius polynomial vs Kelvin input) |
| Author / Claude (v1) | TSM thin-water stability guard (depth ramp + `dTdt_max_per_hour` rate cap) |
| Author / Claude (v1) | Multi-cell-safe debug print removal |
| v3 synthesis | Reconciliation of two chunking implementations |
| v3 synthesis | Resolution of Richardson `-1` factor TODO |
| v3 synthesis | Resolution of `flux_sediment / 0.5` TODO |
| v3 synthesis | A coherent design where wet-mask, hotstart, and per-process substepping coexist |
| v3 synthesis | Test infrastructure that validates v3 against both v1 and v2 outputs |

The combined contribution is a meaningful version increment, not a bug-fix release.

---

## 3. Goals and Non-Goals (umbrella)

### Goals

1. Single coherent codebase that supersedes both v1 and v2 once mature.
2. Backward-compatible YAML configuration: existing v2 configs run on v3 without modification.
3. Test infrastructure that validates v3 against both v1 and v2 outputs.
4. Clear retirement path for v2 and v1.
5. Each module preserves the v2 `Process` class abstraction; no new architectural framework is introduced.

### Non-Goals

1. v3 will not introduce a new architectural framework. The class-based v2 framework remains the architectural baseline.
2. v3 will not break v2's public API. Process names, configuration keys, and registry names remain compatible.
3. v3 will not address NSM2 features (multi-pool organic matter, alkalinity, methane/sulfide, silica). Those are scoped for separate work after NSM1 is complete in v3.
4. v3 will not retire v1 immediately. v1 stays in the codebase as a deprecated reference for one release cycle, then is removed.
5. v3 will not require a different Python environment or dependency set than v2.

Module-specific goals and non-goals are documented in each module's spec.

---

## 4. Architectural Approach

### Thin-overlay strategy

v3 starts as a thin overlay over v2. Components that have not been modified relative to v2 are imported from v2 at the package level; components that have been modified are implemented in v3 directly.

```
src/clearwater_modules_v3/
├── __init__.py                    # re-exports from v2 + v3-specific
├── README.md                      # the framing in Section 2 of this spec
├── model.py                       # v3 Model with kernel optimization + wet-mask + hotstart + chunking
├── config/
│   ├── __init__.py                # re-exports init_from_file from v2 initially
│   └── init.py                    # v3 init_from_file (added wet-mask, hotstart kwargs)
├── processes/
│   ├── __init__.py                # imports unmodified from v2; overrides v3-native
│   ├── base.py                    # re-export from v2 (Process, ProcessFactory)
│   ├── temperature.py             # v3-native (the merged TSM)
│   ├── riverine.py                # re-export from v2 initially
│   ├── nitrogen.py                # eventually v3-native; overlay from v2 initially
│   ├── floating_algae.py          # eventually v3-native; overlay from v2 initially
│   ├── benthic_algae.py           # eventually v3-native; overlay from v2 initially
│   └── (additional NSM1 processes added during NSM1 v3 work)
└── utils/
    ├── __init__.py
    ├── constants.py               # re-export from v2
    ├── conversions.py             # re-export from v2
    ├── reaeration.py              # new in v3 (NSM1 dependency)
    ├── sediment.py                # new in v3 (NSM1 dependency)
    ├── light.py                   # new in v3 (NSM1 dependency)
    └── partitioning.py            # new in v3 (NSM1 dependency)
```

As more components migrate from v2-overlay to v3-native, the imports in `processes/__init__.py` are progressively replaced. When every component is v3-native, the overlay imports are deleted in a final cleanup PR and v3 stands alone.

### Why thin overlay rather than full copy

- Day-one code duplication is minimized.
- The diff between v2 and v3 in version control reflects only the merge work, not file moves.
- LimnoTech's continued v2 work (if any) before v3 is mature can be pulled forward into v3 by replacing the overlay imports incrementally.
- Reduces the merge-conflict surface if v2 and v3 evolve in parallel during the transition.

### Package registration

v3 is registered as a separate editable Python package via the modules repo's `pyproject.toml`. The pixi `dev` environment installs both v2 and v3 editable, so both are importable in the same shell.

The existing `pyproject.toml` entry already installs both `clearwater_modules` (v1) and `clearwater_modules_v2` from the same package (the `path = "."` entry covers both subpackages). Adding `clearwater_modules_v3` requires no new pixi entry — it is part of the same package.

### Integrator pattern (applies to all Process classes in v3)

The v2 NSM1 multiplicative-integrator bug was a symptom of an unclear contract for what `Process.run` should do. v3 establishes the contract explicitly:

1. Each `Process.run` reads its state variables from the registry at the current time.
2. Each `Process.run` computes net rate of change for each state variable, with units of `[state] / second` (additively combining sources and sinks).
3. Each `Process.run` applies the rate to the state via Forward Euler: `state_new = state_old + rate * self.time_step.total_seconds()`.
4. Each `Process.run` writes the updated state back to the registry via `set_at_time`.
5. Negative-state guards (e.g., concentration cannot go below zero) are applied with `xr.where(state_new < 0, 0, state_new)`.

This is the corrected version of what v2 NSM1 attempted. The integrator is the same shape across all Process classes; only the rate computation differs.

---

## 5. Migration Strategy

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

### Migration tables

Each module spec includes a migration table mapping v1/v2 imports to v3 imports for that module's downstream users.

---

## 6. Coordination with LimnoTech (umbrella)

### Communications plan

1. v3 development happens privately in the streaming repo until ready to share. The current LimnoTech-visible state on `memory-refactor-pytestUpdate` includes the original v3 design spec from 2026-05-04 (commit `b21e973`); LimnoTech has been notified of the v3 plan via that commit.
2. When v3 1.0.0 ships, the PR back to `memory-refactor-pytestUpdate` includes all the design specs, the implementation, and the test results.
3. Anthony's pending email response to the coordination questions (sent 2026-05-04) may inform timing; v3 development proceeds in parallel and adapts as his reply arrives.

### What requires LimnoTech input vs. what can proceed without

Module specs document the items that benefit from LimnoTech input. Items that can proceed without input are: pure ports of v1 fixes, mechanical translation of v1 kinetics into v3 framework, test infrastructure, and shared physics primitives.

---

## 7. Risks (umbrella)

| Risk | Probability | Impact | Mitigation |
|---|---|---|---|
| Permanent three-way split: v1, v2, v3 all maintained indefinitely | Low (with explicit retirement plan) | High (compounds maintenance burden) | Section 5 retirement plan; commit to v2 freeze when v3 0.1.0 ships |
| LimnoTech objects to v3 as a directory | Low (consistent with their own v2 precedent) | High (requires re-planning) | Frame as a working merge proposal in the conversation with Anthony; v3 is reversible |
| TSM and NSM1 sessions diverge in framework conventions | Medium | Medium (rework when integrating) | Both sessions follow this architecture spec; periodic check-ins by author across both branches; final integration on a `v3-integration` branch |
| Streaming-repo pixi env diverges from modules-repo pixi env | Low | Low | `pixi.lock` ensures reproducibility; both envs install from same `pyproject.toml` |

Module-specific risks are documented in each module's spec.

---

## Appendix A: Pixi / install verification

After scaffolding the v3 package:

```bash
cd /Users/todd/GitHub/ecohydrology/ClearWater-modules-streaming
pixi run -e dev python -c "import clearwater_modules_v3; print(clearwater_modules_v3.__file__)"
```

Expected: prints the path to `src/clearwater_modules_v3/__init__.py`.

---

## Appendix B: Cross-document map

| Topic | Lives in |
|---|---|
| Quick start, env setup, branch conventions | This spec, Section 1 |
| Why v3 exists, what it combines | This spec, Section 2 |
| Goals/non-goals that apply to all v3 | This spec, Section 3 |
| Module-specific goals/non-goals | Each module spec |
| Thin-overlay package layout | This spec, Section 4 |
| Integrator-pattern contract | This spec, Section 4 |
| v2/v1 retirement plan | This spec, Section 5 |
| Module-specific migration tables | Each module spec |
| LimnoTech communications | This spec, Section 6 |
| Module-specific LimnoTech coordination items | Each module spec |
| Umbrella risks | This spec, Section 7 |
| Module-specific risks | Each module spec |
| Component inventory per module | Each module spec |
| Phased implementation plan per module | Each module spec |
| Test infrastructure per module | Each module spec |
| Performance targets per module | Each module spec |
| Open questions per module | Each module spec |
