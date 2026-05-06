# SSM — Improvement and v3 Migration Plan

**Status:** Draft
**Date:** 2026-05-06
**Scope:** `src/clearwater_modules_v2/processes/sediment/` (the prototype Sediment Simulation Module) and its migration to a v3-native Process under `src/clearwater_modules_v3/processes/sediment/`.

**Read this with the umbrella spec.** The v3 framework contract (`Process`, `ProcessFactory`, kernel-precomputed schedule, registry-level wet-mask, hotstart ordering, retirement plan) is documented in `clearwater_modules_v3_architecture_specification.md`. The original SSM design document (clean-room port from EFDC SEDZLJ + SAND2008-5621) is `ssm_design_spec.md`. This document covers SSM-specific defects and migration only; it does not re-derive the SEDZLJ algorithm.

---

## 1. Executive summary

SSM is a research-prototype Sediment Simulation Module that closes the sediment gap so the ClearWater stack does not depend on HEC-RAS Sediment for bed-state and suspended-sediment dynamics. It is a clean-room Python port of EFDC SEDZLJ on the ClearWater-Riverine mesh, with seven pluggable bedload closures, a Sanford-Maa consolidation model, and a vegetation-cohesion feedback API consumed by ESM. The current v2 prototype is well-architected and scientifically literate, but a code review surfaced four CRITICAL defects (notably: `ssm_bed_elevation` is allocated but never written, so ESM's mode-6 scour mortality reads zeros) and several MAJOR numerical/coverage gaps (no CFL guard on the bedload upwind step; bedload mass double-counted into the suspended `_source` field when the standalone bedload solver is on; the riverine source-injection contract has no consumer). This plan fixes those defects in seven phases, then migrates the cleaned-up module to a v3-native Process. Estimated effort: roughly 2-3 engineer-weeks for Phases 1-6 (the science/correctness work); 1-2 weeks for Phase 7 (mechanical v3 migration once the time-dim collapse from Phase 1 has landed).

---

## 2. Background

### 2.1 What SSM is and why it exists

SSM is a faithful clean-room port of EFDC SEDZLJ (Jones & Lick 2001; Ziegler 2002; SAND2008-5621) onto the ClearWater-Riverine unstructured mesh. It owns:

- A multi-class, multi-layer cohesive + non-cohesive bed state (`g/cm²` mass per layer per class, per-layer τ_crit, layer activation flags, layer age for consolidation).
- An armoring + active-layer reorganization step that mirrors SEDZLJ's `s_sedzlj.f90`.
- Seven peer-reviewed bedload closures behind a `BedloadTransportFunction` Protocol (van Rijn 1984, Wilcock-Crowe 2003, Yang, Wu 2000, Engelund-Hansen 1967, Toffaleti 1968) plus a default; an explicit upwind face-flux step on the bedload mass field.
- Erosion + deposition (Krone-Partheniades for cohesive, SEDflume / Power-Law for non-cohesive; Gessler 1967 probability of deposition).
- A Sanford-Maa age-based consolidation model behind a Protocol contract.
- A vegetation-cohesion feedback API (`erosion.apply_vegetation_cohesion`) that ESM uses to push composite τ_crit back into SSM.

The prototype lives at `src/clearwater_modules_v2/processes/sediment/`:

```
ssm.py            orchestrator (Process subclass; run() and init_process())
contracts.py      single source of truth for variable names, dims, units, dtypes
bed.py            BedState xarray wrapper; layer reorganization; bed-elevation update
bedload.py        BedloadTransportFunction Protocol + 7 closures; upwind step
erosion.py        Power-law / SEDflume erosion; vegetation-cohesion API
deposition.py     Krone-Partheniades + Gessler probability of deposition
settling.py       Cheng 1997 settling velocity
shear.py          τ_b drivers (current-only, RAS-imported, wave-current stub)
armoring.py       D50 mass-weighted average; τ_crit interpolation
consolidation.py  Sanford-Maa age model behind Protocol contract
classes.py        SedimentClassRegistry (n_class fixed at construction)
coupling.py       Riverine bridge helpers
io/               SEDflume + CSV/YAML loaders
```

A v3 of the ClearWater modules has landed at `src/clearwater_modules_v3/processes/`. v3 introduces a kernel-precomputed schedule, registry-level wet-mask, and an opt-in hotstart contract. Sediment has not yet been migrated to v3-native and currently re-exports v2 unchanged — TSM, NSM1, and the riverine constituent solver have moved over, but SSM has not.

### 2.2 How ESM consumes SSM

ESM (`src/esm/`) reads five SSM bed-state fields each step for excess-shear scour mortality (mode 6) and vegetation biostabilization:

- `ssm_bed_change` — per-step Δz (m); used for burial mortality.
- `ssm_bed_elevation` — absolute bed surface (m); reasoning about absolute z.
- `ssm_bed_d50_surface` — surface-layer D50 (μm); biostabilization grain-size logic.
- `ssm_bed_shear_stress` — applied τ_b (Pa); the excess-shear forcing for mode 6.
- `ssm_bed_critical_shear_stress` — surface τ_crit (Pa); the threshold for mode 6.

Read site: `ClearWater-modules-phase2-ESM-streaming/src/esm/io/clearwater_interface.py:383-397`. ESM treats two of the five fields (`bed_elevation`, `d50_surface`) as nominally available even though, today, they are populated with zeros by SSM (see C1 and C2). ESM-side mode-6 wiring is at `src/esm/model.py:811-889` and biostabilization at `src/esm/processes/biostabilization.py`.

### 2.3 Verdict

SSM is well-architected and scientifically literate. The contracts module, the pluggable bedload Protocol, the opt-in consolidation Protocol, the mass-conservation invariant in `reorganize_active_layer`, and the clean-room license/provenance docstrings are all good design and should be preserved. SSM is **not** production-ready: the four CRITICAL defects below mean that two of ESM's five expected inputs are zeros, the riverine source-injection loop is open, and per-layer τ_crit is silently zeroed for any time step beyond t=0 under the current time-indexed bed-state schema. Once those are fixed, the v3 migration is mechanical: the prototype already lives on an xarray Dataset, `init_process(model, registry)` and `run(time, registry)` already match the v3 signatures, and there are no hidden globals or scalar-to-per-cell conversions. The C3 fix (collapse time dim on bed state) **is** the v3 migration shape change, which is why it leads the phasing.

---

## 3. Findings

Severity follows the v3 review convention (`design/clearwater_modules_v3_review_findings.md`): CRITICAL = correctness or coupling defect that produces wrong outputs in default configurations; MAJOR = numerical robustness or coverage gap that can produce wrong outputs in plausible configurations; MINOR = housekeeping that does not affect correctness today.

### 3.1 CRITICAL

**C1. `ssm_bed_elevation` is allocated but never written; ESM consumes zeros.**

- `src/clearwater_modules_v2/processes/sediment/contracts.py:74` declares `VAR_BED_ELEVATION = "ssm_bed_elevation"` and `contracts.py:159-160` adds it to `BED_STATE_SPECS`.
- `src/clearwater_modules_v2/processes/sediment/bed.py:630-712` (`update_bed_elevation`) writes `VAR_BED_LAYER_THICKNESS`, `VAR_BED_TOTAL_THICKNESS`, `VAR_BED_CHANGE`, and `VAR_BED_CUMULATIVE_CHANGE` — but never `VAR_BED_ELEVATION`.
- ESM's `clearwater_interface.py:385` reads `ssm_bed_elevation` directly. Mode-6 scour mortality and any logic that reasons about absolute bed elevation is currently consuming zeros.
- **Fix.** In `update_bed_elevation`, compute `bed_elev = bed_datum + cumulative_bed_change` and write into `VAR_BED_ELEVATION`. The bed-datum reference needs a config knob: typical practice is to read the RAS cell minimum elevation (`min_face_elev` if available) and add the SSM bed thickness. See open question Q1.

**C2. `ssm_bed_d50_surface` is computed but never written back to the mesh.**

- `src/clearwater_modules_v2/processes/sediment/ssm.py:759` computes `d50_surface` via `armoring.compute_d50_avg`; the value is consumed locally but never written to `mesh[VAR_BED_D50_SURFACE]`.
- ESM's `clearwater_interface.py:386` reads `ssm_bed_d50_surface` for biostabilization and grain-size logic. It receives zeros.
- **Fix.** Add a write at the end of `SSM.run` using the same time-resolving pattern as `_write_tau_diagnostics` (`ssm.py:1187`).

**C3. Per-layer `taucor` (and other layer state) is set at t=0 only; never propagated forward.**

- `src/clearwater_modules_v2/processes/sediment/bed.py:327` writes per-layer `tau_crit` at t=0.
- `BedState.set_layer_taucrit_at` exists but is invoked only in the conservation tests' carry-forward helper at `tests/sediment/test_conservation.py:223-226` (current line 224).
- `SSM.run` reads `bed.layer_taucrit_at(time)` but never writes it. Same gap for `layer_mass`, `class_fraction`, `layer_active`, `layer_age`.
- For label-based time stepping or any integer-index step beyond 0, the value at later time slots is zero, making `tau_slln` zero and tripping `tau > tau_slln` for every cell.
- **Fix (preferred).** Collapse all bed-state DataArrays to a single time slot — no `time` dim on bed state, mirroring the v3 registry-backed pattern where state is overwritten in place. This *also* dramatically reduces memory and aligns with v3.
- **Fix (alternative).** Have `SSM.run` explicitly carry the previous time slot's state forward at the start of each step.

**C4. Riverine source-injection contract is one-sided.**

- `src/clearwater_modules_v2/processes/sediment/ssm.py:1024-1028` writes `mesh[f"{cls.suspended_var}_source"]` as g/cm² per step.
- No documented unit attrs, no sign convention attrs, no consumer code in `riverine.py`.
- The conservation tests close the loop by reading the `_source` field directly (`tests/sediment/test_conservation.py:198-209`), but they do not exercise Riverine's actual transport step.
- Until Riverine's constituent solver reads these fields, the suspended-sediment loop is open.
- **Fix.** Implement Riverine-side ingestion as a constituent source: `kg/m³/s = g/cm²/step × cell_area_cm² × 0.001 / volume_m³ / dt`. Document units explicitly in the variable's `attrs`.

### 3.2 MAJOR

**M1. `T_act` formula uses cgs density divisor 10000 without a clear unit-conversion note.**

- `src/clearwater_modules_v2/processes/sediment/bed.py:434`: `t_act = tactm * d50 * factor * bd1 / 10000.0`. The `/10000` is the μm→cm conversion folded into the formula.
- **Fix.** Pull the constant into `contracts.py` as `_UM_TO_CM = 1.0e-4`; rewrite as `t_act = tactm * (d50 * _UM_TO_CM) * factor * bd1`; add a comment.

**M2. Bedload upwind step has no CFL guard.**

- `src/clearwater_modules_v2/processes/sediment/bedload.py:649-664` performs explicit upwind face-flux update. No CFL check on `dt × |u_eff_max| / Δx`.
- For Albany-scale meshes (~4-8 m cells) and typical bedload velocities, CFL > 1 is plausible at hourly time steps. The closed-domain mass conservation test passes only because `np.maximum(new_cbl, 0.0)` (line 669) clips negatives.
- **Fix.** Compute `cfl = dt * np.max(np.abs(u_edge)) / np.min(edge_length_cm)` per class; sub-cycle when `cfl > 1`, or warn. Emit a CFL diagnostic field.

**M3. Negative concentrations from explicit upwind are clipped, not flagged.**

- `src/clearwater_modules_v2/processes/sediment/bedload.py:669` clips `new_cbl = np.maximum(new_cbl, 0.0)` silently.
- **Fix.** Track per-step `clipped_mass = np.sum(np.maximum(-new_cbl_pre_clip, 0.0))` and write to a `ssm_bedload_clip_diagnostic` field, or warn above a threshold.

**M4. Suspended concentration read assumes RAS bottom-cell, not depth-averaged.**

- `src/clearwater_modules_v2/processes/sediment/ssm.py:1107-1135` uses depth-averaged Riverine constituent concentration as a proxy for near-bed concentration in Krone/Gessler deposition.
- **Fix.** Document the assumption in the docstring. (No behavior change required.)

**M5. Bedload mass is double-counted with suspended source when bedload solver is enabled.**

- `src/clearwater_modules_v2/processes/sediment/ssm.py:1024-1028` injects `net_per_class = erosion_per_class - deposition_arr` into the suspended-class source field for *all* classes including sand fraction in saltation.
- The `BedloadStandaloneExplicit` step at `bedload.py:550` (call site `ssm.py:1016`) does not subtract its own erosion.
- The closed-domain conservation test runs with `bedload_solver="off"` (`tests/sediment/test_conservation.py:153`), so the bedload-on path is uncovered. **This is a real mass-conservation defect under bedload-on.**
- **Fix.** Partition `erosion_per_class` for sand classes between bedload and suspended fractions (use van Rijn suspension criterion `u*/w_s` to split). Only inject the suspended fraction into the `_source` field. Alternatively, exclude bedload-eligible classes from `_source` injection entirely with documentation. See open question Q2.

**M6. `bed_change` from total-thickness delta loses internal-reorganization signal.**

- `src/clearwater_modules_v2/processes/sediment/bed.py:687`: `bed_change = (total_thickness - prev_total)`. Internal layer-mass redistributions register as zero change, which is fine for net-mass reporting but blunts the ESM biostabilization signal.
- **Fix.** Document as a deliberate simplification, or compute `bed_change` from the layer-1 mass delta only.

### 3.3 MINOR

| ID | Location | Defect | Fix |
|---|---|---|---|
| m1 | `ssm.py:222` | `Process.variables = []`; `Process.validate()` cannot catch missing Riverine constituents at config time | Populate `variables` after `init_process` (post-class-resolution) or override `validate` |
| m2 | `ssm.py:84-108` and `ssm.py:146-156` | Per-cell Python loop in `_surface_layer_index`, `_surface_class_fraction`, `_surface_taucrit_pa` | Vectorize with `np.argmax(layer_active != LAYER_ABSENT, axis=1)` |
| m3 | `contracts.py:138-169` (`BED_STATE_SPECS`) | Bed-state DataArrays are float32 with float64 compute → conversion every read; round-trip drift could approach `_MASS_CONSERVATION_TOL = 1e-5` over thousands of steps | Switch `BED_STATE_SPECS` to float64 |
| m4 | `ssm.py:1007` | SSM `time_step` independence from Riverine cadence is undocumented | Explicitly require equality, or refactor bedload solver to internal time axis |
| m5 | `deposition.py:32` (`_TAU_FLOOR_PA = 1.0e-12`) | Floor is below numerical noise | Raise to 1e-9 or 1e-8 |
| m6 | `erosion.py:401-403` (`PowerLawErosionModel`) | Constructor does not validate `n_top > 0` | Add `__post_init__` validation |
| m7 | `settling.py:62` (`cheng_1997_settling_velocity`) | Rejects `D50 = 0`; Stokes-regime accuracy below ~0.1 μm undocumented | Add docstring note on accuracy floor |
| m8 | `bed.py:705-711` (`update_bed_elevation`) | Advances `layer_age` by `dt` only when `dt > 0` | Tighten invariant — assert dt ≥ 0 and document the no-op guard |

### 3.4 What to preserve through v3 migration

The following are good design and must not be broken by the cleanup or migration:

1. **`contracts.py` as the single source of truth** for variable names, dims, units, dtypes, Zarr chunks. v3's registry should consume `BED_STATE_SPECS` directly.
2. **The pluggable `BedloadTransportFunction` Protocol** (`bedload.py:939`) with seven peer-reviewed formulas (`bedload.py:974-1615`).
3. **The opt-in Sanford-Maa consolidation model** with `__post_init__` validation and a clean Protocol contract (`consolidation.py`).
4. **Clean-room license/provenance docstrings** in every algorithmic module.
5. **The mass-conservation invariant assert** in `reorganize_active_layer` (`bed.py:567-574`).
6. **The `bind_mesh` test entry point** that decouples SSM from full Model construction.
7. **`BED_STATE_SPECS` schema with explicit `VarSpec` records** — adopt this schema shape for all v3 process schemas (`contracts.py:125-169`).
8. **Vegetation-cohesion feedback API** (`erosion.apply_vegetation_cohesion`) — sign convention is correct on the ESM side.

---

## 4. v3 Migration Assessment

### 4.1 What maps cleanly

The prototype maps cleanly onto v3's `Process` pattern because:

- State already lives on an xarray Dataset (the mesh).
- `init_process(model, registry)` and `run(time, registry)` signatures already match v3's contract (`clearwater_modules_v3/processes/base.py:1-80`).
- No hidden globals; no scalar-to-per-cell conversions are needed.
- The Process subclass already takes its `time_step` as a constructor argument (Process pattern, not class attribute).
- Initial conditions flow through `init_process` rather than through `__init__` side effects.
- `BED_STATE_SPECS` already describes every owned variable with explicit dims, dtype, units, and role — this is the same shape the v3 registry expects for a typed write.

### 4.2 Friction points

- **Time-dimensioned bed-state arrays must collapse to a single slot.** v3 state is overwritten in place via the registry. This is the C3 fix and must precede any other work.
- **`from_hotstart` is not implemented.** v3's base spec at `clearwater_modules_v3/processes/base.py:30-65` makes hotstart opt-in (each process either defines `from_hotstart` or omits it; the Model invokes it conditionally after `init_process` per the M5 ordering contract). SSM bed state is per-cell-per-layer per-class — substantial state to checkpoint, but well-defined. This is "do later," not blocker.
- **The `_source` injection mechanism needs a Riverine consumer.** This is the C4 fix. v3 should formalize it as a typed registry write rather than an undocumented mesh-scratch field — i.e., add a `SUSPENDED_SOURCE_SPECS` to `contracts.py` mirroring `BED_STATE_SPECS`.

### 4.3 What changes shape under v3

The v2 prototype uses the `clearwater_modules_v2.processes.base.Process` ABC (`base.py:14`); v3 currently re-exports v2 `Process` unchanged from `clearwater_modules_v3/processes/base.py:79`. So the import switch is single-line; the *behavioral* differences are:

- v3 exposes the kernel-precomputed schedule and registry wet-mask; SSM's `run` should accept the wet-mask as registry-supplied rather than scanning every cell.
- v3 makes hotstart opt-in: SSM can defer hotstart entirely.
- v3 writes go through the registry, not through direct mesh DataArray mutation. SSM's `_write_per_face_source`, `_write_per_class_flux`, `_write_tau_diagnostics`, and the new C1/C2 writes should be refactored to registry calls.

---

## 5. Phasing

Each phase has scope, files touched, acceptance criteria, effort, and dependencies. Phases 1-6 are sequenced because each layers on the previous one's machinery; Phase 7 is gated on Phases 1-3 having landed.

### Phase 1 — Collapse time dim on bed state (= the v3 migration shape change)

**Scope.** Drop `DIM_TIME` from every `BED_STATE_SPECS` entry whose role is `bed_state`. Update `BedState` accessors (`layer_mass_at`, `layer_taucrit_at`, etc.) to ignore the time argument. Update the conservation tests' carry-forward helper to no-op. This is foundational because C1, C2, and m3 all touch the same machinery, and this *is* the shape change required for v3.

**Files touched.**
- `contracts.py` — drop `DIM_TIME` from bed_state `VarSpec.dims`; bump `BED_STATE_SPECS` accordingly.
- `bed.py` — flatten `set_layer_*_at`/`layer_*_at` to direct DataArray mutation; drop `_assign_time` for bed_state vars; remove the `idx == 0` guard in `update_bed_elevation` and replace with a stored `prev_total_thickness` attribute on `BedState`.
- `ssm.py` — drop `_resolve_time_index` calls for bed-state writes; keep them for diagnostics that remain time-dimensioned.
- `tests/sediment/test_conservation.py:223-226` (and equivalents at lines 335 and 403) — remove the carry-forward helper.

**Acceptance criteria.**
- All 190 currently-passing sediment tests still pass after the schema change.
- `bed.layer_taucrit_at(time)` returns the per-layer τ_crit from the *current* bed state regardless of `time`.
- A new test (`test_bed_state_persists_across_steps`) confirms that running 10 SSM steps and reading per-layer τ_crit at each step returns nonzero values consistent with the SEDflume initialization.

**Effort.** Medium (1-2 days). The change is mechanical but touches every bed-state read/write site.

**Dependencies.** None.

### Phase 2 — Write `bed_elevation` and `d50_surface` to the mesh (C1, C2)

**Scope.** Tiny patches that immediately unblock ESM mode-6 testing.

**Files touched.**
- `bed.py:630-712` (`update_bed_elevation`) — compute and write `VAR_BED_ELEVATION = bed_datum + cumulative_bed_change`. Add `bed_datum` to `BedState`'s constructor with a sensible default (zero) and a config knob.
- `ssm.py:759` site — after `compute_d50_avg`, add `_write_per_face(mesh, VAR_BED_D50_SURFACE, time, d50_surface)` mirroring the `_write_tau_diagnostics` pattern at `ssm.py:1187`.
- `ssm.py` constructor — accept `bed_datum` parameter (scalar or per-cell) with default sourced from RAS HDF or a zero array.

**Acceptance criteria.**
- A new round-trip test (`tests/sediment/test_esm_bed_state_roundtrip.py`) confirms ESM's `clearwater_interface.SedimentState.bed_elevation` reads back nonzero values matching `bed_datum + cumulative_bed_change` after one SSM step.
- Same test asserts `bed_d50_surface` matches the value `armoring.compute_d50_avg` computes for the surface layer.
- Existing tests still pass.

**Effort.** Small (half-day).

**Dependencies.** Phase 1 (the writes need to target single-slot DataArrays).

### Phase 3 — Riverine source ingestion (C4)

**Scope.** Implement the Riverine consumer for `f"{cls.suspended_var}_source"`. Document units in `BED_STATE_SPECS` attrs (or in a new `SUSPENDED_SOURCE_SPECS`). This requires Riverine-side work and is the gate for closed-loop suspended-sediment transport.

**Files touched.**
- `src/clearwater_modules_v2/processes/riverine.py` — read `f"{constituent}_source"` (g/cm²/step) before each constituent solve; convert to kg/m³/s using cell area and volume; add to the constituent right-hand-side.
- `contracts.py` — add `SUSPENDED_SOURCE_SPECS` with explicit `attrs` carrying units, sign convention ("+ = bed → water-column"), and the cadence requirement (one write per SSM step).
- `ssm.py:1024-1028` site — set `attrs` on the source field.

**Acceptance criteria.**
- A new closed-loop conservation test (`tests/sediment/test_riverine_loop_conservation.py`) runs SSM + Riverine on a closed domain and asserts total mass (bed + water column) is conserved to within `_MASS_CONSERVATION_TOL` over 100 steps.
- Existing 190 tests still pass.

**Effort.** Medium (2-3 days). The Riverine integration is the friction.

**Dependencies.** Phase 1 (single-slot bed state simplifies the source-write contract).

### Phase 4 — Bedload double-counting (M5)

**Scope.** Resolve the bedload/suspended mass-conservation defect. Write a bedload-on conservation test first; it will fail and force a clean partition.

**Files touched.**
- `tests/sediment/test_conservation.py` — new test class that runs the same closed-domain scenario with `bedload_solver="standalone"` and asserts (bed + bedload + suspended) total mass conservation.
- `ssm.py:1024-1028` — partition `erosion_per_class[s]` for bedload-eligible sand classes by van Rijn's `u*/w_s` suspension criterion; inject only the suspended fraction into the `_source` field. (Or, alternative: exclude bedload-eligible classes entirely with a clear documentation note. Decide via Q2.)
- `bedload.py` — if partitioning is chosen, accept the "fraction-to-bedload" array from SSM rather than re-computing it inside `BedloadStandaloneExplicit.step`.

**Acceptance criteria.**
- New bedload-on conservation test passes.
- `bedload_solver="off"` test still passes (no regression on the cohesive-only path).
- Documentation note in `ssm.py` docstring explains the partitioning rule and cites van Rijn's suspension criterion.

**Effort.** Medium (3-5 days). The partitioning rule is the design call; the implementation is small once decided.

**Dependencies.** Phase 3 (closed-loop test infrastructure).

### Phase 5 — CFL diagnostic and bedload clip reporting (M2, M3)

**Scope.** Numerical-robustness instrumentation. No behavior change at default config; new diagnostic fields.

**Files touched.**
- `bedload.py:649-669` — compute per-step `cfl = dt × max|u_edge| / min(edge_length_cm)`; warn when `cfl > 1`; emit `ssm_bedload_cfl` diagnostic field. Track `clipped_mass = sum(max(-new_cbl_pre_clip, 0))`; emit `ssm_bedload_clip_diagnostic`.
- `contracts.py` — add `VAR_BEDLOAD_CFL` and `VAR_BEDLOAD_CLIP` to a new `DIAGNOSTIC_SPECS` tuple.

**Acceptance criteria.**
- A new test (`tests/sediment/test_bedload_cfl_diagnostic.py`) constructs a high-velocity scenario, runs one step, and asserts the CFL diagnostic is > 1 and a `UserWarning` was raised.
- Default Albany-scale scenario emits CFL ≤ 1 and zero clipped mass.

**Effort.** Small (1-2 days).

**Dependencies.** Phase 1 (single-slot diagnostic writes simplify the pattern).

### Phase 6 — Vectorize surface-layer scan, populate `variables`, switch to float64 (m1, m2, m3)

**Scope.** Housekeeping for v3 migration. No behavior change; performance and correctness-at-scale only.

**Files touched.**
- `ssm.py:84-108` and `ssm.py:146-156` — replace per-cell loops with `np.argmax(layer_active != LAYER_ABSENT, axis=1)` then gather. Validate against existing tests.
- `ssm.py:222` (`Process.variables`) — populate after `init_process` once classes are known; or override `Process.validate` to walk `registry_classes`.
- `contracts.py:138-169` (`BED_STATE_SPECS`) — flip dtype from `"float32"` to `"float64"` for every bed_state record. Diagnostic-role records can stay float32.
- `deposition.py:32` (m5) — raise `_TAU_FLOOR_PA` to 1e-9.
- `erosion.py:401-403` (m6) — add `__post_init__` for `PowerLawErosionModel`.
- `settling.py:62` (m7) — docstring note.
- `bed.py:705-711` (m8) — assert `dt_seconds >= 0`.

**Acceptance criteria.**
- All existing tests still pass.
- `Process.validate()` raises a clear error when a configured suspended class has no corresponding constituent in the Riverine registry.
- Memory profile of bed-state Dataset shows the expected ~2× growth (float32 → float64).
- A microbenchmark on the Corvallis-Santiam-Albany mesh shows surface-layer scan dropped from O(nface) Python iterations to a single vectorized call.

**Effort.** Small-Medium (2 days).

**Dependencies.** Phases 1-3 (baseline tests must already be green at single-slot bed state).

### Phase 7 — v3-native migration

**Scope.** Wrap the cleaned prototype in the v3 Process pattern. Route writes through the registry rather than direct mesh mutation. Opt into hotstart hooks if needed (deferred per §4.2).

**Files touched.**
- New: `src/clearwater_modules_v3/processes/sediment/` — module package mirroring the v2 layout. Re-export `Process` and `ProcessFactory` from v3 base. Subclass `SSM` from `clearwater_modules_v3.processes.base.Process`.
- New: `src/clearwater_modules_v3/processes/sediment/contracts.py` — typed registry specs derived from `BED_STATE_SPECS`. The registry write API replaces the v2 mesh-mutation helpers (`_write_per_face_source`, `_write_per_class_flux`, `_write_tau_diagnostics`).
- New: `src/clearwater_modules_v3/processes/sediment/__init__.py` — registers the v3 SSM class with `ProcessFactory`.
- Existing: `src/clearwater_modules_v2/processes/sediment/` — retained verbatim (per the v3 retirement plan in the umbrella spec).

**Acceptance criteria.**
- v3 SSM passes the v3 NSM1/TSM-equivalent regression suite on the Sumwere Creek and Corvallis-Santiam-Albany meshes.
- Coupled v3 SSM + v3 Riverine + v3 ESM smoke test runs end-to-end.
- v3 SSM matches v2 SSM outputs to within float64 round-off on a regression mesh.

**Effort.** Medium (1-2 weeks). Mechanical once Phases 1-3 are in.

**Dependencies.** Phases 1-3 must have landed. Phases 4-6 should land for parity; if scheduling forces v3 migration before Phases 4-6, port the defects forward and track them as v3-side work.

### Deferred (not in this plan)

- Wave-current shear closure (currently a stub at `shear.py:398-419`).
- Full Gibson finite-strain consolidation with porosity evolution (`ssm_design_spec.md` §5.10 NG1).
- MMS / analytical-solution convergence study for advection-diffusion of suspended sediment.
- HEC-RAS Sediment or MT3DMS regression comparison.
- Hotstart for v3 SSM (opt-in; defer until first user request).

---

## 6. Test plan additions

| Test file | Asserts |
|---|---|
| `tests/sediment/test_bed_state_persists_across_steps.py` (Phase 1) | `bed.layer_taucrit_at` returns nonzero values at every step after t=0; mass-conservation invariant holds across 10 steps with the single-slot schema. |
| `tests/sediment/test_esm_bed_state_roundtrip.py` (Phase 2) | After one SSM step on a synthetic mesh, ESM's `clearwater_interface.read_sediment_state` returns nonzero `bed_elevation` and `bed_d50_surface` matching SSM's internal values to float64 tolerance. |
| `tests/sediment/test_riverine_loop_conservation.py` (Phase 3) | Closed-loop SSM + Riverine on a closed domain conserves (bed + water-column) total mass to within `_MASS_CONSERVATION_TOL` over 100 steps. |
| `tests/sediment/test_conservation_bedload_on.py` (Phase 4) | Closed-domain conservation test with `bedload_solver="standalone"` enabled — currently uncovered. |
| `tests/sediment/test_bedload_cfl_diagnostic.py` (Phase 5) | High-velocity scenario emits `ssm_bedload_cfl > 1` and a `UserWarning`; default scenario emits zero clipped mass. |
| `tests/sediment/test_v3_regression.py` (Phase 7) | v3 SSM matches v2 SSM outputs to within float64 round-off on a regression mesh. |

Existing coverage as of this writing: 190/191 sediment tests pass (the one failure is a Riverine API drift, not a SSM defect). The existing suite includes closed-domain conservation (`test_conservation.py`), the SAND2008-5621 reference dataset (`test_sand2008_reference.py`), per-module unit tests, and a coupled smoke test (`test_ssm_riverine_smoke.py`).

Deferred test work (not blocking this plan): MMS or analytical-solution convergence for suspended-sediment advection-diffusion; regression against an HEC-RAS Sediment or MT3DMS reference run.

---

## 7. Open questions

**Q1. Bed-datum reference for C1.** Should `bed_datum` come from the RAS HDF `min_face_elev` per cell, from a user-supplied scalar or array, or from a config knob? Recommendation: read `min_face_elev` from RAS HDF when available (`src/esm/io/hecras_reader.py` already loads it); fall back to user-supplied scalar; default to zero with a warning. The `BedState` constructor takes a `bed_datum: float | np.ndarray | None = None` argument and the SSM constructor wires the loader.

**Q2. M5 fix strategy: partition or exclude.** Two options:
- *Partition.* Use van Rijn's suspension criterion `u*/w_s` to split sand-class erosion between bedload and suspended fractions. More physically correct; requires deciding the split in `ssm.run` and passing it to `BedloadStandaloneExplicit.step`.
- *Exclude.* Drop bedload-eligible classes from the `_source` injection entirely with a documentation note. Simpler; loses the "suspension above the threshold" physics but is conservative.

Recommendation: start with *exclude* (simpler, conservative, and the test in Phase 4 will catch any regression); upgrade to *partition* in a follow-up when the suspension criterion is validated on a regression mesh.

**Q3. Phase 7 timing.** Migrate to v3-native in the same effort as the defect fixes, or as a follow-up? Recommendation: defect fixes first (Phases 1-6) on the v2 prototype, then v3 migration as a follow-up. Rationale: Phases 1-3 are the v3 shape change; Phases 4-6 stabilize the prototype before duplicating it; Phase 7 is then mechanical and low-risk.

**Q4. Float64 cost.** Doubling the bed-state Dataset size (m3) is meaningful on Albany-scale meshes (n_face ~ 100k, n_layer = 8, n_class = 7 → 5.6M elements per layer var). Confirm with a memory profile that the doubled Dataset still fits in target deployment envelopes before committing to Phase 6 m3.

**Q5. Hotstart in v3.** Defer or implement in Phase 7? Recommendation: defer. Per `clearwater_modules_v3/processes/base.py:30-65`, hotstart is opt-in; SSM has well-defined state but no user has yet asked for hotstart support. Track as Phase 7+1.
