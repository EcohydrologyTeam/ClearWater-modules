# Phase 0 Close-Out — v3 NSM1 Pattern Alignment

**Date:** 2026-05-13
**Branch:** `streaming`
**Baseline commit:** `186b5c4`
**Spec reference:** `design/clearwater_modules_v3_nsm1_pattern_alignment_specification.md` §6 Phase 0

Phase 0 is complete. This document records the deliverables and resolution status of each Phase 0.x sub-task.

---

## Deliverables

| Sub-task | Deliverable | Location |
|---|---|---|
| 0.1 | Full test suite baseline (805 passed, 2 xfailed, 71.36 s) | `tests/v3/nsm1/baseline/baseline_junit_full_186b5c4.xml`, `baseline_pytest_summary_186b5c4.json`, `baseline_pytest_full_186b5c4.txt` |
| 0.2 | Bit-identical 4,320-substep coupled trajectory NetCDF (20 vars × 5 cells × 4,321 substep indices) | `tests/v3/nsm1/baseline/baseline_coupled_trajectory_186b5c4.nc` |
| 0.3 | Tier 1 conservation baseline (40 passed at `rtol=1e-12, clip_events == {}`) | `tests/v3/nsm1/baseline/baseline_tier1_junit_186b5c4.xml`, `baseline_tier1_summary_186b5c4.json` |
| 0.4 | Dependency pin (329 packages via `pixi list`) | `tests/v3/nsm1/baseline/baseline_pixi_list_186b5c4.txt` |
| 0.5 | Appendix A diff and full registry-diagnostics catalog (70 names) | `design/clearwater_modules_v3_nsm1_appendix_a_diff.md` |
| 0.6 | Q1 / Q2 resolutions in `utils/numerics.py` + 14 new tests | `src/clearwater_modules_v3/utils/numerics.py`, `src/clearwater_modules_v3/model.py`, `tests/v3/test_clip_negative_state_phase0.py` |
| 0.x (bonus) | Per-phase parity-check runner (the §11.2 enforcement script) | `tests/v3/nsm1/baseline/check_baseline_parity.py` |

---

## Q1 — Step-index source resolution

**Decision:** Option (d) — store the substep index on `Diagnostics` itself, updated by `Model`'s substep loop. Picked over option (b) (`init_process` captures `self._step_ref = model.current_step_ref`) because it is *zero* extra surface area on `Process`: existing call sites pass `self.diagnostics` already, and the function reads `diagnostics.current_step` as the step default. No new `init_process` field; no signature changes.

**Implementation:**

- `Diagnostics.current_step: int = -1` (sentinel "no step set" → logged as `step=None` for back-compat with the 1.0.0 convention).
- `Model.__process_loop_full` and `Model.__process_loop_chunked` write `self.diagnostics.current_step = step_index` at the top of every substep iteration (before any `process.run` call so even a `clip_negative_state` invoked from the first Process sees the correct value).
- `clip_negative_state`'s `step=` kwarg is preserved for explicit overrides; when omitted (or `None`), the function reads `diagnostics.current_step`. `-1` is normalized to `None` in the log record.

**Verified by:** `tests/v3/test_clip_negative_state_phase0.py::test_q1_*` (4 tests).

## Q2 — Graceful no-op when diagnostics is None

**Decision:** Make `diagnostics` an optional parameter (`Diagnostics | None = None`). Also harden the function to accept `np.ndarray` and Python scalar inputs in addition to `xr.DataArray`, returning the same container type as the input. This lets call sites drop the `isinstance(... DataArray) and self.diagnostics is not None` guard entirely.

**Implementation:**

- Signature: `clip_negative_state(state: ArrayLike, name: str, diagnostics: Diagnostics | None = None, step: int | None = None)`.
- `xr.DataArray` in → `xr.DataArray` out (coords / dims / attrs / name preserved).
- `np.ndarray` in → `np.ndarray` out.
- Python scalar in → scalar out (preserving the input's numeric type).
- `diagnostics is None`: clipping still happens; counting/logging is skipped.

**Verified by:** `tests/v3/test_clip_negative_state_phase0.py::test_q2_*` and `test_back_compat_*` (10 tests).

---

## Zero-regression invariants verified

| Check | Result |
|---|---|
| Full suite (`pytest tests/`) | 819 passed, 2 xfailed (805 + 14 new Phase 0.6 tests; 0 pre-existing regressions). |
| Phase 0 baseline NetCDF parity (`check_baseline_parity.py`) | OK — bit-identical against `baseline_coupled_trajectory_186b5c4.nc`. |
| Tier 1 conservation tests | 40 passed at `rtol=1e-12, clip_events == {}`. |
| xfailed count | 2 (unchanged from baseline; both `test_changed_kbod_20`). |

The Phase 0 prep changes to `utils/numerics.py` and `model.py` are **structural-only**: they add capability (None-diagnostics graceful path, step-from-Diagnostics default) without changing the behaviour of any existing call site. Phase 1 will be the first phase that actually exercises the new paths by removing the per-Process `isinstance / diagnostics-not-None` guards.

---

## What changes for Phase 1

Phase 1's mechanical alignment pass can now safely call `clip_negative_state(state, name, self.diagnostics)` without any guard branches in:

- `Alkalinity.run` (`alkalinity.py:485-490`)
- `DOX.run` (`dox.py:730-735`)
- `N2.run` (`n2.py:367-370`)
- `Nitrogen._clip` (`nitrogen.py:371-376`) — the entire `_clip` helper can be deleted; replace each call site with a direct `clip_negative_state(state, name, self.diagnostics)`.

Step-index attribution is now automatic — none of those call sites needs to thread `step=` through; the Model already populates `diagnostics.current_step` per substep.
