# Phase 0 Baseline Artifacts — v3 NSM1 Pattern Alignment

Captured 2026-05-13 against commit `186b5c4` ("Add v3 NSM1 pattern alignment specification") on branch `streaming`.

These artifacts are the **gold reference** for the zero-regression contract in `design/clearwater_modules_v3_nsm1_pattern_alignment_specification.md` §11. Every per-Process phase commit (Phase 1 through Phase 10) of the pattern-alignment work must reproduce them bit-identically when no `REGISTRY_DIAGNOSTICS` names are pre-registered.

## Files

| File | Purpose |
|---|---|
| `baseline_coupled_trajectory_186b5c4.nc` | 4,320-substep coupled NSM1 demo trajectory; 20 state/forcing variables × 5 cells × 4,321 substep indices (initial condition + 4,320 substeps). Bit-identical reproducibility verified across two consecutive captures. **Load-bearing for §11.2.** |
| `baseline_junit_full_186b5c4.xml` | Full repo test suite JUnit XML (805 passed, 2 xfailed, 0 failed, 71.36 s wall). |
| `baseline_junit_186b5c4.xml` | `tests/v3` subset JUnit XML (392 passed, 0 xfailed, 11.89 s wall). |
| `baseline_pytest_full_186b5c4.txt` | Verbose pytest text output for the full suite. |
| `baseline_pytest_v3_186b5c4.txt` | Verbose pytest text output for `tests/v3` only. |
| `baseline_pytest_summary_186b5c4.json` | Extracted summary: counts, xfailed list, slow tests >1s, wall times. **Load-bearing for §11.4.** |
| `baseline_tier1_junit_186b5c4.xml` | Tier 1 conservation tests JUnit XML (40 passed, 8.69 s wall). |
| `baseline_tier1_summary_186b5c4.json` | Tier 1 summary: per-test status and timing; asserted contract `rtol=1e-12, clip_events == {}`. **Load-bearing for §11.4.** |
| `baseline_pixi_list_186b5c4.txt` | Full conda+pypi dependency pin (329 packages) via `pixi list`. Re-capture this if dependencies change; do not silently re-baseline. |
| `capture_baseline_trajectory.py` | Script that produced `baseline_coupled_trajectory_*.nc`. Run via `pixi run --environment dev python tests/v3/nsm1/baseline/capture_baseline_trajectory.py <commit>`. |
| `check_baseline_parity.py` | Bit-identical parity check used by every subsequent phase commit. Exit code 0 = match, 1 = mismatch (prints offending variables). Run via `pixi run --environment dev python tests/v3/nsm1/baseline/check_baseline_parity.py`. |

## Reproduction

The baseline was captured in the `dev` pixi environment (only `dev` has the editable installs of `clearwater_modules`, `clearwater_data`, and `clearwater_riverine`):

```sh
pixi run --environment dev python tests/v3/nsm1/baseline/capture_baseline_trajectory.py 186b5c4
pixi run --environment dev pytest tests/ --junitxml=tests/v3/nsm1/baseline/baseline_junit_full_<commit>.xml
pixi run --environment dev pytest tests/v3/nsm1 -k tier1 --junitxml=tests/v3/nsm1/baseline/baseline_tier1_junit_<commit>.xml
```

## What "bit-identical" means here

`check_baseline_parity.py` uses `numpy.array_equal` on the raw `float64` values — `rtol=0, atol=0`, no tolerance whatsoever. A single-bit difference in a single cell at a single substep fails the check. This is the strongest invariant the test suite carries and is the only one that catches the operand-reordering / broadcast-shift class of regressions that motivated the §11.6 refactor-discipline rules.

## When to re-baseline

Only when one of the following changes:

- A dependency version (numpy, xarray, netCDF4, or anything upstream that affects floating-point evaluation order).
- A deliberate kinetics change committed outside the pattern-alignment work.
- A change to the demo's default initial conditions or default parameters that is intentional and reviewed.

Re-baselining is a separate, signed-off commit with its own short hash in the filenames. The old baseline files are not overwritten — they remain in the tree so the history of references is auditable.

## Pathogen warnings during capture

`Pathogen` prints two warnings on the first capture about optional registry variables `Solid` and `ap` not being present. These are the once-only `_get_optional` warn-latches; they are part of the baseline behavior and will be reproduced (silently after the first emission) on every parity run.
