# Phase 0 Baseline Artifacts — v3 NSM1 Pattern Alignment

Originally captured 2026-05-13 against commit `186b5c4` ("Add v3 NSM1 pattern alignment specification") on branch `streaming`.

**Active baseline: `b51df71`** (terminal gold-standard baseline, re-baselined 2026-05-16, incorporating the NSM1-CA-1 and NSM1-SCI-N1 alkalinity kinetics fixes — the only trajectory-perturbing changes in the gold-standard gate set; see "Re-baseline log" below). The `624ed7c` (CA-1 only) and `186b5c4` (pre-fix) artifacts are retained in the tree for auditability and are no longer the active reference.

These artifacts are the **gold reference** for the zero-regression contract in `design/clearwater_modules_v3_nsm1_pattern_alignment_specification.md` §11. Every per-Process phase commit (Phase 1 through Phase 10) of the pattern-alignment work must reproduce them bit-identically when no `REGISTRY_DIAGNOSTICS` names are pre-registered.

## Files

| File | Purpose |
|---|---|
| `baseline_coupled_trajectory_b51df71.nc` | **ACTIVE (terminal).** 4,320-substep coupled NSM1 demo trajectory at the NSM1-CA-1 + NSM1-SCI-N1 fixes. Same shape/contract as `186b5c4`; differs only in `alkalinity`. **Load-bearing for §11.2.** |
| `baseline_coupled_trajectory_624ed7c.nc` | Superseded (CA-1 only; encodes the SCI-N1 4× denitrification-alkalinity defect). Retained for auditability. |
| `baseline_coupled_trajectory_186b5c4.nc` | Superseded (pre-CA-1; encodes the raw-weight alkalinity bug). 4,320-substep coupled NSM1 demo trajectory; 20 state/forcing variables × 5 cells × 4,321 substep indices (initial condition + 4,320 substeps). Bit-identical reproducibility verified across two consecutive captures. Retained for auditability. |
| `baseline_junit_full_b51df71.xml` | **ACTIVE (terminal).** Full repo test suite JUnit XML at the CA-1 + SCI-N1 fixes (995 passed, 2 xfailed, 0 failed, 319.15 s wall; +2 vs `624ed7c` from the new SCI-N1 regression tests). |
| `baseline_tier1_junit_b51df71.xml` | **ACTIVE (terminal).** Tier 1 conservation tests JUnit XML at the CA-1 + SCI-N1 fixes (40 passed). |
| `baseline_junit_full_624ed7c.xml` | Superseded (CA-1 only). Full repo test suite JUnit XML (993 passed, 2 xfailed, 0 failed, 327.77 s wall). |
| `baseline_tier1_junit_624ed7c.xml` | Superseded (CA-1 only). Tier 1 conservation tests JUnit XML (40 passed). |
| `baseline_junit_full_186b5c4.xml` | Superseded (pre-CA-1). Full repo test suite JUnit XML (805 passed, 2 xfailed, 0 failed, 71.36 s wall). |
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

## Re-baseline log

### `b51df71` — 2026-05-16 — NSM1-SCI-N1 (terminal gold-standard baseline)

Trigger: a deliberate kinetics change (gold-standard spec Workstream A2). Commit `b51df71` fixed NSM1-SCI-N1 (MAJOR): the denitrification alkalinity coefficient `r_alkden` was corrected `4/14/1000` → `1/14/1000` eq/mg-N (1 eq alkalinity per mol NO₃-N reduced; CE-QUAL-W2 `water-quality.f90:3157`, Stumm & Morgan), a deliberate divergence from the upstream Fortran/v1 defect.

Scope of change vs `624ed7c`: **only `alkalinity` differs** (21,600/21,600 cell-substeps; the denitrification source term); all other 19 variables bit-identical. Captured under `pixi --environment dev`; parity re-verified under the conda `clearwater` test env.

**Terminal baseline.** NSM1-CA-1 (`624ed7c`) and NSM1-SCI-N1 (`b51df71`) are the only trajectory-perturbing changes in the entire gold-standard gate set — the remaining gate items (C1 carbon-routing, C2 CBOD, C3/C4 DOX docs/guards/off-by-default, C5 doc-only) do not alter the coupled-demo trajectory. Absent a new deliberate kinetics change or dependency bump, `b51df71` is the final gold-standard reference. References in `test_coupled_demo_parity.py` and `check_baseline_parity.py` now point at `b51df71`; `624ed7c` and `186b5c4` artifacts retained unmodified.

### `624ed7c` — 2026-05-16 — NSM1-CA-1 alkalinity kinetics fix

Trigger: a deliberate kinetics change committed outside the pattern-alignment work (gold-standard spec Workstream A1). Commit `624ed7c` fixed NSM1-CA-1 (CRITICAL): the alkalinity algal/benthic coupling now uses the intensive carbon ratios `rca = AWc/AWa`, `rcb = BWc/BWd` instead of the raw stoichiometric weights, correcting a 1000×/100× overstatement of the floating/benthic alkalinity flux.

Scope of change vs `186b5c4` baseline: **only `alkalinity` differs** (21,600/21,600 cell-substeps); all other 19 state/forcing variables are bit-identical. Verified bit-identical across two environments (captured under `pixi --environment dev`; parity re-verified under the conda `clearwater` test env). The `186b5c4` `.nc` and companion XML/JSON/txt artifacts are retained unmodified; references in `test_coupled_demo_parity.py` and `check_baseline_parity.py` now point at `624ed7c`.

## Pathogen warnings during capture

`Pathogen` prints two warnings on the first capture about optional registry variables `Solid` and `ap` not being present. These are the once-only `_get_optional` warn-latches; they are part of the baseline behavior and will be reproduced (silently after the first emission) on every parity run.
