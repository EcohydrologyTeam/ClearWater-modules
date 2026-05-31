"""Bit-identical parity check against the Phase 0 baseline.

Runs the same v3 NSM1 demo (4,320 substeps, fixed RNG seed, no
``REGISTRY_DIAGNOSTICS`` names pre-registered) and asserts the
resulting state-variable trajectory equals the committed Phase 0
baseline NetCDF *bit-identically* across every variable, every cell,
every substep.

This is the runner that enforces clause §11.2 of
``design/clearwater_modules_v3_nsm1_pattern_alignment_specification.md``.
Every per-Process phase commit (Phase 1 through Phase 10) runs this
check; failure rolls the phase back.

Usage:
    pixi run --environment dev python tests/v3/nsm1/baseline/check_baseline_parity.py
        [--baseline tests/v3/nsm1/baseline/baseline_coupled_trajectory_d82a5ed.nc]

Exit code:
    0 — bit-identical match
    1 — mismatch on at least one variable; prints the offending list
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[3]
_SRC = _REPO_ROOT / "src"
if _SRC.is_dir() and str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import numpy as np
import xarray as xr


def _capture_current_trajectory(baseline_attrs: dict) -> xr.Dataset:
    """Re-run the demo with the same configuration as the baseline."""
    from datetime import datetime, timedelta

    from clearwater_modules_v3.examples import build_nsm1_demo

    np.random.seed(int(baseline_attrs["rng_seed"]))
    n_substeps = int(baseline_attrs["n_substeps"])
    start = datetime.fromisoformat(baseline_attrs["start_time_iso"])
    time_step_seconds = int(baseline_attrs["time_step_seconds"])
    expected_n_cells = int(baseline_attrs["n_cells"])

    demo = build_nsm1_demo(time_step=timedelta(seconds=time_step_seconds))
    assert demo.registry.get(sorted(demo.registry.keys())[0]).sizes["cell"] == expected_n_cells, (
        f"baseline n_cells={expected_n_cells}, demo n_cells differs"
    )

    state_names = sorted(demo.registry.keys())
    traj = {
        name: np.empty((n_substeps + 1, expected_n_cells), dtype=np.float64)
        for name in state_names
    }
    for name in state_names:
        traj[name][0] = demo.registry.get(name).values

    t = start
    for i in range(n_substeps):
        demo.step(t)
        t += demo.time_step
        for name in state_names:
            traj[name][i + 1] = demo.registry.get(name).values

    return xr.Dataset(
        data_vars={n: (("substep", "cell"), traj[n]) for n in state_names},
        coords={
            "substep": np.arange(n_substeps + 1, dtype=np.int64),
            "cell": np.arange(expected_n_cells, dtype=np.int64),
        },
    )


def check_parity(baseline_path: Path) -> tuple[bool, list[str]]:
    baseline = xr.open_dataset(baseline_path)
    current = _capture_current_trajectory(baseline.attrs)

    baseline_vars = set(baseline.data_vars)
    current_vars = set(current.data_vars)
    if baseline_vars != current_vars:
        only_baseline = baseline_vars - current_vars
        only_current = current_vars - baseline_vars
        mismatches = []
        if only_baseline:
            mismatches.append(f"variables missing from current run: {sorted(only_baseline)}")
        if only_current:
            mismatches.append(f"new variables not in baseline: {sorted(only_current)}")
        return False, mismatches

    mismatches: list[str] = []
    for name in sorted(baseline_vars):
        b_vals = baseline[name].values
        c_vals = current[name].values
        if b_vals.shape != c_vals.shape:
            mismatches.append(
                f"{name}: shape {c_vals.shape} != baseline {b_vals.shape}"
            )
            continue
        if not np.array_equal(b_vals, c_vals):
            diff_count = int((b_vals != c_vals).sum())
            try:
                max_abs = float(np.nanmax(np.abs(b_vals - c_vals)))
            except (TypeError, ValueError):
                max_abs = float("nan")
            mismatches.append(
                f"{name}: {diff_count} non-identical cells (max|diff|={max_abs:.3e})"
            )

    return (len(mismatches) == 0), mismatches


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--baseline",
        default="tests/v3/nsm1/baseline/baseline_coupled_trajectory_d82a5ed.nc",
        help="Path to the Phase 0 baseline NetCDF",
    )
    args = parser.parse_args()

    baseline_path = Path(args.baseline)
    if not baseline_path.is_file():
        print(f"ERROR: baseline not found at {baseline_path}", file=sys.stderr)
        return 2

    ok, mismatches = check_parity(baseline_path)
    if ok:
        print(f"OK: bit-identical parity vs {baseline_path}")
        return 0

    print(f"FAIL: {len(mismatches)} variable(s) diverge from {baseline_path}", file=sys.stderr)
    for m in mismatches:
        print(f"  - {m}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
