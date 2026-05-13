"""Capture the Phase 0 gold-reference coupled-run trajectory.

Runs the v3 NSM1 5-cell synthetic-mesh demo for N substeps with a fixed
RNG seed and saves the full state-variable trajectory (every substep,
every cell, every state) to a NetCDF file under
``tests/v3/nsm1/baseline/baseline_coupled_trajectory_<commit>.nc``.

This is the load-bearing artifact of the pattern-alignment
specification's §11 zero-regression contract: every subsequent phase
commit must reproduce this NetCDF bit-identically when no
``REGISTRY_DIAGNOSTICS`` names are pre-registered.

Usage:
    pixi run python tests/v3/nsm1/baseline/capture_baseline_trajectory.py <commit-short-hash>

The script writes ``baseline_coupled_trajectory_<commit>.nc`` next to
itself.
"""

from __future__ import annotations

import sys
from datetime import datetime, timedelta
from pathlib import Path

# Match tests/v3/conftest.py: the streaming-repo ``src/`` must be on
# sys.path so we import the in-tree clearwater_modules_v3 and not the
# vendored editable install. Repo root is two parents above the
# baseline directory.
_REPO_ROOT = Path(__file__).resolve().parents[3]
_SRC = _REPO_ROOT / "src"
if _SRC.is_dir() and str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import numpy as np
import xarray as xr

# Spec target: 4,320 substeps. With time_step = 5 min, this is 15 days
# of model time.
N_SUBSTEPS = 4320

# Fixed RNG seed (the demo currently uses deterministic ICs; the seed
# is captured here for future scenarios that might introduce
# stochastic forcing).
RNG_SEED = 20260513


def main(commit: str) -> Path:
    from clearwater_modules_v3.examples import build_nsm1_demo

    np.random.seed(RNG_SEED)
    demo = build_nsm1_demo()
    start = datetime(2026, 1, 1, 0, 0, 0)

    # Step-by-step capture so we get every substep into the trajectory,
    # not just every snapshot_every. The deep-copy snapshot path costs
    # O(n_vars * n_cells) per substep but is acceptable at 5 cells.
    state_names = sorted(demo.registry.keys())
    n_cells = demo.registry.get(state_names[0]).sizes["cell"]

    traj = {name: np.empty((N_SUBSTEPS + 1, n_cells), dtype=np.float64) for name in state_names}
    for name in state_names:
        traj[name][0] = demo.registry.get(name).values

    t = start
    for i in range(N_SUBSTEPS):
        demo.step(t)
        t += demo.time_step
        for name in state_names:
            traj[name][i + 1] = demo.registry.get(name).values

    # Substep coordinate: 0 .. N_SUBSTEPS (inclusive of IC at index 0)
    substep_index = np.arange(N_SUBSTEPS + 1, dtype=np.int64)

    ds = xr.Dataset(
        data_vars={
            name: (("substep", "cell"), traj[name]) for name in state_names
        },
        coords={
            "substep": substep_index,
            "cell": np.arange(n_cells, dtype=np.int64),
        },
        attrs={
            "commit": commit,
            "rng_seed": RNG_SEED,
            "n_substeps": N_SUBSTEPS,
            "n_cells": n_cells,
            "time_step_seconds": int(demo.time_step.total_seconds()),
            "start_time_iso": start.isoformat(),
            "scope": "v3 NSM1 1.0.0 + pattern-alignment baseline",
            "purpose": "Phase 0 gold reference for §11 zero-regression contract",
        },
    )

    out_path = Path(__file__).parent / f"baseline_coupled_trajectory_{commit}.nc"
    # Use NETCDF4 with no compression so bit-identicality is preserved
    # across reads. Compression with shuffle/zlib is not byte-stable
    # across libnetcdf versions.
    ds.to_netcdf(out_path, format="NETCDF4", engine="netcdf4")
    return out_path


if __name__ == "__main__":
    commit = sys.argv[1] if len(sys.argv) > 1 else "unknown"
    out = main(commit)
    print(f"Wrote {out}")
