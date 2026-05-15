"""v1-vs-v3 NSM1 standalone performance benchmark.

Head-to-head per-step kinetics timing (and peak memory) for v1 NSM1
(``clearwater_modules.nsm1.model.NutrientBudget``) versus v3 NSM1
(the 11-Process stack via ``build_nsm1_demo``). No transport, no
orchestrator: the comparison isolates the NSM1 implementations.

This is the version-controlled, reviewed form of the harness used to
produce ``design/clearwater_modules_v3_nsm1_v1_v3_performance_memo.md``.
See that memo for the interpretation and the stated limitations
(NSM1-only; synthetic constant IC; single-run timing; v1 measured in
its production configuration with ``track_dynamic_variables=False``).

v1 is constructed exactly as the coupled pipeline constructs it
(``ClearWater-modules-phase2-ESM-streaming/case_studies/santiam_salem/
scripts/08_run_coupled.py``): same parameter dicts, the five active
states with the other eleven zeroed, ``updateable_static_variables``,
``track_dynamic_variables=False``, ``time_dim="days"``. The synthetic
constant initial conditions match the v3 demo defaults / Santiam-Salem
provenance ICs.

Two modes:

* ``sweep`` — times v1 and v3 across a cell-count range in one
  process and prints a comparison table. Peak RSS is NOT reliable in
  this mode (the two engines' allocations contaminate each other);
  use ``isolated`` for memory.
* ``isolated`` — runs a single engine at a single cell count in this
  process and prints clean peak RSS. Invoke once per engine (in
  separate processes) for an apples-to-apples memory comparison.

Memory note: 2,000,000 cells exceeds the memory available on the
target workstation for v1 (its ``(time_steps+1, n_cells)``
state-history allocation). The default sweep ceiling is therefore
1,000,000 cells. ``--max-cells`` can raise it, with the understanding
that v1 may be killed by the OS at larger sizes (reported as an error
row, not a crash).

Usage:
    pixi run --environment dev python tests/v3/nsm1/baseline/v1_v3_nsm1_benchmark.py
    pixi run --environment dev python tests/v3/nsm1/baseline/v1_v3_nsm1_benchmark.py --mode isolated --engine v3 --n-cells 1000000
    pixi run --environment dev python tests/v3/nsm1/baseline/v1_v3_nsm1_benchmark.py --mode isolated --engine v1 --n-cells 1000000
"""
from __future__ import annotations

import argparse
import gc
import resource
import statistics
import sys
import time
import traceback
from datetime import datetime
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[3]
_SRC = _REPO_ROOT / "src"
if _SRC.is_dir() and str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import numpy as np
import xarray as xr

START = datetime(2026, 1, 1)

# v1 NSM1 state set (NSM_ALL_STATES / NSM_ACTIVE_STATES mirror
# 08_run_coupled.py). Inactive states are zero-initialised.
NSM_ALL_STATES = [
    "Ap", "Ab", "NH4", "NO3", "OrgN", "N2",
    "TIP", "OrgP", "POC", "DOC", "DIC",
    "POM", "CBOD", "DOX", "PX", "Alk",
]
NSM_ACTIVE_STATES = ["Ap", "NH4", "NO3", "TIP", "DOX"]

# Synthetic constant ICs (v3 demo defaults / Santiam-Salem provenance).
IC = {"Ap": 1.6, "NH4": 0.02, "NO3": 0.137, "TIP": 0.029, "DOX": 9.4}
TWC = 17.35
QSOLAR = 200.0

# Parameter dicts mirror 08_run_coupled.py (santiam_salem copy).
GLOBAL_PARAMETERS = {
    "use_NH4": True, "use_NO3": True, "use_TIP": True, "use_DOX": True,
    "use_Algae": True, "use_Balgae": False, "use_OrgN": False,
    "use_OrgP": False, "use_POC": False, "use_DOC": False,
    "use_DIC": False, "use_POM": False, "use_CBOD": False,
    "use_Pathogen": False, "use_Alk": False, "use_N2": False,
    "use_SedFlux": False,
}
GLOBAL_VARS = {
    "vson": 0.01, "vsoc": 0.01, "vsop": 0.01, "vs": 0.01,
    "SOD_20": 0.5, "SOD_theta": 1.047, "vb": 0.01, "fcom": 0.4,
    "kaw_20_user": 0.0, "kah_20_user": 1.0,
    "hydraulic_reaeration_option": 1, "wind_reaeration_option": 1,
    "dt": 1.0, "depth": 1.5, "TwaterC": TWC, "theta": 1.047,
    "velocity": 0.5, "flow": 150.0, "topwidth": 100.0, "slope": 0.0002,
    "shear_velocity": 0.05334, "pressure_mb": 1013.25, "wind_speed": 3.0,
    "q_solar": QSOLAR, "Solid": 1,
    "lambda0": 0.02, "lambda1": 0.0088, "lambda2": 0.054,
    "lambdas": 0.052, "lambdam": 0.0174, "Fr_PAR": 0.47,
}
ALGAE_PARAMETERS = {
    "AWd": 100, "AWc": 40, "AWn": 7.2, "AWp": 1, "AWa": 1000,
    "KL": 10, "KsN": 0.04, "KsP": 0.0012,
    "mu_max_20": 1.0, "kdp_20": 0.15, "krp_20": 0.2,
    "vsap": 0.15, "growth_rate_option": 3, "light_limitation_option": 1,
}

# (n_cells, n_warmup, n_measured). Windows shrink at large N so the
# benchmark stays bounded. Ceiling is 1,000,000 (see module docstring).
DEFAULT_SWEEP = [
    (5, 10, 60),
    (1_000, 10, 40),
    (10_000, 8, 30),
    (100_000, 5, 20),
    (500_000, 3, 10),
    (1_000_000, 3, 8),
]


def _peak_gb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024**3)


def _da(n: int, val: float) -> xr.DataArray:
    return xr.DataArray(
        np.full(n, val, dtype=np.float64),
        dims=["nface"], coords={"nface": np.arange(n)},
    )


def bench_v1(n_cells: int, n_warm: int, n_meas: int) -> float:
    """Mean ms per ``increment_timestep`` at ``n_cells`` (after warmup)."""
    from clearwater_modules.nsm1.model import NutrientBudget

    isv = {
        nm: _da(n_cells, IC.get(nm, 0.0) if nm in NSM_ACTIVE_STATES else 0.0)
        for nm in NSM_ALL_STATES
    }
    reaction = NutrientBudget(
        time_steps=n_warm + n_meas + 2,
        initial_state_values=isv,
        updateable_static_variables=["q_solar", "TwaterC"],
        algae_parameters=ALGAE_PARAMETERS,
        global_parameters=GLOBAL_PARAMETERS,
        global_vars=GLOBAL_VARS,
        track_dynamic_variables=False,
        time_dim="days",
    )

    def _inputs() -> dict:
        d = {nm: _da(n_cells, IC[nm]) for nm in NSM_ACTIVE_STATES}
        d["TwaterC"] = _da(n_cells, TWC)
        d["q_solar"] = _da(n_cells, QSOLAR)
        return d

    for _ in range(n_warm):
        reaction.increment_timestep(_inputs())
    samples = []
    for _ in range(n_meas):
        inp = _inputs()
        a = time.perf_counter()
        reaction.increment_timestep(inp)
        samples.append((time.perf_counter() - a) * 1000.0)
    return statistics.mean(samples)


def bench_v3(n_cells: int, n_warm: int, n_meas: int) -> float:
    """Mean ms per ``demo.step`` at ``n_cells`` (after warmup)."""
    from clearwater_modules_v3.examples import build_nsm1_demo

    demo = build_nsm1_demo(n_cells=n_cells)
    t = START
    for _ in range(n_warm):
        demo.step(t)
    samples = []
    for _ in range(n_meas):
        a = time.perf_counter()
        demo.step(t)
        samples.append((time.perf_counter() - a) * 1000.0)
    return statistics.mean(samples)


def _run_sweep(max_cells: int) -> int:
    sweep = [row for row in DEFAULT_SWEEP if row[0] <= max_cells]
    print("v1-vs-v3 NSM1 standalone speed sweep")
    print("  timing: per-step advance only (excludes construction)")
    print("  peak RSS NOT reliable here — use --mode isolated for memory\n")
    print(f"  {'cells':>9} | {'v1 ms/step':>12} | {'v3 ms/step':>12} | "
          f"{'v1/v3':>6} | faster")
    print("  " + "-" * 60)
    for n_cells, n_warm, n_meas in sweep:
        res = {}
        for eng, fn in (("v1", bench_v1), ("v3", bench_v3)):
            gc.collect()
            try:
                res[eng] = fn(n_cells, n_warm, n_meas)
            except Exception as e:  # noqa: BLE001
                res[eng] = None
                print(f"  [{eng} @ {n_cells:,} cells] ERROR: "
                      f"{type(e).__name__}: {e}", file=sys.stderr)
                traceback.print_exc()
        v1, v3 = res.get("v1"), res.get("v3")
        if v1 and v3:
            ratio = v1 / v3
            faster = ("v1" if ratio < 1 else "v3" if ratio > 1 else "tied")
            print(f"  {n_cells:>9,} | {v1:12.1f} | {v3:12.1f} | "
                  f"{ratio:6.2f} | {faster}")
        else:
            v1s = f"{v1:12.1f}" if v1 else f"{'ERR':>12}"
            v3s = f"{v3:12.1f}" if v3 else f"{'ERR':>12}"
            print(f"  {n_cells:>9,} | {v1s} | {v3s} | {'n/a':>6} | n/a")
    print("\n  v1/v3 > 1 means v1 slower; < 1 means v1 faster.")
    return 0


def _run_isolated(engine: str, n_cells: int) -> int:
    n_warm, n_meas = 3, 8
    gc.collect()
    fn = bench_v1 if engine == "v1" else bench_v3
    try:
        ms = fn(n_cells, n_warm, n_meas)
    except Exception as e:  # noqa: BLE001
        print(f"RESULT engine={engine} n_cells={n_cells} ERROR="
              f"{type(e).__name__}: {e}", file=sys.stderr)
        traceback.print_exc()
        return 1
    print(
        f"RESULT engine={engine} n_cells={n_cells} "
        f"mean_ms_per_step={ms:.1f} "
        f"ms_per_Mcell={ms / (n_cells / 1e6):.1f} "
        f"peak_rss_GB={_peak_gb():.2f}"
    )
    return 0


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mode", choices=("sweep", "isolated"), default="sweep")
    p.add_argument("--engine", choices=("v1", "v3"),
                   help="required for --mode isolated")
    p.add_argument("--n-cells", type=int,
                   help="required for --mode isolated")
    p.add_argument("--max-cells", type=int, default=1_000_000,
                   help="sweep ceiling (default 1,000,000; 2,000,000 "
                        "exceeds target-workstation memory for v1)")
    args = p.parse_args()

    if args.mode == "isolated":
        if not args.engine or not args.n_cells:
            p.error("--mode isolated requires --engine and --n-cells")
        return _run_isolated(args.engine, args.n_cells)
    return _run_sweep(args.max_cells)


if __name__ == "__main__":
    raise SystemExit(main())
