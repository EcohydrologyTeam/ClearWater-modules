"""v1-vs-v3 TSM standalone performance benchmark.

Head-to-head per-step energy-balance timing (and peak memory) for
v1 TSM (``clearwater_modules.tsm.model.EnergyBudget``) versus v3 TSM
(``clearwater_modules_v3.processes.temperature.Temperature``). No
transport, no orchestrator: the comparison isolates the TSM
implementations.

Companion to the v1-vs-v3 NSM1 benchmark
(``v1_v3_nsm1_benchmark.py``) and to
``design/clearwater_modules_v3_nsm1_v1_v3_performance_memo.md``. Same
method and the same stated limitations apply (TSM-only; synthetic
constant IC/forcing; single-run timing; v1 measured in its production
configuration ``track_dynamic_variables=False``).

v1 is constructed exactly as the dedicated v1 TSM driver constructs
it (``ClearWater-modules-phase2-ESM-streaming/case_studies/
corvallis_santiam_albany/scripts/12_run_tsm_hourly.py``): state =
{water_temp_c, volume, surface_area}; the six meteo forcings declared
``updateable_static_variables``; ``use_sed_temp=False``;
``track_dynamic_variables=False``; ``time_dim="hours"``; advanced via
``increment_timestep(inputs)``. v3 uses a ``Temperature`` Process on
an in-memory registry with the same synthetic constant fields,
``use_sediment_temperature=False`` to match v1's ``use_sed_temp``,
advanced via ``run(time, registry)``.

``12_run_tsm_hourly.py`` is cited only as the authoritative recipe
for v1 ``EnergyBudget`` *object construction* (state set, updateable
statics, ``use_sed_temp``, ``track_dynamic_variables``, ``time_dim``);
it is not a validation-provenance source. It is the only dedicated
standalone v1 TSM driver in the repository — Santiam-Salem runs TSM
only in coupled mode (``08_run_coupled*.py``), so no standalone
Santiam-Salem TSM script exists. The construction recipe is
reach-independent (object configuration, not run inputs).

v3 ``Temperature`` skips its first ``run`` call by design
(``__skip_first_time_step``); the warmup loop length accounts for
this so the timed region excludes the skipped step.

Modes (same as the NSM1 benchmark):

* ``sweep`` — times v1 and v3 across a cell-count range in one
  process; peak RSS NOT reliable here, use ``isolated`` for memory.
* ``isolated`` — single engine, single N, own process, clean RSS.

2,000,000 cells exceeds target-workstation memory for v1; default
sweep ceiling is 1,000,000.

Usage:
    pixi run --environment dev python tests/v3/nsm1/baseline/v1_v3_tsm_benchmark.py
    pixi run --environment dev python tests/v3/nsm1/baseline/v1_v3_tsm_benchmark.py --mode isolated --engine v3 --n-cells 1000000
    pixi run --environment dev python tests/v3/nsm1/baseline/v1_v3_tsm_benchmark.py --mode isolated --engine v1 --n-cells 1000000
"""
from __future__ import annotations

import argparse
import gc
import resource
import statistics
import sys
import time
import traceback
from datetime import datetime, timedelta
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[3]
_SRC = _REPO_ROOT / "src"
if _SRC.is_dir() and str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import numpy as np
import xarray as xr

START = datetime(2026, 1, 1)

# Synthetic constant IC / forcing (Santiam-Salem-ish, physical).
T_WATER = 17.35
VOLUME = 1.5
SURF_AREA = 1.0
AIR_T = 20.0
QSOLAR = 400.0
WIND = 3.0
EAIR_MB = 1.0
PRESSURE_MB = 1013.0
CLOUD = 0.1
SED_T = 12.0
SED_THICK = 0.1
DT_DAYS = 1.0 / 24.0  # hourly

# v1 meteo defaults mirror 12_run_tsm_hourly.py.
V1_UPDATEABLE = ["air_temp_c", "q_solar", "wind_speed",
                 "eair_mb", "pressure_mb", "cloudiness"]
V1_METEO = {
    "air_temp_c": AIR_T, "q_solar": QSOLAR, "wind_speed": WIND,
    "eair_mb": EAIR_MB, "pressure_mb": PRESSURE_MB, "cloudiness": CLOUD,
    "sed_temp_c": SED_T, "dt": DT_DAYS,
}

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
    """Mean ms per v1 ``increment_timestep`` (after warmup)."""
    from clearwater_modules.tsm.model import EnergyBudget

    model = EnergyBudget(
        time_steps=n_warm + n_meas + 2,
        initial_state_values={
            "water_temp_c": _da(n_cells, T_WATER),
            "volume": _da(n_cells, VOLUME),
            "surface_area": _da(n_cells, SURF_AREA),
        },
        updateable_static_variables=list(V1_UPDATEABLE),
        meteo_parameters=dict(V1_METEO),
        temp_parameters={"sed_temp_c": SED_T},
        use_sed_temp=False,
        track_dynamic_variables=False,
        time_dim="hours",
    )

    def _inputs() -> dict:
        return {
            "water_temp_c": _da(n_cells, T_WATER),
            "volume": _da(n_cells, VOLUME),
            "surface_area": _da(n_cells, SURF_AREA),
            "air_temp_c": _da(n_cells, AIR_T),
            "q_solar": _da(n_cells, QSOLAR),
            "wind_speed": _da(n_cells, WIND),
            "eair_mb": _da(n_cells, EAIR_MB),
            "pressure_mb": _da(n_cells, PRESSURE_MB),
            "cloudiness": _da(n_cells, CLOUD),
        }

    for _ in range(n_warm):
        model.increment_timestep(_inputs())
    samples = []
    for _ in range(n_meas):
        inp = _inputs()
        a = time.perf_counter()
        model.increment_timestep(inp)
        samples.append((time.perf_counter() - a) * 1000.0)
    return statistics.mean(samples)


def bench_v3(n_cells: int, n_warm: int, n_meas: int) -> float:
    """Mean ms per v3 ``Temperature.run`` (after warmup; the
    design-intended first-call skip is consumed during warmup)."""
    from clearwater_modules_v3.examples.nsm1_demo_setup import InMemoryRegistry
    from clearwater_modules_v3.processes.temperature import Temperature

    reg = InMemoryRegistry()
    fields = {
        "water_temperature": T_WATER,
        "wetted_surface_area": SURF_AREA,
        "volume": VOLUME,
        "cloudiness": CLOUD,
        "air_temperature": AIR_T,
        "solar_radiation": QSOLAR,
        "wind_speed": WIND,
        "atmospheric_pressure": PRESSURE_MB,
        "atmospheric_vapor_pressure": EAIR_MB,
        "sediment_temperature": SED_T,
        "sediment_thickness": SED_THICK,
    }
    for name, val in fields.items():
        reg.register(name, _da(n_cells, val))

    temp = Temperature(
        time_step=timedelta(hours=1),
        use_sediment_temperature=False,
    )

    class _StubModel:  # init_process only does getattr(model,"diagnostics",None)
        pass

    temp.init_process(_StubModel(), reg)

    t = START
    # First run() is a no-op by design; ensure warmup absorbs it.
    for _ in range(max(n_warm, 2)):
        temp.run(t, reg)
    samples = []
    for _ in range(n_meas):
        a = time.perf_counter()
        temp.run(t, reg)
        samples.append((time.perf_counter() - a) * 1000.0)
    return statistics.mean(samples)


def _run_sweep(max_cells: int) -> int:
    sweep = [row for row in DEFAULT_SWEEP if row[0] <= max_cells]
    print("v1-vs-v3 TSM standalone speed sweep")
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
    p.add_argument("--engine", choices=("v1", "v3"))
    p.add_argument("--n-cells", type=int)
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
