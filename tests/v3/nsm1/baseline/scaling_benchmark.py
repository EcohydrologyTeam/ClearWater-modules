"""Cell-count scaling benchmark for v3 NSM1.

Answers the existential question: does v3 NSM1's per-substep cost scale
acceptably from the 5-cell test mesh to the 500K-2M-cell production
regime the project needs?

Sweeps n_cells = 5 / 1K / 10K / 100K / 500K, with step counts adapted
down at large N so the wall time stays bounded. Reports mean / median
ms/step and a normalised ms/step/Mcell so the scaling exponent is
visible.

Usage:
    pixi run --environment dev python tests/v3/nsm1/baseline/scaling_benchmark.py
"""

from __future__ import annotations

import gc
import resource
import statistics
import sys
import time
from datetime import datetime
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[3]
_SRC = _REPO_ROOT / "src"
if _SRC.is_dir() and str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from clearwater_modules_v3.examples import build_nsm1_demo


START = datetime(2026, 1, 1, 0, 0, 0)

# (n_cells, n_warmup, n_measured). Step counts shrink at large N so the
# benchmark finishes in a few minutes total.
SWEEP = [
    (5, 30, 200),
    (1_000, 20, 100),
    (10_000, 15, 60),
    (100_000, 8, 30),
    (500_000, 4, 12),
]


def _bench(n_cells: int, n_warmup: int, n_measured: int) -> dict:
    gc.collect()
    t_build0 = time.perf_counter()
    demo = build_nsm1_demo(n_cells=n_cells)
    t_build1 = time.perf_counter()

    t = START
    for _ in range(n_warmup):
        demo.step(t)
    samples = []
    for _ in range(n_measured):
        a = time.perf_counter()
        demo.step(t)
        b = time.perf_counter()
        samples.append((b - a) * 1000.0)

    peak_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024)
    mean_ms = statistics.mean(samples)
    median_ms = statistics.median(samples)
    p95_ms = sorted(samples)[min(len(samples) - 1, int(0.95 * len(samples)))]
    # Normalised: ms per step per million cells.
    ms_per_mcell = mean_ms / (n_cells / 1_000_000.0)

    result = {
        "n_cells": n_cells,
        "build_s": t_build1 - t_build0,
        "mean_ms": mean_ms,
        "median_ms": median_ms,
        "p95_ms": p95_ms,
        "peak_mb": peak_mb,
        "ms_per_mcell": ms_per_mcell,
    }
    del demo
    gc.collect()
    return result


def main() -> int:
    print(f"v3 NSM1 cell-count scaling benchmark")
    print(f"  sweep: {[s[0] for s in SWEEP]} cells\n")

    rows = []
    for n_cells, n_warm, n_meas in SWEEP:
        r = _bench(n_cells, n_warm, n_meas)
        rows.append(r)
        print(
            f"  {n_cells:>8,} cells | build {r['build_s']:5.2f}s | "
            f"mean {r['mean_ms']:8.1f} ms/step | "
            f"median {r['median_ms']:8.1f} | "
            f"p95 {r['p95_ms']:8.1f} | "
            f"peak {r['peak_mb']:6.0f} MB | "
            f"{r['ms_per_mcell']:8.1f} ms/step/Mcell"
        )

    print("\nExtrapolation to production targets (linear in the")
    print("large-N regime, using the 500K-cell ms/step/Mcell):")
    big = rows[-1]
    ms_per_cell = big["mean_ms"] / big["n_cells"]
    for n_cells in (500_000, 1_000_000, 2_000_000):
        ms_step = ms_per_cell * n_cells
        for label, n_steps in (
            ("3 mo @ 15-min (8,760 steps)", 8_760),
            ("1 yr @ 15-min (35,040 steps)", 35_040),
        ):
            total_s = ms_step / 1000.0 * n_steps
            hrs = total_s / 3600.0
            print(
                f"  {n_cells:>9,} cells, {label:<32} "
                f"~{ms_step:8.0f} ms/step -> ~{hrs:6.1f} hr"
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
