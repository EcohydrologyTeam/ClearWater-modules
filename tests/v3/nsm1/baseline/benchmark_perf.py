"""Phase 10 perf benchmark — pattern-alignment overhead measurement.

Measures per-substep wall time on the 5-cell synthetic mesh in two
modes:

1. **No-subscription** — base pattern-aligned NSM1 with no Appendix A
   diagnostic names pre-registered. The opportunistic-write loop is
   ``n × O(1)`` membership checks per substep with zero registry
   writes. Should be effectively equivalent to the pre-refactor cost.
2. **Full-subscription** — every Appendix A diagnostic name
   pre-registered. The opportunistic-write loop pays
   ``n × set_at_time`` per substep (~80 writes per substep across
   all 11 Processes).

Spec §8 perf budgets (pattern-alignment spec):

- No-subscription: ≤ 5% wall-clock overhead vs pre-refactor baseline.
- Full-subscription: ≤ 15% wall-clock overhead vs no-subscription.

The pre-refactor baseline was 17.6 ms/step on the 5-cell mesh
(measured at Phase 0). This script does NOT enforce the budget — it
prints the measurements so they can be recorded in the Phase 10.A
closeout.

Usage:
    pixi run --environment dev python tests/v3/nsm1/baseline/benchmark_perf.py
"""

from __future__ import annotations

import statistics
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[3]
_SRC = _REPO_ROOT / "src"
if _SRC.is_dir() and str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import numpy as np
import xarray as xr


N_WARMUP_SUBSTEPS = 60
N_MEASURED_SUBSTEPS = 500
START = datetime(2026, 1, 1, 0, 0, 0)
PHASE_0_BASELINE_MS_PER_STEP = 17.6  # documented baseline


def _build_demo_unsubscribed():
    from clearwater_modules_v3.examples import build_nsm1_demo
    return build_nsm1_demo()


def _build_demo_subscribed():
    from clearwater_modules_v3.examples import (
        build_nsm1_demo,
        default_initial_conditions,
    )
    from clearwater_modules_v3.processes import (
        Alkalinity, BenthicAlgae, CBOD, Carbon, DOX, FloatingAlgae,
        N2, Nitrogen, POM, Pathogen, Phosphorus,
    )

    ic = default_initial_conditions()
    reference = ic["pom"]
    all_names = []
    for cls in (
        Carbon, DOX, Nitrogen, FloatingAlgae, BenthicAlgae, Phosphorus,
        POM, CBOD, N2, Pathogen, Alkalinity,
    ):
        all_names.extend(cls.REGISTRY_DIAGNOSTICS)
    for name in all_names:
        if name not in ic:
            ic[name] = xr.zeros_like(reference)
    return build_nsm1_demo(initial_conditions=ic), len(all_names)


def _run(demo, n_warmup: int, n_measured: int) -> float:
    """Returns mean ms/substep over the measurement window."""
    t = START
    # Warmup.
    for _ in range(n_warmup):
        demo.step(t)
        t += demo.time_step

    # Measured window.
    samples: list[float] = []
    for _ in range(n_measured):
        t0 = time.perf_counter()
        demo.step(t)
        t1 = time.perf_counter()
        samples.append((t1 - t0) * 1000.0)  # ms
        t += demo.time_step

    return statistics.mean(samples)


def main() -> int:
    print(f"Phase 10 perf benchmark")
    print(f"  Phase 0 baseline: {PHASE_0_BASELINE_MS_PER_STEP} ms/step")
    print(f"  Warmup substeps:  {N_WARMUP_SUBSTEPS}")
    print(f"  Measured substeps: {N_MEASURED_SUBSTEPS}")
    print()

    # No-subscription.
    demo_unsubscribed = _build_demo_unsubscribed()
    no_sub_ms = _run(
        demo_unsubscribed, N_WARMUP_SUBSTEPS, N_MEASURED_SUBSTEPS
    )
    print(f"No-subscription:    {no_sub_ms:.2f} ms/step")

    # Full-subscription.
    demo_subscribed, n_names = _build_demo_subscribed()
    full_sub_ms = _run(
        demo_subscribed, N_WARMUP_SUBSTEPS, N_MEASURED_SUBSTEPS
    )
    print(f"Full subscription:  {full_sub_ms:.2f} ms/step "
          f"({n_names} Appendix A names registered)")

    print()
    overhead_vs_baseline = (no_sub_ms / PHASE_0_BASELINE_MS_PER_STEP - 1.0) * 100
    overhead_subscribed = (full_sub_ms / no_sub_ms - 1.0) * 100
    print(f"Overhead summary:")
    print(f"  no-sub vs Phase 0 baseline:  {overhead_vs_baseline:+.1f}%  "
          f"(budget: ≤ 5.0%)")
    print(f"  full-sub vs no-sub:          {overhead_subscribed:+.1f}%  "
          f"(budget: ≤ 15.0%)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
