"""Phase 11 profiling — cProfile a coupled NSM1 demo run.

Runs ``build_nsm1_demo()`` + N_PROFILE substeps under cProfile,
writes the .prof to ``profile_perf.prof`` and a text summary to
``profile_perf.txt``. The text summary covers the top-30 hot functions
by both cumulative time and self time.

Usage:
    pixi run --environment dev python tests/v3/nsm1/baseline/profile_perf.py
"""

from __future__ import annotations

import cProfile
import io
import pstats
import sys
from datetime import datetime, timedelta
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[3]
_SRC = _REPO_ROOT / "src"
if _SRC.is_dir() and str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))


# Same warmup + measurement window as the perf benchmark for direct
# comparability.
N_WARMUP = 60
N_PROFILE = 500
START = datetime(2026, 1, 1, 0, 0, 0)


def _do_run() -> None:
    from clearwater_modules_v3.examples import build_nsm1_demo

    demo = build_nsm1_demo()
    t = START
    for _ in range(N_WARMUP):
        demo.step(t)
        t += demo.time_step
    for _ in range(N_PROFILE):
        demo.step(t)
        t += demo.time_step


def main() -> int:
    out_dir = Path(__file__).parent
    prof_path = out_dir / "profile_perf.prof"
    txt_path = out_dir / "profile_perf.txt"

    print(f"Profiling {N_PROFILE} substeps after {N_WARMUP}-substep warmup...")
    profiler = cProfile.Profile()
    profiler.enable()
    _do_run()
    profiler.disable()

    profiler.dump_stats(str(prof_path))
    print(f"  raw stats -> {prof_path}")

    # Build text summary at three angles.
    out = io.StringIO()
    out.write("=" * 80 + "\n")
    out.write(f"v3 NSM1 cProfile summary (N_PROFILE={N_PROFILE})\n")
    out.write("=" * 80 + "\n\n")

    out.write("--- Top 30 by cumulative time ---\n\n")
    stats = pstats.Stats(profiler, stream=out).strip_dirs()
    stats.sort_stats("cumulative")
    stats.print_stats(30)

    out.write("\n\n--- Top 30 by total (self) time ---\n\n")
    stats.sort_stats("tottime")
    stats.print_stats(30)

    out.write("\n\n--- v3 NSM1 functions only, by cumulative time ---\n\n")
    stats.sort_stats("cumulative")
    stats.print_stats("clearwater_modules_v3", 40)

    out.write("\n\n--- Pattern-G hot-spot candidates ---\n\n")
    stats.sort_stats("cumulative")
    stats.print_stats(r"_change_with_components|_rate_with_components|sanitize_rate|clip_negative_state|setattr|set_at_time", 40)

    txt_path.write_text(out.getvalue())
    print(f"  text summary -> {txt_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
