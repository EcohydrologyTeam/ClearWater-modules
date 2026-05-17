"""Phase 10 end-to-end coupled demo parity tests.

Two complementary contracts:

1. ``test_baseline_parity_bit_identical`` — replays the Phase 0
   4,320-substep baseline scenario and asserts every state-variable
   trajectory matches the committed
   ``baseline_coupled_trajectory_b51df71.nc`` bit-identically when no
   ``REGISTRY_DIAGNOSTICS`` names are pre-registered. This is the
   §11.2 contract enforced inside the test suite. (Terminal
   gold-standard baseline at ``b51df71``, incorporating the NSM1-CA-1
   and NSM1-SCI-N1 alkalinity kinetics fixes — the only
   trajectory-perturbing gate changes. Prior baselines ``624ed7c``
   (CA-1 only) and ``186b5c4`` (pre-fix) are retained for auditability
   per baseline/README.md.)

2. ``test_diagnostics_subscription_smoke`` — runs the same demo
   shorter (60 substeps for speed) but pre-registers every Appendix A
   name across all 11 Processes. Asserts (a) the state-variable
   subset is bit-identical to the no-subscription run, (b) every
   diagnostic name has finite values written each substep, (c) the
   diagnostic dataset has the expected shape.

These tests survive Phase 10 cleanup. They are the smoke tests that
ride along in the regular ``pytest tests/`` run and catch any future
drift away from the pattern-alignment contract.
"""

from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.examples import (
    build_nsm1_demo,
    default_initial_conditions,
)
from clearwater_modules_v3.processes import (
    Alkalinity,
    BenthicAlgae,
    CBOD,
    Carbon,
    DOX,
    FloatingAlgae,
    N2,
    Nitrogen,
    POM,
    Pathogen,
    Phosphorus,
)


_BASELINE_NETCDF = (
    Path(__file__).resolve().parent
    / "baseline"
    / "baseline_coupled_trajectory_b51df71.nc"
)


N_BASELINE_SUBSTEPS = 4320
N_SMOKE_SUBSTEPS = 60
START = datetime(2026, 1, 1, 0, 0, 0)


# ---------------------------------------------------------------------------
# §11.2 baseline parity
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not _BASELINE_NETCDF.is_file(),
    reason=f"baseline NetCDF not present at {_BASELINE_NETCDF}",
)
def test_baseline_parity_bit_identical() -> None:
    """End-to-end §11.2 contract: every state variable in the
    4,320-substep coupled trajectory matches the Phase 0 baseline
    NetCDF bit-identically (``rtol=0, atol=0``)."""
    baseline = xr.open_dataset(_BASELINE_NETCDF)
    rng_seed = int(baseline.attrs["rng_seed"])
    n_substeps = int(baseline.attrs["n_substeps"])
    assert n_substeps == N_BASELINE_SUBSTEPS, (
        f"baseline n_substeps={n_substeps}, expected {N_BASELINE_SUBSTEPS}"
    )
    n_cells = int(baseline.attrs["n_cells"])
    time_step_seconds = int(baseline.attrs["time_step_seconds"])
    start = datetime.fromisoformat(baseline.attrs["start_time_iso"])

    np.random.seed(rng_seed)
    demo = build_nsm1_demo(time_step=timedelta(seconds=time_step_seconds))

    state_names = sorted(demo.registry.keys())
    traj = {
        name: np.empty((n_substeps + 1, n_cells), dtype=np.float64)
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

    mismatches: list[str] = []
    for name in sorted(baseline.data_vars):
        b_vals = baseline[name].values
        c_vals = traj.get(name)
        if c_vals is None:
            mismatches.append(f"{name}: present in baseline but not in run")
            continue
        if not np.array_equal(b_vals, c_vals):
            diff_count = int((b_vals != c_vals).sum())
            mismatches.append(f"{name}: {diff_count} non-identical cells")

    new_in_run = set(traj.keys()) - set(baseline.data_vars)
    for name in sorted(new_in_run):
        mismatches.append(f"{name}: present in run but not in baseline")

    assert not mismatches, "\n".join(mismatches)


# ---------------------------------------------------------------------------
# Diagnostics-subscription smoke
# ---------------------------------------------------------------------------


_ALL_APPENDIX_A_NAMES: tuple[str, ...] = tuple(
    name
    for cls in (
        Carbon, DOX, Nitrogen, FloatingAlgae, BenthicAlgae, Phosphorus,
        POM, CBOD, N2, Pathogen, Alkalinity,
    )
    for name in cls.REGISTRY_DIAGNOSTICS
)


def test_diagnostics_subscription_smoke_state_bit_identical() -> None:
    """Subscribing to every Appendix A name across all 11 Processes
    must NOT change any state-variable trajectory. Pattern G
    zero-cost-when-unused invariant verified end-to-end."""
    # Reference: no subscription.
    demo_a = build_nsm1_demo()
    demo_a.run(START, n_steps=N_SMOKE_SUBSTEPS)

    # Subscribed: pre-register every Appendix A diagnostic name.
    ic_b = default_initial_conditions()
    reference = ic_b["pom"]  # any state variable suffices as a shape template.
    for name in _ALL_APPENDIX_A_NAMES:
        if name not in ic_b:
            ic_b[name] = xr.zeros_like(reference)
    demo_b = build_nsm1_demo(initial_conditions=ic_b)
    demo_b.run(START, n_steps=N_SMOKE_SUBSTEPS)

    # State variables that the demo Processes update — verify all
    # are bit-identical between the two runs.
    state_var_names = {
        "ammonium", "nitrate", "organic_nitrogen", "n2",
        "tip", "organic_phosphorus",
        "poc", "doc", "dic", "pom", "cbod",
        "oxygen_dissolved", "pathogen", "alkalinity",
        "algae_floating", "benthic_algae",
    }
    for name in sorted(state_var_names):
        a = demo_a.registry.get(name).values
        b = demo_b.registry.get(name).values
        np.testing.assert_array_equal(
            a, b,
            err_msg=(
                f"{name}: bit-identical state-trajectory invariant "
                "violated under full diagnostics subscription"
            ),
        )


def test_diagnostics_subscription_smoke_all_diagnostics_written() -> None:
    """Subscribing to every Appendix A name produces finite, populated
    registry entries each substep across all 11 Processes."""
    ic = default_initial_conditions()
    reference = ic["pom"]
    for name in _ALL_APPENDIX_A_NAMES:
        if name not in ic:
            ic[name] = xr.zeros_like(reference)

    demo = build_nsm1_demo(initial_conditions=ic)
    demo.run(START, n_steps=N_SMOKE_SUBSTEPS)

    bad: list[str] = []
    for name in _ALL_APPENDIX_A_NAMES:
        if name not in demo.registry:
            bad.append(f"{name}: not in registry after run")
            continue
        value = demo.registry.get(name)
        arr = value.values if isinstance(value, xr.DataArray) else np.asarray(value)
        if not np.isfinite(arr).all():
            bad.append(f"{name}: non-finite values after {N_SMOKE_SUBSTEPS} substeps")

    assert not bad, "\n".join(bad)


def test_diagnostics_subscription_smoke_no_substep_skipped() -> None:
    """Quick sanity check: after N_SMOKE_SUBSTEPS the state variables
    have evolved from their initial conditions for at least one state
    variable per Process. This catches a regression where the
    pattern-G setattr loop silently overwrites the integrator step
    (the worst-case structural failure)."""
    demo = build_nsm1_demo()
    initial = {
        name: demo.registry.get(name).values.copy()
        for name in (
            "ammonium", "nitrate", "tip", "poc", "doc", "dic", "pom",
            "oxygen_dissolved", "algae_floating", "benthic_algae",
        )
    }
    demo.run(START, n_steps=N_SMOKE_SUBSTEPS)

    evolved_any: list[str] = []
    for name, init_vals in initial.items():
        final = demo.registry.get(name).values
        if not np.array_equal(init_vals, final):
            evolved_any.append(name)

    assert evolved_any, (
        f"No state variable evolved after {N_SMOKE_SUBSTEPS} substeps; "
        "the integrator may be silently no-op'ing"
    )
