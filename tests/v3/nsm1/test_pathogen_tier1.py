"""Tier 1 closed-system mass-conservation test for v3 NSM1 Pathogen.

Phase 3.1 (v3 NSM1 design spec, Section 11): asserts that the v3-native
``Pathogen`` Process is conservative (PX invariant) when all loss terms
are zeroed out — natural decay (``kdx_20=0``), light-induced decay
(``apx=0``), and settling (``vx=0``). Under these conditions:

* The PX state must equal its initial value at every cell to roundoff
  (``rtol=1e-12``).
* The clip-with-log diagnostics must remain empty
  (``diagnostics.clip_events == {}``); a clip event under closed-system
  conditions signals either a malformed test or an integrator bug.

This test is the analogue of
``test_tier1_total_n_conservation_closed_system_nitrogen_only`` for the
new Pathogen constituent. Pathogen is independent (no coupling to other
Phase 3 Processes), so the assertion is per-cell PX equality rather
than a derived total-X invariant.
"""

from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np

from clearwater_modules_v3.processes.pathogen import Pathogen
from clearwater_modules_v3.utils.numerics import Diagnostics

from .conftest import InMemoryRegistry


def test_tier1_pathogen_conservation_closed_system_decay_disabled(
    in_memory_registry: InMemoryRegistry,
    closed_system_time_window: tuple[datetime, datetime, timedelta],
) -> None:
    """Closed-system PX conservation when all loss terms are disabled.

    Setup:
    * 5-cell mesh; ``initial_state_5cell`` initial conditions
      (``pathogen`` initialized to ``[1e3, 5e3, 1e4, 5e4, 1e5]``).
    * ``kdx_20=0`` disables natural decay.
    * ``apx=0`` disables light-induced decay (the Process still
      computes PAR and KEXT, but the rate is zeroed by ``apx*...*PX``).
    * ``vx=0`` disables settling.
    * ``solar_radiation`` is registered with a non-zero value so the
      light path is exercised even though ``apx=0`` zeroes the rate;
      this verifies that the light path is numerically stable
      (``KEXT*depth > 0``, no NaN, no clip).

    Expected:
    * ``pathogen_final == pathogen_initial`` per cell (rtol=1e-12).
    * ``diagnostics.clip_events == {}``.
    """
    start, end, time_step = closed_system_time_window

    # Register solar_radiation since it's not in the default
    # initial_state_5cell fixture but Pathogen.run reads it.
    import xarray as xr  # local import; conftest doesn't expose xr
    in_memory_registry.register(
        "solar_radiation",
        xr.DataArray(np.full(5, 300.0, dtype=float), dims="cell"),
    )

    pathogen = Pathogen(
        parameters={
            "kdx_20": 0.0,   # disable natural decay
            "apx": 0.0,      # disable light-induced decay
            "vx": 0.0,       # disable settling
        },
        time_step=time_step,
    )
    # Wire the local diagnostics so we can assert ``clip_events == {}``
    # afterwards. ``Pathogen.__init__`` already creates a fresh
    # Diagnostics on ``self.diagnostics``; we replace it here so the
    # test owns the reference.
    diagnostics = Diagnostics()
    pathogen.diagnostics = diagnostics

    # Snapshot the initial PX state for the per-cell equality check.
    px_initial = in_memory_registry.get_at_time("pathogen", start).copy()

    # Run 100 substeps. Pathogen.run mutates the registry in place via
    # ``set_at_time``.
    current_time = start
    while current_time < end:
        pathogen.run(current_time, in_memory_registry)
        current_time += time_step

    # Tier 1 invariant 1: PX state per-cell invariance under
    # zeroed-loss closed-system conditions.
    px_final = in_memory_registry.get_at_time("pathogen", end)
    np.testing.assert_allclose(
        px_final.values,
        px_initial.values,
        rtol=1e-12,
        err_msg=(
            "Closed-system PX invariance failed under "
            "kdx_20=apx=vx=0. "
            f"initial={px_initial.values!r}, "
            f"final={px_final.values!r}, "
            f"absolute drift={(px_final.values - px_initial.values)!r}"
        ),
    )

    # Tier 1 invariant 2: no clipping under closed-system conditions.
    assert diagnostics.clip_events == {}, (
        f"Clip events fired under closed-system Tier 1 Pathogen "
        f"conditions: {diagnostics.clip_events!r}. The clip log is "
        f"{diagnostics.clip_log!r}."
    )


# ---------------------------------------------------------------------------
# Pathogen instantiation smoke tests
# ---------------------------------------------------------------------------


def test_pathogen_instantiates_with_defaults() -> None:
    """``Pathogen()`` constructs cleanly with no arguments and pulls
    PATHOGEN_DEFAULTS onto the instance."""
    pathogen = Pathogen()
    # Spot-check a few merged DEFAULTS keys.
    assert pathogen.kdx_20 == 0.8
    assert pathogen.kdx_theta == 1.07
    assert pathogen.apx == 1.0
    assert pathogen.vx == 1.0


def test_pathogen_accepts_parameter_override() -> None:
    """``Pathogen(parameters={'kdx_20': 0.5})`` overrides one entry and
    leaves the rest at their DEFAULTS values."""
    pathogen = Pathogen(parameters={"kdx_20": 0.5})
    assert pathogen.kdx_20 == 0.5
    # Others unchanged
    assert pathogen.kdx_theta == 1.07
    assert pathogen.apx == 1.0
    assert pathogen.vx == 1.0
