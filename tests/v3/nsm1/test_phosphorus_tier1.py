"""Tier 1 closed-system mass-conservation tests for v3 NSM1 Phosphorus.

Phase 4 of the v3 NSM1 implementation plan (design spec Section 11
Phase 4, Section 9 Tier 1 contract). Asserts that the v3-native
``Phosphorus`` Process is conservative (total P invariant) when:

* Settling is disabled (``vs=0``, ``vsop=0``).
* Sediment release is disabled (``rpo4_20=0``; the v3 default).
* No algae are present in the model (no ``FloatingAlgae`` /
  ``BenthicAlgae`` couplings).

Under these conditions, the only kinetic pathway is OrgP <-> TIP
hydrolysis (``kop_tc * OrgP``), which is mass-conserving within the
total-P pool: every milligram of P that leaves OrgP enters TIP.

Tier 1 invariants:

* ``total_p = TIP + OrgP`` constant per cell to roundoff
  (``rtol=1e-12``).
* ``diagnostics.clip_events == {}``: a clip event under closed-system
  conditions signals either a malformed test or an integrator bug.

This test is the analogue of
``test_tier1_total_n_conservation_closed_system_nitrogen_only`` for
Phosphorus.
"""

from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np

from clearwater_modules_v3.processes.phosphorus import Phosphorus
from clearwater_modules_v3.utils.numerics import Diagnostics

from .conftest import InMemoryRegistry, total_p


def test_tier1_phosphorus_conservation_closed_system_loss_disabled(
    in_memory_registry: InMemoryRegistry,
    closed_system_time_window: tuple[datetime, datetime, timedelta],
) -> None:
    """Closed-system total-P conservation when settling and sediment
    release are disabled, no algae in the model. The TIP <-> OrgP
    hydrolysis is internal balanced; total P is invariant.

    Setup:
    * 5-cell mesh; ``initial_state_5cell`` initial conditions
      (``tip = [0.10, 0.12, 0.14, 0.16, 0.18]`` mg-P/L,
      ``organic_phosphorus = [0.05, 0.06, 0.07, 0.08, 0.09]`` mg-P/L).
    * ``vs=0`` disables TIP settling.
    * ``vsop=0`` disables OrgP settling.
    * ``rpo4_20=0`` disables sediment P release (v3 default; explicit
      for clarity).
    * No algae are wired into the Process: ``use_floating_algae`` and
      ``use_benthic_algae`` remain False, so all algal-coupling source/
      sink terms collapse to zero.
    * ``kop_20`` left at v3 default (0.1 1/d) so the hydrolysis
      pathway is exercised end-to-end. Hydrolysis is mass-conserving
      within total-P (every mg-P that leaves OrgP enters TIP).

    Expected (Phase 4):
    * Per-cell ``TIP + OrgP`` invariant to roundoff (``rtol=1e-12``).
    * Total-P helper invariant to roundoff (``rtol=1e-12``).
    * ``diagnostics.clip_events == {}``.
    """
    start, end, time_step = closed_system_time_window

    phosphorus = Phosphorus(
        parameters={
            # Disable settling for both states (closed-system).
            "vs": 0.0,
            "vsop": 0.0,
            # Disable sediment release (already 0 by default; explicit).
            "rpo4_20": 0.0,
            # Hydrolysis pathway active at the v3 default kop_20 = 0.1 1/d.
        },
        time_step=time_step,
    )
    # Wire a local Diagnostics so the test owns the reference for the
    # ``clip_events == {}`` assertion. ``Phosphorus.__init__`` already
    # creates a fresh Diagnostics; we replace it here so the test owns it.
    diagnostics = Diagnostics()
    phosphorus.diagnostics = diagnostics

    # No algae in the model: leave use_floating_algae / use_benthic_algae
    # at their __init__ defaults (False). This is the standalone Tier 1
    # mode; ``init_process`` is not called.
    assert phosphorus.use_floating_algae is False
    assert phosphorus.use_benthic_algae is False

    # Snapshot the initial state for per-cell and aggregate comparisons.
    tip_initial = in_memory_registry.get_at_time("tip", start).copy()
    orgp_initial = in_memory_registry.get_at_time(
        "organic_phosphorus", start
    ).copy()
    total_p_initial = float(total_p(in_memory_registry).values)

    # Run 100 substeps. Phosphorus.run mutates the registry in place via
    # ``set_at_time`` for both ``tip`` and ``organic_phosphorus``.
    current_time = start
    while current_time < end:
        phosphorus.run(current_time, in_memory_registry)
        current_time += time_step

    # Tier 1 invariant 1: per-cell total-P (TIP + OrgP) constant to roundoff.
    tip_final = in_memory_registry.get_at_time("tip", end)
    orgp_final = in_memory_registry.get_at_time("organic_phosphorus", end)
    per_cell_initial = tip_initial.values + orgp_initial.values
    per_cell_final = tip_final.values + orgp_final.values
    np.testing.assert_allclose(
        per_cell_final,
        per_cell_initial,
        rtol=1e-12,
        err_msg=(
            "Closed-system per-cell total-P conservation failed. "
            f"initial={per_cell_initial!r}, final={per_cell_final!r}, "
            f"absolute drift={(per_cell_final - per_cell_initial)!r}"
        ),
    )

    # Tier 1 invariant 2: aggregate total-P (helper-derived) constant.
    # Because no algae are in the registry-derived total_p, this is
    # equivalent to summing per-cell (TIP + OrgP).
    total_p_final = float(total_p(in_memory_registry).values)
    np.testing.assert_allclose(
        total_p_final,
        total_p_initial,
        rtol=1e-12,
        err_msg=(
            "Closed-system aggregate total-P conservation failed. "
            f"initial={total_p_initial!r}, final={total_p_final!r}, "
            f"absolute drift={(total_p_final - total_p_initial)!r}"
        ),
    )

    # Tier 1 invariant 3: no clipping under closed-system Tier 1 conditions.
    assert diagnostics.clip_events == {}, (
        f"Clip events fired under closed-system Tier 1 Phosphorus "
        f"conditions: {diagnostics.clip_events!r}. The clip log is "
        f"{diagnostics.clip_log!r}."
    )


# ---------------------------------------------------------------------------
# Phosphorus instantiation smoke tests
# ---------------------------------------------------------------------------


def test_phosphorus_instantiates_with_defaults() -> None:
    """``Phosphorus()`` constructs cleanly with no arguments and pulls
    PHOSPHORUS_DEFAULTS onto the instance."""
    phosphorus = Phosphorus()
    # Spot-check phosphorus DEFAULTS keys.
    assert phosphorus.kop_20 == 0.1
    assert phosphorus.kop_theta == 1.047
    assert phosphorus.rpo4_20 == 0.0
    assert phosphorus.kdpo4 == 0.0
    assert phosphorus.vsop == 0.1
    assert phosphorus.vs == 0.1
    # Spot-check inline partitioning fallback defaults.
    assert phosphorus.use_TIP is True
    assert phosphorus.use_OrgP is True


def test_phosphorus_accepts_parameter_override() -> None:
    """``Phosphorus(parameters={'kop_20': 0.5})`` overrides one entry
    and leaves the rest at their DEFAULTS values."""
    phosphorus = Phosphorus(parameters={"kop_20": 0.5})
    assert phosphorus.kop_20 == 0.5
    # Others unchanged.
    assert phosphorus.kop_theta == 1.047
    assert phosphorus.vsop == 0.1
    assert phosphorus.vs == 0.1
