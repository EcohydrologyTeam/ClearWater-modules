"""Tier 1 regression: v3 NSM1 Processes don't poison state on dry cells.

Addresses the Gemini code-review finding (Phase 9.F follow-up) that
several v3 Process classes divide by ``depth`` (e.g.,
``vsoc / depth * POC``, ``rpo4_tc / depth``, ``vx / depth * PX``,
``SOD_tc / depth``, ``Fb / depth * benthic_term``). At ``depth == 0``
those expressions produce ``inf`` rather than ``NaN``. The
shared ``utils.numerics.sanitize_rate`` helper now catches both
``NaN`` and ``inf`` defensively, even though the orchestration-layer
wet-mask in ``Model`` is the primary line of defense for dry-cell
gating.

These tests run each Process directly with a 5-cell registry where
one cell has ``depth == 0`` and assert that the post-Euler state has
no ``inf`` or ``NaN`` at any cell.
"""
from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

from .conftest import InMemoryRegistry


def _make_dry_cell_registry(in_memory_registry: InMemoryRegistry) -> None:
    """Mutate the registry so cell index 2 has depth = 0; register
    forcings that aren't pre-registered by the standard fixture."""
    depth = in_memory_registry.get("depth")
    new_depth = depth.copy()
    new_depth.values[2] = 0.0
    in_memory_registry.register("depth", new_depth)
    # Register solar_radiation for Pathogen and any other process that
    # reads it; the standard Tier 1 fixture omits it.
    if "solar_radiation" not in in_memory_registry:
        in_memory_registry.register(
            "solar_radiation",
            xr.DataArray(np.full(5, 300.0, dtype=float), dims="cell"),
        )


def _assert_finite_state(state_array: xr.DataArray, name: str) -> None:
    """Assert no inf or NaN in any cell of ``state_array``."""
    finite_mask = np.isfinite(state_array.values)
    assert finite_mask.all(), (
        f"{name}: dry-cell-robustness invariant violated; got "
        f"{state_array.values!r}"
    )


def test_dry_cell_pathogen_does_not_corrupt_state(
    in_memory_registry: InMemoryRegistry,
    closed_system_time_window: tuple[datetime, datetime, timedelta],
) -> None:
    """Pathogen.run with ``depth == 0`` at one cell should not produce
    ``inf`` or ``NaN`` in the post-Euler ``pathogen`` state. Phase 9.F
    follow-up: Pathogen previously had no rate guard at all."""
    from clearwater_modules_v3.processes.pathogen import Pathogen

    _make_dry_cell_registry(in_memory_registry)
    start, _end, time_step = closed_system_time_window

    pathogen = Pathogen(time_step=time_step)
    pathogen.run(start, in_memory_registry)

    _assert_finite_state(
        in_memory_registry.get_at_time("pathogen", start),
        "pathogen",
    )


def test_dry_cell_dox_does_not_corrupt_state(
    in_memory_registry: InMemoryRegistry,
    closed_system_time_window: tuple[datetime, datetime, timedelta],
) -> None:
    """DOX.run with ``depth == 0`` at one cell should not produce
    ``inf`` or ``NaN`` in the post-Euler ``oxygen_dissolved`` state."""
    from clearwater_modules_v3.processes.dox import DOX

    _make_dry_cell_registry(in_memory_registry)
    start, _end, time_step = closed_system_time_window

    dox = DOX(time_step=time_step)
    dox.run(start, in_memory_registry)

    _assert_finite_state(
        in_memory_registry.get_at_time("oxygen_dissolved", start),
        "oxygen_dissolved",
    )


def test_dry_cell_phosphorus_does_not_corrupt_state(
    in_memory_registry: InMemoryRegistry,
    closed_system_time_window: tuple[datetime, datetime, timedelta],
) -> None:
    """Phosphorus.run with ``depth == 0`` at one cell should not produce
    ``inf`` or ``NaN`` in the post-Euler TIP / OrgP states."""
    from clearwater_modules_v3.processes.phosphorus import Phosphorus

    _make_dry_cell_registry(in_memory_registry)
    start, _end, time_step = closed_system_time_window

    phos = Phosphorus(time_step=time_step)
    phos.run(start, in_memory_registry)

    _assert_finite_state(
        in_memory_registry.get_at_time("tip", start),
        "tip",
    )
    _assert_finite_state(
        in_memory_registry.get_at_time("organic_phosphorus", start),
        "organic_phosphorus",
    )


def test_dry_cell_n2_does_not_corrupt_state(
    in_memory_registry: InMemoryRegistry,
    closed_system_time_window: tuple[datetime, datetime, timedelta],
) -> None:
    """N2.run with ``depth == 0`` at one cell should not produce
    ``inf`` or ``NaN`` in the post-Euler ``n2`` state."""
    from clearwater_modules_v3.processes.n2 import N2

    _make_dry_cell_registry(in_memory_registry)
    start, _end, time_step = closed_system_time_window

    n2 = N2(time_step=time_step)
    n2.run(start, in_memory_registry)

    _assert_finite_state(
        in_memory_registry.get_at_time("n2", start),
        "n2",
    )
