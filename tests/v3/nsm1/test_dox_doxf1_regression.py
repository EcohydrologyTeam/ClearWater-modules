"""NSM1-DOX-F1 (MINOR/MAJOR-doc) regression: the freshwater
DO-saturation assumption is explicit, and brackish input cannot pass
silently.

Gold-standard spec Workstream C3 (audit C6).

``dox_sat_apha`` is the fresh-water APHA saturation (no salinity
correction — exact for fresh water, matches v1). The salinity-corrected
form is a documented NSM2 deferral. ``DOX.run`` must (a) leave the
freshwater result numerically unchanged, and (b) emit a one-time
warning if a nonzero ``salinity`` is present in the registry.
"""
from __future__ import annotations

import logging
from datetime import datetime, timedelta

import numpy as np
import xarray as xr

from clearwater_modules_v3.processes.dox import DOX

from .conftest import InMemoryRegistry


DT = timedelta(minutes=5)


def _registry(with_salinity: float | None = None) -> InMemoryRegistry:
    reg = InMemoryRegistry()
    one = lambda v: xr.DataArray(np.array([v]), dims="cell")
    reg.register("oxygen_dissolved", one(7.5))
    reg.register("water_temperature", one(22.0))
    reg.register("depth", one(1.5))
    if with_salinity is not None:
        reg.register("salinity", one(with_salinity))
    return reg


def test_doxf1_freshwater_result_unchanged_and_no_warning(caplog):
    """No salinity (or zero salinity) -> no warning, and dox_sat is the
    fresh-water APHA value (the guard never touches the math)."""
    t = datetime(2026, 5, 16)

    d_none = DOX(time_step=DT)
    with caplog.at_level(logging.WARNING):
        d_none.run(t, _registry(with_salinity=None))
    sat_no_salinity = float(np.asarray(d_none.dox_sat)[0])
    assert not any(
        "salinity" in r.message.lower() for r in caplog.records
    ), "freshwater run must not warn about salinity"

    caplog.clear()
    d_zero = DOX(time_step=DT)
    with caplog.at_level(logging.WARNING):
        d_zero.run(t, _registry(with_salinity=0.0))
    sat_zero_salinity = float(np.asarray(d_zero.dox_sat)[0])
    assert not any(
        "salinity" in r.message.lower() for r in caplog.records
    ), "zero-salinity run must not warn"

    # No numeric change: presence of a zero-salinity variable does not
    # alter the fresh-water saturation.
    np.testing.assert_allclose(sat_zero_salinity, sat_no_salinity, rtol=1e-12)


def test_doxf1_brackish_warns_once_and_no_numeric_change(caplog):
    """Nonzero salinity -> exactly one warning across multiple substeps,
    and dox_sat is still the (unchanged) fresh-water value (the salinity
    correction is a documented deferral, not silently applied)."""
    t0 = datetime(2026, 5, 16)
    t1 = datetime(2026, 5, 16, 0, 5)
    t2 = datetime(2026, 5, 16, 0, 10)

    fresh_ref = DOX(time_step=DT)
    fresh_ref.run(t0, _registry(with_salinity=None))
    fresh_sat = float(np.asarray(fresh_ref.dox_sat)[0])

    d = DOX(time_step=DT)
    with caplog.at_level(logging.WARNING):
        d.run(t0, _registry(with_salinity=35.0))
        d.run(t1, _registry(with_salinity=35.0))
        d.run(t2, _registry(with_salinity=35.0))

    salinity_warnings = [
        r for r in caplog.records if "salinity" in r.message.lower()
    ]
    assert len(salinity_warnings) == 1, (
        f"expected exactly one salinity warning across 3 substeps, "
        f"got {len(salinity_warnings)}"
    )
    assert "NSM1-DOX-F1" in salinity_warnings[0].message or "freshwater" \
        in salinity_warnings[0].message.lower()

    # Brackish input does NOT silently change the math: dox_sat equals
    # the fresh-water value (salinity correction is deferred, audit C6).
    brackish_sat = float(np.asarray(d.dox_sat)[0])
    np.testing.assert_allclose(brackish_sat, fresh_sat, rtol=1e-12)
