"""NSM1-DOX-F2 regression: the silent-zero atmospheric-reaeration path
warns, and an opt-in CE-QUAL-W2-MINKL-style floor is available without
changing default behaviour.

Gold-standard spec Workstream C4 (findings DOX-F2; review_SUMMARY §3).

With ``hydraulic_reaeration_option == 1`` + ``kah_20_user == 0.0`` (and
the wind path off), atmospheric reaeration is silently zero. v3
preserves v1/Fortran parity by default (no implicit floor) but must
(a) warn once, and (b) offer the opt-in ``min_reaeration_ka`` floor
(default 0.0 = OFF).

Non-shared-path contract (spec Section 1(4)): the floor value is an
**independently hardcoded** literal, not read from the process.
"""
from __future__ import annotations

import logging
from datetime import datetime, timedelta

import numpy as np
import xarray as xr

from clearwater_modules_v3.processes.dox import DOX

from .conftest import InMemoryRegistry


DT = timedelta(minutes=5)
FLOOR_LITERAL = 2.0  # independently hardcoded opt-in floor (1/d)

# Silent-zero hydraulic + wind config.
_SILENT_ZERO = {
    "hydraulic_reaeration_option": 1,
    "kah_20_user": 0.0,
    "wind_reaeration_option": 1,
    "kaw_20_user": 0.0,
}


def _registry() -> InMemoryRegistry:
    reg = InMemoryRegistry()
    one = lambda v: xr.DataArray(np.array([v]), dims="cell")
    reg.register("oxygen_dissolved", one(6.0))      # undersaturated
    reg.register("water_temperature", one(22.0))
    reg.register("depth", one(1.5))
    return reg


def test_doxf2_default_floor_off_and_silent_zero_warns_once(caplog):
    """Default ``min_reaeration_ka`` is OFF (0.0); the silent-zero path
    produces zero atmospheric reaeration (v1 parity preserved) and warns
    exactly once across substeps."""
    d = DOX(parameters=dict(_SILENT_ZERO), time_step=DT)
    assert d.min_reaeration_ka == 0.0  # default OFF (parity-preserving)

    t0 = datetime(2026, 5, 16)
    with caplog.at_level(logging.WARNING):
        for k in range(3):
            d.run(t0 + k * DT, _registry())

    f2_warnings = [
        r for r in caplog.records if "DOX-F2" in r.message
    ]
    assert len(f2_warnings) == 1, (
        f"expected exactly one DOX-F2 warning, got {len(f2_warnings)}"
    )
    # Parity preserved: no implicit floor -> zero atmospheric reaeration.
    np.testing.assert_allclose(
        np.asarray(d.atm_reaeration_rate), 0.0, atol=1e-12
    )


def test_doxf2_optin_floor_applies():
    """With ``min_reaeration_ka`` > 0 the silent-zero ka is floored, so
    atmospheric reaeration is nonzero (= floor * (DOsat - DOX))."""
    t0 = datetime(2026, 5, 16)

    no_floor = DOX(parameters=dict(_SILENT_ZERO), time_step=DT)
    no_floor.run(t0, _registry())
    assert float(np.asarray(no_floor.atm_reaeration_rate)[0]) == 0.0

    floored = DOX(
        parameters={**_SILENT_ZERO, "min_reaeration_ka": FLOOR_LITERAL},
        time_step=DT,
    )
    floored.run(t0, _registry())

    sat = float(np.asarray(floored.dox_sat)[0])
    dox = 6.0
    expected_atm = FLOOR_LITERAL * (sat - dox)  # ka_tc floored to 2.0
    np.testing.assert_allclose(
        np.asarray(floored.atm_reaeration_rate), expected_atm, rtol=1e-9
    )
    # The opt-in floor strictly increased atmospheric reaeration vs the
    # default (zero) silent-zero path.
    assert float(np.asarray(floored.atm_reaeration_rate)[0]) > 0.0


def test_doxf2_default_options_do_not_warn(caplog):
    """The shipped defaults (hydraulic_reaeration_option=5) are not the
    silent-zero path: no DOX-F2 warning."""
    d = DOX(time_step=DT)  # all defaults
    t0 = datetime(2026, 5, 16)
    with caplog.at_level(logging.WARNING):
        d.run(t0, _registry())
    assert not any("DOX-F2" in r.message for r in caplog.records)
