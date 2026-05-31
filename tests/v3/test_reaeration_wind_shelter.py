"""Wind-shelter consumption in the oxygen-reaeration path.

Guards `design/clearwater_modules_v3_dox_wind_shelter_fix.md`: the
wind-driven reaeration velocity ``utils.reaeration.kaw_20`` now accepts a
``wind_shelter`` coefficient and applies it to the raw wind BEFORE the
height rescale — the same composition order the TSM wind function uses
(``raw * shelter * height_factor``; CE-QUAL-W2 ``w2_4_unix.f90:480``).
``DOX.run`` reads the optional per-cell ``wind_shelter_coefficient``
forcing (the same one Temperature consumes) and threads it through, so a
sheltered cell gets reduced wind for gas transfer as well as heat
exchange. Default ``1.0`` (no forcing registered) preserves prior output.
"""
from __future__ import annotations

from datetime import datetime

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.utils.reaeration import kaw_20
from clearwater_modules_v3.processes.dox import DOX
from clearwater_modules_v3.examples import InMemoryRegistry


def _da(value, n=2):
    return xr.DataArray(np.full(n, value, dtype=float), dims=["cell"])


# ---------------------------------------------------------------------------
# kaw_20 unit tests (design spec 4.1)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("option", [2, 3, 5, 7])
def test_shelter_default_is_backcompat(option):
    """wind_shelter=1.0 (default) reproduces the pre-change output."""
    u = _da(4.0)
    base = kaw_20(_da(0.0), u, option)
    with_default = kaw_20(_da(0.0), u, option, wind_shelter=1.0)
    np.testing.assert_array_equal(
        np.asarray(base.values), np.asarray(with_default.values)
    )


@pytest.mark.parametrize("option", [2, 3, 5, 7])
def test_shelter_equivalent_to_premultiplying_wind(option):
    """kaw_20(U, shelter=s) == kaw_20(s*U, shelter=1): the shelter is
    mathematically equivalent to scaling the raw wind (holds for the
    piecewise options too, since Uw10 is identical)."""
    u = _da(6.0)
    s = 0.5
    sheltered = kaw_20(_da(0.0), u, option, wind_shelter=s)
    premultiplied = kaw_20(_da(0.0), u * s, option, wind_shelter=1.0)
    np.testing.assert_allclose(
        np.asarray(sheltered.values), np.asarray(premultiplied.values),
        rtol=1e-12, atol=0.0,
    )


def test_shelter_per_cell_array_broadcasts():
    """A per-cell wind_shelter DataArray broadcasts over wind_speed."""
    u = xr.DataArray([5.0, 5.0, 5.0], dims=["cell"])
    shelter = xr.DataArray([1.0, 0.5, 0.0], dims=["cell"])
    out = kaw_20(xr.zeros_like(u), u, 5, wind_shelter=shelter)
    out = np.asarray(out.values)
    # cell 0 (shelter 1) = unsheltered; cell 1 (0.5) = reduced; cell 2 (0) = 0 wind.
    full = float(np.asarray(kaw_20(_da(0.0, 1), _da(5.0, 1), 5).values)[0])
    np.testing.assert_allclose(out[0], full, rtol=1e-12)
    assert out[1] < out[0]
    np.testing.assert_allclose(
        out[2], float(np.asarray(kaw_20(_da(0.0, 1), _da(0.0, 1), 5).values)[0]),
        rtol=1e-12,
    )


def test_shelter_applied_before_option3_threshold():
    """Option 3 is piecewise at Uw10 = 3.5 m/s. At wind_input_height=10
    (identity rescale), wind=4.0 with shelter=0.5 -> Uw10=2.0 -> the
    low-wind branch, matching an unsheltered 2.0 m/s wind."""
    sheltered = kaw_20(
        _da(0.0), _da(4.0), 3, wind_input_height=10.0, wind_shelter=0.5
    )
    low_branch = kaw_20(
        _da(0.0), _da(2.0), 3, wind_input_height=10.0, wind_shelter=1.0
    )
    np.testing.assert_allclose(
        np.asarray(sheltered.values), np.asarray(low_branch.values),
        rtol=1e-12, atol=0.0,
    )
    # And it must differ from the high-wind branch (unsheltered 4.0).
    high_branch = kaw_20(_da(0.0), _da(4.0), 3, wind_input_height=10.0)
    assert not np.allclose(
        np.asarray(sheltered.values), np.asarray(high_branch.values)
    )


# ---------------------------------------------------------------------------
# DOX integration (design spec 4.2)
# ---------------------------------------------------------------------------


def _dox_wind_only():
    """DOX with only the wind reaeration path active (hydraulic off), so
    atm_reaeration_rate isolates the sheltered wind term."""
    return DOX(
        parameters={
            "hydraulic_reaeration_option": 1, "kah_20_user": 0.0,
            "wind_reaeration_option": 5, "kaw_20_user": 0.0,  # Wanninkhof
        },
    )


def _dox_registry(*, wind_speed, shelter=None, n=2):
    reg = InMemoryRegistry()
    reg.register("oxygen_dissolved", _da(5.0, n))   # below saturation
    reg.register("water_temperature", _da(25.0, n))
    reg.register("depth", _da(1.0, n))
    reg.register("ammonium", _da(0.1, n))
    reg.register("wind_speed", _da(wind_speed, n))
    if shelter is not None:
        reg.register("wind_shelter_coefficient", _da(shelter, n))
    return reg


def _reaeration_rate(proc, reg):
    proc.run(datetime(2026, 1, 1), reg)
    return np.asarray(proc.atm_reaeration_rate)


def test_dox_registered_shelter_equals_reduced_wind():
    """Registering wind_shelter_coefficient=s gives the same DOX
    reaeration as an unsheltered run at s * wind_speed; and a shelter
    actually reduces reaeration vs the unsheltered full-wind run."""
    rate_sheltered = _reaeration_rate(
        _dox_wind_only(), _dox_registry(wind_speed=6.0, shelter=0.5)
    )
    rate_reduced_wind = _reaeration_rate(
        _dox_wind_only(), _dox_registry(wind_speed=3.0, shelter=None)
    )
    rate_full_wind = _reaeration_rate(
        _dox_wind_only(), _dox_registry(wind_speed=6.0, shelter=None)
    )

    # shelter=0.5 at U=6 == no shelter at U=3
    np.testing.assert_allclose(
        rate_sheltered, rate_reduced_wind, rtol=1e-12, atol=0.0
    )
    # shelter reduced reaeration relative to full wind (nonzero effect)
    assert np.all(np.abs(rate_sheltered) < np.abs(rate_full_wind))


def test_dox_no_shelter_registered_is_backcompat():
    """No wind_shelter_coefficient registered -> shelter defaults to 1.0,
    output identical to a run that explicitly registers shelter=1.0."""
    rate_absent = _reaeration_rate(
        _dox_wind_only(), _dox_registry(wind_speed=6.0, shelter=None)
    )
    rate_unity = _reaeration_rate(
        _dox_wind_only(), _dox_registry(wind_speed=6.0, shelter=1.0)
    )
    np.testing.assert_array_equal(rate_absent, rate_unity)


# ---------------------------------------------------------------------------
# N2 and Carbon wind-reaeration paths (design spec 3.3 follow-up)
# ---------------------------------------------------------------------------

from clearwater_modules_v3.processes.n2 import N2
from clearwater_modules_v3.processes.carbon import Carbon

_WIND_ONLY = {
    "hydraulic_reaeration_option": 1, "kah_20_user": 0.0,
    "wind_reaeration_option": 5, "kaw_20_user": 0.0,  # Wanninkhof
}


def _n2_registry(*, wind_speed, shelter=None, n=2):
    reg = InMemoryRegistry()
    reg.register("n2", _da(10.0, n))            # below N2 saturation
    reg.register("water_temperature", _da(25.0, n))
    reg.register("depth", _da(1.0, n))
    reg.register("atmospheric_pressure", _da(1013.0, n))
    reg.register("wind_speed", _da(wind_speed, n))
    if shelter is not None:
        reg.register("wind_shelter_coefficient", _da(shelter, n))
    return reg


def _carbon_registry(*, wind_speed, shelter=None, n=2):
    reg = InMemoryRegistry()
    reg.register("poc", _da(4.0, n))
    reg.register("doc", _da(2.0, n))
    reg.register("dic", _da(5.0, n))
    reg.register("water_temperature", _da(25.0, n))
    reg.register("depth", _da(1.0, n))
    reg.register("oxygen_dissolved", _da(8.0, n))
    reg.register("wind_speed", _da(wind_speed, n))
    if shelter is not None:
        reg.register("wind_shelter_coefficient", _da(shelter, n))
    return reg


def test_n2_registered_shelter_equals_reduced_wind():
    """N2 atmospheric exchange: shelter=0.5 at U=6 == no shelter at U=3."""
    def rate(ws, sh):
        p = N2(parameters=dict(_WIND_ONLY))
        p.run(datetime(2026, 1, 1), _n2_registry(wind_speed=ws, shelter=sh))
        return np.asarray(p.n2_atm_exchange_rate)

    np.testing.assert_allclose(
        rate(6.0, 0.5), rate(3.0, None), rtol=1e-12, atol=0.0
    )
    assert np.all(np.abs(rate(6.0, 0.5)) < np.abs(rate(6.0, None)))


def test_carbon_registered_shelter_equals_reduced_wind():
    """Carbon DIC reaeration (via _ka_tc -> kaw_20): shelter=0.5 at U=6 ==
    no shelter at U=3. Also exercises the registry wind_speed now threaded
    to _ka_tc."""
    def rate(ws, sh):
        p = Carbon(parameters=dict(_WIND_ONLY))
        p.run(datetime(2026, 1, 1), _carbon_registry(wind_speed=ws, shelter=sh))
        return np.asarray(p.dic_atm_exchange_rate)

    np.testing.assert_allclose(
        rate(6.0, 0.5), rate(3.0, None), rtol=1e-12, atol=0.0
    )
    assert np.all(np.abs(rate(6.0, 0.5)) < np.abs(rate(6.0, None)))
