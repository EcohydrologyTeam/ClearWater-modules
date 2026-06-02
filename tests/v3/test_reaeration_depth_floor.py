"""Opt-in depth floor in the oxygen-reaeration path.

Guards the ``min_reaeration_depth`` parameter added to
``utils.reaeration.kah_20`` and ``utils.reaeration.ka_tc`` (and threaded
through DOX / N2 / Carbon via the shared ``parameters.dox`` reaeration
menu).

Motivation: the reaeration coefficient carries inverse-depth terms
(``depth**-1.85`` in ``kah_20``; ``kaw_tc / depth`` in ``ka_tc``). A
coupled HEC-RAS-2D transport run delivers a physically faithful per-cell
mean depth that is legitimately ~1e-6 m (or exactly 0) at a wetting
front. Ungated, those cells drive ``ka_tc`` to ~1e11/d and a
Forward-Euler blow-up in the N2 / DIC atmospheric-exchange term.

``min_reaeration_depth`` clamps ``depth = max(depth, min_reaeration_depth)``
BEFORE the inverse-depth terms evaluate. Default ``0.0`` is OFF and
byte-identical to the prior unfloored expression — it mirrors the
existing opt-in CE-QUAL-W2-style ``min_reaeration_ka`` (a floor on the
result) and the TSM ``q_net_depth_skip_threshold`` guards. This test
proves (a) default-off back-compat, (b) the floor clamps a tiny / zero
depth, and (c) the parameter threads through a DOX process run.
"""
from __future__ import annotations

from datetime import datetime

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.utils.reaeration import kah_20, ka_tc
from clearwater_modules_v3.processes.dox import DOX
from clearwater_modules_v3.examples import InMemoryRegistry


def _da(value, n=2):
    return xr.DataArray(np.full(n, value, dtype=float), dims=["cell"])


def _kah(depth, *, min_depth=0.0, option=5, velocity=0.1):
    """``kah_20`` with the non-depth forcings held at benign values so the
    depth-piecewise option-5 (Cover/Owens) branch isolates the depth term."""
    return kah_20(
        kah_20_user=_da(0.0),
        hydraulic_reaeration_option=option,
        velocity=_da(velocity),
        depth=_da(depth),
        flow=_da(1.0),
        topwidth=_da(5.0),
        slope=_da(0.001),
        shear_velocity=_da(0.1),
        min_depth=min_depth,
    )


def _katc_wind_only(depth, *, min_depth=0.0, kaw=2.0):
    """``ka_tc`` with the hydraulic component zeroed and theta=1 at T=20,
    so the result reduces to the wind term ``kaw / depth``."""
    return ka_tc(
        kah_20=_da(0.0),
        kaw_20=_da(kaw),
        kah_theta=_da(1.024),
        kaw_theta=_da(1.024),
        T_water_C=_da(20.0),
        depth=_da(depth),
        min_depth=min_depth,
    )


# ---------------------------------------------------------------------------
# kah_20 (hydraulic inverse-depth term)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("depth", [0.05, 0.5, 1.0, 3.0])
def test_kah_20_default_off_is_backcompat(depth):
    """min_depth=0.0 (default) reproduces the unfloored output exactly."""
    base = _kah(depth)
    explicit_off = _kah(depth, min_depth=0.0)
    np.testing.assert_array_equal(
        np.asarray(base.values), np.asarray(explicit_off.values)
    )


def test_kah_20_no_op_when_depth_above_floor():
    """A floor below the actual depth leaves the result untouched."""
    deep = _kah(1.0, min_depth=0.0)
    floored = _kah(1.0, min_depth=0.01)
    np.testing.assert_array_equal(
        np.asarray(deep.values), np.asarray(floored.values)
    )


def test_kah_20_floor_clamps_tiny_depth():
    """At a sub-physical depth the floor makes kah_20 equal the value AT
    the floor depth, and that is vastly smaller than the unfloored spike."""
    floored = np.asarray(_kah(1e-6, min_depth=0.01).values)
    at_floor = np.asarray(_kah(0.01, min_depth=0.0).values)
    unfloored = np.asarray(_kah(1e-6, min_depth=0.0).values)

    np.testing.assert_allclose(floored, at_floor, rtol=1e-12, atol=0.0)
    # The unfloored spike is enormous (depth**-1.85 at 1e-6 m).
    assert np.all(unfloored > 1e9)
    assert np.all(floored < unfloored / 1e6)


# ---------------------------------------------------------------------------
# ka_tc (wind inverse-depth term)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("depth", [0.05, 0.5, 1.0, 3.0])
def test_ka_tc_default_off_is_backcompat(depth):
    base = _katc_wind_only(depth)
    explicit_off = _katc_wind_only(depth, min_depth=0.0)
    np.testing.assert_array_equal(
        np.asarray(base.values), np.asarray(explicit_off.values)
    )


def test_ka_tc_floor_clamps_tiny_depth():
    floored = np.asarray(_katc_wind_only(1e-6, min_depth=0.01).values)
    at_floor = np.asarray(_katc_wind_only(0.01, min_depth=0.0).values)
    unfloored = np.asarray(_katc_wind_only(1e-6, min_depth=0.0).values)

    np.testing.assert_allclose(floored, at_floor, rtol=1e-12, atol=0.0)
    # kaw=2 m/d at depth 0.01 m -> 200/d; unfloored at 1e-6 m -> 2e6/d.
    np.testing.assert_allclose(floored, 200.0, rtol=1e-9)
    assert np.all(unfloored > 1e6)


def test_ka_tc_floor_handles_zero_depth():
    """depth == 0 yields inf/nan unfloored; the floor makes it finite."""
    unfloored = np.asarray(_katc_wind_only(0.0, min_depth=0.0).values)
    floored = np.asarray(_katc_wind_only(0.0, min_depth=0.01).values)
    assert not np.all(np.isfinite(unfloored))
    assert np.all(np.isfinite(floored))
    np.testing.assert_allclose(floored, 200.0, rtol=1e-9)


# ---------------------------------------------------------------------------
# DOX process integration (parameter threads through the reaeration menu)
# ---------------------------------------------------------------------------


def _dox_hydraulic_only(min_reaeration_depth=0.0):
    """DOX with only the hydraulic reaeration path active (option 5,
    wind off), so atm_reaeration_rate isolates the depth-driven kah_20."""
    return DOX(
        parameters={
            "hydraulic_reaeration_option": 5,
            "wind_reaeration_option": 1, "kaw_20_user": 0.0,  # wind off
            "min_reaeration_depth": min_reaeration_depth,
        },
    )


def _dox_registry(depth, n=2):
    reg = InMemoryRegistry()
    reg.register("oxygen_dissolved", _da(5.0, n))   # below saturation
    reg.register("water_temperature", _da(25.0, n))
    reg.register("depth", _da(depth, n))
    reg.register("ammonium", _da(0.1, n))
    reg.register("wind_speed", _da(0.0, n))
    return reg


def _atm_rate(proc, reg):
    proc.run(datetime(2026, 1, 1), reg)
    return np.asarray(proc.atm_reaeration_rate, dtype=float)


def test_dox_min_reaeration_depth_default_off_backcompat():
    """No / zero min_reaeration_depth reproduces the prior (unfloored) run."""
    rate_absent = _atm_rate(DOX(parameters={
        "hydraulic_reaeration_option": 5,
        "wind_reaeration_option": 1, "kaw_20_user": 0.0,
    }), _dox_registry(1e-6))
    rate_explicit_off = _atm_rate(_dox_hydraulic_only(0.0), _dox_registry(1e-6))
    np.testing.assert_array_equal(rate_absent, rate_explicit_off)


def test_dox_min_reaeration_depth_prevents_blowup():
    """min_reaeration_depth clamps the sub-physical newly-wet depth: the
    floored reaeration matches a run AT the floor depth and is orders of
    magnitude below the unfloored spike."""
    rate_unfloored = _atm_rate(_dox_hydraulic_only(0.0), _dox_registry(1e-6))
    rate_floored = _atm_rate(_dox_hydraulic_only(0.01), _dox_registry(1e-6))
    rate_at_floor = _atm_rate(_dox_hydraulic_only(0.0), _dox_registry(0.01))

    np.testing.assert_allclose(rate_floored, rate_at_floor, rtol=1e-9, atol=0.0)
    assert np.all(np.abs(rate_unfloored) > 1e6)
    assert np.all(np.abs(rate_floored) < np.abs(rate_unfloored) / 1e6)
