"""Parity tests: v3 N2 sub-rate methods/helpers vs v1 nsm1.processes helpers.

Each test exercises one of the v3 N2 stateless helpers (``khn2_tc``,
``pwv``, ``n2sat_henry``) or the ``N2.run`` cached attributes
(``n2_sat``, ``n2_atm_exchange_rate``) against the equivalent v1
reference function.

v1 reference: ``clearwater_modules.nsm1.processes`` ``KHN2_tc``,
``pwv``, ``N2sat``, ``dN2dt``.

Notes on small numerical conventions:

* v1 ``celsius_to_kelvin`` uses ``+273.16`` whereas v3
  ``n2._kelvin`` uses ``+273.15``. To isolate the kinetics formulas
  from this 0.01 K offset, the parity tests below pass matched Kelvin
  values directly into both v1 and v3 (the stateless-helper tests use
  the v1 +273.16 convention, and the run-cached tests construct the
  v1 reference with +273.15 to match what v3 used internally).
* v1 ``N2sat`` uses the literal ``0.000986923 mb -> atm`` conversion
  factor; v3 uses ``1.0 / 1013.25``. Both agree to ~7 significant
  figures, comfortably within ``rtol=1e-6``.

The integrator branch in ``N2.run`` (Forward Euler, clip-with-log,
state set_at_time, derived TDG) is exercised in
``tests/v3/nsm1/test_n2_tier1.py``. This module covers the kinetic
sub-rate forms.

Synthetic mesh: 5-cell numpy arrays.
"""
from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules.nsm1 import processes as v1
from clearwater_modules_v3.processes.n2 import (
    N2,
    _kelvin,
    khn2_tc,
    n2sat_henry,
    pwv,
)


class _StubRegistry:
    """Minimal stand-in for VariableRegistry used by N2.run."""

    def __init__(self) -> None:
        self._data: dict[str, xr.DataArray] = {}

    def register(self, name: str, value: xr.DataArray) -> None:
        self._data[name] = value

    def get_at_time(self, name: str, time: datetime) -> xr.DataArray:
        if name not in self._data:
            raise KeyError(name)
        return self._data[name]

    def set_at_time(self, name: str, time: datetime, value: xr.DataArray) -> None:
        self._data[name] = value

    def __contains__(self, name: str) -> bool:
        return name in self._data


@pytest.fixture(scope="module")
def water_temp_c_5cell():
    return xr.DataArray(np.array([15.0, 18.0, 20.0, 22.0, 25.0]))


@pytest.fixture(scope="module")
def water_temp_k_5cell(water_temp_c_5cell):
    """Kelvin temperatures using the v1 (+273.16) convention so the
    stateless-helper parity tests can pass matched K values to both
    v1 and v3 without the 0.01 K offset between v3 (+273.15) and v1."""
    return water_temp_c_5cell + 273.16


@pytest.fixture(scope="module")
def n2_5cell():
    return xr.DataArray(np.array([10.0, 10.5, 11.0, 11.5, 12.0]))


@pytest.fixture(scope="module")
def depth_5cell():
    return xr.DataArray(np.array([1.0, 1.5, 2.0, 2.5, 3.0]))


@pytest.fixture(scope="function")
def time_zero() -> datetime:
    return datetime(2026, 1, 1)


def test_khn2_tc_matches_v1_KHN2_tc(water_temp_k_5cell):
    """v3 ``khn2_tc`` matches v1 ``KHN2_tc`` at matched Kelvin values."""
    v3_value = khn2_tc(water_temp_k_5cell)
    v1_value = v1.KHN2_tc(water_temp_k_5cell)

    np.testing.assert_allclose(
        np.asarray(v3_value), np.asarray(v1_value), rtol=1e-12
    )


def test_pwv_matches_v1_pwv(water_temp_k_5cell):
    """v3 ``pwv`` (water-vapor partial pressure) matches v1 ``pwv``."""
    v3_value = pwv(water_temp_k_5cell)
    v1_value = v1.pwv(water_temp_k_5cell)

    np.testing.assert_allclose(
        np.asarray(v3_value), np.asarray(v1_value), rtol=1e-12
    )


def test_n2sat_henry_matches_v1_N2sat(water_temp_k_5cell):
    """v3 ``n2sat_henry`` matches v1 ``N2sat`` at matched inputs."""
    pressure_mb = 1013.25
    khn2 = v1.KHN2_tc(water_temp_k_5cell)
    pwv_atm = v1.pwv(water_temp_k_5cell)

    v3_value = n2sat_henry(khn2, pressure_mb, pwv_atm)
    v1_value = v1.N2sat(KHN2_tc=khn2, pressure_mb=pressure_mb, pwv=pwv_atm)

    np.testing.assert_allclose(
        np.asarray(v3_value), np.asarray(v1_value), rtol=1e-6
    )


def test_n2_atmospheric_exchange_matches_v1_dN2dt(
    water_temp_c_5cell, n2_5cell, depth_5cell, time_zero
):
    """v3 cached ``n2_atm_exchange_rate`` matches v1 ``dN2dt = 1.034 *
    ka_tc * (N2sat - N2)`` for a fixed user-supplied ``ka_tc``.

    Setup: pin the reaeration menu to user-supplied (=1) with a
    non-zero ``kah_20_user`` so ``ka_tc`` is deterministic and
    matches both sides. Wind reaeration is zeroed. ``kah_theta=1.0``
    and ``kaw_theta=1.0`` make ``ka_tc == kah_20`` independent of
    temperature for the parity assertion.
    """
    ka_user = 1.5  # 1/d
    n2 = N2(
        parameters={
            "kah_20_user": ka_user,
            "kaw_20_user": 0.0,
            "kah_theta": 1.0,   # disable the temperature correction so
            "kaw_theta": 1.0,   # ka_tc reduces to the user-supplied value
            "hydraulic_reaeration_option": 1,
            "wind_reaeration_option": 1,
            "pressure_mb": 1013.25,
        },
        time_step=timedelta(minutes=5),
    )

    registry = _StubRegistry()
    registry.register("n2", n2_5cell.copy())
    registry.register("water_temperature", water_temp_c_5cell)
    registry.register("depth", depth_5cell)

    # Run a single substep so the cached attributes are populated.
    n2.run(time_zero, registry)

    # v1 reference. v3 ``N2.run`` converts deg C -> K with +273.15
    # whereas v1 ``celsius_to_kelvin`` uses +273.16. Pass +273.15 here
    # so the comparison isolates the kinetic formulas from that
    # convention difference. (The Henry's-law / Van't Hoff formulas
    # are identical between v1 and v3.)
    twater_k = water_temp_c_5cell + 273.15
    khn2 = v1.KHN2_tc(twater_k)
    pwv_atm = v1.pwv(twater_k)
    n2sat_v1 = v1.N2sat(KHN2_tc=khn2, pressure_mb=1013.25, pwv=pwv_atm)
    # ka_tc with kah_theta=1, kaw_theta=1 reduces to ka_user.
    v1_dn2dt = v1.dN2dt(ka_tc=ka_user, N2sat=n2sat_v1, N2=n2_5cell)

    # rtol=1e-5 accommodates the residual ~3e-7 relative difference
    # between v1's literal ``0.000986923`` mb->atm conversion and v3's
    # ``1.0/1013.25``, amplified slightly by the (N2sat - N2)
    # subtraction when N2sat ~ N2 (saturation ~14 vs N2 ~10-12).
    np.testing.assert_allclose(
        np.asarray(n2.n2_atm_exchange_rate),
        np.asarray(v1_dn2dt),
        rtol=1e-5,
    )


def test_n2_sat_cached_matches_v1_N2sat_at_run_time(
    water_temp_c_5cell, n2_5cell, depth_5cell, time_zero
):
    """v3 cached ``n2_sat`` matches v1 ``N2sat`` at matched inputs.

    Cross-checks the saturation cache populated by ``N2.run`` against
    the v1 reference, accepting the 0.01 K offset between v3
    (273.15) and v1 (273.16) as within the rtol=1e-6 tolerance for
    typical surface temperatures.
    """
    n2 = N2(
        parameters={
            "kah_20_user": 0.0,
            "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1,
            "wind_reaeration_option": 1,
            "pressure_mb": 1013.25,
        },
        time_step=timedelta(minutes=5),
    )

    registry = _StubRegistry()
    registry.register("n2", n2_5cell.copy())
    registry.register("water_temperature", water_temp_c_5cell)
    registry.register("depth", depth_5cell)

    n2.run(time_zero, registry)

    # Pass v1 the same +273.15 K convention v3 uses internally so the
    # comparison isolates the saturation formula (see notes on the
    # 0.01 K offset between v1 ``celsius_to_kelvin`` and v3 ``_kelvin``).
    twater_k = water_temp_c_5cell + 273.15
    khn2 = v1.KHN2_tc(twater_k)
    pwv_atm = v1.pwv(twater_k)
    n2sat_v1 = v1.N2sat(KHN2_tc=khn2, pressure_mb=1013.25, pwv=pwv_atm)

    np.testing.assert_allclose(
        np.asarray(n2.n2_sat), np.asarray(n2sat_v1), rtol=1e-6
    )
