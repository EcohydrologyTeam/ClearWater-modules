"""v3 N2 kinetic regression against frozen v1 reference values.

Migration from ``tests/test_5_n2_calculations_v2.py``. v1 references
captured 2026-05-10 and frozen as numpy literals so v1 source can be
retired.

Notes on conventions:
- v1 ``celsius_to_kelvin`` uses ``+273.16``; v3 uses ``+273.15``. The
  stateless-helper tests use Kelvin values from the v1 convention so
  the kinetic formulas match exactly; the run-cached tests use the v3
  convention (+273.15) and the v1 reference was captured accordingly.
- v1 ``N2sat`` uses literal ``0.000986923 mb -> atm`` factor; v3 uses
  ``1/1013.25``. Both agree to ~7 significant figures (rtol=1e-6).
"""
from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.n2 import (
    N2,
    khn2_tc,
    n2sat_henry,
    pwv,
)


V1_KHN2_TC_REFERENCE = np.array([
    0.000756069696483947,
    0.000721729216680316,
    0.0007000764404082914,
    0.0006793536535552897,
    0.0006499049524863102,
])

V1_PWV_REFERENCE = np.array([
    0.016838151070052498,
    0.020379788328724724,
    0.02308854207367629,
    0.026107032428089552,
    0.03128045141547398,
])

V1_N2SAT_HENRY_REFERENCE = np.array([
    16.44265152168544,
    15.639289765216636,
    15.128144651820115,
    14.634980473573602,
    13.926209099075631,
])

V1_DN2DT_REFERENCE = np.array([
    9.996826788783263,
    7.975075842296073,
    6.406645240289634,
    4.866113106528651,
    2.991124664293515,
])

V1_N2_SAT_AT_273_15_REFERENCE = np.array([
    16.445407342864772,
    15.641892870597081,
    15.130654571431098,
    14.637403679257673,
    13.928513645579313,
])


class _StubRegistry:
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
    v3_value = khn2_tc(water_temp_k_5cell)
    np.testing.assert_allclose(
        np.asarray(v3_value), V1_KHN2_TC_REFERENCE, rtol=1e-12
    )


def test_pwv_matches_v1_pwv(water_temp_k_5cell):
    v3_value = pwv(water_temp_k_5cell)
    np.testing.assert_allclose(
        np.asarray(v3_value), V1_PWV_REFERENCE, rtol=1e-12
    )


def test_n2sat_henry_matches_v1_N2sat(water_temp_k_5cell):
    pressure_mb = 1013.25
    khn2 = khn2_tc(water_temp_k_5cell)
    pwv_atm = pwv(water_temp_k_5cell)
    v3_value = n2sat_henry(khn2, pressure_mb, pwv_atm)
    np.testing.assert_allclose(
        np.asarray(v3_value), V1_N2SAT_HENRY_REFERENCE, rtol=1e-6
    )


def test_n2_atmospheric_exchange_matches_v1_dN2dt(
    water_temp_c_5cell, n2_5cell, depth_5cell, time_zero
):
    """v3 cached ``n2_atm_exchange_rate`` matches frozen v1 ``dN2dt``
    reference (rtol=1e-5 to accommodate the literal mb->atm conversion
    factor difference between v1 and v3, amplified by ``N2sat - N2``)."""
    ka_user = 1.5
    n2 = N2(
        parameters={
            "kah_20_user": ka_user,
            "kaw_20_user": 0.0,
            "kah_theta": 1.0,
            "kaw_theta": 1.0,
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

    np.testing.assert_allclose(
        np.asarray(n2.n2_atm_exchange_rate),
        V1_DN2DT_REFERENCE,
        rtol=1e-5,
    )


def test_n2_sat_cached_matches_v1_N2sat_at_run_time(
    water_temp_c_5cell, n2_5cell, depth_5cell, time_zero
):
    """v3 cached ``n2_sat`` matches frozen v1 ``N2sat`` reference
    (computed via v1 at +273.15 to match v3's internal convention)."""
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

    np.testing.assert_allclose(
        np.asarray(n2.n2_sat),
        V1_N2_SAT_AT_273_15_REFERENCE,
        rtol=1e-6,
    )
