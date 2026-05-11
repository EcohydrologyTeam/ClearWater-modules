"""v3 FloatingAlgae kinetic regression against frozen v1 reference values.

Migration from ``tests/test_5_floating_algae_calculations_v2.py``.
"""
from datetime import timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.floating_algae import FloatingAlgae


V1_AP_RESPIRATION_REFERENCE = np.array([
    1.5896319653691853,
    3.6489392077605642,
    6.0,
    8.769671999999998,
    12.581528577500064,
])

V1_AP_DEATH_REFERENCE = np.array([
    1.192223974026889,
    2.736704405820423,
    4.5,
    6.577253999999998,
    9.43614643312505,
])

V1_AP_SETTLING_REFERENCE = np.array([
    3.0,
    3.0,
    2.9999999999999996,
    3.0,
    2.5,
])

V1_FP_REFERENCE = np.array([
    0.9541984732824428,
    0.9668508287292817,
    0.9765625,
    0.9803921568627451,
    0.9842519685039369,
])

V1_FN_REFERENCE = np.array([
    0.963302752293578,
    0.9813084112149533,
    0.9874608150470219,
    0.9905660377358491,
    0.9925093632958801,
])

V1_AP_GROWTH_REFERENCE = np.array([
    0.2384447948053778,
    1.313618114793803,
    4.199999999999999,
    8.139351825,
    15.097834293000078,
])


@pytest.fixture(scope="module")
def algae_5cell():
    return xr.DataArray(np.array([10.0, 20.0, 30.0, 40.0, 50.0]))


@pytest.fixture(scope="module")
def water_temp_5cell():
    return xr.DataArray(np.array([15.0, 18.0, 20.0, 22.0, 25.0]))


@pytest.fixture(scope="module")
def depth_5cell():
    return xr.DataArray(np.array([0.5, 1.0, 1.5, 2.0, 3.0]))


@pytest.fixture(scope="function")
def nh4_5cell():
    return xr.DataArray(np.array([0.05, 0.10, 0.15, 0.20, 0.30]))


@pytest.fixture(scope="function")
def no3_5cell():
    return xr.DataArray(np.array([1.0, 2.0, 3.0, 4.0, 5.0]))


@pytest.fixture(scope="function")
def tip_5cell():
    return xr.DataArray(np.array([0.05, 0.07, 0.10, 0.12, 0.15]))


@pytest.fixture(scope="function")
def fa_instance():
    inst = FloatingAlgae(
        time_step=timedelta(minutes=5),
        settling_velocity=0.15,
        repiration_rate=0.2,
        repiration_rate_correction_factor=1.047,
        death_rate=0.15,
        death_rate_correction_factor=1.047,
        growth_rate_option=1,
        growth_rate_max=1.0,
        growth_rate_correction=1.047,
        phosphorus_michaelis_menton_constant=0.0012,
        nitrogen_michaelis_menton_constant=0.04,
        light_limitation_option=1,
        light_limitation_constant=10.0,
        light_attenuation_coefficient=1.0,
    )
    inst.use_nitrate = True
    inst.use_ammonium = True
    inst.use_phosphate = True
    return inst


def test_rate_respiration_matches_v1_ApRespiration(
    fa_instance, algae_5cell, water_temp_5cell
):
    v3_rate = fa_instance.rate_respiration(algae_5cell, water_temp_5cell)
    np.testing.assert_allclose(
        v3_rate.values, V1_AP_RESPIRATION_REFERENCE, rtol=1e-6
    )


def test_rate_death_matches_v1_ApDeath(
    fa_instance, algae_5cell, water_temp_5cell
):
    v3_rate = fa_instance.rate_death(algae_5cell, water_temp_5cell)
    np.testing.assert_allclose(
        v3_rate.values, V1_AP_DEATH_REFERENCE, rtol=1e-6
    )


def test_rate_settling_matches_v1_ApSettling(
    fa_instance, algae_5cell, depth_5cell
):
    v3_rate = fa_instance.rate_settling(algae_5cell, depth_5cell)
    np.testing.assert_allclose(
        v3_rate.values, V1_AP_SETTLING_REFERENCE, rtol=1e-6
    )


def test_limit_phosphorus_matches_v1_FP(fa_instance, tip_5cell):
    fdp = 0.5
    v3_limit = fa_instance.limit_phosphorus(
        concentration=tip_5cell, fraction_dissolved=fdp
    )
    np.testing.assert_allclose(
        np.asarray(v3_limit), V1_FP_REFERENCE, rtol=1e-6
    )


def test_limit_nitrogen_matches_v1_FN(fa_instance, nh4_5cell, no3_5cell):
    """v3 ``limit_nitrogen`` matches frozen v1 ``FN`` reference.

    Pass deep copies into v3 to defend against the legacy in-place
    mutation of the nitrate input (the original parity test flagged
    this as a known v2-side bug; v3 may still have the same issue
    pending review).
    """
    v3_limit = fa_instance.limit_nitrogen(
        nitrate=no3_5cell.copy(deep=True), ammonium=nh4_5cell.copy(deep=True)
    )
    np.testing.assert_allclose(
        np.asarray(v3_limit), V1_FN_REFERENCE, rtol=1e-6
    )


def test_rate_growth_multiplicative_matches_v1_mu_times_Ap(
    fa_instance, algae_5cell, water_temp_5cell
):
    limit_p = xr.DataArray(np.array([0.30, 0.40, 0.50, 0.55, 0.60]))
    limit_n = xr.DataArray(np.array([0.50, 0.60, 0.70, 0.75, 0.80]))
    limit_l = xr.DataArray(np.array([0.20, 0.30, 0.40, 0.45, 0.50]))

    v3_rate = fa_instance.rate_growth(
        algae_5cell, water_temp_5cell, limit_p, limit_n, limit_l
    )
    np.testing.assert_allclose(
        v3_rate.values, V1_AP_GROWTH_REFERENCE, rtol=1e-6
    )
