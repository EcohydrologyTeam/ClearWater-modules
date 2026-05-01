"""Parity tests: v2 BenthicAlgae sub-rate methods vs v1 nsm1.processes helpers.

Each test constructs a v2 BenthicAlgae instance, calls one of its rate
sub-methods directly, and compares to the equivalent v1 helper-function
output computed with the same inputs.

Scope: only the v1-equivalent sub-rate computations on BenthicAlgae are
exercised (rate_growth multiplicative branch, rate_respiration, rate_death
inherited from FloatingAlgae, limit_density). The v2 ``run()`` method
itself is intentionally NOT exercised here because it has unrelated open
issues that need design input.

Synthetic mesh: 5-cell numpy arrays, single time step.
"""
from datetime import timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules.nsm1 import processes as v1
from clearwater_modules_v2.processes.benthic_algae import BenthicAlgae


@pytest.fixture(scope="module")
def benthic_5cell():
    return xr.DataArray(np.array([5.0, 10.0, 15.0, 20.0, 25.0]))


@pytest.fixture(scope="module")
def water_temp_5cell():
    return xr.DataArray(np.array([15.0, 18.0, 20.0, 22.0, 25.0]))


@pytest.fixture(scope="module")
def depth_5cell():
    return xr.DataArray(np.array([0.5, 1.0, 1.5, 2.0, 3.0]))


@pytest.fixture(scope="function")
def ba_instance():
    """BenthicAlgae instance with v1-aligned defaults."""
    inst = BenthicAlgae(
        time_step=timedelta(minutes=5),
        settling_velocity=0.0,
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
        density_michaelis_menton_constant=10.0,
    )
    inst.use_nitrate = True
    inst.use_ammonium = True
    inst.use_phosphate = True
    return inst


def test_rate_respiration_matches_v1_AbRespiration(ba_instance, benthic_5cell, water_temp_5cell):
    """v2 rate_respiration (inherited) == v1 AbRespiration(krb_tc, Ab)."""
    v2_rate = ba_instance.rate_respiration(benthic_5cell, water_temp_5cell)

    krb_tc = v1.krb_tc(0.2, water_temp_5cell, 1.047)
    v1_rate = v1.AbRespiration(krb_tc, benthic_5cell)

    np.testing.assert_allclose(v2_rate.values, v1_rate.values, rtol=1e-6)


def test_rate_death_matches_v1_AbDeath(ba_instance, benthic_5cell, water_temp_5cell):
    """v2 rate_death (inherited) == v1 AbDeath(kdb_tc, Ab)."""
    v2_rate = ba_instance.rate_death(benthic_5cell, water_temp_5cell)

    kdb_tc = v1.kdb_tc(0.15, water_temp_5cell, 1.047)
    v1_rate = v1.AbDeath(kdb_tc, benthic_5cell)

    np.testing.assert_allclose(v2_rate.values, v1_rate.values, rtol=1e-6)


def test_limit_density_matches_v1_FSb(ba_instance, benthic_5cell):
    """v2 limit_density(Ab) == v1 FSb(Ab, Ksb)."""
    v2_limit = ba_instance.limit_density(algae=benthic_5cell)
    v1_limit = v1.FSb(Ab=benthic_5cell, Ksb=10.0)

    np.testing.assert_allclose(v2_limit.values, np.asarray(v1_limit), rtol=1e-6)


def test_rate_growth_multiplicative_matches_v1_mub_times_Ab(
    ba_instance, benthic_5cell, water_temp_5cell
):
    """v2 rate_growth (option=1) == v1 mub_max_tc * FLb * FPb * FNb * FSb * Ab."""
    # Synthetic limit factors so the test does not depend on the v2 limit_light
    # implementation (which differs in parenthesization from v1 FLb).
    limit_p = xr.DataArray(np.array([0.30, 0.40, 0.50, 0.55, 0.60]))
    limit_n = xr.DataArray(np.array([0.50, 0.60, 0.70, 0.75, 0.80]))
    limit_l = xr.DataArray(np.array([0.20, 0.30, 0.40, 0.45, 0.50]))
    limit_s = xr.DataArray(np.array([0.50, 0.40, 0.30, 0.25, 0.20]))

    v2_rate = ba_instance.rate_growth(
        benthic_5cell, water_temp_5cell, limit_p, limit_n, limit_l, limit_s
    )

    mub_max_tc_val = v1.mub_max_tc(1.0, water_temp_5cell, 1.047)
    v1_mub = v1.mub(
        mub_max_tc=mub_max_tc_val,
        b_growth_rate_option=1,
        FLb=limit_l,
        FPb=limit_p,
        FNb=limit_n,
        FSb=limit_s,
    )
    v1_rate = v1.AbGrowth(v1_mub, benthic_5cell)

    np.testing.assert_allclose(v2_rate.values, np.asarray(v1_rate), rtol=1e-6)
