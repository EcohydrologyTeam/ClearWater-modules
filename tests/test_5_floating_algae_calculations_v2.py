"""Parity tests: v2 FloatingAlgae sub-rate methods vs v1 nsm1.processes helpers.

Each test constructs a v2 FloatingAlgae instance, calls one of its rate
sub-methods directly, and compares to the equivalent v1 helper-function
output computed with the same inputs.

Scope: only the v1-equivalent sub-rate computations on FloatingAlgae are
exercised (rate_growth multiplicative branch, rate_respiration, rate_death,
rate_settling, limit_phosphorus, limit_nitrogen). The v2 ``run()`` method
itself is intentionally NOT exercised here because it has unrelated open
issues that need design input (broken update equation, broken NaN check
patterns elsewhere). Those tests will land once the design is settled.

Synthetic mesh: 5-cell numpy arrays, single time step.
"""
from datetime import timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules.nsm1 import processes as v1
from clearwater_modules_v3.processes.floating_algae import FloatingAlgae


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
    # Function-scoped to defend against in-place mutation in v2 limit_nitrogen
    # (v2 currently does n_concentration = nitrate; n_concentration += ammonium,
    # which aliases the nitrate fixture if scope="module").
    return xr.DataArray(np.array([0.05, 0.10, 0.15, 0.20, 0.30]))


@pytest.fixture(scope="function")
def no3_5cell():
    return xr.DataArray(np.array([1.0, 2.0, 3.0, 4.0, 5.0]))


@pytest.fixture(scope="function")
def tip_5cell():
    return xr.DataArray(np.array([0.05, 0.07, 0.10, 0.12, 0.15]))


@pytest.fixture(scope="function")
def fa_instance():
    """FloatingAlgae instance with v1-aligned defaults."""
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
    # init_process is normally invoked by Model; mimic the use-flag setup directly
    inst.use_nitrate = True
    inst.use_ammonium = True
    inst.use_phosphate = True
    return inst


def test_rate_respiration_matches_v1_ApRespiration(fa_instance, algae_5cell, water_temp_5cell):
    """v2 rate_respiration == v1 ApRespiration(krp_tc, Ap)."""
    v2_rate = fa_instance.rate_respiration(algae_5cell, water_temp_5cell)

    krp_tc = v1.krp_tc(water_temp_5cell, 0.2, 1.047)
    v1_rate = v1.ApRespiration(krp_tc, algae_5cell)

    np.testing.assert_allclose(v2_rate.values, v1_rate.values, rtol=1e-6)


def test_rate_death_matches_v1_ApDeath(fa_instance, algae_5cell, water_temp_5cell):
    """v2 rate_death == v1 ApDeath(kdp_tc, Ap)."""
    v2_rate = fa_instance.rate_death(algae_5cell, water_temp_5cell)

    kdp_tc = v1.kdp_tc(water_temp_5cell, 0.15, 1.047)
    v1_rate = v1.ApDeath(kdp_tc, algae_5cell)

    np.testing.assert_allclose(v2_rate.values, v1_rate.values, rtol=1e-6)


def test_rate_settling_matches_v1_ApSettling(fa_instance, algae_5cell, depth_5cell):
    """v2 rate_settling == v1 ApSettling(vsap, Ap, depth)."""
    v2_rate = fa_instance.rate_settling(algae_5cell, depth_5cell)
    v1_rate = v1.ApSettling(0.15, algae_5cell, depth_5cell)

    np.testing.assert_allclose(v2_rate.values, v1_rate.values, rtol=1e-6)


def test_limit_phosphorus_matches_v1_FP(fa_instance, tip_5cell):
    """v2 limit_phosphorus(TIP, fdp) == v1 FP(fdp, TIP, use_TIP=True, KsP)."""
    fdp = 0.5  # matches the hardcoded TODO in v2 floating_algae.run
    v2_limit = fa_instance.limit_phosphorus(concentration=tip_5cell, fraction_dissolved=fdp)
    v1_limit = v1.FP(fdp=fdp, TIP=tip_5cell, use_TIP=True, KsP=0.0012)

    np.testing.assert_allclose(np.asarray(v2_limit), np.asarray(v1_limit), rtol=1e-6)


def test_limit_nitrogen_matches_v1_FN(fa_instance, nh4_5cell, no3_5cell):
    """v2 limit_nitrogen(NO3, NH4) == v1 FN(use_NH4, use_NO3, NH4, NO3, KsN).

    NOTE: v2 limit_nitrogen aliases its nitrate input via ``n_concentration = nitrate``
    and then does ``n_concentration += ammonium``, which mutates the caller's
    nitrate DataArray in place. We pass deep copies into the v2 call so that
    the v1 reference call below sees the original values. This is a known
    additional bug; do not remove the .copy() until that bug is fixed.
    """
    v2_limit = fa_instance.limit_nitrogen(
        nitrate=no3_5cell.copy(deep=True), ammonium=nh4_5cell.copy(deep=True)
    )
    v1_limit = v1.FN(use_NH4=True, use_NO3=True, NH4=nh4_5cell, NO3=no3_5cell, KsN=0.04)

    np.testing.assert_allclose(np.asarray(v2_limit), np.asarray(v1_limit), rtol=1e-6)


def test_rate_growth_multiplicative_matches_v1_mu_times_Ap(
    fa_instance, algae_5cell, water_temp_5cell
):
    """v2 rate_growth (option=1, multiplicative) == v1 mu_max_tc * FL * FP * FN * Ap."""
    # Use synthetic limit factors so the test does not depend on limit_light
    # (v2 limit_light option 1 has a different parenthesization than v1 FL).
    limit_p = xr.DataArray(np.array([0.30, 0.40, 0.50, 0.55, 0.60]))
    limit_n = xr.DataArray(np.array([0.50, 0.60, 0.70, 0.75, 0.80]))
    limit_l = xr.DataArray(np.array([0.20, 0.30, 0.40, 0.45, 0.50]))

    v2_rate = fa_instance.rate_growth(
        algae_5cell, water_temp_5cell, limit_p, limit_n, limit_l
    )

    mu_max_tc_val = v1.mu_max_tc(water_temp_5cell, 1.0, 1.047)
    v1_mu = v1.mu(
        mu_max_tc=mu_max_tc_val,
        growth_rate_option=1,
        FL=limit_l,
        FP=limit_p,
        FN=limit_n,
    )
    v1_rate = v1.ApGrowth(v1_mu, algae_5cell)

    np.testing.assert_allclose(v2_rate.values, np.asarray(v1_rate), rtol=1e-6)
