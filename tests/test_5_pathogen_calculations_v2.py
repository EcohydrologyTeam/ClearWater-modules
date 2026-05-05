"""Parity tests: v3 Pathogen sub-rate methods vs v1 nsm1.processes helpers.

Each test constructs a v3 ``Pathogen`` instance, calls one of its rate
sub-methods directly, and compares to the equivalent v1 helper-function
output computed with the same inputs.

Scope: the v1-equivalent sub-rate computations on Pathogen
(``_rate_natural_decay``, ``_rate_light_decay``, ``_rate_settling``)
plus the lumped ``rate`` aggregator. The integrator branch in
``Pathogen.run`` (Forward Euler + clip-with-log + registry write) is
exercised in ``tests/v3/nsm1/test_pathogen_tier1.py``.

v1 reference: ``clearwater_modules.nsm1.processes`` ``kdx_tc``,
``PathogenDeath``, ``PathogenDecay``, ``PathogenSettling``, ``dPXdt``.

v3 deviation note (``_rate_light_decay``): v1 ``PathogenDecay`` uses
raw ``q_solar`` (W/m^2 incident). v3 ``_rate_light_decay`` uses
``PAR(q_solar, Fr_PAR) = q_solar * Fr_PAR`` so the effective rate is
scaled by ``Fr_PAR`` relative to v1. The Phase 3.1 docstring documents
this as an intentional deviation absorbed into the calibration target
``apx``. The light-decay parity test below pins ``Fr_PAR=1.0`` to make
the v3 formula exactly equivalent to v1 at the kinetics level. With
``Fr_PAR=1.0`` v3 should match v1 to roundoff.

Synthetic mesh: 5-cell numpy arrays, single time step.
"""
from datetime import timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules.nsm1 import processes as v1
from clearwater_modules_v3.processes.pathogen import Pathogen
from clearwater_modules_v3.utils.light import L


@pytest.fixture(scope="module")
def px_5cell():
    return xr.DataArray(np.array([1.0e3, 5.0e3, 1.0e4, 5.0e4, 1.0e5]))


@pytest.fixture(scope="module")
def water_temp_5cell():
    return xr.DataArray(np.array([15.0, 18.0, 20.0, 22.0, 25.0]))


@pytest.fixture(scope="module")
def depth_5cell():
    return xr.DataArray(np.array([0.5, 1.0, 1.5, 2.0, 3.0]))


@pytest.fixture(scope="module")
def q_solar_5cell():
    return xr.DataArray(np.array([200.0, 250.0, 300.0, 350.0, 400.0]))


@pytest.fixture(scope="module")
def solid_5cell():
    return xr.DataArray(np.array([10.0, 12.0, 15.0, 18.0, 20.0]))


@pytest.fixture(scope="module")
def poc_5cell():
    return xr.DataArray(np.array([1.0, 1.2, 1.4, 1.6, 1.8]))


@pytest.fixture(scope="module")
def ap_5cell():
    return xr.DataArray(np.array([5.0, 6.0, 7.0, 8.0, 10.0]))


@pytest.fixture(scope="function")
def pathogen_instance():
    """Pathogen instance with v1-aligned defaults and Fr_PAR=1.0.

    Fr_PAR=1.0 makes ``_rate_light_decay`` exactly equivalent to v1's
    ``PathogenDecay`` (which uses raw ``q_solar`` rather than PAR).
    """
    return Pathogen(
        parameters={
            "kdx_20": 0.8,
            "kdx_theta": 1.07,
            "apx": 1.0,
            "vx": 1.0,
            "Fr_PAR": 1.0,  # remove the PAR scaling so v3 form == v1 form
        },
        time_step=timedelta(minutes=5),
    )


def test_pathogen_natural_decay_matches_v1_PathogenDeath(
    pathogen_instance, px_5cell, water_temp_5cell
):
    """v3 ``_rate_natural_decay`` == v1 ``PathogenDeath(kdx_tc, PX)``."""
    v3_rate = pathogen_instance._rate_natural_decay(px_5cell, water_temp_5cell)

    kdx_tc = v1.kdx_tc(water_temp_5cell, 0.8, 1.07)
    v1_rate = v1.PathogenDeath(kdx_tc, px_5cell)

    np.testing.assert_allclose(np.asarray(v3_rate), np.asarray(v1_rate), rtol=1e-6)


def test_pathogen_light_decay_matches_v1_PathogenDecay(
    pathogen_instance,
    px_5cell,
    depth_5cell,
    q_solar_5cell,
    solid_5cell,
    poc_5cell,
    ap_5cell,
):
    """v3 ``_rate_light_decay`` == v1 ``PathogenDecay(apx, q_solar, L, depth, PX)``.

    With ``Fr_PAR=1.0`` (set in the fixture) the v3 effective surface
    irradiance equals ``q_solar``, matching v1 exactly at the kinetics
    level.
    """
    v3_rate = pathogen_instance._rate_light_decay(
        px=px_5cell,
        depth=depth_5cell,
        q_solar=q_solar_5cell,
        solid=solid_5cell,
        poc=poc_5cell,
        ap=ap_5cell,
    )

    # Compute v1 KEXT via the v3 utility (same Beer-Lambert formula as the
    # v1 shared L). Inputs match the Pathogen instance attributes.
    kext = L(
        lambda0=pathogen_instance.lambda0,
        lambda1=pathogen_instance.lambda1,
        lambda2=pathogen_instance.lambda2,
        lambdas=pathogen_instance.lambdas,
        lambdam=pathogen_instance.lambdam,
        Solid=solid_5cell,
        POC=poc_5cell,
        fcom=pathogen_instance.fcom,
        Ap=ap_5cell,
        use_Algae=pathogen_instance.use_Algae,
        use_POC=pathogen_instance.use_POC,
    )
    v1_rate = v1.PathogenDecay(
        apx=1.0,
        q_solar=q_solar_5cell,
        L=kext,
        depth=depth_5cell,
        PX=px_5cell,
    )

    np.testing.assert_allclose(np.asarray(v3_rate), np.asarray(v1_rate), rtol=1e-6)


def test_pathogen_settling_matches_v1_PathogenSettling(
    pathogen_instance, px_5cell, depth_5cell
):
    """v3 ``_rate_settling`` == v1 ``PathogenSettling(vx, depth, PX)``."""
    v3_rate = pathogen_instance._rate_settling(px_5cell, depth_5cell)
    v1_rate = v1.PathogenSettling(vx=1.0, depth=depth_5cell, PX=px_5cell)

    np.testing.assert_allclose(np.asarray(v3_rate), np.asarray(v1_rate), rtol=1e-6)


def test_pathogen_total_rate_matches_v1_dPXdt(
    pathogen_instance,
    px_5cell,
    water_temp_5cell,
    depth_5cell,
    q_solar_5cell,
    solid_5cell,
    poc_5cell,
    ap_5cell,
):
    """v3 ``rate`` == v1 ``dPXdt(PathogenDeath, PathogenDecay, PathogenSettling)``.

    v3 ``rate`` returns the *signed* rate of change ``dPX/dt`` (the
    sum of the three negative loss terms). v1 ``dPXdt`` is the same
    quantity computed by negating the three positive loss helpers.
    """
    v3_rate = pathogen_instance.rate(
        px=px_5cell,
        water_temperature=water_temp_5cell,
        depth=depth_5cell,
        q_solar=q_solar_5cell,
        solid=solid_5cell,
        poc=poc_5cell,
        ap=ap_5cell,
    )

    kdx_tc = v1.kdx_tc(water_temp_5cell, 0.8, 1.07)
    death = v1.PathogenDeath(kdx_tc, px_5cell)
    kext = L(
        lambda0=pathogen_instance.lambda0,
        lambda1=pathogen_instance.lambda1,
        lambda2=pathogen_instance.lambda2,
        lambdas=pathogen_instance.lambdas,
        lambdam=pathogen_instance.lambdam,
        Solid=solid_5cell,
        POC=poc_5cell,
        fcom=pathogen_instance.fcom,
        Ap=ap_5cell,
        use_Algae=pathogen_instance.use_Algae,
        use_POC=pathogen_instance.use_POC,
    )
    decay = v1.PathogenDecay(
        apx=1.0,
        q_solar=q_solar_5cell,
        L=kext,
        depth=depth_5cell,
        PX=px_5cell,
    )
    settling = v1.PathogenSettling(vx=1.0, depth=depth_5cell, PX=px_5cell)
    v1_rate = v1.dPXdt(
        PathogenDeath=death,
        PathogenDecay=decay,
        PathogenSettling=settling,
    )

    np.testing.assert_allclose(np.asarray(v3_rate), np.asarray(v1_rate), rtol=1e-6)
