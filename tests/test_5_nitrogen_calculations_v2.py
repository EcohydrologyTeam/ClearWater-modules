"""Parity tests: v2 Nitrogen sub-rate methods vs v1 nsm1.processes helpers.

Each test constructs a v2 Nitrogen instance, calls one of its rate
sub-methods directly, and compares to the equivalent v1 helper-function
output computed with the same inputs.

This file is the most important of the three v1-parity tests because the
Nitrogen process had 3 of the 4 LimnoTech bugs (missing set_at_time for
NH4 and NO3, time_step_frequency typo, broken NaN check). These tests
lock in the post-fix sub-rate behavior.

Scope: only the v1-equivalent sub-rate computations are exercised
(ammonium_nitrification, ammonium_from_bed, ammonium_decay_nitrate,
nitrate_denitrification, nitrate_bed_denitrification, nitrification_inhibition).
The full ``run()`` time-stepping is intentionally NOT exercised here because
the update equation has unrelated open issues that need design input.

Synthetic mesh: 5-cell numpy arrays, single time step.
"""
from datetime import timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules.nsm1 import processes as v1
from clearwater_modules_v2.processes.nitrogen import Nitrogen


@pytest.fixture(scope="function")
def nh4_5cell():
    # Function-scoped: some v2 sub-rate methods alias inputs (e.g., FloatingAlgae.limit_nitrogen
    # uses an in-place += that mutates its input). Function scope keeps tests independent.
    return xr.DataArray(np.array([0.05, 0.10, 0.15, 0.20, 0.30]))


@pytest.fixture(scope="function")
def no3_5cell():
    return xr.DataArray(np.array([1.0, 2.0, 3.0, 4.0, 5.0]))


@pytest.fixture(scope="function")
def water_temp_5cell():
    return xr.DataArray(np.array([15.0, 18.0, 20.0, 22.0, 25.0]))


@pytest.fixture(scope="function")
def depth_5cell():
    return xr.DataArray(np.array([0.5, 1.0, 1.5, 2.0, 3.0]))


@pytest.fixture(scope="function")
def dox_5cell():
    return xr.DataArray(np.array([4.0, 6.0, 7.0, 8.0, 10.0]))


@pytest.fixture(scope="function")
def n_instance():
    """Nitrogen instance with v1-aligned defaults."""
    inst = Nitrogen(
        time_step=timedelta(minutes=5),
        denitrification_rate=0.002,
        denitrification_theta=1.045,
        nitrification_rate=0.1,
        nitrification_theta=1.083,
        sediment_denitrification_rate=0.1,
        sediment_denitrification_theta=1.045,
        sediment_ammonium_release_rate=0.1,
        sediment_ammonium_release_theta=1.074,
        ammonium_decay_rate=0.0,
        ammonium_decay_theta=1.0,
        nitrification_oxygen_inhibition_factor=0.6,
    )
    inst.use_nitrate = True
    inst.use_ammonium = True
    inst.use_floating_algae = False
    inst.use_benthic_algae = False
    return inst


def test_nitrification_inhibition_matches_v1_NitrificationInhibition(n_instance, dox_5cell):
    """v2 nitrification_inhibition(DOX) == v1 NitrificationInhibition(use_DOX, KNR, DOX)."""
    v2_val = n_instance.nitrification_inhibition(dox_5cell)
    v1_val = v1.NitrificationInhibition(use_DOX=True, KNR=0.6, DOX=dox_5cell)

    np.testing.assert_allclose(v2_val.values, np.asarray(v1_val), rtol=1e-6)


def test_ammonium_nitrification_matches_v1_NH4_Nitrification(
    n_instance, nh4_5cell, water_temp_5cell, dox_5cell
):
    """v2 ammonium_nitrification == v1 NH4_Nitrification(NitrificationInhibition, NH4, knit_tc, use_NH4)."""
    v2_rate = n_instance.ammonium_nitrification(nh4_5cell, water_temp_5cell, dox_5cell)

    knit_tc_val = v1.knit_tc(water_temp_5cell, 0.1, 1.083)
    inhib = v1.NitrificationInhibition(use_DOX=True, KNR=0.6, DOX=dox_5cell)
    v1_rate = v1.NH4_Nitrification(
        NitrificationInhibition=inhib,
        NH4=nh4_5cell,
        knit_tc=knit_tc_val,
        use_NH4=True,
    )

    np.testing.assert_allclose(v2_rate.values, np.asarray(v1_rate), rtol=1e-6)


def test_ammonium_from_bed_matches_v1_NH4fromBed(n_instance, depth_5cell, water_temp_5cell):
    """v2 ammonium_from_bed == v1 NH4fromBed(depth, rnh4_tc)."""
    v2_rate = n_instance.ammonium_from_bed(depth=depth_5cell, temperature=water_temp_5cell)

    rnh4_tc_val = v1.rnh4_tc(water_temp_5cell, 0.1, 1.074)
    v1_rate = v1.NH4fromBed(depth=depth_5cell, rnh4_tc=rnh4_tc_val)

    np.testing.assert_allclose(v2_rate.values, v1_rate.values, rtol=1e-6)


def test_nitrate_bed_denitrification_matches_v1_NO3_BedDenit(
    n_instance, depth_5cell, no3_5cell, water_temp_5cell
):
    """v2 nitrate_bed_denitrification == v1 NO3_BedDenit(depth, vno3_tc, NO3)."""
    v2_rate = n_instance.nitrate_bed_denitrification(
        depth=depth_5cell, nitrate=no3_5cell, temperature=water_temp_5cell
    )

    vno3_tc_val = v1.vno3_tc(water_temp_5cell, 0.1, 1.045)
    v1_rate = v1.NO3_BedDenit(depth=depth_5cell, vno3_tc=vno3_tc_val, NO3=no3_5cell)

    np.testing.assert_allclose(v2_rate.values, v1_rate.values, rtol=1e-6)


def test_nitrate_denitrification_water_column_matches_v1_NO3_Denit(
    n_instance, no3_5cell, water_temp_5cell, dox_5cell
):
    """v2 nitrate_denitrification == v1 NO3_Denit(use_DOX, DOX, KsOxdn, kdnit_tc, NO3)
    (when DOX values are well-defined and yield no NaN denominator).
    """
    half_sat_o2 = 1.0
    v2_rate = n_instance.nitrate_denitrification(
        dissolved_oxygen=dox_5cell,
        half_saturation_oxygen=half_sat_o2,
        nitrate=no3_5cell,
        temperature=water_temp_5cell,
    )

    kdnit_tc_val = v1.kdnit_tc(water_temp_5cell, 0.002, 1.045)
    v1_rate = v1.NO3_Denit(
        use_DOX=True,
        DOX=dox_5cell,
        KsOxdn=half_sat_o2,
        kdnit_tc=kdnit_tc_val,
        NO3=no3_5cell,
    )

    np.testing.assert_allclose(np.asarray(v2_rate), np.asarray(v1_rate), rtol=1e-6)


def test_change_ammonium_no_algae_drops_to_decay_minus_nitrification_plus_bed(
    n_instance, nh4_5cell, no3_5cell, water_temp_5cell, depth_5cell, dox_5cell
):
    """When use_floating_algae=use_benthic_algae=False, change_ammonium reduces to:
    ammonium_decay_nitrate - ammonium_nitrification + ammonium_from_bed.
    Verifies the post-fix isnull NaN replacement (bug #4) does not perturb finite inputs.
    """
    n_instance.use_floating_algae = False
    n_instance.use_benthic_algae = False

    v2_rate = n_instance.change_ammonium(
        nitrate=no3_5cell,
        ammonium=nh4_5cell,
        temperature=water_temp_5cell,
        depth=depth_5cell,
        oxygen_dissolved=dox_5cell,
    )

    expected = (
        n_instance.ammonium_decay_nitrate(nh4_5cell, water_temp_5cell)
        - n_instance.ammonium_nitrification(nh4_5cell, water_temp_5cell, dox_5cell)
        + n_instance.ammonium_from_bed(depth=depth_5cell, temperature=water_temp_5cell)
    )

    np.testing.assert_allclose(v2_rate.values, expected.values, rtol=1e-6)
    # And no NaNs survive the isnull replacement
    assert not np.any(np.isnan(v2_rate.values))


def test_change_nitrate_no_algae_drops_to_nitrification_minus_denit_minus_beddenit(
    n_instance, nh4_5cell, no3_5cell, water_temp_5cell, depth_5cell, dox_5cell
):
    """When use_floating_algae=use_benthic_algae=False, change_nitrate reduces to:
    ammonium_nitrification - nitrate_denitrification(half_sat=1) - nitrate_bed_denitrification.
    Verifies the post-fix isnull NaN replacement and that the time_step_frequency
    typo (bug #3) is no longer reachable from change_nitrate.
    """
    n_instance.use_floating_algae = False
    n_instance.use_benthic_algae = False

    v2_rate = n_instance.change_nitrate(
        nitrate=no3_5cell,
        ammonium=nh4_5cell,
        temperature=water_temp_5cell,
        depth=depth_5cell,
        oxygen_dissolved=dox_5cell,
    )

    expected = (
        n_instance.ammonium_nitrification(nh4_5cell, water_temp_5cell, dox_5cell)
        - n_instance.nitrate_denitrification(dox_5cell, 1, no3_5cell, water_temp_5cell)
        - n_instance.nitrate_bed_denitrification(depth_5cell, no3_5cell, water_temp_5cell)
    )

    np.testing.assert_allclose(np.asarray(v2_rate), np.asarray(expected), rtol=1e-6)
    assert not np.any(np.isnan(np.asarray(v2_rate)))
