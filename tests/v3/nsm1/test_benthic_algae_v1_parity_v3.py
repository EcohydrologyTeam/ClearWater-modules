"""v3 BenthicAlgae kinetic regression against frozen v1 reference values.

Migration from ``tests/test_5_benthic_algae_calculations_v2.py``.
"""
from datetime import timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.benthic_algae import BenthicAlgae


V1_AB_RESPIRATION_REFERENCE = np.array([
    0.7948159826845926,
    1.8244696038802821,
    3.0,
    4.384835999999999,
    6.290764288750032,
])

V1_AB_DEATH_REFERENCE = np.array([
    0.5961119870134445,
    1.3683522029102115,
    2.25,
    3.288626999999999,
    4.718073216562525,
])

V1_FSB_REFERENCE = np.array([
    0.6666666666666667,
    0.5,
    0.4,
    0.33333333333333337,
    0.2857142857142857,
])

V1_AB_GROWTH_REFERENCE = np.array([
    0.05961119870134445,
    0.2627236229587606,
    0.6299999999999999,
    1.017418978125,
    1.509783429300008,
])


@pytest.fixture(scope="module")
def benthic_5cell():
    return xr.DataArray(np.array([5.0, 10.0, 15.0, 20.0, 25.0]))


@pytest.fixture(scope="module")
def water_temp_5cell():
    return xr.DataArray(np.array([15.0, 18.0, 20.0, 22.0, 25.0]))


@pytest.fixture(scope="function")
def ba_instance():
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


def test_rate_respiration_matches_v1_AbRespiration(
    ba_instance, benthic_5cell, water_temp_5cell
):
    v3_rate = ba_instance.rate_respiration(benthic_5cell, water_temp_5cell)
    np.testing.assert_allclose(
        v3_rate.values, V1_AB_RESPIRATION_REFERENCE, rtol=1e-6
    )


def test_rate_death_matches_v1_AbDeath(
    ba_instance, benthic_5cell, water_temp_5cell
):
    v3_rate = ba_instance.rate_death(benthic_5cell, water_temp_5cell)
    np.testing.assert_allclose(
        v3_rate.values, V1_AB_DEATH_REFERENCE, rtol=1e-6
    )


def test_limit_density_matches_v1_FSb(ba_instance, benthic_5cell):
    v3_limit = ba_instance.limit_density(algae=benthic_5cell)
    np.testing.assert_allclose(v3_limit.values, V1_FSB_REFERENCE, rtol=1e-6)


def test_rate_growth_multiplicative_matches_v1_mub_times_Ab(
    ba_instance, benthic_5cell, water_temp_5cell
):
    limit_p = xr.DataArray(np.array([0.30, 0.40, 0.50, 0.55, 0.60]))
    limit_n = xr.DataArray(np.array([0.50, 0.60, 0.70, 0.75, 0.80]))
    limit_l = xr.DataArray(np.array([0.20, 0.30, 0.40, 0.45, 0.50]))
    limit_s = xr.DataArray(np.array([0.50, 0.40, 0.30, 0.25, 0.20]))

    v3_rate = ba_instance.rate_growth(
        benthic_5cell, water_temp_5cell, limit_p, limit_n, limit_l, limit_s
    )
    np.testing.assert_allclose(
        v3_rate.values, V1_AB_GROWTH_REFERENCE, rtol=1e-6
    )


def test_phase9e_bwa_harmonized_to_wasp7_canonical():
    """Phase 9.E follow-up: BWa = 1000 gives rab = BWa/BWd = 10
    mg-Chla/g-DW (WASP7 canonical benthic Chla:DW). Pin so any future
    change requires explicit reconciliation against the WASP7 reference."""
    from clearwater_modules_v3.parameters.balgae import (
        DEFAULTS as BALGAE_DEFAULTS,
    )
    from clearwater_modules_v3.parameters.algae import (
        DEFAULTS as ALGAE_DEFAULTS,
    )

    bwa = BALGAE_DEFAULTS["BWa"]
    bwd = BALGAE_DEFAULTS["BWd"]

    assert bwa == 1000.0
    assert bwd == 100.0

    # WASP7 canonical Chla:DW = 10 mg/g.
    np.testing.assert_allclose(bwa / bwd, 10.0, rtol=1e-12)

    # v3 benthic and floating algae share Chla:DW basis.
    awa = ALGAE_DEFAULTS["AWa"]
    awd = ALGAE_DEFAULTS["AWd"]
    np.testing.assert_allclose(awa / awd, bwa / bwd, rtol=1e-12)

    # Confirm v3 deliberately differs from v1 (3500) and Fortran (5000).
    assert bwa < 3500.0
    assert bwa < 5000.0
