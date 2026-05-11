"""v3 Nitrogen kinetic regression against frozen v1 reference values.

Migration from ``tests/test_5_nitrogen_calculations_v2.py``.
"""
from datetime import timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.nitrogen import Nitrogen


V1_NITR_INHIB_REFERENCE = np.array([
    0.9092820467105875,
    0.9726762775527075,
    0.9850044231795223,
    0.99177025295098,
    0.9975212478233336,
])

V1_NH4_NITR_REFERENCE = np.array([
    0.003051591085018409,
    0.008292995138949277,
    0.014775066347692836,
    0.023264728404268437,
    0.04458468401833538,
])

V1_NH4_FROM_BED_REFERENCE = np.array([
    0.13996150018542117,
    0.08669447825529096,
    0.06666666666666667,
    0.05767380000000001,
    0.04763214639622081,
])

V1_NO3_BED_DENIT_REFERENCE = np.array([
    0.1604902093001368,
    0.18314599024747605,
    0.20000000000000004,
    0.21840499999999996,
    0.2076969896088541,
])

V1_NO3_DENIT_REFERENCE = np.array([
    0.00032098041860027354,
    0.0005232742578499317,
    0.00075,
    0.0009706888888888891,
    0.0011328926705937502,
])

V1_ORGN_SETTLING_REFERENCE = np.array([0.05, 0.05, 0.05, 0.05, 0.05])


@pytest.fixture(scope="function")
def nh4_5cell():
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
    inst = Nitrogen(
        time_step=timedelta(minutes=5),
        denitrification_rate=0.002, denitrification_theta=1.045,
        nitrification_rate=0.1, nitrification_theta=1.083,
        sediment_denitrification_rate=0.1, sediment_denitrification_theta=1.045,
        sediment_ammonium_release_rate=0.1, sediment_ammonium_release_theta=1.074,
        ammonium_decay_rate=0.0, ammonium_decay_theta=1.0,
        nitrification_oxygen_inhibition_factor=0.6,
    )
    inst.use_nitrate = True
    inst.use_ammonium = True
    inst.use_floating_algae = False
    inst.use_benthic_algae = False
    return inst


# ---------------------------------------------------------------------------
# Frozen v1 parity tests
# ---------------------------------------------------------------------------


def test_nitrification_inhibition_matches_v1(n_instance, dox_5cell):
    v3_val = n_instance.nitrification_inhibition(dox_5cell)
    np.testing.assert_allclose(
        v3_val.values, V1_NITR_INHIB_REFERENCE, rtol=1e-6
    )


def test_ammonium_nitrification_matches_v1(
    n_instance, nh4_5cell, water_temp_5cell, dox_5cell
):
    v3_rate = n_instance.ammonium_nitrification(nh4_5cell, water_temp_5cell, dox_5cell)
    np.testing.assert_allclose(
        v3_rate.values, V1_NH4_NITR_REFERENCE, rtol=1e-6
    )


def test_ammonium_from_bed_matches_v1(n_instance, depth_5cell, water_temp_5cell):
    v3_rate = n_instance.ammonium_from_bed(depth=depth_5cell, temperature=water_temp_5cell)
    np.testing.assert_allclose(
        v3_rate.values, V1_NH4_FROM_BED_REFERENCE, rtol=1e-6
    )


def test_nitrate_bed_denitrification_matches_v1(
    n_instance, depth_5cell, no3_5cell, water_temp_5cell
):
    v3_rate = n_instance.nitrate_bed_denitrification(
        depth=depth_5cell, nitrate=no3_5cell, temperature=water_temp_5cell
    )
    np.testing.assert_allclose(
        v3_rate.values, V1_NO3_BED_DENIT_REFERENCE, rtol=1e-6
    )


def test_nitrate_denitrification_water_column_matches_v1(
    n_instance, no3_5cell, water_temp_5cell, dox_5cell
):
    half_sat_o2 = 1.0
    v3_rate = n_instance.nitrate_denitrification(
        dissolved_oxygen=dox_5cell,
        half_saturation_oxygen=half_sat_o2,
        nitrate=no3_5cell,
        temperature=water_temp_5cell,
    )
    np.testing.assert_allclose(
        np.asarray(v3_rate), V1_NO3_DENIT_REFERENCE, rtol=1e-6
    )


# ---------------------------------------------------------------------------
# Self-contained / analytical tests (no v1 dependency)
# ---------------------------------------------------------------------------


def test_change_ammonium_no_algae_drops_to_minus_nitrification_plus_bed(
    n_instance, nh4_5cell, no3_5cell, water_temp_5cell, depth_5cell, dox_5cell
):
    """Phase 9.A.2 N2: change_ammonium without algae = -nitrification + bed."""
    n_instance.use_floating_algae = False
    n_instance.use_benthic_algae = False

    v3_rate = n_instance.change_ammonium(
        nitrate=no3_5cell, ammonium=nh4_5cell, temperature=water_temp_5cell,
        depth=depth_5cell, oxygen_dissolved=dox_5cell,
    )

    expected = (
        - n_instance.ammonium_nitrification(nh4_5cell, water_temp_5cell, dox_5cell)
        + n_instance.ammonium_from_bed(depth=depth_5cell, temperature=water_temp_5cell)
    )
    np.testing.assert_allclose(v3_rate.values, expected.values, rtol=1e-6)
    assert not np.any(np.isnan(v3_rate.values))


def test_default_nitrogen_uses_v3_defaults_for_nitrification(
    nh4_5cell, water_temp_5cell, dox_5cell
):
    """Phase 9.A.2 N1: default Nitrogen() reads knit_20=0.1, knit_theta=1.083, KNR=0.6."""
    inst = Nitrogen()
    inst.use_nitrate = True
    inst.use_ammonium = True
    inst.use_floating_algae = False
    inst.use_benthic_algae = False

    assert inst.knit_20 == 0.1
    assert inst.knit_theta == 1.083
    assert inst.KNR == 0.6

    inhib = inst.nitrification_inhibition(dox_5cell)
    expected_inhib = 1.0 - np.exp(-0.6 * dox_5cell.values)
    np.testing.assert_allclose(np.asarray(inhib), expected_inhib, rtol=1e-12)

    rate = inst.ammonium_nitrification(nh4_5cell, water_temp_5cell, dox_5cell)
    knit_tc = 0.1 * 1.083 ** (water_temp_5cell.values - 20.0)
    expected = nh4_5cell.values * knit_tc * expected_inhib
    np.testing.assert_allclose(np.asarray(rate), expected, rtol=1e-6)


def test_default_nitrogen_uses_v3_defaults_for_denitrification(
    no3_5cell, water_temp_5cell, dox_5cell
):
    """Phase 9.A.2 N10 / Phase 9.E: default kdnit_20=0.002, kdnit_theta=1.045."""
    inst = Nitrogen()
    inst.use_nitrate = True
    inst.use_ammonium = True
    inst.use_floating_algae = False
    inst.use_benthic_algae = False

    assert inst.kdnit_20 == 0.002
    assert inst.kdnit_theta == 1.045

    rate = inst.nitrate_denitrification(
        dissolved_oxygen=dox_5cell, half_saturation_oxygen=inst.KsOxdn,
        nitrate=no3_5cell, temperature=water_temp_5cell,
    )
    kdnit_tc = 0.002 * 1.045 ** (water_temp_5cell.values - 20.0)
    expected = (
        no3_5cell.values * kdnit_tc
        * (1.0 - dox_5cell.values / (dox_5cell.values + inst.KsOxdn))
    )
    np.testing.assert_allclose(np.asarray(rate), expected, rtol=1e-6)


def test_default_nitrogen_sediment_rates_zero_at_v3_defaults(
    depth_5cell, water_temp_5cell, no3_5cell
):
    """Phase 9.A.2 N4, N11: default rnh4_20=0, vno3_20=0 (matches v1/Fortran)."""
    inst = Nitrogen()
    inst.use_nitrate = True
    inst.use_ammonium = True

    assert inst.rnh4_20 == 0.0
    assert inst.vno3_20 == 0.0

    nh4_from_bed = inst.ammonium_from_bed(depth=depth_5cell, temperature=water_temp_5cell)
    np.testing.assert_allclose(np.asarray(nh4_from_bed), 0.0)

    no3_bed_denit = inst.nitrate_bed_denitrification(
        depth=depth_5cell, nitrate=no3_5cell, temperature=water_temp_5cell
    )
    np.testing.assert_allclose(np.asarray(no3_bed_denit), 0.0)


def test_phantom_ammonium_decay_term_dropped_from_change_ammonium(
    nh4_5cell, no3_5cell, water_temp_5cell, depth_5cell, dox_5cell
):
    """Phase 9.A.2 N2: phantom ammonium_decay_nitrate term dropped from NH4 budget."""
    inst = Nitrogen(ammonium_decay_rate=99.0)
    inst.use_nitrate = True
    inst.use_ammonium = True
    inst.use_floating_algae = False
    inst.use_benthic_algae = False

    rate = inst.change_ammonium(
        nitrate=no3_5cell, ammonium=nh4_5cell, temperature=water_temp_5cell,
        depth=depth_5cell, oxygen_dissolved=dox_5cell,
    )
    expected = -inst.ammonium_nitrification(nh4_5cell, water_temp_5cell, dox_5cell)
    np.testing.assert_allclose(np.asarray(rate), np.asarray(expected), rtol=1e-6)


def test_legacy_kwargs_still_override_defaults():
    """Phase 9.A.2 wiring contract: legacy v2 kwargs override defaults and sync."""
    inst = Nitrogen(
        nitrification_rate=0.5, nitrification_theta=1.1,
        denitrification_rate=0.05,
        sediment_ammonium_release_rate=0.2,
        sediment_denitrification_rate=0.3,
        nitrification_oxygen_inhibition_factor=0.7,
    )
    assert inst.nitrification_rate == 0.5
    assert inst.knit_20 == 0.5
    assert inst.nitrification_theta == 1.1
    assert inst.knit_theta == 1.1
    assert inst.denitrification_rate == 0.05
    assert inst.kdnit_20 == 0.05
    assert inst.sediment_ammonium_release_rate == 0.2
    assert inst.rnh4_20 == 0.2
    assert inst.sediment_denitrification_rate == 0.3
    assert inst.vno3_20 == 0.3
    assert inst.nitrification_oxygen_inhibition_factor == 0.7
    assert inst.KNR == 0.7


def test_nitrate_uptake_floating_algae_uses_dynamic_split(no3_5cell, nh4_5cell):
    """Phase 9.A.2 N12: NH4 + NO3 algal-uptake paths sum to rna * algal_growth_rate."""
    inst = Nitrogen()
    inst.use_floating_algae = True
    inst.use_benthic_algae = False
    inst.use_nitrate = True
    inst.use_ammonium = True

    class FakeFloatingAlgae:
        AWn = 7.2
        AWa = 1000.0
        algal_growth_rate = xr.DataArray(np.array([0.5, 0.6, 0.7, 0.8, 1.0]))
        algal_nh4_uptake_fraction = xr.DataArray(np.array([0.2, 0.3, 0.5, 0.7, 0.9]))

        def ammonium_growth(self):
            rna = self.AWn / self.AWa
            return self.algal_nh4_uptake_fraction * rna * self.algal_growth_rate

    fake = FakeFloatingAlgae()
    inst.floating_algae_process = fake

    nh4_uptake = inst.ammonium_floating_growth()
    no3_uptake = inst.nitrate_uptake_floating_algae(
        nitrate=no3_5cell, ammonium=nh4_5cell,
        algea_growth_rate=fake.algal_growth_rate,
    )
    rna = fake.AWn / fake.AWa
    np.testing.assert_allclose(
        np.asarray(nh4_uptake + no3_uptake),
        np.asarray(rna * fake.algal_growth_rate),
        rtol=1e-12,
    )


def test_nitrate_uptake_benthic_algae_uses_dynamic_split_and_correct_units(
    no3_5cell, nh4_5cell, depth_5cell
):
    """Phase 9.A.2 N13: benthic NO3 uptake uses rnb=BWn/BWd, /depth, Fb, dynamic split."""
    inst = Nitrogen()
    inst.use_floating_algae = False
    inst.use_benthic_algae = True
    inst.use_nitrate = True
    inst.use_ammonium = True

    class FakeBenthicAlgae:
        BWn = 7.2
        BWd = 100.0
        Fb = 0.9
        balgae_growth_rate = xr.DataArray(np.array([0.5, 0.6, 0.7, 0.8, 1.0]))
        balgae_nh4_uptake_fraction = xr.DataArray(np.array([0.2, 0.3, 0.5, 0.7, 0.9]))

    fake = FakeBenthicAlgae()
    inst.benthic_algae_process = fake

    no3_uptake = inst.nitrate_uptake_benthic_algae(
        nitrate=no3_5cell, ammonium=nh4_5cell,
        algea_growth_rate=fake.balgae_growth_rate, depth=depth_5cell,
    )
    rnb = fake.BWn / fake.BWd
    expected = (
        (1.0 - fake.balgae_nh4_uptake_fraction) * rnb * fake.Fb
        * fake.balgae_growth_rate / depth_5cell
    )
    np.testing.assert_allclose(np.asarray(no3_uptake), np.asarray(expected), rtol=1e-12)


def test_change_nitrate_no_algae_drops_to_nitrification_minus_denit_minus_beddenit(
    n_instance, nh4_5cell, no3_5cell, water_temp_5cell, depth_5cell, dox_5cell
):
    """Phase 2.B Bug #9: change_nitrate uses self.KsOxdn (not hardcoded 1)."""
    n_instance.use_floating_algae = False
    n_instance.use_benthic_algae = False

    v3_rate = n_instance.change_nitrate(
        nitrate=no3_5cell, ammonium=nh4_5cell, temperature=water_temp_5cell,
        depth=depth_5cell, oxygen_dissolved=dox_5cell,
    )
    expected = (
        n_instance.ammonium_nitrification(nh4_5cell, water_temp_5cell, dox_5cell)
        - n_instance.nitrate_denitrification(
            dox_5cell, n_instance.KsOxdn, no3_5cell, water_temp_5cell
        )
        - n_instance.nitrate_bed_denitrification(depth_5cell, no3_5cell, water_temp_5cell)
    )
    np.testing.assert_allclose(np.asarray(v3_rate), np.asarray(expected), rtol=1e-6)
    assert not np.any(np.isnan(np.asarray(v3_rate)))


# ---------------------------------------------------------------------------
# Phase 9.E theta-transposition fix pin tests
# ---------------------------------------------------------------------------


def test_phase9e_kon_theta_matches_fortran():
    assert Nitrogen().kon_theta == 1.047


def test_phase9e_rnh4_theta_matches_fortran():
    assert Nitrogen().rnh4_theta == 1.074


def test_phase9e_kdnit_theta_matches_fortran():
    assert Nitrogen().kdnit_theta == 1.045


def test_phase9e_vno3_theta_matches_fortran():
    assert Nitrogen().vno3_theta == 1.08


def test_phase9e_nitrogen_theta_pairs_consistent_with_phosphorus():
    """Phase 9.E: kon/rnh4 nitrogen pair mirrors kop/rpo4 phosphorus pair."""
    from clearwater_modules_v3.parameters.phosphorus import (
        DEFAULTS as PHOSPHORUS_DEFAULTS,
    )
    from clearwater_modules_v3.parameters.nitrogen import (
        DEFAULTS as NITROGEN_DEFAULTS,
    )
    assert NITROGEN_DEFAULTS["kon_theta"] == PHOSPHORUS_DEFAULTS["kop_theta"]
    assert NITROGEN_DEFAULTS["rnh4_theta"] == PHOSPHORUS_DEFAULTS["rpo4_theta"]


def test_phase9e_orgn_settling_matches_v1_no_arrhenius(
    n_instance, water_temp_5cell, depth_5cell
):
    """Phase 9.E: v3 OrgN settling is temperature-invariant (raw vson, no Arrhenius)."""
    n_instance.vson_20 = 0.05
    n_instance.use_OrgN = True
    organic_nitrogen = xr.DataArray(np.array([0.5, 1.0, 1.5, 2.0, 3.0]))

    rate_at_5C = n_instance.organic_nitrogen_settling(
        organic_nitrogen=organic_nitrogen,
        temperature=xr.DataArray(np.full(5, 5.0)),
        depth=depth_5cell,
    )
    rate_at_30C = n_instance.organic_nitrogen_settling(
        organic_nitrogen=organic_nitrogen,
        temperature=xr.DataArray(np.full(5, 30.0)),
        depth=depth_5cell,
    )
    np.testing.assert_allclose(
        np.asarray(rate_at_5C), np.asarray(rate_at_30C), rtol=1e-12
    )

    # v3 == v1 OrgN_Settling reference exactly: vson_20 / depth * OrgN.
    # With vson_20=0.05, depth=[0.5..3.0], OrgN=[0.5..3.0]:
    # rate = 0.05 / depth * OrgN. With depth and OrgN scaled together, value
    # stays at 0.05 per cell (matches the frozen reference).
    np.testing.assert_allclose(
        np.asarray(rate_at_5C), V1_ORGN_SETTLING_REFERENCE, rtol=1e-12
    )


def test_phase9e_vson_theta_removed_from_defaults():
    """Phase 9.E: vson_theta removed from NITROGEN_DEFAULTS."""
    from clearwater_modules_v3.parameters.nitrogen import (
        DEFAULTS as NITROGEN_DEFAULTS,
    )
    assert "vson_theta" not in NITROGEN_DEFAULTS
    assert "vson_20" in NITROGEN_DEFAULTS
