"""v3 DOX kinetic regression against frozen v1 reference values.

Migration from ``tests/test_5_dox_calculations_v2.py``.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.dox import DOX, dox_sat_apha

from tests.v3.nsm1.conftest import InMemoryRegistry


V1_DOX_SAT_REFERENCE = np.array([
    10.08395641037939,
    9.467092544178955,
    9.09251445406875,
    8.743797364039978,
    8.263537028982684,
])

V1_ATM_REAERATION_REFERENCE = np.array([
    5.403638790095212,
    3.3064771119870713,
    2.09251445406875,
    0.7799280647955834,
    -1.955083497304062,
])

V1_DOX_NITRIFICATION_REFERENCE = np.array([
    0.01395013067436987,
    0.03791083492091098,
    0.06754316044659582,
    0.10635304413379856,
    0.20381569836953317,
])

V1_SOD_WITH_ATTEN_REFERENCE = np.array([
    1.1956130765856912,
    0.762854091440777,
    0.5833333333333334,
    0.49937777777777786,
    0.4055229023030304,
])

V1_DOX_AP_GROWTH_REFERENCE = np.array([
    0.06621383647798741,
    0.07559245283018867,
    0.08593710691823898,
    0.0956377358490566,
    0.1131069182389937,
])


@pytest.fixture(scope="function")
def water_temp_5cell() -> xr.DataArray:
    return xr.DataArray(np.array([15.0, 18.0, 20.0, 22.0, 25.0]), dims="cell")


@pytest.fixture(scope="function")
def depth_5cell() -> xr.DataArray:
    return xr.DataArray(np.array([0.5, 1.0, 1.5, 2.0, 3.0]), dims="cell")


@pytest.fixture(scope="function")
def dox_5cell() -> xr.DataArray:
    return xr.DataArray(np.array([4.0, 6.0, 7.0, 8.0, 10.0]), dims="cell")


@pytest.fixture(scope="function")
def nh4_5cell() -> xr.DataArray:
    return xr.DataArray(np.array([0.05, 0.10, 0.15, 0.20, 0.30]), dims="cell")


def _build_registry(dox, water_temp, depth, ammonium=None) -> InMemoryRegistry:
    registry = InMemoryRegistry()
    registry.register("oxygen_dissolved", dox.copy())
    registry.register("water_temperature", water_temp.copy())
    registry.register("depth", depth.copy())
    if ammonium is not None:
        registry.register("ammonium", ammonium.copy())
    return registry


@dataclass
class _MockNitrogen:
    nitrification_flux_rate: xr.DataArray
    denitrification_flux_rate: float = 0.0


@dataclass
class _MockFloatingAlgae:
    algal_growth_rate: xr.DataArray
    algal_respiration_rate: xr.DataArray
    algal_nh4_uptake_fraction: xr.DataArray | float


def test_o2sat_full_apha_matches_v1(water_temp_5cell):
    """v3 ``dox_sat_apha`` matches frozen v1 ``DOX_sat`` reference."""
    pressure_mb = 1013.25
    v3_value = dox_sat_apha(water_temp_5cell, pressure_mb)
    np.testing.assert_allclose(
        np.asarray(v3_value), V1_DOX_SAT_REFERENCE, rtol=1e-6
    )


def test_atmospheric_reaeration_matches_v1(
    water_temp_5cell, depth_5cell, dox_5cell
):
    """v3 ``atm_reaeration_rate`` matches frozen v1 reference."""
    kah_20_user = 1.0
    pressure_mb = 1013.25

    dox_proc = DOX(
        parameters={
            "kah_20_user": kah_20_user, "kah_theta": 1.024,
            "kaw_20_user": 0.0, "kaw_theta": 1.024,
            "hydraulic_reaeration_option": 1, "wind_reaeration_option": 1,
            "pressure_mb": pressure_mb, "SOD_20": 0.0,
        },
        time_step=timedelta(minutes=5),
    )
    registry = _build_registry(dox_5cell, water_temp_5cell, depth_5cell)
    dox_proc.run(datetime(2026, 1, 1), registry)

    np.testing.assert_allclose(
        np.asarray(dox_proc.atm_reaeration_rate),
        V1_ATM_REAERATION_REFERENCE,
        rtol=1e-6,
    )


def test_dox_nitrification_sink_matches_v1(
    water_temp_5cell, depth_5cell, dox_5cell, nh4_5cell
):
    """v3 ``dox_nitrification_rate`` matches frozen v1 reference."""
    KNR = 0.6
    knit_20 = 0.1
    knit_theta = 1.083
    ron = 2.0 * 32.0 / 14.0

    from clearwater_modules_v3.utils.conversions import arrhenius_correction
    knit_tc = arrhenius_correction(water_temp_5cell, knit_20, knit_theta)
    nitrification_flux = (
        (1.0 - np.exp(-KNR * dox_5cell)) * knit_tc * nh4_5cell
    )
    mock_nitrogen = _MockNitrogen(
        nitrification_flux_rate=nitrification_flux,
        denitrification_flux_rate=0.0,
    )

    dox_proc = DOX(
        parameters={
            "ron": ron, "KNR": KNR,
            "knit_20": knit_20, "knit_theta": knit_theta,
            "kah_20_user": 0.0, "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1, "wind_reaeration_option": 1,
            "SOD_20": 0.0,
        },
        time_step=timedelta(minutes=5),
    )
    dox_proc.use_nitrogen = True
    dox_proc.nitrogen_process = mock_nitrogen
    dox_proc.use_NH4 = True

    registry = _build_registry(
        dox_5cell, water_temp_5cell, depth_5cell, nh4_5cell
    )
    dox_proc.run(datetime(2026, 1, 1), registry)

    np.testing.assert_allclose(
        np.asarray(dox_proc.dox_nitrification_rate),
        V1_DOX_NITRIFICATION_REFERENCE,
        rtol=1e-6,
    )


def test_dox_sod_sink_matches_v1_with_attenuation(
    water_temp_5cell, depth_5cell, dox_5cell
):
    """v3 ``dox_sod_rate`` matches frozen v1 reference (with attenuation)."""
    SOD_20 = 1.0
    dox_proc = DOX(
        parameters={
            "SOD_20": SOD_20, "SOD_theta": 1.060, "KsSOD": 1.0,
            "use_DOX": True,
            "kah_20_user": 0.0, "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1, "wind_reaeration_option": 1,
        },
        time_step=timedelta(minutes=5),
    )
    registry = _build_registry(dox_5cell, water_temp_5cell, depth_5cell)
    dox_proc.run(datetime(2026, 1, 1), registry)

    np.testing.assert_allclose(
        np.asarray(dox_proc.dox_sod_rate),
        V1_SOD_WITH_ATTEN_REFERENCE,
        rtol=1e-6,
    )


def test_dox_algal_photosynthesis_source_matches_fortran_anchored(
    water_temp_5cell, depth_5cell, dox_5cell, nh4_5cell
):
    """v3 floating-algae photosynthesis O2 source matches Fortran-anchored
    expected value and frozen v1 reference."""
    AWc = 40.0
    AWa = 1000.0
    rca = AWc / AWa
    roc = 32.0 / 12.0

    algal_growth = xr.DataArray(np.array([0.5, 0.6, 0.7, 0.8, 1.0]), dims="cell")
    nh4_uptake_fr = xr.DataArray(np.array([0.2, 0.4, 0.5, 0.6, 0.8]), dims="cell")
    mock_algae = _MockFloatingAlgae(
        algal_growth_rate=algal_growth,
        algal_respiration_rate=xr.zeros_like(algal_growth),
        algal_nh4_uptake_fraction=nh4_uptake_fr,
    )

    dox_proc = DOX(
        parameters={
            "AWc": AWc, "AWa": AWa, "roc": roc,
            "kah_20_user": 0.0, "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1, "wind_reaeration_option": 1,
            "SOD_20": 0.0,
        },
        time_step=timedelta(minutes=5),
    )
    dox_proc.use_floating_algae = True
    dox_proc.use_Algae = True
    dox_proc.floating_algae_process = mock_algae

    registry = _build_registry(
        dox_5cell, water_temp_5cell, depth_5cell, nh4_5cell
    )
    dox_initial = registry.get_at_time("oxygen_dissolved", datetime(2026, 1, 1)).copy()
    dox_proc.run(datetime(2026, 1, 1), registry)
    dox_final = registry.get_at_time("oxygen_dissolved", datetime(2026, 1, 1))

    dt_days = timedelta(minutes=5).total_seconds() / 86400.0
    v3_dox_rate = (dox_final - dox_initial) / dt_days

    expected_rate = (
        algal_growth * rca * roc
        * (138.0 / 106.0 - 32.0 / 106.0 * nh4_uptake_fr)
    )
    np.testing.assert_allclose(
        np.asarray(v3_dox_rate), np.asarray(expected_rate), rtol=1e-6
    )

    # Cross-check against frozen v1 DOX_ApGrowth reference (computed with
    # the correct rca = AWc/AWa).
    np.testing.assert_allclose(
        np.asarray(v3_dox_rate), V1_DOX_AP_GROWTH_REFERENCE, rtol=1e-6
    )


def test_dox_algal_photosynthesis_uses_correct_rca(
    water_temp_5cell, depth_5cell, dox_5cell, nh4_5cell
):
    """Audit C1 (Phase 9.B): v3 uses rca = 0.04 mg-C/ug-Chla, not raw AWc."""
    rca = 40.0 / 1000.0
    roc = 32.0 / 12.0
    ap_growth = xr.DataArray(np.full(5, 0.5), dims="cell")
    nh4_fr = xr.DataArray(np.full(5, 1.0), dims="cell")
    mock_algae = _MockFloatingAlgae(
        algal_growth_rate=ap_growth,
        algal_respiration_rate=xr.zeros_like(ap_growth),
        algal_nh4_uptake_fraction=nh4_fr,
    )

    dox_proc = DOX(
        parameters={
            "kah_20_user": 0.0, "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1, "wind_reaeration_option": 1,
            "SOD_20": 0.0,
        },
        time_step=timedelta(minutes=5),
    )
    dox_proc.use_floating_algae = True
    dox_proc.use_Algae = True
    dox_proc.floating_algae_process = mock_algae

    registry = _build_registry(dox_5cell, water_temp_5cell, depth_5cell, nh4_5cell)
    dox_initial = registry.get_at_time("oxygen_dissolved", datetime(2026, 1, 1)).copy()
    dox_proc.run(datetime(2026, 1, 1), registry)
    dox_final = registry.get_at_time("oxygen_dissolved", datetime(2026, 1, 1))

    dt_days = timedelta(minutes=5).total_seconds() / 86400.0
    rate = (dox_final - dox_initial) / dt_days

    expected_rate_corrected = 0.5 * rca * roc * 1.0
    np.testing.assert_allclose(
        np.asarray(rate), np.full(5, expected_rate_corrected), rtol=1e-6
    )

    raw_AWc_rate = 0.5 * 40.0 * roc * 1.0
    assert not np.allclose(np.asarray(rate), np.full(5, raw_AWc_rate), rtol=1e-3)


def test_dox_sod_attenuates_at_low_dox(water_temp_5cell, depth_5cell):
    """Audit C2: SOD attenuates to 50% at DOX=KsSOD, approaches 0 at hypoxia."""
    SOD_20 = 1.0
    SOD_theta = 1.060
    KsSOD = 1.0
    dox_at_ksod = xr.DataArray(np.full(5, KsSOD), dims="cell")

    dox_proc = DOX(
        parameters={
            "SOD_20": SOD_20, "SOD_theta": SOD_theta, "KsSOD": KsSOD,
            "use_DOX": True,
            "kah_20_user": 0.0, "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1, "wind_reaeration_option": 1,
        },
        time_step=timedelta(minutes=5),
    )
    registry = _build_registry(dox_at_ksod, water_temp_5cell, depth_5cell)
    dox_proc.run(datetime(2026, 1, 1), registry)
    attenuated_sod = np.asarray(dox_proc.dox_sod_rate)

    dox_proc_unattenuated = DOX(
        parameters={
            "SOD_20": SOD_20, "SOD_theta": SOD_theta, "KsSOD": KsSOD,
            "use_DOX": False,
            "kah_20_user": 0.0, "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1, "wind_reaeration_option": 1,
        },
        time_step=timedelta(minutes=5),
    )
    registry2 = _build_registry(dox_at_ksod.copy(), water_temp_5cell, depth_5cell)
    dox_proc_unattenuated.run(datetime(2026, 1, 1), registry2)
    unattenuated_sod = np.asarray(dox_proc_unattenuated.dox_sod_rate)

    np.testing.assert_allclose(attenuated_sod, 0.5 * unattenuated_sod, rtol=1e-12)

    dox_zero = xr.DataArray(np.full(5, 1e-6), dims="cell")
    dox_proc_low = DOX(
        parameters={
            "SOD_20": SOD_20, "SOD_theta": SOD_theta, "KsSOD": KsSOD,
            "use_DOX": True,
            "kah_20_user": 0.0, "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1, "wind_reaeration_option": 1,
        },
        time_step=timedelta(minutes=5),
    )
    registry3 = _build_registry(dox_zero, water_temp_5cell, depth_5cell)
    dox_proc_low.run(datetime(2026, 1, 1), registry3)
    sod_at_low_dox = np.asarray(dox_proc_low.dox_sod_rate)
    assert np.all(sod_at_low_dox < 1e-5 * unattenuated_sod)


def test_phase9e_sod_20_value_pinned():
    """Phase 9.E: SOD_20 = 1.0 (conservative midpoint of Chapra 1997)."""
    from clearwater_modules_v3.parameters.dox import DEFAULTS as DOX_DEFAULTS

    sod_20 = DOX_DEFAULTS["SOD_20"]
    sod_theta = DOX_DEFAULTS["SOD_theta"]
    assert sod_20 == 1.0
    assert sod_theta == 1.060
    assert sod_20 > 0.2
    assert 0.2 <= sod_20 <= 3.0


def test_phase9e_default_hydraulic_reaeration_option_is_5():
    """Phase 9.E: hydraulic_reaeration_option default = 5 (Cover 1976)."""
    from clearwater_modules_v3.parameters.dox import DEFAULTS as DOX_DEFAULTS
    assert DOX_DEFAULTS["hydraulic_reaeration_option"] == 5
    assert DOX_DEFAULTS["kah_20_user"] == 0.0


def test_phase9e_default_dox_reaeration_uses_hydraulics():
    """Phase 9.E: default reaeration must use stream hydraulics."""
    from clearwater_modules_v3.utils.reaeration import kah_20
    from clearwater_modules_v3.parameters.dox import DEFAULTS as DOX_DEFAULTS

    velocity = xr.DataArray(np.array([0.3, 0.5, 1.0]), dims="cell")
    depth = xr.DataArray(np.array([2.0, 1.0, 0.5]), dims="cell")
    flow = xr.DataArray(np.array([1.0, 1.0, 1.0]), dims="cell")
    topwidth = xr.DataArray(np.array([10.0, 10.0, 10.0]), dims="cell")
    slope = xr.DataArray(np.array([0.001, 0.001, 0.001]), dims="cell")
    shear_v = xr.DataArray(np.array([0.05, 0.05, 0.05]), dims="cell")

    kah_value = kah_20(
        kah_20_user=xr.DataArray(np.array([DOX_DEFAULTS["kah_20_user"]] * 3), dims="cell"),
        hydraulic_reaeration_option=xr.DataArray(
            np.array([DOX_DEFAULTS["hydraulic_reaeration_option"]] * 3), dims="cell"
        ),
        velocity=velocity,
        depth=depth,
        flow=flow,
        topwidth=topwidth,
        slope=slope,
        shear_velocity=shear_v,
    )
    kah_arr = np.asarray(kah_value)

    assert np.all(kah_arr > 0)
    assert not np.allclose(kah_arr, kah_arr[0])
