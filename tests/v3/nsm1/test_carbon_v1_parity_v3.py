"""v3 Carbon kinetic regression against frozen v1 reference values.

Migration from ``tests/test_5_carbon_calculations_v2.py``.

Notes on Phase 9.B / 9.E v3-vs-v1 deviations:
- POC hydrolysis: v3 no longer applies the DOX-Monod factor (Phase
  9.B audit C4); v1 (no factor) and v3 (no factor) now agree.
- DIC CO2 reaeration: v3 multiplies the Henry's-law-equilibrium term
  by 12000 to convert mol-C/L -> mg-C/L (Phase 9.E unit
  reconciliation); v1 leaves the unit mismatch. v3 rate = v1 rate +
  delta where delta is the 12000-1 correction on the equilibrium
  contribution.
- DIC algal respiration: v3 omits the legacy ``/12000``; v3 rate =
  v1_legacy * 12000 (mg-C/L/d vs mol-C/L/d).
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.carbon import Carbon, henrys_k_co2
from clearwater_modules_v3.utils.conversions import arrhenius_correction

from tests.v3.nsm1.conftest import InMemoryRegistry


V1_POC_HYDROLYSIS_REFERENCE = np.array([
    0.003974079913422963,
    0.0054734088116408465,
    0.006999999999999999,
    0.008769671999999997,
    0.011323375719750058,
])

V1_POC_SETTLING_REFERENCE = np.array([
    1.0,
    0.6,
    0.4666666666666666,
    0.4,
    0.3,
])

V1_DOC_OXIDATION_REFERENCE = np.array([
    0.012717055722953483,
    0.01720214197944266,
    0.021,
    0.025334607999999995,
    0.0320257091063638,
])

V3_CO2_REAERATION_EXPECTED_REFERENCE = np.array([
    -0.6477866148499153,
    -0.7996433443549138,
    -0.940821979012813,
    -1.0929443026674583,
    -1.2914019368515635,
])

V1_CO2_REAERATION_LEGACY_REFERENCE = np.array([
    -0.819774347877671,
    -0.9682514815981748,
    -1.1075861018315845,
    -1.2581725722300223,
    -1.4548742357951416,
])

V1_DIC_ALGAL_RESP_LEGACY_REFERENCE = np.array([
    1.6666666666666667e-06,
    2e-06,
    2.333333333333333e-06,
    2.666666666666667e-06,
    3.3333333333333333e-06,
])

V1_HENRYS_K_REFERENCE = np.array([
    0.04565115048616903,
    0.041680499541667006,
    0.03931489273426968,
    0.037148297328475896,
    0.03422936627957086,
])


@pytest.fixture(scope="function")
def water_temp_5cell() -> xr.DataArray:
    return xr.DataArray(np.array([15.0, 18.0, 20.0, 22.0, 25.0]), dims="cell")


@pytest.fixture(scope="function")
def depth_5cell() -> xr.DataArray:
    return xr.DataArray(np.array([0.5, 1.0, 1.5, 2.0, 3.0]), dims="cell")


@pytest.fixture(scope="function")
def poc_5cell() -> xr.DataArray:
    return xr.DataArray(np.array([1.0, 1.2, 1.4, 1.6, 1.8]), dims="cell")


@pytest.fixture(scope="function")
def doc_5cell() -> xr.DataArray:
    return xr.DataArray(np.array([2.0, 2.2, 2.4, 2.6, 2.8]), dims="cell")


@pytest.fixture(scope="function")
def dic_5cell() -> xr.DataArray:
    return xr.DataArray(np.array([5.0, 5.5, 6.0, 6.5, 7.0]), dims="cell")


@pytest.fixture(scope="function")
def dox_5cell() -> xr.DataArray:
    return xr.DataArray(np.array([4.0, 6.0, 7.0, 8.0, 10.0]), dims="cell")


def _build_registry(poc, doc, dic, water_temp, depth, dox=None) -> InMemoryRegistry:
    registry = InMemoryRegistry()
    registry.register("poc", poc.copy())
    registry.register("doc", doc.copy())
    registry.register("dic", dic.copy())
    registry.register("water_temperature", water_temp.copy())
    registry.register("depth", depth.copy())
    if dox is not None:
        registry.register("oxygen_dissolved", dox.copy())
    return registry


@dataclass
class _MockFloatingAlgae:
    algal_growth_rate: xr.DataArray
    algal_respiration_rate: xr.DataArray
    algal_poc_from_mortality_rate: float = 0.0
    algal_doc_from_mortality_rate: float = 0.0


@dataclass
class _MockCBOD:
    cbod_oxidation_rate: xr.DataArray | float


def test_poc_hydrolysis_matches_v1(
    water_temp_5cell, depth_5cell, poc_5cell, doc_5cell, dic_5cell, dox_5cell
):
    carbon = Carbon(
        parameters={
            "kpoc_20": 0.005, "kpoc_theta": 1.047,
            "KsOxmc": 1.0,
            "vsoc": 0.0, "kdoc_20": 0.0, "JDIC": 0.0,
            "kah_20_user": 0.0, "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1, "wind_reaeration_option": 1,
        },
        time_step=timedelta(minutes=5),
    )
    registry = _build_registry(
        poc_5cell, doc_5cell, dic_5cell, water_temp_5cell, depth_5cell, dox_5cell
    )
    carbon.run(datetime(2026, 1, 1), registry)
    np.testing.assert_allclose(
        np.asarray(carbon.poc_hydrolysis_rate),
        V1_POC_HYDROLYSIS_REFERENCE,
        rtol=1e-6,
    )


def test_poc_settling_matches_v1(
    water_temp_5cell, depth_5cell, poc_5cell, doc_5cell, dic_5cell, dox_5cell
):
    carbon = Carbon(
        parameters={
            "vsoc": 0.5, "kpoc_20": 0.0, "kdoc_20": 0.0, "JDIC": 0.0,
            "kah_20_user": 0.0, "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1, "wind_reaeration_option": 1,
        },
        time_step=timedelta(minutes=5),
    )
    registry = _build_registry(
        poc_5cell, doc_5cell, dic_5cell, water_temp_5cell, depth_5cell, dox_5cell
    )
    poc_initial = registry.get_at_time("poc", datetime(2026, 1, 1)).copy()
    carbon.run(datetime(2026, 1, 1), registry)
    poc_final = registry.get_at_time("poc", datetime(2026, 1, 1))

    dt_days = timedelta(minutes=5).total_seconds() / 86400.0
    v3_settling_rate = -(poc_final - poc_initial) / dt_days

    np.testing.assert_allclose(
        np.asarray(v3_settling_rate),
        V1_POC_SETTLING_REFERENCE,
        rtol=1e-6,
    )


def test_doc_oxidation_matches_v1(
    water_temp_5cell, depth_5cell, poc_5cell, doc_5cell, dic_5cell, dox_5cell
):
    carbon = Carbon(
        parameters={
            "kdoc_20": 0.01, "kdoc_theta": 1.047, "KsOxmc": 1.0,
            "kpoc_20": 0.0, "vsoc": 0.0, "JDIC": 0.0,
            "kah_20_user": 0.0, "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1, "wind_reaeration_option": 1,
        },
        time_step=timedelta(minutes=5),
    )
    registry = _build_registry(
        poc_5cell, doc_5cell, dic_5cell, water_temp_5cell, depth_5cell, dox_5cell
    )
    carbon.run(datetime(2026, 1, 1), registry)

    np.testing.assert_allclose(
        np.asarray(carbon.doc_dic_oxidation_rate),
        V1_DOC_OXIDATION_REFERENCE,
        rtol=1e-6,
    )


def test_dic_co2_reaeration_matches_v1(
    water_temp_5cell, depth_5cell, poc_5cell, doc_5cell, dic_5cell, dox_5cell
):
    """v3 DIC CO2 reaeration rate (mg-C/L/d, Phase 9.E unit-corrected)
    matches frozen v3-expected reference. Also asserts the v3 vs v1
    legacy delta is the (12000-1) correction on the equilibrium term."""
    kah_20_user = 1.0
    carbon = Carbon(
        parameters={
            "pCO2": 383.0, "FCO2": 0.2,
            "kah_20_user": kah_20_user, "kah_theta": 1.024,
            "kaw_20_user": 0.0, "kaw_theta": 1.024,
            "hydraulic_reaeration_option": 1, "wind_reaeration_option": 1,
            "kdoc_20": 0.0, "kpoc_20": 0.0, "vsoc": 0.0, "JDIC": 0.0,
        },
        time_step=timedelta(minutes=5),
    )
    registry = _build_registry(
        poc_5cell, doc_5cell, dic_5cell, water_temp_5cell, depth_5cell, dox_5cell
    )
    dic_initial = registry.get_at_time("dic", datetime(2026, 1, 1)).copy()
    carbon.run(datetime(2026, 1, 1), registry)
    dic_final = registry.get_at_time("dic", datetime(2026, 1, 1))

    dt_days = timedelta(minutes=5).total_seconds() / 86400.0
    v3_dic_rate = (dic_final - dic_initial) / dt_days

    np.testing.assert_allclose(
        np.asarray(v3_dic_rate),
        V3_CO2_REAERATION_EXPECTED_REFERENCE,
        rtol=1e-6,
    )

    # Cross-check the v3 vs v1-legacy delta relationship.
    ka_tc = arrhenius_correction(water_temp_5cell, kah_20_user, 1.024)
    kh = henrys_k_co2(water_temp_5cell)
    delta = 0.923 * ka_tc * kh * 383.0 / 1.0e6 * (12000.0 - 1.0)
    np.testing.assert_allclose(
        np.asarray(v3_dic_rate - V1_CO2_REAERATION_LEGACY_REFERENCE),
        np.asarray(delta),
        rtol=1e-6,
    )


def test_dic_algal_respiration_source_matches_fortran_anchored(
    water_temp_5cell, depth_5cell, poc_5cell, doc_5cell, dic_5cell, dox_5cell
):
    """v3 DIC algal-respiration coupling = rca * ApRespiration (mg-C/L/d).
    Phase 9.E removed the legacy ``/ 12000`` so v3 = v1_legacy * 12000."""
    AWc = 40.0
    AWa = 1000.0
    rca = AWc / AWa

    algal_resp = xr.DataArray(np.array([0.5, 0.6, 0.7, 0.8, 1.0]), dims="cell")
    mock_algae = _MockFloatingAlgae(
        algal_growth_rate=xr.zeros_like(algal_resp),
        algal_respiration_rate=algal_resp,
    )

    carbon = Carbon(
        parameters={
            "AWc": AWc, "AWa": AWa,
            "kdoc_20": 0.0, "kpoc_20": 0.0, "vsoc": 0.0, "JDIC": 0.0,
            "kah_20_user": 0.0, "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1, "wind_reaeration_option": 1,
        },
        time_step=timedelta(minutes=5),
    )
    carbon.use_floating_algae = True
    carbon.use_Algae = True
    carbon.floating_algae_process = mock_algae

    registry = _build_registry(
        poc_5cell, doc_5cell, dic_5cell, water_temp_5cell, depth_5cell, dox_5cell
    )
    dic_initial = registry.get_at_time("dic", datetime(2026, 1, 1)).copy()
    carbon.run(datetime(2026, 1, 1), registry)
    dic_final = registry.get_at_time("dic", datetime(2026, 1, 1))

    dt_days = timedelta(minutes=5).total_seconds() / 86400.0
    v3_dic_rate = (dic_final - dic_initial) / dt_days

    expected_rate = rca * algal_resp
    np.testing.assert_allclose(
        np.asarray(v3_dic_rate), np.asarray(expected_rate), rtol=1e-6
    )

    # v3 = v1_legacy * 12000 (Phase 9.E unit promotion).
    np.testing.assert_allclose(
        np.asarray(v3_dic_rate),
        V1_DIC_ALGAL_RESP_LEGACY_REFERENCE * 12000.0,
        rtol=1e-6,
    )


def test_dic_algal_growth_uses_correct_rca(
    water_temp_5cell, depth_5cell, poc_5cell, doc_5cell, dic_5cell, dox_5cell
):
    """Audit C1 (Phase 9.B): v3 dDIC/dt from algal growth uses
    ``rca = 0.04`` mg-C/ug-Chla, not the raw ``AWc = 40``."""
    AWc = 40.0
    AWa = 1000.0
    rca = AWc / AWa

    algal_growth = xr.DataArray(np.array([1.0, 2.0, 3.0, 4.0, 5.0]), dims="cell")
    mock_algae = _MockFloatingAlgae(
        algal_growth_rate=algal_growth,
        algal_respiration_rate=xr.zeros_like(algal_growth),
    )

    carbon = Carbon(
        parameters={
            "AWc": AWc, "AWa": AWa,
            "kdoc_20": 0.0, "kpoc_20": 0.0, "vsoc": 0.0, "JDIC": 0.0,
            "kah_20_user": 0.0, "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1, "wind_reaeration_option": 1,
        },
        time_step=timedelta(minutes=5),
    )
    carbon.use_floating_algae = True
    carbon.use_Algae = True
    carbon.floating_algae_process = mock_algae

    registry = _build_registry(
        poc_5cell, doc_5cell, dic_5cell, water_temp_5cell, depth_5cell, dox_5cell
    )
    dic_initial = registry.get_at_time("dic", datetime(2026, 1, 1)).copy()
    carbon.run(datetime(2026, 1, 1), registry)
    dic_final = registry.get_at_time("dic", datetime(2026, 1, 1))

    dt_days = timedelta(minutes=5).total_seconds() / 86400.0
    v3_dic_rate = (dic_final - dic_initial) / dt_days

    expected_rate = -rca * algal_growth
    np.testing.assert_allclose(
        np.asarray(v3_dic_rate), np.asarray(expected_rate), rtol=1e-6
    )

    # Negative regression: confirm NOT using raw AWc.
    raw_rate = -AWc * algal_growth
    assert not np.allclose(
        np.asarray(v3_dic_rate), np.asarray(raw_rate), rtol=1e-3
    )


def test_dic_includes_cbod_oxidation(
    water_temp_5cell, depth_5cell, poc_5cell, doc_5cell, dic_5cell, dox_5cell
):
    """Audit C3 (Phase 9.B): v3 dDIC/dt includes the source
    ``cbod_oxidation_rate / roc`` (Phase 9.E mg-C/L/d convention)."""
    roc = 32.0 / 12.0
    cbod_ox_rate = xr.DataArray(
        np.array([0.10, 0.20, 0.30, 0.40, 0.50]), dims="cell"
    )
    mock_cbod = _MockCBOD(cbod_oxidation_rate=cbod_ox_rate)

    carbon = Carbon(
        parameters={
            "roc": roc,
            "kdoc_20": 0.0, "kpoc_20": 0.0, "vsoc": 0.0, "JDIC": 0.0,
            "kah_20_user": 0.0, "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1, "wind_reaeration_option": 1,
        },
        time_step=timedelta(minutes=5),
    )
    carbon.use_cbod = True
    carbon.cbod_process = mock_cbod

    registry = _build_registry(
        poc_5cell, doc_5cell, dic_5cell, water_temp_5cell, depth_5cell, dox_5cell
    )
    dic_initial = registry.get_at_time("dic", datetime(2026, 1, 1)).copy()
    carbon.run(datetime(2026, 1, 1), registry)
    dic_final = registry.get_at_time("dic", datetime(2026, 1, 1))

    dt_days = timedelta(minutes=5).total_seconds() / 86400.0
    v3_dic_rate = (dic_final - dic_initial) / dt_days

    expected_rate = cbod_ox_rate / roc
    np.testing.assert_allclose(
        np.asarray(v3_dic_rate), np.asarray(expected_rate), rtol=1e-6
    )

    assert np.all(np.asarray(v3_dic_rate) > 0)


def test_poc_hydrolysis_no_longer_dox_attenuated(
    water_temp_5cell, depth_5cell, poc_5cell, doc_5cell, dic_5cell
):
    """Audit C4 (Phase 9.B): POC hydrolysis is independent of DOX."""
    kpoc_20 = 0.005
    kpoc_theta = 1.047
    KsOxmc = 1.0

    def _run_with_dox(dox_value):
        dox_array = xr.DataArray(np.full(5, dox_value), dims="cell")
        carbon = Carbon(
            parameters={
                "kpoc_20": kpoc_20, "kpoc_theta": kpoc_theta, "KsOxmc": KsOxmc,
                "vsoc": 0.0, "kdoc_20": 0.0, "JDIC": 0.0,
                "kah_20_user": 0.0, "kaw_20_user": 0.0,
                "hydraulic_reaeration_option": 1, "wind_reaeration_option": 1,
            },
            time_step=timedelta(minutes=5),
        )
        registry = _build_registry(
            poc_5cell, doc_5cell, dic_5cell, water_temp_5cell, depth_5cell, dox_array
        )
        carbon.run(datetime(2026, 1, 1), registry)
        return carbon.poc_hydrolysis_rate

    rate_low_dox = _run_with_dox(0.5)
    rate_high_dox = _run_with_dox(20.0)

    np.testing.assert_allclose(
        np.asarray(rate_low_dox), np.asarray(rate_high_dox), rtol=1e-12
    )


def test_henrys_k_co2_matches_v1(water_temp_5cell):
    """v3 ``henrys_k_co2`` matches frozen v1 ``Henrys_k`` reference."""
    v3_value = henrys_k_co2(water_temp_5cell)
    np.testing.assert_allclose(
        np.asarray(v3_value), V1_HENRYS_K_REFERENCE, rtol=1e-6
    )
