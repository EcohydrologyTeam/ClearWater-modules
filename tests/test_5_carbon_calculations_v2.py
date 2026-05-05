"""Parity tests: v3 Carbon sub-rate methods vs v1 nsm1.processes helpers.

Phase 7.B of the v3 NSM1 implementation plan. Mirrors the established
Phase 2 pattern at ``tests/test_5_floating_algae_calculations_v2.py``.

Each test constructs a v3 Carbon instance (POC + DOC + DIC), drives it
through one ``Process.run`` substep against an in-memory registry, and
compares either a cached step-scoped rate (``poc_hydrolysis_rate``,
``doc_dic_oxidation_rate``) or the back-calculated dDIC/dt to the
equivalent v1 helper-function output computed with the same inputs.

Scope: v1-equivalent ``dPOCdt`` / ``dDOCdt`` / ``dDICdt`` sub-terms:

* POC hydrolysis (``kpoc_tc * POC * DOX_attenuation``) -- v3 deviation
  documented inline (v1 uses ``kpoc_tc * POC`` with no DOX coupling).
* POC settling (``vsoc / depth * POC``)
* DOC oxidation (``kdoc_tc * DOC * DOX_attenuation``)
* DIC CO2 atmospheric reaeration
  (``0.923 * ka_tc * (KH * pCO2 / 1e6 - FCO2 * DIC)``)
* DIC algal respiration coupling (``ApRespiration * rca / 12000``)

Synthetic mesh: 5-cell xarray DataArrays, single time step.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules.nsm1 import processes as v1
from clearwater_modules_v3.processes.carbon import Carbon, henrys_k_co2

from tests.v3.nsm1.conftest import InMemoryRegistry


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="function")
def water_temp_5cell() -> xr.DataArray:
    return xr.DataArray(
        np.array([15.0, 18.0, 20.0, 22.0, 25.0]), dims="cell"
    )


@pytest.fixture(scope="function")
def depth_5cell() -> xr.DataArray:
    return xr.DataArray(
        np.array([0.5, 1.0, 1.5, 2.0, 3.0]), dims="cell"
    )


@pytest.fixture(scope="function")
def poc_5cell() -> xr.DataArray:
    return xr.DataArray(
        np.array([1.0, 1.2, 1.4, 1.6, 1.8]), dims="cell"
    )


@pytest.fixture(scope="function")
def doc_5cell() -> xr.DataArray:
    return xr.DataArray(
        np.array([2.0, 2.2, 2.4, 2.6, 2.8]), dims="cell"
    )


@pytest.fixture(scope="function")
def dic_5cell() -> xr.DataArray:
    return xr.DataArray(
        np.array([5.0, 5.5, 6.0, 6.5, 7.0]), dims="cell"
    )


@pytest.fixture(scope="function")
def dox_5cell() -> xr.DataArray:
    return xr.DataArray(
        np.array([4.0, 6.0, 7.0, 8.0, 10.0]), dims="cell"
    )


def _build_registry(
    poc: xr.DataArray,
    doc: xr.DataArray,
    dic: xr.DataArray,
    water_temp: xr.DataArray,
    depth: xr.DataArray,
    dox: xr.DataArray | None = None,
) -> InMemoryRegistry:
    """Wire the state variables that Carbon.run reads."""
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
    """Stand-in FloatingAlgae for the DIC algal-coupling parity test.

    Carbon._floating_algae_respiration_rate / _growth_rate read the
    cached attributes on this object via getattr; this mock exposes
    exactly those attributes plus the C:Chla stoichiometric ratio
    (Carbon reads AWc from its own DEFAULTS, not from the algae
    process; AWc here is informational only).
    """
    algal_growth_rate: xr.DataArray
    algal_respiration_rate: xr.DataArray
    algal_poc_from_mortality_rate: float = 0.0
    algal_doc_from_mortality_rate: float = 0.0


# ---------------------------------------------------------------------------
# Parity tests
# ---------------------------------------------------------------------------


def test_poc_hydrolysis_matches_v1(
    water_temp_5cell, depth_5cell, poc_5cell, doc_5cell, dic_5cell, dox_5cell
):
    """v3 ``poc_hydrolysis_rate`` cache after run == v1 ``POC_hydrolysis``.

    v1 formula: ``kpoc_tc * POC`` (line 2455-2465; NO DOX-Monod factor).
    v3 formula: ``kpoc_tc * POC * DOX_attenuation`` where
    ``DOX_attenuation = DOX / (KsOxmc + DOX)``.

    DEVIATION (documented): v3 multiplies POC hydrolysis by the same
    DOX-Monod factor used for DOC->DIC oxidation, while v1 applies that
    factor only to DOC->DIC oxidation. To make the parity check robust,
    this test sets ``KsOxmc = 0`` so that
    ``DOX_attenuation = DOX / DOX = 1`` exactly, making the v3 formula
    collapse to v1's ``kpoc_tc * POC``.
    """
    carbon = Carbon(
        parameters={
            "kpoc_20": 0.005,
            "kpoc_theta": 1.047,
            # Force DOX_attenuation == 1 so v3's POC-hydrolysis term
            # collapses to v1's (kpoc_tc * POC). This is the documented
            # deviation: v3 added DOX-Monod to POC hydrolysis where v1
            # had none.
            "KsOxmc": 0.0,
            # Disable other terms by zeroing their coefficients.
            "vsoc": 0.0,
            "kdoc_20": 0.0,
            "JDIC": 0.0,
            "kah_20_user": 0.0,
            "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1,
            "wind_reaeration_option": 1,
        },
        time_step=timedelta(minutes=5),
    )
    registry = _build_registry(
        poc_5cell, doc_5cell, dic_5cell,
        water_temp_5cell, depth_5cell, dox_5cell,
    )
    carbon.run(datetime(2026, 1, 1), registry)

    kpoc_tc = v1.kpoc_tc(water_temp_5cell, 0.005, 1.047)
    v1_rate = v1.POC_hydrolysis(kpoc_tc, poc_5cell)

    np.testing.assert_allclose(
        np.asarray(carbon.poc_hydrolysis_rate),
        np.asarray(v1_rate),
        rtol=1e-6,
    )


def test_poc_settling_matches_v1(
    water_temp_5cell, depth_5cell, poc_5cell, doc_5cell, dic_5cell, dox_5cell
):
    """v3 POC settling sub-term == v1 ``POC_settling``.

    v1 formula: ``vsoc / depth * POC`` (line 2469-2481).
    v3 formula: same (carbon.py line 373, ``poc_settling = self.vsoc /
    depth * poc``). v3 does not cache this individually, so this test
    isolates POC settling by zeroing every other POC term and reading
    the post-run POC state to recover the integrated rate via
    ``-(poc_new - poc_old) / dt_days``.
    """
    carbon = Carbon(
        parameters={
            "vsoc": 0.5,                # m/d; non-zero
            "kpoc_20": 0.0,             # disable POC hydrolysis
            "kdoc_20": 0.0,             # disable DOC oxidation (no POC effect anyway)
            "JDIC": 0.0,
            "kah_20_user": 0.0,
            "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1,
            "wind_reaeration_option": 1,
        },
        time_step=timedelta(minutes=5),
    )
    registry = _build_registry(
        poc_5cell, doc_5cell, dic_5cell,
        water_temp_5cell, depth_5cell, dox_5cell,
    )
    poc_initial = registry.get_at_time("poc", datetime(2026, 1, 1)).copy()
    carbon.run(datetime(2026, 1, 1), registry)
    poc_final = registry.get_at_time("poc", datetime(2026, 1, 1))

    dt_days = timedelta(minutes=5).total_seconds() / 86400.0
    # POC integrator is forward Euler; with only the settling sink active,
    # poc_new = poc_old - settling_rate * dt_days, so settling_rate ==
    # -(poc_new - poc_old) / dt_days.
    v3_settling_rate = -(poc_final - poc_initial) / dt_days

    v1_rate = v1.POC_settling(0.5, depth_5cell, poc_5cell)

    np.testing.assert_allclose(
        np.asarray(v3_settling_rate),
        np.asarray(v1_rate),
        rtol=1e-6,
    )


def test_doc_oxidation_matches_v1(
    water_temp_5cell, depth_5cell, poc_5cell, doc_5cell, dic_5cell, dox_5cell
):
    """v3 ``doc_dic_oxidation_rate`` cache after run == v1 ``DOC_DIC_oxidation``.

    v1 formula: ``DOX / (KsOxmc + DOX) * kdoc_tc * DOC`` (line 2629-2647,
    when ``use_DOX==True``).
    v3 formula: ``kdoc_tc * DOC * DOX_attenuation`` (carbon.py line 385),
    same expression.

    Both should match exactly with matched parameters.
    """
    carbon = Carbon(
        parameters={
            "kdoc_20": 0.01,
            "kdoc_theta": 1.047,
            "KsOxmc": 1.0,
            # Disable other terms so doc_dic_oxidation_rate is isolated.
            "kpoc_20": 0.0,
            "vsoc": 0.0,
            "JDIC": 0.0,
            "kah_20_user": 0.0,
            "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1,
            "wind_reaeration_option": 1,
        },
        time_step=timedelta(minutes=5),
    )
    registry = _build_registry(
        poc_5cell, doc_5cell, dic_5cell,
        water_temp_5cell, depth_5cell, dox_5cell,
    )
    carbon.run(datetime(2026, 1, 1), registry)

    kdoc_tc = v1.kdoc_tc(water_temp_5cell, 0.01, 1.047)
    v1_rate = v1.DOC_DIC_oxidation(
        DOX=dox_5cell,
        KsOxmc=1.0,
        kdoc_tc=kdoc_tc,
        DOC=doc_5cell,
        use_DOX=True,
    )

    np.testing.assert_allclose(
        np.asarray(carbon.doc_dic_oxidation_rate),
        np.asarray(v1_rate),
        rtol=1e-6,
    )


def test_dic_co2_reaeration_matches_v1(
    water_temp_5cell, depth_5cell, poc_5cell, doc_5cell, dic_5cell, dox_5cell
):
    """v3 DIC CO2 atmospheric reaeration sub-term == v1 ``Atmospheric_CO2_reaeration``.

    v1 formula: ``0.923 * ka_tc * (K_H * pCO2 / 1e6 - FCO2 * DIC)`` (line
    2698-2714, where ``K_H = Henrys_k(TwaterC)``).
    v3 formula: same form (carbon.py line 418-421). v3 imports its own
    ``henrys_k_co2`` (port of v1 ``Henrys_k``) which is parity-checked
    here by reusing the v1 ``Henrys_k`` for the reference.

    The parity test isolates the CO2 reaeration term by zeroing every
    other DIC source/sink, then back-calculates the integrated DIC rate
    from the post-run DIC state.
    """
    # Use a non-zero user-defined hydraulic reaeration so ka_tc != 0.
    kah_20_user = 1.0    # 1/d
    carbon = Carbon(
        parameters={
            # CO2 reaeration knobs.
            "pCO2": 383.0,
            "FCO2": 0.2,
            "kah_20_user": kah_20_user,
            "kah_theta": 1.024,
            "kaw_20_user": 0.0,         # zero out wind contribution
            "kaw_theta": 1.024,
            "hydraulic_reaeration_option": 1,
            "wind_reaeration_option": 1,
            # Disable every other DIC source/sink.
            "kdoc_20": 0.0,
            "kpoc_20": 0.0,
            "vsoc": 0.0,
            "JDIC": 0.0,
        },
        time_step=timedelta(minutes=5),
    )
    registry = _build_registry(
        poc_5cell, doc_5cell, dic_5cell,
        water_temp_5cell, depth_5cell, dox_5cell,
    )
    dic_initial = registry.get_at_time("dic", datetime(2026, 1, 1)).copy()
    carbon.run(datetime(2026, 1, 1), registry)
    dic_final = registry.get_at_time("dic", datetime(2026, 1, 1))

    dt_days = timedelta(minutes=5).total_seconds() / 86400.0
    v3_dic_rate = (dic_final - dic_initial) / dt_days

    # v1 reference: ka_tc combined hydraulic + wind, but with wind
    # zeroed, ka_tc == arrhenius_correction(T, kah_20, kah_theta).
    v1_ka_tc = v1.arrhenius_correction(water_temp_5cell, kah_20_user, 1.024)
    v1_kh = v1.Henrys_k(water_temp_5cell)
    v1_co2_reaer = v1.Atmospheric_CO2_reaeration(
        ka_tc=v1_ka_tc,
        K_H=v1_kh,
        pCO2=383.0,
        FCO2=0.2,
        DIC=dic_5cell,
    )

    np.testing.assert_allclose(
        np.asarray(v3_dic_rate),
        np.asarray(v1_co2_reaer),
        rtol=1e-6,
    )


def test_dic_algal_respiration_source_matches_v1(
    water_temp_5cell, depth_5cell, poc_5cell, doc_5cell, dic_5cell, dox_5cell
):
    """v3 DIC algal-respiration coupling sub-term == v1 ``DIC_algal_respiration``.

    v1 formula: ``ApRespiration * rca / 12000`` (line 2717-2731), where
    ``rca = AWc / AWa = AWc`` in v3's per-Chla convention.
    v3 formula: ``algae_respiration * AWc / 12000`` (carbon.py line 429),
    same expression. v3 reads ``algal_respiration_rate`` from a
    FloatingAlgae sibling Process; the parity test wires a mock with a
    fixed ``algal_respiration_rate`` value.

    To isolate the DIC algal-respiration source, the test zeros every
    other DIC source/sink (including algal photosynthesis by setting
    ``algal_growth_rate = 0`` on the mock).
    """
    # AWc default: 40.0 mg-C/ug-Chla.
    AWc = 40.0
    # Synthetic algal respiration rate (ug-Chla/L/d).
    algal_resp = xr.DataArray(
        np.array([0.5, 0.6, 0.7, 0.8, 1.0]), dims="cell"
    )
    mock_algae = _MockFloatingAlgae(
        algal_growth_rate=xr.zeros_like(algal_resp),
        algal_respiration_rate=algal_resp,
    )

    carbon = Carbon(
        parameters={
            "AWc": AWc,
            # Disable every other DIC source/sink.
            "kdoc_20": 0.0,
            "kpoc_20": 0.0,
            "vsoc": 0.0,
            "JDIC": 0.0,
            "kah_20_user": 0.0,
            "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1,
            "wind_reaeration_option": 1,
            # ka_tc -> 0 means CO2 reaeration -> 0; FCO2 doesn't matter.
        },
        time_step=timedelta(minutes=5),
    )
    # Wire the mock FloatingAlgae directly (init_process is normally invoked
    # by Model; mimic the use-flag setup directly here).
    carbon.use_floating_algae = True
    carbon.use_Algae = True
    carbon.floating_algae_process = mock_algae

    registry = _build_registry(
        poc_5cell, doc_5cell, dic_5cell,
        water_temp_5cell, depth_5cell, dox_5cell,
    )
    dic_initial = registry.get_at_time("dic", datetime(2026, 1, 1)).copy()
    carbon.run(datetime(2026, 1, 1), registry)
    dic_final = registry.get_at_time("dic", datetime(2026, 1, 1))

    dt_days = timedelta(minutes=5).total_seconds() / 86400.0
    v3_dic_rate = (dic_final - dic_initial) / dt_days

    # v1 reference: DIC_algal_respiration with use_Algae=True.
    v1_rate = v1.DIC_algal_respiration(
        ApRespiration=algal_resp, rca=AWc, use_Algae=True
    )

    np.testing.assert_allclose(
        np.asarray(v3_dic_rate),
        np.asarray(v1_rate),
        rtol=1e-6,
    )


def test_henrys_k_co2_matches_v1(water_temp_5cell):
    """v3 ``henrys_k_co2`` (carbon.py:124) == v1 ``Henrys_k`` (line 2687-2695).

    Spot-check the Henry's law constant utility v3 ports inline. v3 uses
    the same formula:
    ``10**(2385.73 / Tk + 0.0152642 * Tk - 14.0184)`` with ``Tk =
    TwaterC + 273.15``.
    """
    v3_value = henrys_k_co2(water_temp_5cell)
    v1_value = v1.Henrys_k(water_temp_5cell)

    np.testing.assert_allclose(
        np.asarray(v3_value), np.asarray(v1_value), rtol=1e-6
    )
