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
    v3 formula (post Phase 9.B audit fix C4): ``kpoc_tc * POC`` (no
    DOX-Monod factor). v3 previously multiplied by
    ``DOX / (KsOxmc + DOX)``; the Phase 9.B audit removed that factor
    because neither Fortran (``modCarbon.f90:170``) nor v1 attenuate POC
    hydrolysis by DOX. With non-zero ``KsOxmc`` and finite DOX, v3 and
    v1 should now agree directly.
    """
    carbon = Carbon(
        parameters={
            "kpoc_20": 0.005,
            "kpoc_theta": 1.047,
            # KsOxmc is now irrelevant for POC hydrolysis (no Monod
            # attenuation applied) but still gates DOC->DIC oxidation.
            # Use the default-typical value to confirm the fix decoupled
            # POC hydrolysis from the Monod factor.
            "KsOxmc": 1.0,
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
    """v3 DIC CO2 atmospheric reaeration in mg-C/L/d.

    v1 formula: ``0.923 * ka_tc * (K_H * pCO2 / 1e6 - FCO2 * DIC)``
    (line 2698-2714). The Henry's-law term ``K_H * pCO2 / 1e6`` is in
    mol-C/L; ``FCO2 * DIC`` (with v1 DIC stored in mg-C/L) is in
    mg-C/L. v1's formula has unit-mixed inputs producing a numerically
    inconsistent rate.

    Phase 9.E unit reconciliation: v3 converts the Henry's-law term
    from mol-C/L to mg-C/L by multiplying by 12000 (= 12 g/mol-C *
    1000 mg/g) so both terms inside the parentheses are in mg-C/L
    and the resulting rate is mg-C/L/d. v3 = v1 * 12000 in the
    Henry's-law-equilibrium contribution; v3 = v1 in the
    ``-0.923 * ka_tc * FCO2 * DIC`` outflow contribution. The full
    expected is computed manually below.
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

    # v3 expected (mg-C/L/d): 0.923 * ka_tc * (KH * pCO2/1e6 * 12000 - FCO2 * DIC)
    # where ka_tc combined hydraulic+wind, but with wind zeroed,
    # ka_tc == arrhenius_correction(T, kah_20, kah_theta).
    v1_ka_tc = v1.arrhenius_correction(water_temp_5cell, kah_20_user, 1.024)
    v1_kh = v1.Henrys_k(water_temp_5cell)
    expected_v3_rate = (
        0.923
        * v1_ka_tc
        * (v1_kh * 383.0 / 1.0e6 * 12000.0 - 0.2 * dic_5cell)
    )

    np.testing.assert_allclose(
        np.asarray(v3_dic_rate),
        np.asarray(expected_v3_rate),
        rtol=1e-6,
    )

    # Cross-check: v3 = v1 * 12000 in the Henry's-law-equilibrium term
    # alone. Compute v1's mol-C/L/d form and confirm the 12000x relation
    # for the equilibrium-driven contribution.
    v1_co2_reaer_legacy = v1.Atmospheric_CO2_reaeration(
        ka_tc=v1_ka_tc,
        K_H=v1_kh,
        pCO2=383.0,
        FCO2=0.2,
        DIC=dic_5cell,
    )
    # v1: 0.923 * ka_tc * (KH * pCO2/1e6 - FCO2 * DIC)
    # v3: 0.923 * ka_tc * (KH * pCO2/1e6 * 12000 - FCO2 * DIC)
    # v3 - v1 = 0.923 * ka_tc * (KH * pCO2/1e6) * (12000 - 1)
    delta = (
        0.923 * v1_ka_tc * v1_kh * 383.0 / 1.0e6 * (12000.0 - 1.0)
    )
    np.testing.assert_allclose(
        np.asarray(v3_dic_rate - v1_co2_reaer_legacy),
        np.asarray(delta),
        rtol=1e-6,
    )


def test_dic_algal_respiration_source_matches_fortran_anchored(
    water_temp_5cell, depth_5cell, poc_5cell, doc_5cell, dic_5cell, dox_5cell
):
    """v3 DIC algal-respiration coupling sub-term in mg-C/L/d.

    Fortran formula (``modCarbon.f90:247``):
        ``ApRespiration_DIC = rca * ApRespiration / 12000.0``
    where ``rca = AWc / AWa = 40 / 1000 = 0.04 mg-C/ug-Chla``.

    The Fortran ``/ 12000.0`` produces a rate in mol-C/L/d, but Fortran
    labels DIC's units as mg-C/L (modMain.f90:301). The mismatch was
    flagged in Phase 9.B audit C9 and corrected in Phase 9.E: v3 now
    integrates dDIC/dt in mg-C/L/d throughout, so this term is
    ``rca * ApRespiration`` (no ``/ 12000``).

    Phase 9.B audit C1 (rca = AWc/AWa, not raw AWc): still applies; v3
    now derives the ratio at run time. The earlier parity test masked
    the C1 bug by passing ``rca = AWc = 40`` to the v1 reference (same
    wrong number on both sides).
    """
    AWc = 40.0      # mg-C raw weight
    AWa = 1000.0    # ug-Chla algal unit
    rca = AWc / AWa  # 0.04 mg-C / ug-Chla

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
            "AWa": AWa,
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

    # Phase 9.E expected: rate = rca * ApRespiration (mg-C/L/d). The
    # Fortran ``/ 12000`` is dropped because v3 integrates dDIC/dt in
    # mg-C/L/d directly.
    expected_rate = rca * algal_resp

    np.testing.assert_allclose(
        np.asarray(v3_dic_rate),
        np.asarray(expected_rate),
        rtol=1e-6,
    )

    # v1 helper still divides by 12000 (legacy mol-C/L/d form). v3
    # produces mg-C/L/d, which equals v1's value * 12000. Confirm the
    # 12000x relation explicitly so the divergence is documented in the
    # test rather than masked.
    v1_rate_legacy = v1.DIC_algal_respiration(
        ApRespiration=algal_resp, rca=rca, use_Algae=True
    )
    np.testing.assert_allclose(
        np.asarray(v3_dic_rate),
        np.asarray(v1_rate_legacy) * 12000.0,
        rtol=1e-6,
    )


# ---------------------------------------------------------------------------
# Audit-anchored regression tests (Phase 9.B)
# ---------------------------------------------------------------------------


@dataclass
class _MockCBOD:
    """Stand-in CBOD Process exposing ``cbod_oxidation_rate``."""
    cbod_oxidation_rate: xr.DataArray | float


def test_dic_algal_growth_uses_correct_rca(
    water_temp_5cell, depth_5cell, poc_5cell, doc_5cell, dic_5cell, dox_5cell
):
    """Audit C1 (Phase 9.B): v3 dDIC/dt from algal *growth* uses
    ``rca = 0.04`` mg-C/ug-Chla, not the raw ``AWc = 40``.

    Default-instantiated Carbon + a mock FloatingAlgae with non-zero
    ``algal_growth_rate``. After the Phase 9.E unit reconciliation
    (DIC budget integrated in mg-C/L/d throughout, no ``/ 12000``),
    the DIC sink magnitude is ``-rca * ApGrowth`` (mg-C/L/d), which
    is 1000x smaller than the prior v3 (raw-AWc) magnitude but 12000x
    larger than the legacy Fortran/v1 mol-C/L/d form.
    """
    AWc = 40.0
    AWa = 1000.0
    rca = AWc / AWa  # 0.04
    algal_growth = xr.DataArray(
        np.array([1.0, 2.0, 3.0, 4.0, 5.0]), dims="cell"
    )
    mock_algae = _MockFloatingAlgae(
        algal_growth_rate=algal_growth,
        algal_respiration_rate=xr.zeros_like(algal_growth),
    )

    carbon = Carbon(
        parameters={
            "AWc": AWc,
            "AWa": AWa,
            # Isolate algal photosynthesis sink.
            "kdoc_20": 0.0,
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

    # Algal growth is a DIC *sink*; expected dDIC/dt is negative.
    # Phase 9.E: mg-C/L/d (no ``/ 12000``).
    expected_rate = -rca * algal_growth
    np.testing.assert_allclose(
        np.asarray(v3_dic_rate),
        np.asarray(expected_rate),
        rtol=1e-6,
    )

    # Negative regression: confirm v3 is NOT using the raw AWc (which
    # would yield 1000x larger magnitude).
    raw_rate = -AWc * algal_growth
    assert not np.allclose(
        np.asarray(v3_dic_rate), np.asarray(raw_rate), rtol=1e-3
    )


def test_dic_includes_cbod_oxidation(
    water_temp_5cell, depth_5cell, poc_5cell, doc_5cell, dic_5cell, dox_5cell
):
    """Audit C3 (Phase 9.B): with CBOD wired and ``cbod_oxidation_rate``
    non-zero, v3 dDIC/dt includes the source ``cbod_oxidation_rate /
    roc`` (Fortran ``modCarbon.f90:262-269``, Phase 9.E mg-C/L/d
    convention without the legacy ``/ 12000`` mol-C/L/d conversion).

    Prior v3 omitted the CBOD->DIC source entirely; Phase 9.B added it
    via the ``self.cbod_process.cbod_oxidation_rate`` cache. Phase 9.E
    further removed the legacy ``/ 12000`` to align units.
    """
    roc = 32.0 / 12.0
    cbod_ox_rate = xr.DataArray(
        np.array([0.10, 0.20, 0.30, 0.40, 0.50]), dims="cell"
    )  # mg-O2/L/d
    mock_cbod = _MockCBOD(cbod_oxidation_rate=cbod_ox_rate)

    carbon = Carbon(
        parameters={
            "roc": roc,
            # Isolate CBOD source.
            "kdoc_20": 0.0,
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
    carbon.use_cbod = True
    carbon.cbod_process = mock_cbod

    registry = _build_registry(
        poc_5cell, doc_5cell, dic_5cell,
        water_temp_5cell, depth_5cell, dox_5cell,
    )
    dic_initial = registry.get_at_time("dic", datetime(2026, 1, 1)).copy()
    carbon.run(datetime(2026, 1, 1), registry)
    dic_final = registry.get_at_time("dic", datetime(2026, 1, 1))

    dt_days = timedelta(minutes=5).total_seconds() / 86400.0
    v3_dic_rate = (dic_final - dic_initial) / dt_days

    # Phase 9.E expected source: cbod_ox / roc (mg-C/L/d, no /12000).
    expected_rate = cbod_ox_rate / roc

    np.testing.assert_allclose(
        np.asarray(v3_dic_rate),
        np.asarray(expected_rate),
        rtol=1e-6,
    )

    # Non-zero magnitude check (regression: prior v3 was 0).
    assert np.all(np.asarray(v3_dic_rate) > 0)


def test_poc_hydrolysis_no_longer_dox_attenuated(
    water_temp_5cell, depth_5cell, poc_5cell, doc_5cell, dic_5cell
):
    """Audit C4 (Phase 9.B): POC hydrolysis is independent of DOX.

    Run Carbon twice with different DOX values; with the C4 fix the
    cached ``poc_hydrolysis_rate`` is identical across both runs because
    the DOX-Monod factor was removed. With the prior buggy form, low
    DOX would attenuate the rate to 50% at DOX = KsOxmc.
    """
    kpoc_20 = 0.005
    kpoc_theta = 1.047
    KsOxmc = 1.0  # so DOX/(DOX+KsOxmc) varies meaningfully across runs.

    def _run_with_dox(dox_value):
        dox_array = xr.DataArray(
            np.full(5, dox_value), dims="cell"
        )
        carbon = Carbon(
            parameters={
                "kpoc_20": kpoc_20,
                "kpoc_theta": kpoc_theta,
                "KsOxmc": KsOxmc,
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
            water_temp_5cell, depth_5cell, dox_array,
        )
        carbon.run(datetime(2026, 1, 1), registry)
        return carbon.poc_hydrolysis_rate

    rate_low_dox = _run_with_dox(0.5)   # Monod factor would be 1/3.
    rate_high_dox = _run_with_dox(20.0)  # Monod factor ~ 1.

    np.testing.assert_allclose(
        np.asarray(rate_low_dox),
        np.asarray(rate_high_dox),
        rtol=1e-12,
        err_msg=(
            "POC hydrolysis rate must be independent of DOX (Phase 9.B "
            "audit C4 fix); v3 should match Fortran modCarbon.f90:170 "
            "and v1 processes.py:2455 form ``kpoc_tc * POC``."
        ),
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
