"""Parity tests: v3 DOX sub-rate methods vs v1 nsm1.processes helpers.

Phase 7.B of the v3 NSM1 implementation plan. Mirrors the established
Phase 2 pattern at ``tests/test_5_floating_algae_calculations_v2.py``.

Each test constructs a v3 DOX instance, drives it through one
``Process.run`` substep against an in-memory registry, and compares
either a cached step-scoped rate (``dox_sat``, ``atm_reaeration_rate``,
``dox_nitrification_rate``, ``dox_sod_rate``) or a back-calculated
``dDOX/dt`` sub-term to the equivalent v1 helper-function output
computed with the same inputs.

Scope: v1-equivalent ``dDOXdt`` sub-terms:

* O2 saturation (full APHA / Benson-Krause + pressure + water-vapor)
* Atmospheric reaeration (``ka_tc * (O2sat - DOX)``)
* Nitrification O2 sink (``ron * nitrification_flux_rate``; v3 reads
  this from a Nitrogen sibling Process per Integration Item 1)
* SOD O2 sink (``SOD_tc / depth``; v3 uses pure Arrhenius SOD_tc with
  no DOX-Monod attenuation)
* Floating-algae photosynthesis O2 source (Redfield 138/106 - 32/106 *
  fNH4 stoichiometric factor)

Synthetic mesh: 5-cell xarray DataArrays, single time step.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules.nsm1 import processes as v1
from clearwater_modules_v3.processes.dox import DOX, dox_sat_apha

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
def dox_5cell() -> xr.DataArray:
    return xr.DataArray(
        np.array([4.0, 6.0, 7.0, 8.0, 10.0]), dims="cell"
    )


@pytest.fixture(scope="function")
def nh4_5cell() -> xr.DataArray:
    return xr.DataArray(
        np.array([0.05, 0.10, 0.15, 0.20, 0.30]), dims="cell"
    )


def _build_registry(
    dox: xr.DataArray,
    water_temp: xr.DataArray,
    depth: xr.DataArray,
    ammonium: xr.DataArray | None = None,
) -> InMemoryRegistry:
    """Wire the state variables that DOX.run reads."""
    registry = InMemoryRegistry()
    registry.register("oxygen_dissolved", dox.copy())
    registry.register("water_temperature", water_temp.copy())
    registry.register("depth", depth.copy())
    if ammonium is not None:
        registry.register("ammonium", ammonium.copy())
    return registry


@dataclass
class _MockNitrogen:
    """Stand-in Nitrogen for the DOX nitrification-coupling parity test.

    DOX._nitrification_flux reads ``nitrification_flux_rate`` (mg-N/L/d,
    positive magnitude) from the Nitrogen Process per Integration Item 1
    (registry rate-variable convention, spec resolved Q10).
    """
    nitrification_flux_rate: xr.DataArray
    denitrification_flux_rate: float = 0.0


@dataclass
class _MockFloatingAlgae:
    """Stand-in FloatingAlgae for the DOX algal-photosynthesis parity test."""
    algal_growth_rate: xr.DataArray
    algal_respiration_rate: xr.DataArray
    algal_nh4_uptake_fraction: xr.DataArray | float


# ---------------------------------------------------------------------------
# Parity tests
# ---------------------------------------------------------------------------


def test_o2sat_full_apha_matches_v1(water_temp_5cell):
    """v3 ``dox_sat_apha`` (full Benson-Krause + pressure + water-vapor)
    == v1 ``DOX_sat`` (line 2901-2923).

    Both implement the APHA formulation:
        DOX_sat_uncorrected = Benson-Krause polynomial in 1/T_K
        DOX_sat = DOX_sat_uncorrected * P_atm
                  * (1 - p_wv/P_atm) * (1 - alpha * P_atm)
                  / ((1 - p_wv) * (1 - alpha))
    """
    pressure_mb = 1013.25       # sea-level standard

    v3_value = dox_sat_apha(water_temp_5cell, pressure_mb)

    # v1 takes TwaterK + p_wv + alpha precomputed.
    t_water_k = water_temp_5cell + 273.15
    v1_pwv = v1.pwv(t_water_k)
    v1_alpha = v1.DOs_atm_alpha(water_temp_5cell)
    v1_value = v1.DOX_sat(
        TwaterK=t_water_k,
        pressure_mb=pressure_mb,
        pwv=v1_pwv,
        DOs_atm_alpha=v1_alpha,
    )

    np.testing.assert_allclose(
        np.asarray(v3_value), np.asarray(v1_value), rtol=1e-6
    )


def test_atmospheric_reaeration_matches_v1(
    water_temp_5cell, depth_5cell, dox_5cell
):
    """v3 ``atm_reaeration_rate`` cache after run == v1 ``Atm_O2_reaeration``.

    v1 formula: ``ka_tc * (DOX_sat - DOX)`` (line 2927-2939).
    v3 formula: same (dox.py ``_atm_reaeration_flux`` line 391-398).
    """
    kah_20_user = 1.0    # 1/d; user-defined hydraulic reaeration
    pressure_mb = 1013.25

    dox_proc = DOX(
        parameters={
            "kah_20_user": kah_20_user,
            "kah_theta": 1.024,
            "kaw_20_user": 0.0,         # disable wind contribution
            "kaw_theta": 1.024,
            "hydraulic_reaeration_option": 1,
            "wind_reaeration_option": 1,
            "pressure_mb": pressure_mb,
            # Disable every other DOX source/sink.
            "SOD_20": 0.0,
        },
        time_step=timedelta(minutes=5),
    )
    # Standalone-mode: leave coupling flags False so algae/Nitrogen/Carbon/
    # CBOD contributions are zero, and run() exercises only atm reaeration
    # + SOD (SOD_20=0 zeroes that as well).
    registry = _build_registry(dox_5cell, water_temp_5cell, depth_5cell)
    dox_proc.run(datetime(2026, 1, 1), registry)

    # v1 reference: ka_tc with wind=0 simplifies to arrhenius_correction(
    # T, kah_20, kah_theta).
    v1_ka_tc = v1.arrhenius_correction(water_temp_5cell, kah_20_user, 1.024)
    t_water_k = water_temp_5cell + 273.15
    v1_pwv = v1.pwv(t_water_k)
    v1_alpha = v1.DOs_atm_alpha(water_temp_5cell)
    v1_dox_sat = v1.DOX_sat(t_water_k, pressure_mb, v1_pwv, v1_alpha)
    v1_atm_reaer = v1.Atm_O2_reaeration(
        ka_tc=v1_ka_tc, DOX_sat=v1_dox_sat, DOX=dox_5cell
    )

    np.testing.assert_allclose(
        np.asarray(dox_proc.atm_reaeration_rate),
        np.asarray(v1_atm_reaer),
        rtol=1e-6,
    )


def test_dox_nitrification_sink_matches_v1(
    water_temp_5cell, depth_5cell, dox_5cell, nh4_5cell
):
    """v3 ``dox_nitrification_rate`` cache after run == v1 net nitrification O2 demand.

    v1 ``DOX_Nitrification`` (line 2980-2999):
        (1 - exp(-KNR * DOX)) * ron * knit_tc * NH4

    v3 reads from a Nitrogen sibling Process: ``ron * nitrification_flux_rate``
    where ``nitrification_flux_rate`` is the upstream Nitrogen.run-cached
    ``ammonium * knit_tc * (1 - exp(-KNR * DOX))`` (per Integration Item 1).

    With matched parameters and a mock Nitrogen pre-loaded with the same
    ``(1 - exp(-KNR*DOX)) * knit_tc * NH4`` flux, the two should agree
    to floating-point precision.
    """
    KNR = 0.6           # 1/(mg-O2/L)
    knit_20 = 0.1       # 1/d
    knit_theta = 1.083
    ron = 2.0 * 32.0 / 14.0     # mg-O2/mg-N

    # Pre-compute the v1 nitrification flux at our test point and stuff
    # it onto the mock Nitrogen. v3 DOX multiplies by ``ron`` only.
    v1_knit_tc = v1.arrhenius_correction(water_temp_5cell, knit_20, knit_theta)
    nitrification_flux = (
        (1.0 - np.exp(-KNR * dox_5cell)) * v1_knit_tc * nh4_5cell
    )
    mock_nitrogen = _MockNitrogen(
        nitrification_flux_rate=nitrification_flux,
        denitrification_flux_rate=0.0,
    )

    dox_proc = DOX(
        parameters={
            "ron": ron,
            "KNR": KNR,
            "knit_20": knit_20,
            "knit_theta": knit_theta,
            # Disable every other DOX source/sink.
            "kah_20_user": 0.0,
            "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1,
            "wind_reaeration_option": 1,
            "SOD_20": 0.0,
        },
        time_step=timedelta(minutes=5),
    )
    # Mimic init_process: enable Nitrogen coupling.
    dox_proc.use_nitrogen = True
    dox_proc.nitrogen_process = mock_nitrogen
    # use_NH4 is True by default from GLOBAL_PARAM_DEFAULTS.
    dox_proc.use_NH4 = True

    registry = _build_registry(
        dox_5cell, water_temp_5cell, depth_5cell, nh4_5cell
    )
    dox_proc.run(datetime(2026, 1, 1), registry)

    # v1 reference: DOX_Nitrification computes the full inline form,
    # which numerically equals ron * nitrification_flux when KNR/knit/NH4
    # are the same.
    v1_dox_nitr = v1.DOX_Nitrification(
        KNR=KNR,
        DOX=dox_5cell,
        ron=ron,
        knit_tc=v1_knit_tc,
        NH4=nh4_5cell,
        use_NH4=True,
    )

    np.testing.assert_allclose(
        np.asarray(dox_proc.dox_nitrification_rate),
        np.asarray(v1_dox_nitr),
        rtol=1e-6,
    )


def test_dox_sod_sink_matches_v1_with_attenuation(
    water_temp_5cell, depth_5cell, dox_5cell
):
    """v3 ``dox_sod_rate`` cache after run == v1 ``DOX_SOD`` *with*
    DOX-Monod attenuation.

    Phase 9.B audit fix C2: v3 now applies the Fortran
    ``modGlobalParam.f90:254`` DOX-Monod attenuation
    (``SOD *= DOX / (DOX + KsSOD)``) inside ``_sod_flux`` (the
    primitive ``utils.sediment.SOD_tc`` remains pure Arrhenius for
    architectural reasons). The reference v1 call is now made with
    ``use_DOX=True`` so v1 applies the same attenuation; v3 should
    agree.
    """
    SOD_20 = 1.0
    SOD_theta = 1.060
    KsSOD = 1.0

    dox_proc = DOX(
        parameters={
            "SOD_20": SOD_20,
            "SOD_theta": SOD_theta,
            "KsSOD": KsSOD,
            "use_DOX": True,
            # Disable every other DOX source/sink.
            "kah_20_user": 0.0,
            "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1,
            "wind_reaeration_option": 1,
        },
        time_step=timedelta(minutes=5),
    )
    registry = _build_registry(dox_5cell, water_temp_5cell, depth_5cell)
    dox_proc.run(datetime(2026, 1, 1), registry)

    # v1 reference: SOD_tc with use_DOX=True applies the same attenuation.
    v1_sod_tc = v1.SOD_tc(
        SOD_20=SOD_20,
        TwaterC=water_temp_5cell,
        SOD_theta=SOD_theta,
        DOX=dox_5cell,
        KsSOD=KsSOD,
        use_DOX=True,
    )
    v1_sod = v1.DOX_SOD(depth=depth_5cell, SOD_tc=v1_sod_tc)

    np.testing.assert_allclose(
        np.asarray(dox_proc.dox_sod_rate),
        np.asarray(v1_sod),
        rtol=1e-6,
    )


def test_dox_algal_photosynthesis_source_matches_fortran_anchored(
    water_temp_5cell, depth_5cell, dox_5cell, nh4_5cell
):
    """v3 floating-algae photosynthesis O2 source == Fortran-anchored
    expected value with ``rca = AWc / AWa``.

    Fortran formula (``modDOX.f90:135``):
        ``O2_ApGrowth = (138/106 - 32/106 * fNH4) * roc * rca * ApGrowth``
    where ``rca = AWc / AWa`` (mg-C/ug-Chla). v1 derives the same.

    Phase 9.B audit C1: prior v3 used ``self.AWc = 40`` (raw weight)
    instead of the derived ratio (0.04 mg-C/ug-Chla), inflating the
    photosynthesis O2 source by 1000x. The previous parity test masked
    this by passing ``rca = AWc = 40`` to the v1 reference (same wrong
    value on both sides). The fix derives ``rca = AWc / AWa`` and this
    test asserts against the Fortran-anchored magnitude.
    """
    AWc = 40.0      # mg-C raw weight
    AWa = 1000.0    # ug-Chla algal unit
    rca = AWc / AWa  # 0.04 mg-C / ug-Chla (Fortran modAlgae.f90)
    roc = 32.0 / 12.0

    # Synthetic floating-algae rates.
    algal_growth = xr.DataArray(
        np.array([0.5, 0.6, 0.7, 0.8, 1.0]), dims="cell"
    )
    nh4_uptake_fr = xr.DataArray(
        np.array([0.2, 0.4, 0.5, 0.6, 0.8]), dims="cell"
    )
    mock_algae = _MockFloatingAlgae(
        algal_growth_rate=algal_growth,
        algal_respiration_rate=xr.zeros_like(algal_growth),
        algal_nh4_uptake_fraction=nh4_uptake_fr,
    )

    dox_proc = DOX(
        parameters={
            "AWc": AWc,
            "AWa": AWa,
            "roc": roc,
            # Disable every other DOX source/sink.
            "kah_20_user": 0.0,
            "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1,
            "wind_reaeration_option": 1,
            "SOD_20": 0.0,
        },
        time_step=timedelta(minutes=5),
    )
    # Wire the mock FloatingAlgae directly.
    dox_proc.use_floating_algae = True
    dox_proc.use_Algae = True
    dox_proc.floating_algae_process = mock_algae

    registry = _build_registry(
        dox_5cell, water_temp_5cell, depth_5cell, nh4_5cell
    )
    dox_initial = registry.get_at_time(
        "oxygen_dissolved", datetime(2026, 1, 1)
    ).copy()
    dox_proc.run(datetime(2026, 1, 1), registry)
    dox_final = registry.get_at_time(
        "oxygen_dissolved", datetime(2026, 1, 1)
    )

    dt_days = timedelta(minutes=5).total_seconds() / 86400.0
    v3_dox_rate = (dox_final - dox_initial) / dt_days

    # Fortran-anchored expected value.
    expected_rate = (
        algal_growth
        * rca
        * roc
        * (138.0 / 106.0 - 32.0 / 106.0 * nh4_uptake_fr)
    )

    np.testing.assert_allclose(
        np.asarray(v3_dox_rate),
        np.asarray(expected_rate),
        rtol=1e-6,
    )

    # Confirm v1 helper agrees when fed the *correct* rca.
    v1_dox_apg = v1.DOX_ApGrowth(
        ApGrowth=algal_growth,
        rca=rca,
        roc=roc,
        ApUptakeFr_NH4=nh4_uptake_fr,
        use_Algae=True,
    )
    np.testing.assert_allclose(
        np.asarray(v3_dox_rate),
        np.asarray(v1_dox_apg),
        rtol=1e-6,
    )


# ---------------------------------------------------------------------------
# Audit-anchored regression tests (Phase 9.B)
# ---------------------------------------------------------------------------


def test_dox_algal_photosynthesis_uses_correct_rca(
    water_temp_5cell, depth_5cell, dox_5cell, nh4_5cell
):
    """Audit C1 (Phase 9.B): default-instantiated DOX with default
    FloatingAlgae stoichiometry produces an O2 source ~1000x smaller
    than the prior raw-AWc magnitude.

    Setup: ``ApGrowth = 0.5 ug-Chla/L/d``, all-NH4 uptake (factor = 1).
    Fortran-anchored expected = ``0.5 * 0.04 * (32/12) * 1.0 = 0.0533
    mg-O2/L/d``. The prior raw-AWc form would give ~53.3 mg-O2/L/d.
    """
    rca = 40.0 / 1000.0  # 0.04
    roc = 32.0 / 12.0
    ap_growth = xr.DataArray(np.full(5, 0.5), dims="cell")
    nh4_fr = xr.DataArray(np.full(5, 1.0), dims="cell")  # factor = 1.0
    mock_algae = _MockFloatingAlgae(
        algal_growth_rate=ap_growth,
        algal_respiration_rate=xr.zeros_like(ap_growth),
        algal_nh4_uptake_fraction=nh4_fr,
    )

    dox_proc = DOX(
        parameters={
            "kah_20_user": 0.0,
            "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1,
            "wind_reaeration_option": 1,
            "SOD_20": 0.0,
        },
        time_step=timedelta(minutes=5),
    )
    dox_proc.use_floating_algae = True
    dox_proc.use_Algae = True
    dox_proc.floating_algae_process = mock_algae

    registry = _build_registry(dox_5cell, water_temp_5cell, depth_5cell, nh4_5cell)
    dox_initial = registry.get_at_time(
        "oxygen_dissolved", datetime(2026, 1, 1)
    ).copy()
    dox_proc.run(datetime(2026, 1, 1), registry)
    dox_final = registry.get_at_time("oxygen_dissolved", datetime(2026, 1, 1))

    dt_days = timedelta(minutes=5).total_seconds() / 86400.0
    rate = (dox_final - dox_initial) / dt_days

    expected_rate_corrected = 0.5 * rca * roc * 1.0  # ~0.0533
    np.testing.assert_allclose(
        np.asarray(rate),
        np.full(5, expected_rate_corrected),
        rtol=1e-6,
    )

    # Negative regression: v3 must NOT produce the raw-AWc magnitude.
    raw_AWc_rate = 0.5 * 40.0 * roc * 1.0
    assert not np.allclose(np.asarray(rate), np.full(5, raw_AWc_rate), rtol=1e-3)


def test_dox_sod_attenuates_at_low_dox(
    water_temp_5cell, depth_5cell
):
    """Audit C2 (Phase 9.B): SOD effective rate equals half the unattenuated
    Arrhenius value when ``DOX = KsSOD``.

    Fortran (``modGlobalParam.f90:254``): ``SOD_tc *= DOX / (DOX + KsSOD)``
    when ``use_DOX``. At ``DOX = KsSOD``, the factor is 0.5.
    """
    SOD_20 = 1.0
    SOD_theta = 1.060
    KsSOD = 1.0

    # DOX exactly at KsSOD -> Monod factor = 0.5.
    dox_at_ksod = xr.DataArray(np.full(5, KsSOD), dims="cell")

    dox_proc = DOX(
        parameters={
            "SOD_20": SOD_20,
            "SOD_theta": SOD_theta,
            "KsSOD": KsSOD,
            "use_DOX": True,
            "kah_20_user": 0.0,
            "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1,
            "wind_reaeration_option": 1,
        },
        time_step=timedelta(minutes=5),
    )
    registry = _build_registry(dox_at_ksod, water_temp_5cell, depth_5cell)
    dox_proc.run(datetime(2026, 1, 1), registry)
    attenuated_sod = np.asarray(dox_proc.dox_sod_rate)

    # Compute the unattenuated rate by disabling use_DOX.
    dox_proc_unattenuated = DOX(
        parameters={
            "SOD_20": SOD_20,
            "SOD_theta": SOD_theta,
            "KsSOD": KsSOD,
            "use_DOX": False,
            "kah_20_user": 0.0,
            "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1,
            "wind_reaeration_option": 1,
        },
        time_step=timedelta(minutes=5),
    )
    registry2 = _build_registry(
        dox_at_ksod.copy(), water_temp_5cell, depth_5cell
    )
    dox_proc_unattenuated.run(datetime(2026, 1, 1), registry2)
    unattenuated_sod = np.asarray(dox_proc_unattenuated.dox_sod_rate)

    # At DOX = KsSOD the attenuated rate is exactly half.
    np.testing.assert_allclose(
        attenuated_sod, 0.5 * unattenuated_sod, rtol=1e-12
    )

    # And under hypoxia (DOX -> 0) the SOD sink approaches zero.
    dox_zero = xr.DataArray(np.full(5, 1e-6), dims="cell")
    dox_proc_low = DOX(
        parameters={
            "SOD_20": SOD_20,
            "SOD_theta": SOD_theta,
            "KsSOD": KsSOD,
            "use_DOX": True,
            "kah_20_user": 0.0,
            "kaw_20_user": 0.0,
            "hydraulic_reaeration_option": 1,
            "wind_reaeration_option": 1,
        },
        time_step=timedelta(minutes=5),
    )
    registry3 = _build_registry(dox_zero, water_temp_5cell, depth_5cell)
    dox_proc_low.run(datetime(2026, 1, 1), registry3)
    sod_at_low_dox = np.asarray(dox_proc_low.dox_sod_rate)
    # Should be ~ DOX/KsSOD * unattenuated ~ 1e-6 * unattenuated.
    assert np.all(sod_at_low_dox < 1e-5 * unattenuated_sod)
