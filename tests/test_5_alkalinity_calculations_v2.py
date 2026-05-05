"""Parity tests: v3 Alkalinity sub-rate methods vs v1 nsm1.processes helpers.

Phase 7.B of the v3 NSM1 implementation plan. Mirrors the established
Phase 2 pattern at ``tests/test_5_floating_algae_calculations_v2.py``.

Each test constructs a v3 Alkalinity instance, drives it through one
``Process.run`` substep against an in-memory registry, and compares
each cached sub-rate (``alk_nitrification_rate``,
``alk_denitrification_rate``, ``alk_algal_growth_rate``,
``alk_algal_respiration_rate``) to the equivalent v1 helper-function
output computed with the same inputs.

Scope: v1-equivalent ``dAlkdt`` sub-terms:

* Nitrification Alk sink (``r_alkn * nitrification_flux_rate * 50000``)
* Denitrification Alk source (``r_alkden * denitrification_flux_rate * 50000``)
* Algal-growth Alk coupling
  (``(r_alkaa * fNH4 - r_alkan * (1 - fNH4)) * ApGrowth * AWc * 50000``)
* Algal-respiration Alk source (``r_alkaa * ApRespiration * AWc * 50000``)

Known v1 vs v3 difference:

The v1 ``Alk_nitrification`` / ``Alk_denitrification`` functions apply
the DOX-Monod / oxygen-inhibition factor locally inside the function;
the v3 Alkalinity Process consumes Nitrogen's pre-cached
``nitrification_flux_rate`` / ``denitrification_flux_rate`` (Phase 2.B
Item 1), where the same Monod/inhibition factor is already baked in
upstream by Nitrogen.run. The parity test handles this by passing a
Nitrogen mock whose ``*_flux_rate`` already includes the Monod factor,
matching what v1 computes locally.

Synthetic mesh: 5-cell xarray DataArrays, single time step.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules.nsm1 import processes as v1
from clearwater_modules_v3.processes.alkalinity import Alkalinity

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
def alk_5cell() -> xr.DataArray:
    return xr.DataArray(
        np.array([100.0, 105.0, 110.0, 115.0, 120.0]), dims="cell"
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


@pytest.fixture(scope="function")
def no3_5cell() -> xr.DataArray:
    return xr.DataArray(
        np.array([1.0, 2.0, 3.0, 4.0, 5.0]), dims="cell"
    )


def _build_registry(
    alk: xr.DataArray,
    water_temp: xr.DataArray,
    depth: xr.DataArray,
) -> InMemoryRegistry:
    """Wire the state variables that Alkalinity.run reads."""
    registry = InMemoryRegistry()
    registry.register("alkalinity", alk.copy())
    registry.register("water_temperature", water_temp.copy())
    registry.register("depth", depth.copy())
    return registry


@dataclass
class _MockNitrogen:
    """Stand-in Nitrogen process exposing pre-cached flux rates.

    Alkalinity reads ``nitrification_flux_rate`` and
    ``denitrification_flux_rate`` (Phase 2.B Item 1 contract;
    DOX-Monod / O2-inhibition factor already baked in upstream).
    """
    nitrification_flux_rate: xr.DataArray
    denitrification_flux_rate: xr.DataArray


@dataclass
class _MockFloatingAlgae:
    """Stand-in FloatingAlgae for the Alkalinity algal-coupling parity tests."""
    algal_growth_rate: xr.DataArray
    algal_respiration_rate: xr.DataArray
    algal_nh4_uptake_fraction: xr.DataArray | float


# ---------------------------------------------------------------------------
# Parity tests
# ---------------------------------------------------------------------------


def test_alk_nitrification_sink_matches_v1(
    water_temp_5cell, depth_5cell, alk_5cell, dox_5cell, nh4_5cell
):
    """v3 ``alk_nitrification_rate`` cache after run == v1 ``Alk_nitrification``.

    v1 formula (line 3284-3319, when use_NH4 and use_DOX):
        r_alkn * (1 - exp(-KNR * DOX)) * knit_tc * NH4 * 50000

    v3 formula (alkalinity.py line 251-274):
        r_alkn * nitrification_flux_rate * 50000

    KNOWN DIFFERENCE: v1 applies the (1 - exp(-KNR*DOX)) DOX-Monod
    factor locally inside ``Alk_nitrification``. v3 has it baked into
    Nitrogen's ``nitrification_flux_rate`` cache (Phase 2.B Item 1).
    With matched DOX, KNR, knit_tc, NH4, the two should agree to
    floating-point precision -- the parity test wires a Nitrogen mock
    whose ``nitrification_flux_rate == (1 - exp(-KNR*DOX)) * knit_tc *
    NH4`` (matching the v1 inline form).
    """
    KNR = 0.6
    knit_20 = 0.1
    knit_theta = 1.083

    # Pre-compute the v1-equivalent flux that Nitrogen.run would have
    # cached: (1 - exp(-KNR*DOX)) * knit_tc * NH4. This is the form v1
    # ``Alk_nitrification`` evaluates inline.
    v1_knit_tc = v1.arrhenius_correction(water_temp_5cell, knit_20, knit_theta)
    nitrification_flux = (
        (1.0 - np.exp(-KNR * dox_5cell)) * v1_knit_tc * nh4_5cell
    )
    mock_nitrogen = _MockNitrogen(
        nitrification_flux_rate=nitrification_flux,
        denitrification_flux_rate=xr.zeros_like(nitrification_flux),
    )

    alk_proc = Alkalinity(time_step=timedelta(minutes=5))
    alk_proc.use_nitrogen = True
    alk_proc.nitrogen_process = mock_nitrogen
    # use_NH4 / use_NO3 default True from GLOBAL_PARAM_DEFAULTS.

    registry = _build_registry(alk_5cell, water_temp_5cell, depth_5cell)
    alk_proc.run(datetime(2026, 1, 1), registry)

    # v1 reference: full inline form.
    v1_alk_nitr = v1.Alk_nitrification(
        DOX=dox_5cell,
        NH4=nh4_5cell,
        knit_tc=v1_knit_tc,
        KNR=KNR,
        r_alkn=alk_proc.r_alkn,
        use_NH4=True,
        use_DOX=True,
    )

    np.testing.assert_allclose(
        np.asarray(alk_proc.alk_nitrification_rate),
        np.asarray(v1_alk_nitr),
        rtol=1e-6,
    )


def test_alk_denitrification_source_matches_v1(
    water_temp_5cell, depth_5cell, alk_5cell, dox_5cell, no3_5cell
):
    """v3 ``alk_denitrification_rate`` cache after run == v1 ``Alk_denitrification``.

    v1 formula (line 3246-3281, when use_NO3 and use_DOX):
        r_alkden * (1 - DOX/(DOX+KsOxdn)) * kdnit_tc * NO3 * 50000

    v3 formula (alkalinity.py line 276-299):
        r_alkden * denitrification_flux_rate * 50000

    Same DOX-Monod inhibition convention as nitrification: v1 applies
    the factor inline; v3 has it baked into Nitrogen's
    ``denitrification_flux_rate`` upstream. Parity is achieved by
    wiring a Nitrogen mock pre-loaded with the v1-form flux.
    """
    KsOxdn = 0.1
    kdnit_20 = 0.002
    kdnit_theta = 1.045

    v1_kdnit_tc = v1.arrhenius_correction(water_temp_5cell, kdnit_20, kdnit_theta)
    denitrification_flux = (
        (1.0 - dox_5cell / (dox_5cell + KsOxdn)) * v1_kdnit_tc * no3_5cell
    )
    mock_nitrogen = _MockNitrogen(
        nitrification_flux_rate=xr.zeros_like(denitrification_flux),
        denitrification_flux_rate=denitrification_flux,
    )

    alk_proc = Alkalinity(time_step=timedelta(minutes=5))
    alk_proc.use_nitrogen = True
    alk_proc.nitrogen_process = mock_nitrogen

    registry = _build_registry(alk_5cell, water_temp_5cell, depth_5cell)
    alk_proc.run(datetime(2026, 1, 1), registry)

    # v1 reference: full inline form.
    v1_alk_denit = v1.Alk_denitrification(
        DOX=dox_5cell,
        NO3=no3_5cell,
        kdnit_tc=v1_kdnit_tc,
        KsOxdn=KsOxdn,
        r_alkden=alk_proc.r_alkden,
        use_NO3=True,
        use_DOX=True,
    )

    np.testing.assert_allclose(
        np.asarray(alk_proc.alk_denitrification_rate),
        np.asarray(v1_alk_denit),
        rtol=1e-6,
    )


def test_alk_algal_growth_term_matches_v1(
    water_temp_5cell, depth_5cell, alk_5cell
):
    """v3 ``alk_algal_growth_rate`` cache after run == v1 ``Alk_algal_growth``.

    v1 formula (line 3322-3342):
        (r_alkaa * fNH4 - r_alkan * (1 - fNH4)) * ApGrowth * rca * 50000

    v3 formula (alkalinity.py line 301-343): identical (with rca == AWc).

    The sign convention in both v1 and v3 is "Alk SINK contribution":
    the term is subtracted from dAlk/dt. With all-NH4 uptake (fNH4==1)
    the term is positive (NH4 uptake consumes alk -> net sink). With
    all-NO3 uptake (fNH4==0) the term is negative (NO3 uptake produces
    alk -> net source).
    """
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

    alk_proc = Alkalinity(time_step=timedelta(minutes=5))
    alk_proc.use_floating_algae = True
    alk_proc.floating_algae_process = mock_algae

    registry = _build_registry(alk_5cell, water_temp_5cell, depth_5cell)
    alk_proc.run(datetime(2026, 1, 1), registry)

    # v1 reference: Alk_algal_growth.
    v1_alk_growth = v1.Alk_algal_growth(
        ApGrowth=algal_growth,
        r_alkaa=alk_proc.r_alkaa,
        r_alkan=alk_proc.r_alkan,
        ApUptakeFr_NH4=nh4_uptake_fr,
        rca=alk_proc.AWc,
        use_Algae=True,
    )

    np.testing.assert_allclose(
        np.asarray(alk_proc.alk_algal_growth_rate),
        np.asarray(v1_alk_growth),
        rtol=1e-6,
    )


def test_alk_algal_respiration_source_matches_v1(
    water_temp_5cell, depth_5cell, alk_5cell
):
    """v3 ``alk_algal_respiration_rate`` cache after run == v1 ``Alk_algal_respiration``.

    v1 formula (line 3345-3361):
        ApRespiration * r_alkaa * rca * 50000

    v3 formula (alkalinity.py line 345-360): identical (with rca == AWc).

    Algal respiration always produces alk (DIC release pathway), so the
    sign convention is "+ source" in both v1 and v3.
    """
    algal_respiration = xr.DataArray(
        np.array([0.5, 0.6, 0.7, 0.8, 1.0]), dims="cell"
    )
    mock_algae = _MockFloatingAlgae(
        algal_growth_rate=xr.zeros_like(algal_respiration),
        algal_respiration_rate=algal_respiration,
        algal_nh4_uptake_fraction=0.5,    # arbitrary; respiration term doesn't use it
    )

    alk_proc = Alkalinity(time_step=timedelta(minutes=5))
    alk_proc.use_floating_algae = True
    alk_proc.floating_algae_process = mock_algae

    registry = _build_registry(alk_5cell, water_temp_5cell, depth_5cell)
    alk_proc.run(datetime(2026, 1, 1), registry)

    # v1 reference: Alk_algal_respiration.
    v1_alk_resp = v1.Alk_algal_respiration(
        ApRespiration=algal_respiration,
        r_alkaa=alk_proc.r_alkaa,
        rca=alk_proc.AWc,
        use_Algae=True,
    )

    np.testing.assert_allclose(
        np.asarray(alk_proc.alk_algal_respiration_rate),
        np.asarray(v1_alk_resp),
        rtol=1e-6,
    )
