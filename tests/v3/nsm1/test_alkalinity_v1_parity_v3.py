"""v3 Alkalinity kinetic regression against frozen v1 reference values.

Migration history: originally ``tests/test_5_alkalinity_calculations_v2.py``
imported v1 directly to compute the reference inline. This file
freezes those v1 reference values as numpy literals so the test runs
without requiring ``src/clearwater_modules/`` (v1) source.

NSM1-CA-1 re-derivation (2026-05-16, gold-standard spec A1): the
``V1_ALGAL_GROWTH_REFERENCE`` / ``V1_ALGAL_RESP_REFERENCE`` literals
were originally captured by calling v1's ``Alk_algal_growth`` /
``Alk_algal_respiration`` with ``rca`` bound to the *raw* stoichiometric
weight ``AWc`` (=40) -- the same defective convention v3 carried in
``alkalinity.py``. Because both sides shared that wrong input, the
parity test passed while masking a 1000x error (the structural reason
NSM1-CA-1 escaped; see gold-standard spec Section 1(4)). The literals
below are re-derived by calling the **same v1 functions** with the
verified-correct *intensive* ratio ``rca = AWc/AWa`` (=0.04 mg-C/ug-Chla,
v1 ``processes.py:337-347``). They are now a faithful v1 capture and
this is a genuine non-shared-path v1<->v3 parity check.

Scope: v1-equivalent ``dAlkdt`` sub-terms (nitrification sink,
denitrification source, algal-growth coupling, algal-respiration source).
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.alkalinity import Alkalinity
from clearwater_modules_v3.utils.conversions import arrhenius_correction

from tests.v3.nsm1.conftest import InMemoryRegistry


# ---------------------------------------------------------------------------
# Frozen v1 reference values (captured 2026-05-10 from
# clearwater_modules.nsm1.processes against the 5-cell fixtures below)
# ---------------------------------------------------------------------------

V1_NITR_SINK_REFERENCE = np.array([
    0.021797079178702924,
    0.059235679563923394,
    0.10553618819780596,
    0.16617663145906025,
    0.3184620287023956,
])

# NSM1-SCI-N1 (gold-standard spec A2): these are the *upstream-defect*
# v1/Fortran values, captured with the stoichiometrically-wrong
# ``r_alkden = 4/14/1000`` that Fortran (``modAlkalinity.f90:54``), v1,
# and pre-fix v3 all shared (the canonical "wrong at all stages" case,
# invisible to v1<->v3 parity). v3 deliberately diverges: denitrification
# produces 1 eq alkalinity per mol NO3-N (CE-QUAL-W2
# ``water-quality.f90:3157``; Stumm & Morgan), so corrected v3 =
# this reference / 4. Retained verbatim as the auditable divergence
# baseline; ``test_alk_denitrification_source_matches_v1`` asserts the
# /4 corrected value.
V1_DENIT_SOURCE_REFERENCE = np.array([
    0.0005591993355405442,
    0.0008578266522129975,
    0.001207243460764586,
    0.00154077601410934,
    0.0017626335751812243,
])

# SCI-N1 intentional divergence factor: r_alkden 4/14/1000 -> 1/14/1000.
SCI_N1_DENIT_DIVERGENCE = 4.0

# Re-derived 2026-05-16 by calling v1 ``Alk_algal_growth`` with the
# correct intensive ``rca = AWc/AWa = 40/1000 = 0.04`` (was captured at
# raw ``AWc = 40``; exactly 1000x too large). NSM1-CA-1.
V1_ALGAL_GROWTH_REFERENCE = np.array([
    -0.00911949685534591,
    -0.00490566037735849,
    -0.002201257861635219,
    0.001509433962264152,
    0.011949685534591196,
])

# Re-derived 2026-05-16 from v1 ``Alk_algal_respiration`` at intensive
# ``rca = 0.04`` (was raw ``AWc = 40``; exactly 1000x too large).
V1_ALGAL_RESP_REFERENCE = np.array([
    0.0110062893081761,
    0.01320754716981132,
    0.01540880503144654,
    0.01761006289308176,
    0.0220125786163522,
])


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
    registry = InMemoryRegistry()
    registry.register("alkalinity", alk.copy())
    registry.register("water_temperature", water_temp.copy())
    registry.register("depth", depth.copy())
    return registry


@dataclass
class _MockNitrogen:
    nitrification_flux_rate: xr.DataArray
    denitrification_flux_rate: xr.DataArray


@dataclass
class _MockFloatingAlgae:
    algal_growth_rate: xr.DataArray
    algal_respiration_rate: xr.DataArray
    algal_nh4_uptake_fraction: xr.DataArray | float


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_alk_nitrification_sink_matches_v1(
    water_temp_5cell, depth_5cell, alk_5cell, dox_5cell, nh4_5cell
):
    """v3 ``alk_nitrification_rate`` matches frozen v1 reference."""
    KNR = 0.6
    knit_tc = arrhenius_correction(water_temp_5cell, 0.1, 1.083)
    nitrification_flux = (
        (1.0 - np.exp(-KNR * dox_5cell)) * knit_tc * nh4_5cell
    )
    mock_nitrogen = _MockNitrogen(
        nitrification_flux_rate=nitrification_flux,
        denitrification_flux_rate=xr.zeros_like(nitrification_flux),
    )

    alk_proc = Alkalinity(time_step=timedelta(minutes=5))
    alk_proc.use_nitrogen = True
    alk_proc.nitrogen_process = mock_nitrogen

    registry = _build_registry(alk_5cell, water_temp_5cell, depth_5cell)
    alk_proc.run(datetime(2026, 1, 1), registry)

    np.testing.assert_allclose(
        np.asarray(alk_proc.alk_nitrification_rate),
        V1_NITR_SINK_REFERENCE,
        rtol=1e-6,
    )


def test_alk_denitrification_source_matches_v1(
    water_temp_5cell, depth_5cell, alk_5cell, dox_5cell, no3_5cell
):
    """v3 ``alk_denitrification_rate`` is the SCI-N1-corrected value:
    exactly ``V1_DENIT_SOURCE_REFERENCE / 4``. v3 deliberately diverges
    from the v1/Fortran upstream defect (``r_alkden`` 4/14/1000 ->
    1/14/1000); see the reference block above and spec A2.
    """
    KsOxdn = 0.1
    kdnit_tc = arrhenius_correction(water_temp_5cell, 0.002, 1.045)
    denitrification_flux = (
        (1.0 - dox_5cell / (dox_5cell + KsOxdn)) * kdnit_tc * no3_5cell
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

    np.testing.assert_allclose(
        np.asarray(alk_proc.alk_denitrification_rate),
        V1_DENIT_SOURCE_REFERENCE / SCI_N1_DENIT_DIVERGENCE,
        rtol=1e-6,
    )


def test_alk_algal_growth_term_matches_v1(
    water_temp_5cell, depth_5cell, alk_5cell
):
    """v3 ``alk_algal_growth_rate`` matches frozen v1 reference."""
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

    np.testing.assert_allclose(
        np.asarray(alk_proc.alk_algal_growth_rate),
        V1_ALGAL_GROWTH_REFERENCE,
        rtol=1e-6,
    )


def test_alk_algal_respiration_source_matches_v1(
    water_temp_5cell, depth_5cell, alk_5cell
):
    """v3 ``alk_algal_respiration_rate`` matches frozen v1 reference."""
    algal_respiration = xr.DataArray(
        np.array([0.5, 0.6, 0.7, 0.8, 1.0]), dims="cell"
    )
    mock_algae = _MockFloatingAlgae(
        algal_growth_rate=xr.zeros_like(algal_respiration),
        algal_respiration_rate=algal_respiration,
        algal_nh4_uptake_fraction=0.5,
    )

    alk_proc = Alkalinity(time_step=timedelta(minutes=5))
    alk_proc.use_floating_algae = True
    alk_proc.floating_algae_process = mock_algae

    registry = _build_registry(alk_5cell, water_temp_5cell, depth_5cell)
    alk_proc.run(datetime(2026, 1, 1), registry)

    np.testing.assert_allclose(
        np.asarray(alk_proc.alk_algal_respiration_rate),
        V1_ALGAL_RESP_REFERENCE,
        rtol=1e-6,
    )
