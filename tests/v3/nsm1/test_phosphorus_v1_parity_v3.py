"""v3 Phosphorus kinetic regression against frozen v1 reference values.

Migration from ``tests/test_5_phosphorus_calculations_v2.py``.

Note: tests `test_tip_settling_matches_fortran_anchored` and
`test_fdp_unit_factor_dimensionally_correct` use Fortran-anchored
formulas (not v1) as the canonical reference, because Phase 9.B
identified a dimensional bug in v1's ``shared.processes.fdp``. The
Fortran-anchored formula is inlined below (no v1 import needed).
"""
from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.phosphorus import Phosphorus
from clearwater_modules_v3.utils.partitioning import fdp as v3_fdp

from tests.v3.nsm1.conftest import InMemoryRegistry


V1_ORGP_HYDROLYSIS_REFERENCE = np.array([
    0.003974079913422963,
    0.0054734088116408465,
    0.007000000000000001,
    0.008769672,
    0.011323375719750058,
])

V1_ORGP_SETTLING_REFERENCE = np.array([
    0.010000000000000002,
    0.006,
    0.004666666666666667,
    0.004,
    0.003,
])

V1_TIP_SETTLING_FORTRAN_FDP_REFERENCE = np.array([
    9.950248756218638e-05,
    5.970149253731183e-05,
    4.643449419568699e-05,
    3.980099502487455e-05,
    2.9850746268655914e-05,
])


def _fortran_anchored_fdp(use_TIP, Solid, kdpo4):
    """Dimensionally correct fdp matching Fortran ``modGlobalParam.f90:228``."""
    return xr.where(use_TIP, 1.0 / (1.0 + kdpo4 * Solid * 1.0e-6), 0.0)


def _v1_shared_fdp_buggy(use_TIP, Solid, kdpo4):
    """Inline reproduction of v1 ``shared.processes.fdp`` (the buggy form,
    inverted unit factor)."""
    return xr.where(use_TIP, 1.0 / (1.0 + kdpo4 * Solid / 0.000001), 0.0)


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
def tip_5cell() -> xr.DataArray:
    return xr.DataArray(
        np.array([0.10, 0.12, 0.14, 0.16, 0.18]), dims="cell"
    )


@pytest.fixture(scope="function")
def orgp_5cell() -> xr.DataArray:
    return xr.DataArray(
        np.array([0.05, 0.06, 0.07, 0.08, 0.09]), dims="cell"
    )


def _build_registry(tip, orgp, water_temp, depth) -> InMemoryRegistry:
    registry = InMemoryRegistry()
    registry.register("tip", tip.copy())
    registry.register("organic_phosphorus", orgp.copy())
    registry.register("water_temperature", water_temp.copy())
    registry.register("depth", depth.copy())
    return registry


def test_orgp_hydrolysis_matches_v1(
    water_temp_5cell, depth_5cell, tip_5cell, orgp_5cell
):
    phosphorus = Phosphorus(
        parameters={
            "kop_20": 0.1,
            "kop_theta": 1.047,
            "vs": 0.0,
            "vsop": 0.0,
            "rpo4_20": 0.0,
        },
        time_step=timedelta(minutes=5),
    )
    registry = _build_registry(
        tip_5cell, orgp_5cell, water_temp_5cell, depth_5cell
    )
    phosphorus.run(datetime(2026, 1, 1), registry)

    np.testing.assert_allclose(
        np.asarray(phosphorus.orgp_to_tip_hydrolysis_rate),
        V1_ORGP_HYDROLYSIS_REFERENCE,
        rtol=1e-6,
    )


def test_orgp_settling_matches_v1(
    water_temp_5cell, depth_5cell, tip_5cell, orgp_5cell
):
    phosphorus = Phosphorus(
        parameters={
            "vsop": 0.1,
            "kop_20": 0.0,
            "vs": 0.0,
            "rpo4_20": 0.0,
        },
        time_step=timedelta(minutes=5),
    )
    registry = _build_registry(
        tip_5cell, orgp_5cell, water_temp_5cell, depth_5cell
    )
    phosphorus.run(datetime(2026, 1, 1), registry)

    np.testing.assert_allclose(
        np.asarray(phosphorus.orgp_settling_rate),
        V1_ORGP_SETTLING_REFERENCE,
        rtol=1e-6,
    )


def test_tip_settling_matches_fortran_anchored(
    water_temp_5cell, depth_5cell, tip_5cell, orgp_5cell
):
    """v3 TIP settling matches Fortran-anchored ``TIP_Settling`` with the
    dimensionally correct ``fdp``."""
    kdpo4 = 1000.0
    solid = 5.0
    phosphorus = Phosphorus(
        parameters={
            "vs": 0.1,
            "kdpo4": kdpo4,
            "Solid": solid,
            "use_TIP": True,
            "kop_20": 0.0,
            "vsop": 0.0,
            "rpo4_20": 0.0,
        },
        time_step=timedelta(minutes=5),
    )
    registry = _build_registry(
        tip_5cell, orgp_5cell, water_temp_5cell, depth_5cell
    )
    phosphorus.run(datetime(2026, 1, 1), registry)

    np.testing.assert_allclose(
        np.asarray(phosphorus.tip_settling_rate),
        V1_TIP_SETTLING_FORTRAN_FDP_REFERENCE,
        rtol=1e-6,
    )


def test_tip_partitioning_fdp_matches_fortran():
    """v3 ``fdp`` matches Fortran ``modGlobalParam.f90:228`` form."""
    use_TIP = True
    solid = xr.DataArray(np.array([1.0, 2.5, 5.0, 10.0, 25.0]), dims="cell")
    kdpo4 = 100.0

    v3_value = v3_fdp(use_TIP=use_TIP, Solid=solid, kdpo4=kdpo4)
    fortran_value = _fortran_anchored_fdp(
        use_TIP=use_TIP, Solid=solid, kdpo4=kdpo4
    )

    np.testing.assert_allclose(
        np.asarray(v3_value), np.asarray(fortran_value), rtol=1e-6
    )


def test_fdp_unit_factor_dimensionally_correct():
    """Audit C5 (Phase 9.B): v3 ``fdp`` matches the Fortran dimensionally
    correct form, not the v1 inverted form."""
    use_TIP = True
    kdpo4 = 0.001
    solid = 10.0

    v3_value = v3_fdp(
        use_TIP=use_TIP, Solid=xr.DataArray([solid]), kdpo4=kdpo4
    )
    fortran_value = _fortran_anchored_fdp(
        use_TIP=use_TIP, Solid=xr.DataArray([solid]), kdpo4=kdpo4
    )
    buggy_v1_value = _v1_shared_fdp_buggy(
        use_TIP=use_TIP, Solid=xr.DataArray([solid]), kdpo4=kdpo4
    )

    expected_corrected = 1.0 / (1.0 + 0.001 * 10.0 * 1.0e-6)
    np.testing.assert_allclose(
        float(v3_value.values[0]), expected_corrected, rtol=1e-12
    )
    np.testing.assert_allclose(
        float(v3_value.values[0]), float(fortran_value.values[0]),
        rtol=1e-12,
    )
    expected_buggy = 1.0 / (1.0 + 0.001 * 10.0 / 1.0e-6)
    np.testing.assert_allclose(
        float(buggy_v1_value.values[0]), expected_buggy, rtol=1e-6
    )
    assert float(v3_value.values[0]) > 0.99
    assert float(buggy_v1_value.values[0]) < 1e-3


def test_phase9e_vsop_consistent_with_vsap():
    """Phase 9.E: ``vsop = 0.1`` m/d consistent with ``vsap = 0.15`` m/d."""
    from clearwater_modules_v3.parameters.phosphorus import (
        DEFAULTS as PHOSPHORUS_DEFAULTS,
    )
    from clearwater_modules_v3.parameters.algae import (
        DEFAULTS as ALGAE_DEFAULTS,
    )

    vsop = PHOSPHORUS_DEFAULTS["vsop"]
    vsap = ALGAE_DEFAULTS["vsap"]

    assert vsop == 0.1
    assert vsap == 0.15
    assert vsap / 10.0 <= vsop <= vsap * 10.0
    assert vsop > 0.05
