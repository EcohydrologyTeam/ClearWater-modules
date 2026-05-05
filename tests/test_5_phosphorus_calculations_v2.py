"""Parity tests: v3 Phosphorus sub-rate methods vs v1 nsm1.processes helpers.

Phase 7.B of the v3 NSM1 implementation plan. Mirrors the established
Phase 2 pattern at ``tests/test_5_floating_algae_calculations_v2.py``.

Each test constructs a v3 Phosphorus instance, drives it through one
``Process.run`` substep against an in-memory registry, and compares the
resulting cached step-scoped rate (or its constitutive sub-term) to the
equivalent v1 helper-function output computed with the same inputs.

Scope: v1-equivalent ``dTIPdt`` / ``dOrgPdt`` sub-terms:

* OrgP -> TIP hydrolysis (``kop_tc * OrgP``)
* OrgP settling (``vsop / depth * OrgP``)
* TIP settling (``vs / depth * (1 - fdp) * TIP``)
* TIP partitioning fdp utility (v3 ``utils.partitioning.fdp`` ==
  v1 ``shared.processes.fdp``)

Synthetic mesh: 5-cell xarray DataArrays, single time step.
"""
from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules.nsm1 import processes as v1
from clearwater_modules_v3.processes.phosphorus import Phosphorus
from clearwater_modules_v3.utils.partitioning import fdp as v3_fdp

from tests.v3.nsm1.conftest import InMemoryRegistry


def _v1_shared_fdp(use_TIP, Solid, kdpo4):
    """Inline reproduction of v1 ``shared.processes.fdp``.

    The actual ``clearwater_modules.shared.processes.fdp`` is decorated
    with ``@numba.njit``, which cannot accept ``xarray.DataArray`` inputs
    — calling it from this parity test under xr.DataArray inputs raises
    a numba TypingError. The formula itself is trivial and stable; we
    inline it here so the parity check exercises the same arithmetic
    without the numba decorator.

    v1 source (clearwater_modules/shared/processes.py:257-271):
        return xr.where(use_TIP, 1 / (1 + kdpo4 * Solid / 0.000001), 0)
    """
    return xr.where(use_TIP, 1.0 / (1.0 + kdpo4 * Solid / 0.000001), 0.0)


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
def tip_5cell() -> xr.DataArray:
    return xr.DataArray(
        np.array([0.10, 0.12, 0.14, 0.16, 0.18]), dims="cell"
    )


@pytest.fixture(scope="function")
def orgp_5cell() -> xr.DataArray:
    return xr.DataArray(
        np.array([0.05, 0.06, 0.07, 0.08, 0.09]), dims="cell"
    )


def _build_registry(
    tip: xr.DataArray,
    orgp: xr.DataArray,
    water_temp: xr.DataArray,
    depth: xr.DataArray,
) -> InMemoryRegistry:
    """Wire the four state variables that Phosphorus.run reads."""
    registry = InMemoryRegistry()
    registry.register("tip", tip.copy())
    registry.register("organic_phosphorus", orgp.copy())
    registry.register("water_temperature", water_temp.copy())
    registry.register("depth", depth.copy())
    return registry


# ---------------------------------------------------------------------------
# Parity tests
# ---------------------------------------------------------------------------


def test_orgp_hydrolysis_matches_v1(
    water_temp_5cell, depth_5cell, tip_5cell, orgp_5cell
):
    """v3 ``orgp_to_tip_hydrolysis_rate`` cache after run == v1 ``OrgP_DIP_decay``.

    v1 formula: ``kop_tc(T, kop_20, kop_theta) * OrgP`` (when use_OrgP).
    v3 caches the same expression under ``self.orgp_to_tip_hydrolysis_rate``
    after ``run`` (phosphorus.py line 287).
    """
    phosphorus = Phosphorus(
        parameters={
            "kop_20": 0.1,
            "kop_theta": 1.047,
            # Disable other terms to keep the sub-rate test isolated.
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

    # v1 reference.
    kop_tc = v1.kop_tc(water_temp_5cell, 0.1, 1.047)
    v1_rate = v1.OrgP_DIP_decay(kop_tc, orgp_5cell, use_OrgP=True)

    np.testing.assert_allclose(
        np.asarray(phosphorus.orgp_to_tip_hydrolysis_rate),
        np.asarray(v1_rate),
        rtol=1e-6,
    )


def test_orgp_settling_matches_v1(
    water_temp_5cell, depth_5cell, tip_5cell, orgp_5cell
):
    """v3 ``orgp_settling_rate`` cache after run == v1 ``OrgP_Settling``.

    v1 formula: ``vsop / depth * OrgP`` (line 1882-1895).
    v3 caches the same expression under ``self.orgp_settling_rate``
    (phosphorus.py line 301).
    """
    phosphorus = Phosphorus(
        parameters={
            "vsop": 0.1,
            # Disable hydrolysis & TIP settling so OrgP settling is isolated.
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

    v1_rate = v1.OrgP_Settling(0.1, depth_5cell, orgp_5cell)

    np.testing.assert_allclose(
        np.asarray(phosphorus.orgp_settling_rate),
        np.asarray(v1_rate),
        rtol=1e-6,
    )


def test_tip_settling_matches_v1(
    water_temp_5cell, depth_5cell, tip_5cell, orgp_5cell
):
    """v3 ``tip_settling_rate`` cache after run == v1 ``TIP_Settling`` with v3 ``fdp``.

    v1 formula: ``vs / depth * (1 - fdp) * TIP`` (line 1973-1988).
    v3 caches the same expression under ``self.tip_settling_rate``
    (phosphorus.py line 294). The fdp here comes from the v3 partitioning
    utility (``fdp(use_TIP, Solid, kdpo4)``); we pass a non-zero
    ``kdpo4 > 0`` so ``fdp < 1`` and the test exercises the full term.
    """
    kdpo4 = 1000.0     # L/kg; non-zero to force fdp < 1
    solid = 5.0        # mg/L; non-zero to force fdp < 1
    phosphorus = Phosphorus(
        parameters={
            "vs": 0.1,
            "kdpo4": kdpo4,
            "Solid": solid,
            "use_TIP": True,
            # Disable other terms so TIP settling is isolated.
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

    # v1 fdp formula (shared module): 1 / (1 + kdpo4 * Solid / 1e-6).
    v1_fdp = _v1_shared_fdp(use_TIP=True, Solid=solid, kdpo4=kdpo4)
    v1_rate = v1.TIP_Settling(0.1, depth_5cell, v1_fdp, tip_5cell)

    np.testing.assert_allclose(
        np.asarray(phosphorus.tip_settling_rate),
        np.asarray(v1_rate),
        rtol=1e-6,
    )


def test_tip_partitioning_fdp_matches_v1():
    """v3 ``utils.partitioning.fdp`` == v1 ``shared.processes.fdp``.

    v1 signature: ``fdp(use_TIP, Solid, kdpo4)``; v3 utility uses the
    same signature. Both implement
    ``xr.where(use_TIP, 1 / (1 + kdpo4 * Solid / 1e-6), 0)``.

    Note: there are TWO v1 ``fdp`` definitions (Phase 0.2 audit):
    ``nsm1.processes.fdp`` is a degenerate stub returning 1.0 when
    use_TIP, while ``shared.processes.fdp`` is the proper formula. The
    v3 utility ports the ``shared`` formula, so we parity-check against
    that.
    """
    use_TIP = True
    solid = xr.DataArray(np.array([1.0, 2.5, 5.0, 10.0, 25.0]), dims="cell")
    kdpo4 = 100.0

    v3_value = v3_fdp(use_TIP=use_TIP, Solid=solid, kdpo4=kdpo4)
    v1_value = _v1_shared_fdp(use_TIP=use_TIP, Solid=solid, kdpo4=kdpo4)

    np.testing.assert_allclose(
        np.asarray(v3_value), np.asarray(v1_value), rtol=1e-6
    )
