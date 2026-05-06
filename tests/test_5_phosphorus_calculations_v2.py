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


def _v1_shared_fdp_buggy(use_TIP, Solid, kdpo4):
    """Inline reproduction of v1 ``shared.processes.fdp`` (the buggy form).

    The actual ``clearwater_modules.shared.processes.fdp`` is decorated
    with ``@numba.njit``, which cannot accept ``xarray.DataArray`` inputs
    — calling it from this parity test under xr.DataArray inputs raises
    a numba TypingError. The formula itself is trivial and stable; we
    inline it here so the parity check exercises the same arithmetic
    without the numba decorator.

    v1 source (clearwater_modules/shared/processes.py:257-271):
        return xr.where(use_TIP, 1 / (1 + kdpo4 * Solid / 0.000001), 0)

    Note (Phase 9.B audit): this v1 form is dimensionally inverted;
    ``kdpo4 * Solid`` is in mg/kg and must be multiplied by 1e-6 kg/mg
    (equivalently, divided by 1e6 mg/kg) to be dimensionless. The v1
    inline ``/0.000001`` is the wrong direction; Fortran
    (``modGlobalParam.f90:228``) writes ``/ 1.0E6`` which is correct.
    v3 follows Fortran. This helper is retained only for negative
    regression tests (asserting v3 deliberately diverges from v1).
    """
    return xr.where(use_TIP, 1.0 / (1.0 + kdpo4 * Solid / 0.000001), 0.0)


def _fortran_anchored_fdp(use_TIP, Solid, kdpo4):
    """Dimensionally correct fdp matching Fortran ``modGlobalParam.f90:228``.

    ``fdp = 1 / (1 + kdpo4 [L/kg] * Solid [mg/L] * 1e-6 [kg/mg])``. v3
    ports this form; this helper is used as the parity reference.
    """
    return xr.where(use_TIP, 1.0 / (1.0 + kdpo4 * Solid * 1.0e-6), 0.0)


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


def test_tip_settling_matches_fortran_anchored(
    water_temp_5cell, depth_5cell, tip_5cell, orgp_5cell
):
    """v3 ``tip_settling_rate`` cache after run == Fortran-anchored
    ``TIP_Settling`` with the dimensionally correct ``fdp``.

    v1 formula: ``vs / depth * (1 - fdp) * TIP`` (line 1973-1988); v3
    caches the same expression under ``self.tip_settling_rate``. The v3
    ``fdp`` utility was corrected in Phase 9.B from v1's
    ``/ 0.000001`` (== multiply by 1e6, dimensionally inverted) to
    Fortran's ``* 1e-6`` (correct). This test asserts v3's ``fdp``
    output drives the TIP settling term toward the Fortran-anchored
    expectation, not the v1 inverted form.
    """
    kdpo4 = 1000.0     # L/kg
    solid = 5.0        # mg/L
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

    # Fortran-anchored fdp (matches v3 partitioning utility).
    fortran_fdp = _fortran_anchored_fdp(
        use_TIP=True, Solid=solid, kdpo4=kdpo4
    )
    v1_rate_with_fortran_fdp = v1.TIP_Settling(
        0.1, depth_5cell, fortran_fdp, tip_5cell
    )

    np.testing.assert_allclose(
        np.asarray(phosphorus.tip_settling_rate),
        np.asarray(v1_rate_with_fortran_fdp),
        rtol=1e-6,
    )


def test_tip_partitioning_fdp_matches_fortran():
    """v3 ``utils.partitioning.fdp`` == Fortran ``modGlobalParam.f90:228`` form.

    Phase 9.B audit fix: v3 corrects the unit factor from v1's inverted
    ``/ 0.000001`` (multiplies dimensional product by 1e6) to Fortran's
    ``* 1e-6`` (correct dimensionless conversion). This test compares
    v3 against the Fortran-anchored reference, not the buggy v1 inline.
    """
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


# ---------------------------------------------------------------------------
# Audit-anchored regression tests (Phase 9.B)
# ---------------------------------------------------------------------------


def test_fdp_unit_factor_dimensionally_correct():
    """Audit C5 (Phase 9.B): v3 ``fdp`` matches the Fortran dimensionally
    correct form, not the v1 inverted form.

    Setup: ``kdpo4 = 0.001`` L/kg, ``Solid = 10`` mg/L. The dimensional
    product ``kdpo4 * Solid * 1e-6 = 1e-8``, so corrected ``fdp ~= 1``
    (essentially all P stays dissolved at low solids load). The v1 form
    yields ``1 / (1 + 1e4) ~= 1e-4`` (essentially all P sorbs to
    nonexistent solids, physically nonsensical).
    """
    use_TIP = True
    kdpo4 = 0.001       # L/kg
    solid = 10.0        # mg/L

    v3_value = v3_fdp(
        use_TIP=use_TIP,
        Solid=xr.DataArray([solid]),
        kdpo4=kdpo4,
    )
    fortran_value = _fortran_anchored_fdp(
        use_TIP=use_TIP,
        Solid=xr.DataArray([solid]),
        kdpo4=kdpo4,
    )
    buggy_v1_value = _v1_shared_fdp_buggy(
        use_TIP=use_TIP,
        Solid=xr.DataArray([solid]),
        kdpo4=kdpo4,
    )

    # v3 should match Fortran ~ 1.0 (all dissolved at low solids).
    expected_corrected = 1.0 / (1.0 + 0.001 * 10.0 * 1.0e-6)  # ~ 1 - 1e-8
    np.testing.assert_allclose(
        float(v3_value.values[0]), expected_corrected, rtol=1e-12
    )
    np.testing.assert_allclose(
        float(v3_value.values[0]), float(fortran_value.values[0]),
        rtol=1e-12,
    )
    # And v3 should NOT match the buggy v1 form (which yields ~ 1e-4).
    expected_buggy = 1.0 / (1.0 + 0.001 * 10.0 / 1.0e-6)  # ~ 1e-4
    np.testing.assert_allclose(
        float(buggy_v1_value.values[0]), expected_buggy, rtol=1e-6
    )
    assert float(v3_value.values[0]) > 0.99
    assert float(buggy_v1_value.values[0]) < 1e-3


# ---------------------------------------------------------------------------
# Phase 9.E regression: vsop physical-consistency with vsap
# ---------------------------------------------------------------------------
# OrgP in NSM1 originates predominantly from dead-algae detritus (via the
# algal mortality routing ``algal_orgp_from_mortality_rate`` and the
# benthic ``balgae_orgp_from_mortality_rate``). The algae from which OrgP
# derives have a universal NSM1 settling velocity ``vsap = 0.15`` m/d
# (Fortran/v1/v3 all agree). For internal physical consistency, OrgP
# detritus inherited from the same algal pool should settle at a comparable
# rate. Phase 9.E pins ``vsop = 0.1`` m/d on this physical-consistency
# basis. Fortran's 0.01 m/d is 15x slower than the algae from which the
# OrgP derives and is implausible as a representative default. See
# parameter_defaults_corrections.md Section 1.1.

def test_phase9e_vsop_consistent_with_vsap():
    """Phase 9.E: ``vsop = 0.1`` m/d (mid-range of literature 0.01-1.0)
    consistent with ``vsap = 0.15`` m/d algal settling. OrgP detritus
    derived from dead algae should settle on the same order of magnitude
    as the algae itself; Fortran's 0.01 m/d (15x slower) is physically
    inconsistent with the algal-detritus origin."""
    from clearwater_modules_v3.parameters.phosphorus import (
        DEFAULTS as PHOSPHORUS_DEFAULTS,
    )
    from clearwater_modules_v3.parameters.algae import (
        DEFAULTS as ALGAE_DEFAULTS,
    )

    vsop = PHOSPHORUS_DEFAULTS["vsop"]
    vsap = ALGAE_DEFAULTS["vsap"]

    # Pin the canonical Phase 9.E value.
    assert vsop == 0.1, (
        f"vsop should be 0.1 m/d (Phase 9.E physical-consistency choice); "
        f"got {vsop}"
    )
    # Pin the universal NSM1 algal settling default for cross-check.
    assert vsap == 0.15

    # Physical-consistency check: vsop within an order of magnitude of
    # vsap. Fortran's 0.01 m/d (15x slower than vsap) would fail this.
    assert vsap / 10.0 <= vsop <= vsap * 10.0, (
        f"vsop = {vsop} m/d should be within an order of magnitude of "
        f"vsap = {vsap} m/d (algal detritus inheriting OrgP)"
    )
    # And confirm vsop > 0.05 (i.e., not at Fortran's lower-bound 0.01).
    assert vsop > 0.05, (
        f"vsop = {vsop} m/d should not be at Fortran's 0.01 lower bound; "
        f"v3 deliberately chose a value consistent with algal detritus."
    )
