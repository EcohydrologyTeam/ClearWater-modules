"""Standalone unit tests for v3 ``Temperature.mixing_ratio_air``.

Purpose: pin the air mixing ratio formula's behavior, especially the
edge-case guard introduced by the C4 fix described in
``design/clearwater_modules_v3_review_findings.md``.

Background: the dimensionless mass mixing ratio of water vapor in air is

    r = 0.622 * e_air / (P_air - e_air)

where ``e_air`` is atmospheric water-vapor pressure and ``P_air`` is
atmospheric pressure (both in millibars). The denominator vanishes when
``e_air == P_air`` (saturation at total pressure) and goes negative when
``e_air > P_air``, which can occur from data-entry errors, mis-scaled
forcing fields, or sensor noise near saturation.

The original v3 guard caught only ``denom == 0.0`` (an exact-equality
on floats: a measure-zero set). The C4 fix changes the guard to
``denom <= 0.0`` so the negative-denominator case also returns 0.0
rather than producing a negative mixing ratio that propagates through
``density_air``'s ``(1 + r) / (1 + 1.61 r)`` factor and yields
sign-flipped or near-singular air densities, poisoning every flux that
depends on ``density_air``.

These tests pin both the normal and edge-case behavior so that any
regression of the guard fails with a specific, readable error.
"""
from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.temperature import Temperature


@pytest.fixture(scope="module")
def temperature_module() -> Temperature:
    """Stub Temperature instance.

    The mixing-ratio method does not depend on wind or sediment
    parameters, so any valid stub values for the constructor are fine.
    """
    return Temperature(0.3, 1.5, 3.0)


def _r(temp: Temperature, e_air: float, P_air: float) -> float:
    """Helper: invoke the mixing-ratio formula on scalar inputs and
    return the float result.

    Note the kernel signature is
    ``mixing_ratio_air(atmospheric_vapor_pressure, atmospheric_pressure)``
    -- vapor pressure first, total pressure second.
    """
    return float(temp.mixing_ratio_air(e_air, P_air))


def test_mixing_ratio_normal_case(temperature_module: Temperature):
    """Normal case: ``e_air = 20 mb``, ``P_air = 1013 mb`` returns the
    standard formula value to floating-point precision.

    This is the dominant operating regime; any breakage here would mean
    every air-density-derived flux is wrong even before the edge-case
    guard fires.
    """
    e_air = 20.0
    P_air = 1013.0
    expected = 0.622 * e_air / (P_air - e_air)
    result = _r(temperature_module, e_air, P_air)
    assert result == pytest.approx(expected, rel=1e-12, abs=1e-15), (
        f"mixing_ratio_air({e_air}, {P_air}) = {result}; "
        f"expected {expected} (= 0.622 * e / (P - e))"
    )


def test_mixing_ratio_equality_case(temperature_module: Temperature):
    """Equality case: ``e_air == P_air`` returns 0.0.

    This case was already handled by the original guard
    (``denom == 0.0``); pin it so a future refactor that removes the
    equality branch (e.g. by switching to ``denom < 0.0``) is caught.
    """
    e_air = 1013.0
    P_air = 1013.0
    result = _r(temperature_module, e_air, P_air)
    assert result == 0.0, (
        f"mixing_ratio_air({e_air}, {P_air}) = {result}; "
        "expected 0.0 (equality-case guard)."
    )


def test_mixing_ratio_negative_denominator(temperature_module: Temperature):
    """C4 guard: ``e_air > P_air`` returns 0.0.

    Without the C4 fix, ``denom = P - e = -87 mb`` is negative and the
    formula returns ``0.622 * 1100 / (-87) ~= -7.86``, a negative
    mixing ratio that poisons ``density_air``'s
    ``(1 + r) / (1 + 1.61 r)`` factor. With the C4 fix, the guard
    catches this and returns 0.0.

    This is the primary regression guard for the C4 finding.
    """
    e_air = 1100.0
    P_air = 1013.0
    result = _r(temperature_module, e_air, P_air)
    assert result == 0.0, (
        f"mixing_ratio_air({e_air}, {P_air}) = {result}; "
        "expected 0.0 (C4 negative-denominator guard). "
        "If the result is negative (~ -7.86), the C4 fix has regressed "
        "and the guard is back to ``denom == 0.0``."
    )
    # Defensive: also verify the result is not negative, which is the
    # specific failure mode the C4 fix prevents.
    assert result >= 0.0, (
        f"mixing_ratio_air returned a negative value ({result}); "
        "the C4 guard must prevent this."
    )


def test_mixing_ratio_vectorized_mixed_case(temperature_module: Temperature):
    """Vectorized mixed case: a DataArray with cells in normal /
    equality / e>P regimes returns ``[normal, 0.0, 0.0]``.

    This exercises the ``xr.where`` path for both branches of the guard
    simultaneously and confirms the array-valued kernel matches the
    scalar-valued behavior cell-by-cell.
    """
    e_air = xr.DataArray(np.array([20.0, 1013.0, 1100.0]))
    P_air = xr.DataArray(np.array([1013.0, 1013.0, 1013.0]))

    expected_normal = 0.622 * 20.0 / (1013.0 - 20.0)
    expected = np.array([expected_normal, 0.0, 0.0])

    out = temperature_module.mixing_ratio_air(e_air, P_air)
    out_np = np.asarray(out)

    assert out_np.shape == (3,), (
        f"vectorized output has shape {out_np.shape}; expected (3,)"
    )
    np.testing.assert_allclose(
        out_np,
        expected,
        rtol=1e-12,
        atol=1e-15,
        err_msg=(
            "vectorized mixing_ratio_air result does not match "
            "expected [normal, 0.0, 0.0]; check that the C4 ``<= 0.0`` "
            "guard is applied in both inner and outer ``xr.where``."
        ),
    )

    # Also pin the cell-wise values explicitly so a regression has a
    # readable failure message.
    assert out_np[0] == pytest.approx(expected_normal, rel=1e-12, abs=1e-15)
    assert out_np[1] == 0.0
    assert out_np[2] == 0.0
    assert out_np[2] >= 0.0, (
        "C4 guard regression: the e>P cell returned a negative value."
    )
