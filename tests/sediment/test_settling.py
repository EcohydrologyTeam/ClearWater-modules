"""Unit tests for SSM settling-velocity helpers (Cheng 1997).

References
----------
- Cheng, N.-S. (1997). "Simplified settling velocity formula for sediment
  particle." J. Hydraul. Eng. 123(2), 149-152.
- SAND2008-5621 eq. 3 (Thanh, Grace & James 2008).

Locked-in expected values (computed by hand at default constants
``g = 980 cm/s^2``, ``nu = 0.01 cm^2/s``, ``s_s = 2.65``,
``rho_w = 1.0 g/cm^3``):

    factor = ((s_s - 1) * g / nu^2)^(1/3)
           = (1.65 * 980 / 1e-4)^(1/3)
           = 16_170_000^(1/3)
           = 252.83...

For D50 = 125 um -> D = 0.0125 cm:
    d*    = 0.0125 * 252.83 = 3.1604
    inner = sqrt(25 + 1.2 * 3.1604^2) - 5
          = sqrt(25 + 11.984) - 5
          = sqrt(36.984) - 5
          = 6.0815 - 5
          = 1.0815
    w_s   = (0.01 / 0.0125) * 1.0815^1.5
          = 0.8 * 1.1247
          = 0.8997 cm/s

For D50 = 250 um -> D = 0.025 cm:
    d*    = 6.3208
    inner = sqrt(25 + 47.929) - 5 = 8.5398 - 5 = 3.5398
    w_s   = (0.01 / 0.025) * 3.5398^1.5
          = 0.4 * 6.6595
          = 2.664 cm/s

For D50 = 100 um -> D = 0.01 cm:
    d*    = 2.5283
    inner = sqrt(25 + 7.6710) - 5 = 5.7159 - 5 = 0.7159
    w_s   = 1.0 * 0.7159^1.5 = 0.6058 cm/s
"""

from __future__ import annotations

import numpy as np
import pytest

from clearwater_modules_v2.processes.sediment.settling import (
    cheng_1997_settling_velocity,
    resolve_settling_velocities,
)
from clearwater_modules_v2.processes.sediment.classes import (
    SedimentClass,
    SedimentClassRegistry,
)


# ---------------------------------------------------------------------------
# cheng_1997_settling_velocity — scalar and vector behaviour
# ---------------------------------------------------------------------------


def test_cheng_scalar_d50_125um_matches_hand_calc():
    """D50 = 125 μm should give w_s ≈ 0.8997 cm/s (see module docstring)."""
    ws = cheng_1997_settling_velocity(125.0)
    assert isinstance(ws, float)
    assert ws == pytest.approx(0.8997, rel=1e-3)


def test_cheng_scalar_d50_250um_matches_hand_calc():
    """D50 = 250 μm should give w_s ≈ 2.664 cm/s; spec sanity gate ~2.7 cm/s."""
    ws = cheng_1997_settling_velocity(250.0)
    assert ws == pytest.approx(2.664, rel=1e-3)
    # Spec sanity gate.
    assert 2.5 < ws < 2.9


def test_cheng_scalar_d50_100um_in_sanity_band():
    """D50 = 100 μm should give w_s ≈ 0.6 cm/s (spec sanity ~0.71, not exact)."""
    ws = cheng_1997_settling_velocity(100.0)
    assert ws == pytest.approx(0.6058, rel=1e-3)
    # Loose sanity gate from spec.
    assert 0.4 < ws < 0.9


def test_cheng_vectorized_over_numpy_array():
    """Vectorized call should match per-element scalar calls."""
    d50 = np.array([50.0, 100.0, 125.0, 250.0, 500.0])
    ws_vec = cheng_1997_settling_velocity(d50)
    ws_scalar = np.array([cheng_1997_settling_velocity(float(d)) for d in d50])
    assert isinstance(ws_vec, np.ndarray)
    assert ws_vec.shape == d50.shape
    np.testing.assert_allclose(ws_vec, ws_scalar, rtol=1e-12)


def test_cheng_monotonic_increasing_with_d50():
    """w_s should increase monotonically with D50 over the spec range."""
    d50 = np.geomspace(20.0, 2000.0, 25)
    ws = cheng_1997_settling_velocity(d50)
    assert np.all(np.diff(ws) > 0)


def test_cheng_rejects_nonpositive_d50():
    with pytest.raises(ValueError):
        cheng_1997_settling_velocity(0.0)
    with pytest.raises(ValueError):
        cheng_1997_settling_velocity(np.array([100.0, -1.0, 200.0]))


def test_cheng_density_dependence():
    """Heavier particles should settle faster than lighter ones at same D50."""
    ws_quartz = cheng_1997_settling_velocity(125.0, solid_specific_gravity=2.65)
    ws_heavy = cheng_1997_settling_velocity(125.0, solid_specific_gravity=4.0)
    assert ws_heavy > ws_quartz


# ---------------------------------------------------------------------------
# resolve_settling_velocities — registry-driven dispatch
# ---------------------------------------------------------------------------


def test_resolve_uses_user_supplied_when_positive():
    reg = SedimentClassRegistry.from_iterable(
        [
            SedimentClass(label="silt", d50_um=20.0, settling_cm_s=0.05),
            SedimentClass(label="sand", d50_um=250.0, settling_cm_s=3.5),
        ]
    )
    ws = resolve_settling_velocities(reg)
    assert ws.shape == (2,)
    np.testing.assert_allclose(ws, [0.05, 3.5])


def test_resolve_falls_back_to_cheng_when_unset():
    reg = SedimentClassRegistry.from_iterable(
        [
            SedimentClass(label="fine", d50_um=125.0, settling_cm_s=None),
            SedimentClass(label="coarse", d50_um=250.0, settling_cm_s=-1.0),
        ]
    )
    ws = resolve_settling_velocities(reg)
    assert ws.shape == (2,)
    assert ws[0] == pytest.approx(0.8997, rel=1e-3)
    assert ws[1] == pytest.approx(2.664, rel=1e-3)


def test_resolve_mixed_user_and_cheng():
    reg = SedimentClassRegistry.from_iterable(
        [
            SedimentClass(label="user", d50_um=125.0, settling_cm_s=1.234),
            SedimentClass(label="cheng", d50_um=125.0, settling_cm_s=None),
        ]
    )
    ws = resolve_settling_velocities(reg)
    assert ws[0] == pytest.approx(1.234)
    assert ws[1] == pytest.approx(0.8997, rel=1e-3)
