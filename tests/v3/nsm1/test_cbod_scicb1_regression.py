"""NSM1-SCI-CB1 (MAJOR) regression: CBOD settling is a first-order
**rate** (1/d at 20 °C), not a settling velocity.

Gold-standard spec Workstream C2; E3 (match Fortran 1/d, θ=1.024).

Fortran NSM1 ``modCBOD.f90:114`` and QUAL2E apply ``ksbod_tc * CBOD``
with **no depth division**, and use the *settling* Arrhenius
coefficient ``ksbod_theta = 1.024`` (Bowie 1985 / QUAL2E). Pre-fix v3
used ``ksbod_tc / depth * cbod`` (a velocity, m/d) with
``ksbod_theta = 1.047`` (the *oxidation* coefficient) — silent at the
shipped ``ksbod_20 = 0`` default but divergent by ``1/depth`` and a
wrong θ for any calibrated nonzero value.

Non-shared-path contract (spec Section 1(4)): the expected 1/d form is
built from **independently hardcoded** literals (θ = 1.024, the van't
Hoff form, no depth term), NOT by reading the process ``ksbod_theta``
or its settling code.
"""
from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np
import xarray as xr

from clearwater_modules_v3.processes.cbod import CBOD

from .conftest import InMemoryRegistry


KSBOD_20 = 0.2          # 1/d at 20 °C (a calibrated nonzero settling rate)
KSBOD_THETA_LIT = 1.024  # independently hardcoded; the constant under test
KSBOD_THETA_PREFIX = 1.047  # the pre-fix (oxidation) coefficient
T_WATER = 25.0
CBOD_CONC = 12.0
DT = timedelta(minutes=5)


def _registry(depth: float) -> InMemoryRegistry:
    reg = InMemoryRegistry()
    one = lambda v: xr.DataArray(np.array([v]), dims="cell")
    reg.register("cbod", one(CBOD_CONC))
    reg.register("water_temperature", one(T_WATER))
    reg.register("depth", one(depth))
    reg.register("oxygen_dissolved", one(8.0))
    return reg


def test_scicb1_ksbod_theta_default_is_1024():
    c = CBOD(time_step=DT)
    assert c.ksbod_theta == KSBOD_THETA_LIT
    assert c.ksbod_theta != KSBOD_THETA_PREFIX


def test_scicb1_settling_is_first_order_rate_no_depth_division():
    """Nonzero ksbod_20: settling rate equals the Fortran 1/d form and
    is identical at two very different depths (proves no 1/depth)."""
    # Independent expectation (van't Hoff, hardcoded θ=1.024, no depth).
    ksbod_tc = KSBOD_20 * (KSBOD_THETA_LIT ** (T_WATER - 20.0))
    expected_settling = ksbod_tc * CBOD_CONC  # mg-O2/L/d

    c_shallow = CBOD(parameters={"ksbod_20": KSBOD_20}, time_step=DT)
    c_deep = CBOD(parameters={"ksbod_20": KSBOD_20}, time_step=DT)
    t = datetime(2026, 5, 16)
    c_shallow.run(t, _registry(depth=1.0))
    c_deep.run(t, _registry(depth=5.0))

    s_shallow = np.asarray(c_shallow.cbod_settling_rate)
    s_deep = np.asarray(c_deep.cbod_settling_rate)

    # 1/d form: settling is independent of depth.
    np.testing.assert_allclose(s_shallow, s_deep, rtol=1e-12)
    # Matches the independently computed Fortran 1/d form.
    np.testing.assert_allclose(s_shallow, expected_settling, rtol=1e-12)
    # Hard anti-regression: NOT the pre-fix velocity form (which would
    # be 1/depth smaller and depth-dependent).
    prefix_shallow = expected_settling / 1.0
    prefix_deep = expected_settling / 5.0
    assert not np.allclose(s_deep, prefix_deep), (
        "settling still divides by depth (pre-fix m/d velocity form)"
    )
    # (At depth=1.0 the pre-fix and fixed forms coincide; the depth=5.0
    # leg is the discriminating one.)
    assert np.all(np.abs(s_deep - expected_settling) < 1e-12)
