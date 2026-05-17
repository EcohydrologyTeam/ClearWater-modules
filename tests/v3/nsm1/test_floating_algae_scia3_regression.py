"""NSM1-SCI-A3 (MAJOR) regression: FloatingAlgae must drive light
limitation with **PAR**, not total broadband shortwave.

Gold-standard spec Workstream B1.

``solar_radiation`` in the registry is total shortwave (W/m^2). NSM1 v1
applied ``PAR = q_solar * Fr_PAR`` (Fr_PAR=0.47, ``nsm1/constants.py:350``)
upstream of the floating-algae kinetic. Pre-fix v3 dropped that
conversion and passed total shortwave straight into ``limit_light``
against the PAR-scale ``KL`` (``light_limitation_constant``), a v1->v3
regression that under-limits light and over-predicts algal growth
~30-60% wherever light is binding. v3 restores the conversion at the
process boundary.

Non-shared-path contract (spec Section 1(4)): the expected PAR factor
is built from an **independently hardcoded** ``Fr_PAR`` literal (0.47),
NOT by reading ``fa.Fr_PAR`` or importing the parameter DEFAULTS.
"""
from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np
import xarray as xr

from clearwater_modules_v3.processes.floating_algae import FloatingAlgae

from .conftest import InMemoryRegistry


FR_PAR_LITERAL = 0.47  # independently hardcoded; the constant under test
SOLAR_SHORTWAVE = 300.0  # W/m^2 total broadband shortwave
DT = timedelta(minutes=5)


def _fa() -> FloatingAlgae:
    inst = FloatingAlgae(
        time_step=DT,
        growth_rate_option=1,
        growth_rate_max=1.0,
        light_limitation_option=1,
        light_limitation_constant=10.0,
        light_attenuation_coefficient=1.0,
    )
    inst.use_nitrate = True
    inst.use_ammonium = True
    inst.use_phosphate = True
    return inst


def _registry() -> InMemoryRegistry:
    reg = InMemoryRegistry()
    one = lambda v: xr.DataArray(np.array([v]), dims="cell")
    reg.register("algae_floating", one(30.0))
    reg.register("ammonium", one(0.15))
    reg.register("nitrate", one(3.0))
    reg.register("phosphorus_total_inorganic", one(0.10))
    reg.register("depth", one(1.5))
    reg.register("water_temperature", one(20.0))
    reg.register("solar_radiation", one(SOLAR_SHORTWAVE))
    return reg


def test_scia3_light_limitation_uses_par_not_shortwave():
    """The cached ``algal_light_limitation`` after ``run`` equals
    ``limit_light`` evaluated at PAR = SW * 0.47 (independently
    computed), and NOT at raw shortwave."""
    fa = _fa()
    reg = _registry()
    t = datetime(2026, 5, 16)
    fa.run(t, reg)

    algae = reg.get_at_time("algae_floating", t)
    depth = reg.get_at_time("depth", t)

    # Independent "v1-mirror" expectation: v1 applies Fr_PAR upstream,
    # then the FL kinetic sees PAR. Hardcoded 0.47 (not fa.Fr_PAR).
    expected_par = SOLAR_SHORTWAVE * FR_PAR_LITERAL
    expected_limit = fa.limit_light(
        algae=algae, depth=depth, surface_light_intensity=expected_par
    )
    prefix_defect_limit = fa.limit_light(
        algae=algae, depth=depth, surface_light_intensity=SOLAR_SHORTWAVE
    )

    cached = np.asarray(fa.algal_light_limitation)
    np.testing.assert_allclose(
        cached, np.asarray(expected_limit), rtol=1e-12
    )
    # Hard anti-regression: must NOT be the no-Fr_PAR (raw shortwave)
    # value, and PAR must be strictly more light-limiting (smaller
    # half-saturation factor) than raw shortwave.
    assert np.all(np.abs(cached - np.asarray(prefix_defect_limit)) > 1e-9)
    assert np.all(cached < np.asarray(prefix_defect_limit))


def test_scia3_registry_solar_is_total_shortwave_unchanged():
    """The fix converts at the process boundary only; it must not
    mutate the ``solar_radiation`` registry variable (it remains total
    shortwave for any other consumer)."""
    fa = _fa()
    reg = _registry()
    t = datetime(2026, 5, 16)
    fa.run(t, reg)
    np.testing.assert_array_equal(
        np.asarray(reg.get_at_time("solar_radiation", t)),
        np.array([SOLAR_SHORTWAVE]),
    )
