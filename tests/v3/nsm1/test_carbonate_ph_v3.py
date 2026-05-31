"""Carbonate-system pH diagnostic (NSM1-I).

The carbonate solver ``clearwater_modules_v3.utils.carbonate`` is shared
byte-for-byte with the nsm2-and-hab line (NSM2 step S4-3); the first four
tests are the shared solver tests (Tier-5 byte-exact NSM2 constants at
I=0, Tier-3 residual/speciation, D-A-4 graceful failure, Emerson f_NH3).
The remaining tests cover this branch's wiring: ``Alkalinity.run``
computes pH from (alkalinity, dic, water_temperature) at freshwater I=0
as an opportunistic diagnostic, so it is written only when ``"pH"`` is
pre-registered and never changes the alkalinity trajectory (in NSM1-I the
DIC reaeration term uses a constant Fco2, so pH feeds no kinetics).
"""
from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np
import xarray as xr

from clearwater_modules_v3.utils import carbonate as cb
from clearwater_modules_v3.processes.alkalinity import Alkalinity
from clearwater_modules_v3.utils.numerics import Diagnostics
from .conftest import InMemoryRegistry


START = datetime(2026, 1, 1, 0, 0, 0)


def _nsm2_constants(t_c: float):
    tk = t_c + 273.15
    kw = 10.0 ** (-4787.3 / tk - 7.1321 * np.log10(tk)
                  - 0.010365 * tk + 22.80)
    k1 = 10.0 ** (-356.3094 - 0.06091964 * tk + 21834.37 / tk
                  + 126.8339 * np.log10(tk) - 1684915.0 / tk**2)
    k2 = 10.0 ** (-107.8871 - 0.03252849 * tk + 5151.79 / tk
                  + 38.92561 * np.log10(tk) - 563713.9 / tk**2)
    return kw, k1, k2


# ---------------------------------------------------------------------------
# Shared carbonate-solver tests (utils.carbonate; byte-identical to nsm2-and-hab)
# ---------------------------------------------------------------------------


def test_apparent_constants_byte_exact_at_I0():
    """At I=0 the Davies correction vanishes and Kw/K1/K2 are
    byte-identical to the modAlkalinity temperature formulas."""
    for t_c in (5.0, 15.0, 22.0, 30.0):
        kw, k1, k2 = cb.apparent_constants(t_c, 0.0)
        nkw, nk1, nk2 = _nsm2_constants(t_c)
        assert kw == nkw and k1 == nk1 and k2 == nk2, (
            f"I=0 not byte-exact at T={t_c}"
        )
    # I>0 must actually shift the constants.
    kI = cb.apparent_constants(22.0, 0.7)[1]
    assert kI != _nsm2_constants(22.0)[1]


def test_solver_residual_and_speciation():
    kw, k1, k2 = _nsm2_constants(25.0)
    alk = np.array([50.0, 100.0, 200.0])
    dic = np.array([18.0, 24.0, 30.0])
    ph, fb = cb.solve_ph(alk, dic, np.array([kw] * 3), np.array([k1] * 3),
                         np.array([k2] * 3))
    res = cb._residual(ph, dic / 12000.0, alk / 50000.0, kw, k1, k2)
    np.testing.assert_allclose(res, 0.0, atol=1e-8)
    assert not fb.any()
    assert np.all((ph > 4.0) & (ph < 11.0))
    # Higher Alk at fixed DIC -> higher pH.
    ph2, _ = cb.solve_ph(np.array([50.0, 150.0]), np.array([24.0, 24.0]),
                         np.array([kw, kw]), np.array([k1, k1]),
                         np.array([k2, k2]))
    assert ph2[1] > ph2[0]
    # Speciation partitions DIC exactly.
    c0, c1, c2 = cb.speciation(ph, dic, k1, k2)
    np.testing.assert_allclose(c0 + c1 + c2, dic, rtol=1e-12)


def test_da4_graceful_failure_never_raises():
    kw, k1, k2 = _nsm2_constants(20.0)
    n = 5
    # Extreme/garbage inputs must never raise/NaN; pH finite, clamped.
    alk = np.array([1e-6, 0.0, 5e4, -1e3, 1e9])
    dic = np.array([0.0, 1e6, 24.0, 24.0, 0.0])
    ph, fb = cb.solve_ph(alk, dic, np.array([kw] * n), np.array([k1] * n),
                         np.array([k2] * n), prev_ph=np.array([7.3] * n))
    assert np.isfinite(ph).all()
    assert np.all((ph >= 3.0) & (ph <= 13.0))
    assert np.isfinite(fb).all()
    # Strongly-negative alkalinity: no root in [3, 13] -> hold (flagged).
    ph2, fb2 = cb.solve_ph(np.array([-1.0e9]), np.array([24.0]),
                           np.array([kw]), np.array([k1]),
                           np.array([k2]), prev_ph=np.array([7.3]))
    assert np.isfinite(ph2).all()
    assert bool(fb2[0])
    assert 3.0 <= float(ph2[0]) <= 13.0


def test_f_nh3_emerson_pka():
    pka = 0.09018 + 2729.92 / 298.15  # ~9.246 at 25 C; f_NH3(pH=pKa)=0.5
    np.testing.assert_allclose(
        float(cb.f_nh3(np.array([pka]), 25.0)[0]), 0.5, rtol=1e-9
    )
    assert cb.f_nh3(np.array([7.0]), 25.0)[0] < cb.f_nh3(np.array([10.0]), 25.0)[0]


# ---------------------------------------------------------------------------
# Alkalinity.run pH-diagnostic wiring (this branch)
# ---------------------------------------------------------------------------


def _da(value, n=5):
    return xr.DataArray(np.full(n, value, dtype=float), dims=["cell"])


def _registry(*, with_ph=True, with_dic=True):
    reg = InMemoryRegistry()
    reg.register("alkalinity", _da(100.0))
    reg.register("water_temperature", _da(25.0))
    reg.register("depth", _da(1.0))
    if with_dic:
        reg.register("dic", _da(24.0))  # mg-C/L
    if with_ph:
        reg.register("pH", _da(7.0))
    return reg


def _alkalinity():
    alk = Alkalinity(time_step=timedelta(minutes=5))
    alk.diagnostics = Diagnostics()
    return alk


def test_alkalinity_writes_ph_when_registered():
    """With 'pH' and 'dic' registered, run() writes a physically-sensible
    pH (consistent with the standalone solver)."""
    reg = _registry(with_ph=True, with_dic=True)
    _alkalinity().run(START, reg)
    ph = np.asarray(reg.get("pH").values)
    assert np.all((ph > 4.0) & (ph < 11.0))
    # Matches the solver called directly on the same (alk, dic, T).
    kw, k1, k2 = cb.apparent_constants(_da(25.0))
    ph_ref, _ = cb.solve_ph(_da(100.0), _da(24.0), kw, k1, k2)
    np.testing.assert_allclose(ph, np.asarray(ph_ref), rtol=1e-12)


def test_alkalinity_ph_does_not_change_trajectory():
    """Computing pH must not perturb the alkalinity state: alk_new is
    bit-identical whether or not 'pH' is registered (pH is diagnostic)."""
    reg_with = _registry(with_ph=True, with_dic=True)
    reg_without = _registry(with_ph=False, with_dic=True)
    _alkalinity().run(START, reg_with)
    _alkalinity().run(START, reg_without)
    np.testing.assert_array_equal(
        np.asarray(reg_with.get("alkalinity").values),
        np.asarray(reg_without.get("alkalinity").values),
    )
    assert "pH" not in reg_without  # not written when not requested


def test_alkalinity_ph_skipped_without_dic(caplog):
    """'pH' registered but 'dic' absent: pH is not computed (warns once),
    and the run still completes."""
    import logging

    reg = _registry(with_ph=True, with_dic=False)
    with caplog.at_level(logging.WARNING,
                         logger="clearwater_modules_v3.processes.alkalinity"):
        _alkalinity().run(START, reg)
    # pH left at its pre-registered placeholder (not overwritten).
    np.testing.assert_array_equal(
        np.asarray(reg.get("pH").values), np.full(5, 7.0)
    )
    assert any("dic" in r.getMessage() for r in caplog.records)
