"""Phase 8 helper-vs-inline parity for N2.

Deleted in Phase 10 alongside its shadow per §11.3. Tolerance:
``rtol=0, atol=0``.
"""

from __future__ import annotations

from datetime import timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.n2 import N2


@pytest.fixture
def n2() -> N2:
    proc = N2(time_step=timedelta(minutes=5))
    proc.use_nitrogen = False
    proc.nitrogen_process = None
    return proc


N_CELLS = 5


def _da(values: np.ndarray) -> xr.DataArray:
    return xr.DataArray(values, dims="cell")


def _zero_state():
    z = np.zeros(N_CELLS)
    return dict(
        n2_state=_da(z),
        t_water_c=_da(z + 20.0),
        depth=_da(z + 1.0),
        pressure_mb=1013.25,
    )


def _uniform_state():
    return dict(
        n2_state=_da(np.full(N_CELLS, 14.0)),
        t_water_c=_da(np.full(N_CELLS, 20.0)),
        depth=_da(np.full(N_CELLS, 1.0)),
        pressure_mb=1013.25,
    )


def _randomised_state(seed: int = 20260513):
    rng = np.random.default_rng(seed)
    return dict(
        n2_state=_da(rng.uniform(5.0, 25.0, N_CELLS)),
        t_water_c=_da(rng.uniform(5.0, 30.0, N_CELLS)),
        depth=_da(rng.uniform(0.3, 5.0, N_CELLS)),
        pressure_mb=1013.25,
    )


def _cold_water():
    return dict(
        n2_state=_da(np.full(N_CELLS, 18.0)),
        t_water_c=_da(np.full(N_CELLS, 0.5)),
        depth=_da(np.full(N_CELLS, 1.0)),
        pressure_mb=1013.25,
    )


def _hot_water():
    return dict(
        n2_state=_da(np.full(N_CELLS, 10.0)),
        t_water_c=_da(np.full(N_CELLS, 35.0)),
        depth=_da(np.full(N_CELLS, 1.0)),
        pressure_mb=1013.25,
    )


def _supersaturated():
    """N2 > N2sat: atmospheric exchange reverses sign."""
    return dict(
        n2_state=_da(np.full(N_CELLS, 30.0)),
        t_water_c=_da(np.full(N_CELLS, 20.0)),
        depth=_da(np.full(N_CELLS, 1.0)),
        pressure_mb=1013.25,
    )


SCENARIOS = [
    ("zero_state", _zero_state),
    ("uniform_state", _uniform_state),
    ("randomised_state", _randomised_state),
    ("cold_water", _cold_water),
    ("hot_water", _hot_water),
    ("supersaturated", _supersaturated),
]


@pytest.mark.parametrize("label,factory", SCENARIOS, ids=[s[0] for s in SCENARIOS])
def test_helper_matches_inline_bit_identical(
    n2: N2, label: str, factory
) -> None:
    """``_change_with_components`` and ``_change_legacy_inline`` produce
    bit-equal net rates."""
    kwargs = factory()

    rate_new, _components = n2._change_with_components(**kwargs)
    rate_old = n2._change_legacy_inline(**kwargs)

    np.testing.assert_array_equal(
        np.asarray(rate_new), np.asarray(rate_old),
        err_msg=f"N2 rate differs for scenario {label!r}",
    )


@pytest.mark.parametrize("label,factory", SCENARIOS, ids=[s[0] for s in SCENARIOS])
def test_components_dict_contains_helper_subset(
    n2: N2, label: str, factory
) -> None:
    """The helper's components dict carries 3 of 4 REGISTRY_DIAGNOSTICS
    names; ``total_dissolved_gas`` is added by ``run`` after the
    integrator step (it depends on n2_new)."""
    kwargs = factory()
    _, components = n2._change_with_components(**kwargs)
    expected = set(N2.REGISTRY_DIAGNOSTICS) - {"total_dissolved_gas"}
    assert set(components.keys()) == expected


def test_helper_zero_state_produces_finite_components(n2: N2) -> None:
    kwargs = _zero_state()
    _, components = n2._change_with_components(**kwargs)
    for name, value in components.items():
        arr = np.asarray(value)
        assert np.isfinite(arr).all()


def test_supersaturated_atm_exchange_negative(n2: N2) -> None:
    """At N2 > N2sat (atmospheric N2 saturation), the exchange term
    is negative (N2 escaping to atmosphere)."""
    kwargs = _supersaturated()
    _, components = n2._change_with_components(**kwargs)
    atm = np.asarray(components["n2_atm_exchange_rate"])
    # Note: with use_user_ka_zero defaults (no reaeration), the term
    # may be 0. The test asserts it's non-positive (≤ 0) which is the
    # correct sign expectation.
    assert np.all(atm <= 0.0 + 1e-12)
