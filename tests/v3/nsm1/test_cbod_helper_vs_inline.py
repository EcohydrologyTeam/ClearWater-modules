"""Phase 7 helper-vs-inline parity for CBOD.

Asserts that ``CBOD._change_with_components`` and the shadow
``CBOD._change_legacy_inline`` produce bit-identical net rates
(mg-O2/L/d).

Deleted in Phase 10 alongside its shadow per §11.3. Tolerance:
``rtol=0, atol=0``.
"""

from __future__ import annotations

from datetime import timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.cbod import CBOD


@pytest.fixture
def c() -> CBOD:
    proc = CBOD(time_step=timedelta(minutes=5))
    proc.use_DOX = True
    return proc


@pytest.fixture
def c_no_dox() -> CBOD:
    proc = CBOD(time_step=timedelta(minutes=5))
    proc.use_DOX = False
    return proc


N_CELLS = 5


def _da(values: np.ndarray) -> xr.DataArray:
    return xr.DataArray(values, dims="cell")


def _zero_state():
    z = np.zeros(N_CELLS)
    return dict(
        cbod=_da(z),
        water_temperature=_da(z + 20.0),
        depth=_da(z + 1.0),
        dox=_da(z + 8.0),
    )


def _uniform_state():
    return dict(
        cbod=_da(np.full(N_CELLS, 5.0)),
        water_temperature=_da(np.full(N_CELLS, 20.0)),
        depth=_da(np.full(N_CELLS, 1.0)),
        dox=_da(np.full(N_CELLS, 8.0)),
    )


def _randomised_state(seed: int = 20260513):
    rng = np.random.default_rng(seed)
    return dict(
        cbod=_da(rng.uniform(0.5, 30.0, N_CELLS)),
        water_temperature=_da(rng.uniform(5.0, 30.0, N_CELLS)),
        depth=_da(rng.uniform(0.3, 5.0, N_CELLS)),
        dox=_da(rng.uniform(0.5, 14.0, N_CELLS)),
    )


def _hypoxic():
    """DOX → 0: Monod attenuation reduces oxidation toward zero."""
    return dict(
        cbod=_da(np.full(N_CELLS, 5.0)),
        water_temperature=_da(np.full(N_CELLS, 20.0)),
        depth=_da(np.full(N_CELLS, 1.0)),
        dox=_da(np.full(N_CELLS, 0.05)),
    )


def _hot_water():
    return dict(
        cbod=_da(np.full(N_CELLS, 5.0)),
        water_temperature=_da(np.full(N_CELLS, 35.0)),
        depth=_da(np.full(N_CELLS, 1.0)),
        dox=_da(np.full(N_CELLS, 8.0)),
    )


SCENARIOS = [
    ("zero_state", _zero_state),
    ("uniform_state", _uniform_state),
    ("randomised_state", _randomised_state),
    ("hypoxic", _hypoxic),
    ("hot_water", _hot_water),
]


@pytest.mark.parametrize("label,factory", SCENARIOS, ids=[s[0] for s in SCENARIOS])
def test_helper_matches_inline_use_dox(c: CBOD, label: str, factory) -> None:
    kwargs = factory()

    rate_new, _components = c._change_with_components(**kwargs)
    rate_old = c._change_legacy_inline(**kwargs)

    np.testing.assert_array_equal(
        np.asarray(rate_new), np.asarray(rate_old),
        err_msg=f"CBOD rate (use_DOX=True) differs for scenario {label!r}",
    )


@pytest.mark.parametrize("label,factory", SCENARIOS, ids=[s[0] for s in SCENARIOS])
def test_helper_matches_inline_no_dox(c_no_dox: CBOD, label: str, factory) -> None:
    """Same parity, but with use_DOX = False — oxidation is first-order
    in CBOD only (no Monod attenuation)."""
    kwargs = factory()

    rate_new, _components = c_no_dox._change_with_components(**kwargs)
    rate_old = c_no_dox._change_legacy_inline(**kwargs)

    np.testing.assert_array_equal(
        np.asarray(rate_new), np.asarray(rate_old),
        err_msg=f"CBOD rate (use_DOX=False) differs for scenario {label!r}",
    )


@pytest.mark.parametrize("label,factory", SCENARIOS, ids=[s[0] for s in SCENARIOS])
def test_components_dict_contains_all_registry_diagnostics(
    c: CBOD, label: str, factory
) -> None:
    kwargs = factory()
    _, components = c._change_with_components(**kwargs)
    assert set(components.keys()) == set(c.REGISTRY_DIAGNOSTICS)


def test_components_are_positive_magnitudes(c: CBOD) -> None:
    """``cbod_oxidation_rate`` and ``cbod_settling_rate`` are
    positive magnitudes; the integrator applies the negation in the
    net rate sum."""
    kwargs = _uniform_state()
    _, components = c._change_with_components(**kwargs)
    for name in c.REGISTRY_DIAGNOSTICS:
        arr = np.asarray(components[name])
        assert np.all(arr >= 0.0), f"{name} contains negative values"


def test_helper_zero_state_produces_finite_components(c: CBOD) -> None:
    kwargs = _zero_state()
    _, components = c._change_with_components(**kwargs)
    for name, value in components.items():
        arr = np.asarray(value)
        assert np.isfinite(arr).all()


def test_hypoxic_oxidation_attenuates(c: CBOD) -> None:
    """At low DOX, oxidation rate must be smaller than at high DOX
    (Monod attenuation)."""
    normal = _uniform_state()
    _, c_normal = c._change_with_components(**normal)

    low = _uniform_state()
    low["dox"] = _da(np.full(N_CELLS, 0.01))
    _, c_low = c._change_with_components(**low)

    ox_normal = np.asarray(c_normal["cbod_oxidation_rate"])
    ox_low = np.asarray(c_low["cbod_oxidation_rate"])
    assert np.all(ox_low <= ox_normal + 1e-12)
