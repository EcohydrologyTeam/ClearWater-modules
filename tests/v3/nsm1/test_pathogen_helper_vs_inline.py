"""Phase 8 helper-vs-inline parity for Pathogen.

Pathogen uses the rate-form integrator naming convention per spec §10
Q5: the canonical helper is ``_rate_with_components`` (single-state
rate, not ``_change_with_components`` which returns deltas).

Deleted in Phase 10 alongside its shadow per §11.3. Tolerance:
``rtol=0, atol=0``.
"""

from __future__ import annotations

from datetime import timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.pathogen import Pathogen


@pytest.fixture
def path() -> Pathogen:
    return Pathogen(time_step=timedelta(minutes=5))


N_CELLS = 5


def _da(values: np.ndarray) -> xr.DataArray:
    return xr.DataArray(values, dims="cell")


def _zero_state():
    z = np.zeros(N_CELLS)
    return dict(
        px=_da(z),
        water_temperature=_da(z + 20.0),
        depth=_da(z + 1.0),
        q_solar=_da(z + 200.0),
        solid=_da(z),
        poc=_da(z),
        ap=_da(z),
    )


def _uniform_state():
    return dict(
        px=_da(np.full(N_CELLS, 1000.0)),
        water_temperature=_da(np.full(N_CELLS, 20.0)),
        depth=_da(np.full(N_CELLS, 1.0)),
        q_solar=_da(np.full(N_CELLS, 200.0)),
        solid=_da(np.full(N_CELLS, 10.0)),
        poc=_da(np.full(N_CELLS, 1.0)),
        ap=_da(np.full(N_CELLS, 5.0)),
    )


def _randomised_state(seed: int = 20260513):
    rng = np.random.default_rng(seed)
    return dict(
        px=_da(rng.uniform(10.0, 10000.0, N_CELLS)),
        water_temperature=_da(rng.uniform(5.0, 30.0, N_CELLS)),
        depth=_da(rng.uniform(0.3, 5.0, N_CELLS)),
        q_solar=_da(rng.uniform(50.0, 800.0, N_CELLS)),
        solid=_da(rng.uniform(0.0, 50.0, N_CELLS)),
        poc=_da(rng.uniform(0.0, 5.0, N_CELLS)),
        ap=_da(rng.uniform(0.0, 50.0, N_CELLS)),
    )


def _high_light():
    return dict(
        px=_da(np.full(N_CELLS, 1000.0)),
        water_temperature=_da(np.full(N_CELLS, 20.0)),
        depth=_da(np.full(N_CELLS, 1.0)),
        q_solar=_da(np.full(N_CELLS, 1000.0)),
        solid=_da(np.full(N_CELLS, 10.0)),
        poc=_da(np.full(N_CELLS, 1.0)),
        ap=_da(np.full(N_CELLS, 5.0)),
    )


def _hot_water():
    return dict(
        px=_da(np.full(N_CELLS, 1000.0)),
        water_temperature=_da(np.full(N_CELLS, 35.0)),
        depth=_da(np.full(N_CELLS, 1.0)),
        q_solar=_da(np.full(N_CELLS, 200.0)),
        solid=_da(np.full(N_CELLS, 10.0)),
        poc=_da(np.full(N_CELLS, 1.0)),
        ap=_da(np.full(N_CELLS, 5.0)),
    )


def _thin_depth():
    return dict(
        px=_da(np.full(N_CELLS, 1000.0)),
        water_temperature=_da(np.full(N_CELLS, 20.0)),
        depth=_da(np.full(N_CELLS, 0.05)),
        q_solar=_da(np.full(N_CELLS, 200.0)),
        solid=_da(np.full(N_CELLS, 10.0)),
        poc=_da(np.full(N_CELLS, 1.0)),
        ap=_da(np.full(N_CELLS, 5.0)),
    )


SCENARIOS = [
    ("zero_state", _zero_state),
    ("uniform_state", _uniform_state),
    ("randomised_state", _randomised_state),
    ("high_light", _high_light),
    ("hot_water", _hot_water),
    ("thin_depth", _thin_depth),
]


@pytest.mark.parametrize("label,factory", SCENARIOS, ids=[s[0] for s in SCENARIOS])
def test_helper_matches_inline_bit_identical(
    path: Pathogen, label: str, factory
) -> None:
    """``_rate_with_components`` and ``_rate_legacy_inline`` produce
    bit-equal net rates."""
    kwargs = factory()

    rate_new, _components = path._rate_with_components(**kwargs)
    rate_old = path._rate_legacy_inline(**kwargs)

    np.testing.assert_array_equal(
        np.asarray(rate_new), np.asarray(rate_old),
        err_msg=f"Pathogen rate differs for scenario {label!r}",
    )


@pytest.mark.parametrize("label,factory", SCENARIOS, ids=[s[0] for s in SCENARIOS])
def test_components_dict_contains_all_registry_diagnostics(
    path: Pathogen, label: str, factory
) -> None:
    kwargs = factory()
    _, components = path._rate_with_components(**kwargs)
    assert set(components.keys()) == set(path.REGISTRY_DIAGNOSTICS)


def test_components_are_positive_magnitudes(path: Pathogen) -> None:
    """The three sub-rates are positive-magnitude sinks; the integrator
    negates them in the net rate sum."""
    kwargs = _uniform_state()
    _, components = path._rate_with_components(**kwargs)
    for name in path.REGISTRY_DIAGNOSTICS:
        arr = np.asarray(components[name])
        assert np.all(arr >= 0.0), f"{name} contains negative values"


def test_helper_zero_state_produces_finite_components(path: Pathogen) -> None:
    kwargs = _zero_state()
    _, components = path._rate_with_components(**kwargs)
    for name, value in components.items():
        arr = np.asarray(value)
        assert np.isfinite(arr).all(), (
            f"{name} contains non-finite values at zero state"
        )


def test_high_light_increases_light_death(path: Pathogen) -> None:
    """Higher solar radiation → stronger light-induced decay."""
    normal = _uniform_state()
    _, c_normal = path._rate_with_components(**normal)

    high = _uniform_state()
    high["q_solar"] = _da(np.full(N_CELLS, 1000.0))
    _, c_high = path._rate_with_components(**high)

    light_normal = np.asarray(c_normal["pathogen_light_death_rate"])
    light_high = np.asarray(c_high["pathogen_light_death_rate"])
    assert np.all(light_high >= light_normal)
