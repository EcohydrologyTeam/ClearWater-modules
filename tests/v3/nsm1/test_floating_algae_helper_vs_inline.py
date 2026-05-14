"""Phase 5 helper-vs-inline parity for FloatingAlgae.

Asserts that ``FloatingAlgae._change_with_components`` and the shadow
``FloatingAlgae._change_legacy_inline`` produce bit-identical net
rate (ug-Chla/L/d) across a parametrised state/forcing matrix.

Deleted in Phase 10 alongside its shadow per §11.3.

Tolerance: ``rtol=0, atol=0``.
"""

from __future__ import annotations

from datetime import timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.floating_algae import FloatingAlgae


@pytest.fixture
def fa() -> FloatingAlgae:
    f = FloatingAlgae(time_step=timedelta(minutes=5))
    f.use_ammonium = True
    f.use_nitrate = True
    f.use_phosphate = True
    return f


N_CELLS = 5


def _da(values: np.ndarray) -> xr.DataArray:
    return xr.DataArray(values, dims="cell")


def _zero_state():
    z = np.zeros(N_CELLS)
    return dict(
        algae=_da(z),
        ammonium=_da(z),
        nitrate=_da(z),
        phosphorus_total_inorganic=_da(z),
        depth=_da(z + 1.0),
        water_temperature=_da(z + 20.0),
        solar=_da(z + 200.0),
    )


def _uniform_state():
    return dict(
        algae=_da(np.full(N_CELLS, 5.0)),
        ammonium=_da(np.full(N_CELLS, 0.5)),
        nitrate=_da(np.full(N_CELLS, 1.0)),
        phosphorus_total_inorganic=_da(np.full(N_CELLS, 0.05)),
        depth=_da(np.full(N_CELLS, 1.0)),
        water_temperature=_da(np.full(N_CELLS, 20.0)),
        solar=_da(np.full(N_CELLS, 200.0)),
    )


def _randomised_state(seed: int = 20260513):
    rng = np.random.default_rng(seed)
    return dict(
        algae=_da(rng.uniform(0.5, 30.0, N_CELLS)),
        ammonium=_da(rng.uniform(0.0, 2.0, N_CELLS)),
        nitrate=_da(rng.uniform(0.05, 5.0, N_CELLS)),
        phosphorus_total_inorganic=_da(rng.uniform(0.001, 0.5, N_CELLS)),
        depth=_da(rng.uniform(0.3, 5.0, N_CELLS)),
        water_temperature=_da(rng.uniform(5.0, 30.0, N_CELLS)),
        solar=_da(rng.uniform(50.0, 800.0, N_CELLS)),
    )


def _high_light():
    return dict(
        algae=_da(np.full(N_CELLS, 5.0)),
        ammonium=_da(np.full(N_CELLS, 0.5)),
        nitrate=_da(np.full(N_CELLS, 1.0)),
        phosphorus_total_inorganic=_da(np.full(N_CELLS, 0.05)),
        depth=_da(np.full(N_CELLS, 1.0)),
        water_temperature=_da(np.full(N_CELLS, 20.0)),
        solar=_da(np.full(N_CELLS, 1000.0)),
    )


def _low_light():
    return dict(
        algae=_da(np.full(N_CELLS, 5.0)),
        ammonium=_da(np.full(N_CELLS, 0.5)),
        nitrate=_da(np.full(N_CELLS, 1.0)),
        phosphorus_total_inorganic=_da(np.full(N_CELLS, 0.05)),
        depth=_da(np.full(N_CELLS, 1.0)),
        water_temperature=_da(np.full(N_CELLS, 20.0)),
        solar=_da(np.full(N_CELLS, 5.0)),
    )


def _nutrient_limited_n():
    return dict(
        algae=_da(np.full(N_CELLS, 5.0)),
        ammonium=_da(np.full(N_CELLS, 0.001)),
        nitrate=_da(np.full(N_CELLS, 0.001)),
        phosphorus_total_inorganic=_da(np.full(N_CELLS, 1.0)),
        depth=_da(np.full(N_CELLS, 1.0)),
        water_temperature=_da(np.full(N_CELLS, 20.0)),
        solar=_da(np.full(N_CELLS, 200.0)),
    )


def _nutrient_limited_p():
    return dict(
        algae=_da(np.full(N_CELLS, 5.0)),
        ammonium=_da(np.full(N_CELLS, 1.0)),
        nitrate=_da(np.full(N_CELLS, 5.0)),
        phosphorus_total_inorganic=_da(np.full(N_CELLS, 0.0001)),
        depth=_da(np.full(N_CELLS, 1.0)),
        water_temperature=_da(np.full(N_CELLS, 20.0)),
        solar=_da(np.full(N_CELLS, 200.0)),
    )


SCENARIOS = [
    ("zero_state", _zero_state),
    ("uniform_state", _uniform_state),
    ("randomised_state", _randomised_state),
    ("high_light", _high_light),
    ("low_light", _low_light),
    ("nutrient_limited_n", _nutrient_limited_n),
    ("nutrient_limited_p", _nutrient_limited_p),
]


@pytest.mark.parametrize("label,factory", SCENARIOS, ids=[s[0] for s in SCENARIOS])
def test_helper_matches_inline_bit_identical(
    fa: FloatingAlgae, label: str, factory
) -> None:
    """``_change_with_components`` and ``_change_legacy_inline`` must
    produce bit-equal net rates."""
    kwargs = factory()

    rate_new, _components = fa._change_with_components(**kwargs)
    rate_old = fa._change_legacy_inline(**kwargs)

    np.testing.assert_array_equal(
        np.asarray(rate_new), np.asarray(rate_old),
        err_msg=f"FloatingAlgae rate differs for scenario {label!r}",
    )


@pytest.mark.parametrize("label,factory", SCENARIOS, ids=[s[0] for s in SCENARIOS])
def test_components_dict_contains_all_registry_diagnostics(
    fa: FloatingAlgae, label: str, factory
) -> None:
    kwargs = factory()
    _, components = fa._change_with_components(**kwargs)
    assert set(components.keys()) == set(fa.REGISTRY_DIAGNOSTICS)


def test_helper_zero_state_produces_finite_components(fa: FloatingAlgae) -> None:
    kwargs = _zero_state()
    _, components = fa._change_with_components(**kwargs)
    for name, value in components.items():
        arr = np.asarray(value)
        assert np.isfinite(arr).all(), (
            f"{name} contains non-finite values at zero state"
        )


def test_limit_diagnostics_in_unit_interval(fa: FloatingAlgae) -> None:
    """``algal_light_limitation`` and ``algal_nutrient_limitation_*`` are
    limitation factors and must lie in [0, 1] in physically reasonable
    regimes."""
    kwargs = _uniform_state()
    _, components = fa._change_with_components(**kwargs)
    for name in (
        "algal_light_limitation",
        "algal_nutrient_limitation_n",
        "algal_nutrient_limitation_p",
    ):
        arr = np.asarray(components[name])
        assert np.all(arr >= 0.0), f"{name} has negative values"
        assert np.all(arr <= 1.0 + 1e-12), f"{name} exceeds 1.0"


def test_preserved_attribute_caches_match_components(fa: FloatingAlgae) -> None:
    """The seven preserved-name caches that sibling Processes consume
    via getattr (algal_growth_rate, algal_respiration_rate,
    algal_*_from_mortality_rate, algal_pom_from_settling_rate,
    algal_nh4_uptake_fraction) must equal the components dict values
    after _change_with_components runs."""
    kwargs = _uniform_state()
    _, components = fa._change_with_components(**kwargs)

    for name in (
        "algal_growth_rate",
        "algal_respiration_rate",
        "algal_death_rate",
        "algal_settling_rate",
        "algal_orgn_from_mortality_rate",
        "algal_orgp_from_mortality_rate",
        "algal_poc_from_mortality_rate",
        "algal_doc_from_mortality_rate",
        "algal_pom_from_settling_rate",
        "algal_nh4_uptake_fraction",
    ):
        cached = np.asarray(getattr(fa, name))
        from_components = np.asarray(components[name])
        np.testing.assert_array_equal(
            cached, from_components,
            err_msg=f"{name}: self.{name} != components[{name!r}]",
        )
