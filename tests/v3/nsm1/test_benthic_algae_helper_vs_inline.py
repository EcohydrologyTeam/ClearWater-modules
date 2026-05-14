"""Phase 5 helper-vs-inline parity for BenthicAlgae.

Asserts that ``BenthicAlgae._change_with_components`` and the shadow
``BenthicAlgae._change_legacy_inline`` produce bit-identical net rate
(g-D/m^2/d) across a parametrised state/forcing matrix.

Specifically pins the **rate_death dedup contract**: the helper invokes
``rate_death`` once and reuses the cached value, while the shadow
invokes it twice (once via ``rate()`` and once via
``_cache_benthic_mortality_rates``). Because ``rate_death`` is pure,
both paths produce bit-identical outputs.

Deleted in Phase 10 alongside its shadow per §11.3.

Tolerance: ``rtol=0, atol=0``.
"""

from __future__ import annotations

from datetime import timedelta
from unittest.mock import patch

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.benthic_algae import BenthicAlgae


@pytest.fixture
def ba() -> BenthicAlgae:
    b = BenthicAlgae(time_step=timedelta(minutes=5))
    b.use_ammonium = True
    b.use_nitrate = True
    b.use_phosphate = True
    return b


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
        algae=_da(np.full(N_CELLS, 50.0)),
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
        algae=_da(rng.uniform(1.0, 100.0, N_CELLS)),
        ammonium=_da(rng.uniform(0.0, 2.0, N_CELLS)),
        nitrate=_da(rng.uniform(0.05, 5.0, N_CELLS)),
        phosphorus_total_inorganic=_da(rng.uniform(0.001, 0.5, N_CELLS)),
        depth=_da(rng.uniform(0.3, 5.0, N_CELLS)),
        water_temperature=_da(rng.uniform(5.0, 30.0, N_CELLS)),
        solar=_da(rng.uniform(50.0, 800.0, N_CELLS)),
    )


def _high_density():
    """Density-limitation regime: high benthic-algae density triggers
    the limit_density Monod attenuation."""
    return dict(
        algae=_da(np.full(N_CELLS, 500.0)),
        ammonium=_da(np.full(N_CELLS, 0.5)),
        nitrate=_da(np.full(N_CELLS, 1.0)),
        phosphorus_total_inorganic=_da(np.full(N_CELLS, 0.05)),
        depth=_da(np.full(N_CELLS, 1.0)),
        water_temperature=_da(np.full(N_CELLS, 20.0)),
        solar=_da(np.full(N_CELLS, 200.0)),
    )


def _high_temp():
    return dict(
        algae=_da(np.full(N_CELLS, 50.0)),
        ammonium=_da(np.full(N_CELLS, 0.5)),
        nitrate=_da(np.full(N_CELLS, 1.0)),
        phosphorus_total_inorganic=_da(np.full(N_CELLS, 0.05)),
        depth=_da(np.full(N_CELLS, 1.0)),
        water_temperature=_da(np.full(N_CELLS, 35.0)),
        solar=_da(np.full(N_CELLS, 200.0)),
    )


def _thin_depth():
    return dict(
        algae=_da(np.full(N_CELLS, 50.0)),
        ammonium=_da(np.full(N_CELLS, 0.5)),
        nitrate=_da(np.full(N_CELLS, 1.0)),
        phosphorus_total_inorganic=_da(np.full(N_CELLS, 0.05)),
        depth=_da(np.full(N_CELLS, 0.05)),
        water_temperature=_da(np.full(N_CELLS, 20.0)),
        solar=_da(np.full(N_CELLS, 200.0)),
    )


SCENARIOS = [
    ("zero_state", _zero_state),
    ("uniform_state", _uniform_state),
    ("randomised_state", _randomised_state),
    ("high_density", _high_density),
    ("high_temp", _high_temp),
    ("thin_depth", _thin_depth),
]


@pytest.mark.parametrize("label,factory", SCENARIOS, ids=[s[0] for s in SCENARIOS])
def test_helper_matches_inline_bit_identical(
    ba: BenthicAlgae, label: str, factory
) -> None:
    """``_change_with_components`` (single rate_death call) and
    ``_change_legacy_inline`` (two rate_death calls) must produce
    bit-equal net rates."""
    kwargs = factory()

    rate_new, _components = ba._change_with_components(**kwargs)
    rate_old = ba._change_legacy_inline(**kwargs)

    np.testing.assert_array_equal(
        np.asarray(rate_new), np.asarray(rate_old),
        err_msg=f"BenthicAlgae rate differs for scenario {label!r}",
    )


@pytest.mark.parametrize("label,factory", SCENARIOS, ids=[s[0] for s in SCENARIOS])
def test_components_dict_contains_all_registry_diagnostics(
    ba: BenthicAlgae, label: str, factory
) -> None:
    kwargs = factory()
    _, components = ba._change_with_components(**kwargs)
    assert set(components.keys()) == set(ba.REGISTRY_DIAGNOSTICS)


def test_dedup_calls_rate_death_exactly_once(ba: BenthicAlgae) -> None:
    """Phase 5 dedup contract: ``_change_with_components`` invokes
    ``rate_death`` exactly once. Pre-Phase-5 (and the shadow) invoke it
    twice. The values are identical because rate_death is pure, but the
    call count is the actual dedup we want to verify."""
    kwargs = _uniform_state()

    with patch.object(
        ba, "rate_death", wraps=ba.rate_death
    ) as spy:
        ba._change_with_components(**kwargs)
    assert spy.call_count == 1, (
        f"_change_with_components called rate_death {spy.call_count} times; "
        "expected 1 (Phase 5 dedup)"
    )


def test_shadow_calls_rate_death_twice(ba: BenthicAlgae) -> None:
    """Pin the legacy two-call behaviour for the shadow so the dedup
    is meaningful (otherwise this test would silently pass even if the
    shadow were also de-duplicated)."""
    kwargs = _uniform_state()

    with patch.object(
        ba, "rate_death", wraps=ba.rate_death
    ) as spy:
        ba._change_legacy_inline(**kwargs)
    assert spy.call_count == 2, (
        f"_change_legacy_inline called rate_death {spy.call_count} times; "
        "expected 2 (legacy duplicate)"
    )


def test_helper_zero_state_produces_finite_components(ba: BenthicAlgae) -> None:
    kwargs = _zero_state()
    _, components = ba._change_with_components(**kwargs)
    for name, value in components.items():
        arr = np.asarray(value)
        assert np.isfinite(arr).all(), (
            f"{name} contains non-finite values at zero state"
        )


def test_limit_diagnostics_in_unit_interval(ba: BenthicAlgae) -> None:
    kwargs = _uniform_state()
    _, components = ba._change_with_components(**kwargs)
    for name in (
        "balgae_light_limitation",
        "balgae_nutrient_limitation_n",
        "balgae_nutrient_limitation_p",
    ):
        arr = np.asarray(components[name])
        assert np.all(arr >= 0.0), f"{name} has negative values"
        assert np.all(arr <= 1.0 + 1e-12), f"{name} exceeds 1.0"


def test_preserved_attribute_caches_match_components(ba: BenthicAlgae) -> None:
    """The eight preserved-name caches that sibling Processes consume
    must equal the components dict values after _change_with_components
    runs."""
    kwargs = _uniform_state()
    _, components = ba._change_with_components(**kwargs)

    for name in (
        "balgae_growth_rate",
        "balgae_respiration_rate",
        "balgae_death_rate",
        "balgae_orgn_from_mortality_rate",
        "balgae_orgp_from_mortality_rate",
        "balgae_poc_from_mortality_rate",
        "balgae_doc_from_mortality_rate",
        "balgae_nh4_uptake_fraction",
    ):
        cached = np.asarray(getattr(ba, name))
        from_components = np.asarray(components[name])
        np.testing.assert_array_equal(
            cached, from_components,
            err_msg=f"{name}: self.{name} != components[{name!r}]",
        )


def test_pom_from_mortality_rate_still_set(ba: BenthicAlgae) -> None:
    """``balgae_pom_from_mortality_rate`` is consumed by POM but is
    NOT in REGISTRY_DIAGNOSTICS (it's not in Appendix A's BenthicAlgae
    list — see closeout). It must still be set as a side effect of
    the dedup helper for POM consumption."""
    kwargs = _uniform_state()
    ba._change_with_components(**kwargs)
    pom_routing = np.asarray(getattr(ba, "balgae_pom_from_mortality_rate"))
    assert np.isfinite(pom_routing).all()
