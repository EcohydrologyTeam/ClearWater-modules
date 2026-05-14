"""Phase 6 helper-vs-inline parity for Phosphorus.

Asserts that ``Phosphorus._change_with_components`` and the shadow
``Phosphorus._change_legacy_inline`` produce bit-identical
``(dtip_dt, dorgp_dt)`` across a parametrised state/forcing matrix.

Deleted in Phase 10 alongside its shadow per §11.3.

Tolerance: ``rtol=0, atol=0``.
"""

from __future__ import annotations

from datetime import timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.phosphorus import Phosphorus


@pytest.fixture
def p() -> Phosphorus:
    proc = Phosphorus(time_step=timedelta(minutes=5))
    proc.use_TIP = True
    proc.use_OrgP = True
    proc.use_floating_algae = False
    proc.use_benthic_algae = False
    return proc


@pytest.fixture
def p_tip_only() -> Phosphorus:
    proc = Phosphorus(time_step=timedelta(minutes=5))
    proc.use_TIP = True
    proc.use_OrgP = False
    proc.use_floating_algae = False
    proc.use_benthic_algae = False
    return proc


@pytest.fixture
def p_orgp_only() -> Phosphorus:
    proc = Phosphorus(time_step=timedelta(minutes=5))
    proc.use_TIP = False
    proc.use_OrgP = True
    proc.use_floating_algae = False
    proc.use_benthic_algae = False
    return proc


N_CELLS = 5


def _da(values: np.ndarray) -> xr.DataArray:
    return xr.DataArray(values, dims="cell")


def _zero_state():
    z = np.zeros(N_CELLS)
    return dict(
        tip=_da(z),
        orgp=_da(z),
        water_temperature=_da(z + 20.0),
        depth=_da(z + 1.0),
    )


def _uniform_state():
    return dict(
        tip=_da(np.full(N_CELLS, 0.05)),
        orgp=_da(np.full(N_CELLS, 0.1)),
        water_temperature=_da(np.full(N_CELLS, 20.0)),
        depth=_da(np.full(N_CELLS, 1.0)),
    )


def _randomised_state(seed: int = 20260513):
    rng = np.random.default_rng(seed)
    return dict(
        tip=_da(rng.uniform(0.001, 0.5, N_CELLS)),
        orgp=_da(rng.uniform(0.001, 1.0, N_CELLS)),
        water_temperature=_da(rng.uniform(5.0, 30.0, N_CELLS)),
        depth=_da(rng.uniform(0.3, 5.0, N_CELLS)),
    )


def _cold_water():
    return dict(
        tip=_da(np.full(N_CELLS, 0.05)),
        orgp=_da(np.full(N_CELLS, 0.1)),
        water_temperature=_da(np.full(N_CELLS, 0.5)),
        depth=_da(np.full(N_CELLS, 1.0)),
    )


def _hot_water():
    return dict(
        tip=_da(np.full(N_CELLS, 0.05)),
        orgp=_da(np.full(N_CELLS, 0.1)),
        water_temperature=_da(np.full(N_CELLS, 35.0)),
        depth=_da(np.full(N_CELLS, 1.0)),
    )


def _thin_depth():
    return dict(
        tip=_da(np.full(N_CELLS, 0.05)),
        orgp=_da(np.full(N_CELLS, 0.1)),
        water_temperature=_da(np.full(N_CELLS, 20.0)),
        depth=_da(np.full(N_CELLS, 0.05)),
    )


SCENARIOS = [
    ("zero_state", _zero_state),
    ("uniform_state", _uniform_state),
    ("randomised_state", _randomised_state),
    ("cold_water", _cold_water),
    ("hot_water", _hot_water),
    ("thin_depth", _thin_depth),
]


@pytest.mark.parametrize("label,factory", SCENARIOS, ids=[s[0] for s in SCENARIOS])
def test_helper_matches_inline_use_tip_and_orgp(
    p: Phosphorus, label: str, factory
) -> None:
    """``_change_with_components`` and ``_change_legacy_inline`` produce
    bit-equal deltas when both states are active."""
    kwargs = factory()

    dtip_new, dorgp_new, _components = p._change_with_components(**kwargs)
    dtip_old, dorgp_old = p._change_legacy_inline(**kwargs)

    np.testing.assert_array_equal(
        np.asarray(dtip_new), np.asarray(dtip_old),
        err_msg=f"dtip_dt differs for scenario {label!r}",
    )
    np.testing.assert_array_equal(
        np.asarray(dorgp_new), np.asarray(dorgp_old),
        err_msg=f"dorgp_dt differs for scenario {label!r}",
    )


@pytest.mark.parametrize("label,factory", SCENARIOS, ids=[s[0] for s in SCENARIOS])
def test_helper_matches_inline_tip_only(
    p_tip_only: Phosphorus, label: str, factory
) -> None:
    """Same parity, but with use_OrgP = False. dorgp_dt is 0 in both."""
    kwargs = factory()

    dtip_new, dorgp_new, _components = p_tip_only._change_with_components(**kwargs)
    dtip_old, dorgp_old = p_tip_only._change_legacy_inline(**kwargs)

    np.testing.assert_array_equal(np.asarray(dtip_new), np.asarray(dtip_old))
    np.testing.assert_array_equal(np.asarray(dorgp_new), np.asarray(dorgp_old))


@pytest.mark.parametrize("label,factory", SCENARIOS, ids=[s[0] for s in SCENARIOS])
def test_helper_matches_inline_orgp_only(
    p_orgp_only: Phosphorus, label: str, factory
) -> None:
    """Same parity, but with use_TIP = False. dtip_dt is 0 in both."""
    kwargs = factory()

    dtip_new, dorgp_new, _components = p_orgp_only._change_with_components(**kwargs)
    dtip_old, dorgp_old = p_orgp_only._change_legacy_inline(**kwargs)

    np.testing.assert_array_equal(np.asarray(dtip_new), np.asarray(dtip_old))
    np.testing.assert_array_equal(np.asarray(dorgp_new), np.asarray(dorgp_old))


@pytest.mark.parametrize("label,factory", SCENARIOS, ids=[s[0] for s in SCENARIOS])
def test_components_dict_contains_all_registry_diagnostics(
    p: Phosphorus, label: str, factory
) -> None:
    kwargs = factory()
    _, _, components = p._change_with_components(**kwargs)
    assert set(components.keys()) == set(p.REGISTRY_DIAGNOSTICS)


def test_orgp_hydrolysis_alias_matches_legacy_attribute(p: Phosphorus) -> None:
    """``orgp_hydrolysis_rate`` (Appendix A name) and
    ``orgp_to_tip_hydrolysis_rate`` (legacy attribute consumed by
    ``test_phosphorus_v1_parity_v3.py``) must point at the same value
    after _change_with_components runs."""
    kwargs = _uniform_state()
    _, _, components = p._change_with_components(**kwargs)

    new_name = np.asarray(components["orgp_hydrolysis_rate"])
    legacy_attr = np.asarray(p.orgp_to_tip_hydrolysis_rate)
    np.testing.assert_array_equal(
        new_name, legacy_attr,
        err_msg=(
            "orgp_hydrolysis_rate (Appendix A) must alias "
            "self.orgp_to_tip_hydrolysis_rate (legacy)"
        ),
    )


def test_helper_zero_state_produces_finite_components(p: Phosphorus) -> None:
    kwargs = _zero_state()
    _, _, components = p._change_with_components(**kwargs)
    for name, value in components.items():
        arr = np.asarray(value)
        assert np.isfinite(arr).all(), (
            f"{name} contains non-finite values at zero state"
        )
