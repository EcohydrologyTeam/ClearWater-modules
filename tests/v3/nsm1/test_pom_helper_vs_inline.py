"""Phase 7 helper-vs-inline parity for POM.

Asserts that ``POM._change_with_components`` and the shadow
``POM._change_legacy_inline`` produce bit-identical net rates
(mg/L/d). Also pins the Phase 7 cache-relocation: the
``self.pom_doc_source_rate`` attribute (consumed by Carbon via getattr)
is now set inside ``_change_with_components`` rather than inside
``rate()``.

Deleted in Phase 10 alongside its shadow per §11.3. Tolerance:
``rtol=0, atol=0``.
"""

from __future__ import annotations

from datetime import timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.pom import POM


@pytest.fixture
def p() -> POM:
    proc = POM(time_step=timedelta(minutes=5))
    proc.use_POC = False
    proc.use_floating_algae = False
    proc.use_benthic_algae = False
    proc.use_Algae = True
    proc.use_Balgae = True
    proc.floating_algae_process = None
    proc.benthic_algae_process = None
    return proc


N_CELLS = 5


def _da(values: np.ndarray) -> xr.DataArray:
    return xr.DataArray(values, dims="cell")


def _zero_state():
    z = np.zeros(N_CELLS)
    return dict(
        pom=_da(z),
        water_temperature=_da(z + 20.0),
        poc=_da(z),
        depth=_da(z + 1.0),
    )


def _uniform_state():
    return dict(
        pom=_da(np.full(N_CELLS, 2.0)),
        water_temperature=_da(np.full(N_CELLS, 20.0)),
        poc=_da(np.full(N_CELLS, 1.0)),
        depth=_da(np.full(N_CELLS, 1.0)),
    )


def _randomised_state(seed: int = 20260513):
    rng = np.random.default_rng(seed)
    return dict(
        pom=_da(rng.uniform(0.1, 10.0, N_CELLS)),
        water_temperature=_da(rng.uniform(5.0, 30.0, N_CELLS)),
        poc=_da(rng.uniform(0.1, 5.0, N_CELLS)),
        depth=_da(rng.uniform(0.3, 5.0, N_CELLS)),
    )


def _cold_water():
    return dict(
        pom=_da(np.full(N_CELLS, 2.0)),
        water_temperature=_da(np.full(N_CELLS, 0.5)),
        poc=_da(np.full(N_CELLS, 1.0)),
        depth=_da(np.full(N_CELLS, 1.0)),
    )


def _hot_water():
    return dict(
        pom=_da(np.full(N_CELLS, 2.0)),
        water_temperature=_da(np.full(N_CELLS, 35.0)),
        poc=_da(np.full(N_CELLS, 1.0)),
        depth=_da(np.full(N_CELLS, 1.0)),
    )


def _thin_depth():
    return dict(
        pom=_da(np.full(N_CELLS, 2.0)),
        water_temperature=_da(np.full(N_CELLS, 20.0)),
        poc=_da(np.full(N_CELLS, 1.0)),
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
def test_helper_matches_inline_bit_identical(
    p: POM, label: str, factory
) -> None:
    """``_change_with_components`` and ``_change_legacy_inline`` produce
    bit-equal net rates."""
    kwargs = factory()

    rate_new, _components = p._change_with_components(**kwargs)
    rate_old = p._change_legacy_inline(**kwargs)

    np.testing.assert_array_equal(
        np.asarray(rate_new), np.asarray(rate_old),
        err_msg=f"POM rate differs for scenario {label!r}",
    )


@pytest.mark.parametrize("label,factory", SCENARIOS, ids=[s[0] for s in SCENARIOS])
def test_components_dict_contains_all_registry_diagnostics(
    p: POM, label: str, factory
) -> None:
    kwargs = factory()
    _, components = p._change_with_components(**kwargs)
    assert set(components.keys()) == set(p.REGISTRY_DIAGNOSTICS)


def test_phase7_pom_doc_source_rate_set_by_helper(p: POM) -> None:
    """Phase 7 cache-relocation contract: ``self.pom_doc_source_rate``
    (consumed by Carbon via getattr) is set as a side effect of
    ``_change_with_components``, not just inside ``rate()``."""
    kwargs = _uniform_state()

    # Reset the cache to a sentinel so we can verify the helper
    # actually wrote to it.
    p.pom_doc_source_rate = -999.0

    p._change_with_components(**kwargs)

    cached = np.asarray(p.pom_doc_source_rate)
    assert np.all(cached != -999.0), (
        "_change_with_components did not set self.pom_doc_source_rate"
    )
    assert np.isfinite(cached).all()


def test_phase7_pom_doc_source_rate_matches_legacy_value(p: POM) -> None:
    """The Phase 7 cache-relocation must produce the same numerical
    value the pre-Phase-7 ``rate()`` did. Verified by computing both
    paths and comparing the resulting ``pom_doc_source_rate``."""
    kwargs = _uniform_state()

    p.pom_doc_source_rate = -999.0
    p._change_with_components(**kwargs)
    helper_value = np.asarray(p.pom_doc_source_rate).copy()

    p.pom_doc_source_rate = -999.0
    p._change_legacy_inline(**kwargs)
    shadow_value = np.asarray(p.pom_doc_source_rate).copy()

    np.testing.assert_array_equal(
        helper_value, shadow_value,
        err_msg=(
            "Phase 7 cache-relocation broke pom_doc_source_rate parity "
            "vs the pre-Phase-7 path"
        ),
    )


def test_helper_zero_state_produces_finite_components(p: POM) -> None:
    kwargs = _zero_state()
    _, components = p._change_with_components(**kwargs)
    for name, value in components.items():
        arr = np.asarray(value)
        assert np.isfinite(arr).all(), (
            f"{name} contains non-finite values at zero state"
        )
