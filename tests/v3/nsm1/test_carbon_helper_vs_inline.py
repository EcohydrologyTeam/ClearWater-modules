"""Phase 2 helper-vs-inline parity for Carbon.

Asserts that ``Carbon._change_with_components`` and the shadow
``Carbon._change_legacy_inline`` produce bit-identical ``(d_poc, d_doc,
d_dic)`` across a parametrised state/forcing matrix.

This file is **deleted in Phase 10** of
``design/clearwater_modules_v3_nsm1_pattern_alignment_specification.md``
together with the ``_change_legacy_inline`` shadow method, after the
final end-to-end baseline parity passes (§11.3).

Tolerance: ``rtol=0, atol=0`` — bit-identical only. A single-bit
difference fails the test (§11.6 refactor-discipline rules).
"""

from __future__ import annotations

from datetime import timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.carbon import Carbon


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def carbon() -> Carbon:
    """Standalone Carbon Process with v3 DEFAULTS. No sibling Processes
    are wired in — ``use_floating_algae`` etc. default to ``False`` and
    the corresponding source/sink terms collapse to 0. This isolates
    the bare integrator arithmetic and is sufficient for the
    helper-vs-inline parity (which only checks that both methods
    produce the same arithmetic, not that they exercise every
    coupling path)."""
    return Carbon(time_step=timedelta(minutes=5))


N_CELLS = 5


def _da(values: np.ndarray) -> xr.DataArray:
    return xr.DataArray(values, dims="cell")


# Parametrised forcing matrix. Each entry is a (label, factory)
# pair where the factory returns a dict of kwargs for
# ``_change_with_components`` / ``_change_legacy_inline``. The factory
# is called once per test so each test sees a fresh array (avoids
# accidental in-place mutation).
def _zero_state():
    z = np.zeros(N_CELLS)
    return dict(
        poc=_da(z), doc=_da(z), dic=_da(z),
        t_water_c=_da(z + 20.0), depth=_da(z + 1.0), dox=_da(z + 8.0),
    )


def _uniform_state():
    return dict(
        poc=_da(np.full(N_CELLS, 0.5)),
        doc=_da(np.full(N_CELLS, 2.0)),
        dic=_da(np.full(N_CELLS, 10.0)),
        t_water_c=_da(np.full(N_CELLS, 20.0)),
        depth=_da(np.full(N_CELLS, 1.0)),
        dox=_da(np.full(N_CELLS, 8.0)),
    )


def _randomised_state(seed: int = 20260513):
    rng = np.random.default_rng(seed)
    return dict(
        poc=_da(rng.uniform(0.01, 5.0, N_CELLS)),
        doc=_da(rng.uniform(0.5, 15.0, N_CELLS)),
        dic=_da(rng.uniform(5.0, 50.0, N_CELLS)),
        t_water_c=_da(rng.uniform(5.0, 30.0, N_CELLS)),
        depth=_da(rng.uniform(0.3, 5.0, N_CELLS)),
        dox=_da(rng.uniform(2.0, 12.0, N_CELLS)),
    )


def _cold_water():
    return dict(
        poc=_da(np.full(N_CELLS, 1.0)),
        doc=_da(np.full(N_CELLS, 2.0)),
        dic=_da(np.full(N_CELLS, 10.0)),
        t_water_c=_da(np.full(N_CELLS, 0.5)),
        depth=_da(np.full(N_CELLS, 1.0)),
        dox=_da(np.full(N_CELLS, 12.0)),
    )


def _hot_water():
    return dict(
        poc=_da(np.full(N_CELLS, 1.0)),
        doc=_da(np.full(N_CELLS, 2.0)),
        dic=_da(np.full(N_CELLS, 10.0)),
        t_water_c=_da(np.full(N_CELLS, 35.0)),
        depth=_da(np.full(N_CELLS, 1.0)),
        dox=_da(np.full(N_CELLS, 4.0)),
    )


def _thin_depth():
    return dict(
        poc=_da(np.full(N_CELLS, 1.0)),
        doc=_da(np.full(N_CELLS, 2.0)),
        dic=_da(np.full(N_CELLS, 10.0)),
        t_water_c=_da(np.full(N_CELLS, 20.0)),
        depth=_da(np.full(N_CELLS, 0.05)),
        dox=_da(np.full(N_CELLS, 8.0)),
    )


def _low_dox():
    return dict(
        poc=_da(np.full(N_CELLS, 1.0)),
        doc=_da(np.full(N_CELLS, 2.0)),
        dic=_da(np.full(N_CELLS, 10.0)),
        t_water_c=_da(np.full(N_CELLS, 20.0)),
        depth=_da(np.full(N_CELLS, 1.0)),
        dox=_da(np.full(N_CELLS, 0.1)),
    )


SCENARIOS = [
    ("zero_state", _zero_state),
    ("uniform_state", _uniform_state),
    ("randomised_state", _randomised_state),
    ("cold_water", _cold_water),
    ("hot_water", _hot_water),
    ("thin_depth", _thin_depth),
    ("low_dox", _low_dox),
]


@pytest.mark.parametrize("label,factory", SCENARIOS, ids=[s[0] for s in SCENARIOS])
def test_helper_matches_inline_bit_identical(
    carbon: Carbon, label: str, factory
) -> None:
    """``_change_with_components`` and ``_change_legacy_inline`` must
    produce ``(d_poc, d_doc, d_dic)`` that are byte-equal."""
    kwargs = factory()

    d_poc_new, d_doc_new, d_dic_new, _components = (
        carbon._change_with_components(**kwargs)
    )
    d_poc_old, d_doc_old, d_dic_old = carbon._change_legacy_inline(**kwargs)

    np.testing.assert_array_equal(
        np.asarray(d_poc_new), np.asarray(d_poc_old),
        err_msg=f"d_poc differs for scenario {label!r}",
    )
    np.testing.assert_array_equal(
        np.asarray(d_doc_new), np.asarray(d_doc_old),
        err_msg=f"d_doc differs for scenario {label!r}",
    )
    np.testing.assert_array_equal(
        np.asarray(d_dic_new), np.asarray(d_dic_old),
        err_msg=f"d_dic differs for scenario {label!r}",
    )


@pytest.mark.parametrize("label,factory", SCENARIOS, ids=[s[0] for s in SCENARIOS])
def test_components_dict_contains_all_registry_diagnostics(
    carbon: Carbon, label: str, factory
) -> None:
    """The ``components`` dict returned by ``_change_with_components``
    must contain every name in ``Carbon.REGISTRY_DIAGNOSTICS`` and
    nothing else. This pins the pattern G contract: every name in the
    class tuple has a value to write to the registry."""
    kwargs = factory()
    _, _, _, components = carbon._change_with_components(**kwargs)
    assert set(components.keys()) == set(carbon.REGISTRY_DIAGNOSTICS)


def test_helper_zero_state_produces_finite_components(carbon: Carbon) -> None:
    """Zero state should yield finite components (no NaN/inf from
    division-by-zero in any sub-rate). Also confirms the components
    dict survives the zero-state path."""
    kwargs = _zero_state()
    _, _, _, components = carbon._change_with_components(**kwargs)
    for name, value in components.items():
        arr = np.asarray(value)
        assert np.isfinite(arr).all(), (
            f"{name} contains non-finite values at zero state"
        )
