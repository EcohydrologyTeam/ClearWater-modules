"""Phase 9 helper-vs-inline parity for Alkalinity.

Deleted in Phase 10 alongside its shadow per §11.3. Tolerance:
``rtol=0, atol=0``.

Alkalinity's helper signature is unusual: it takes only ``depth`` as a
kwarg because the per-source / per-sink fluxes (nitrification,
denitrification, algal growth / respiration) are read from sibling-
process caches via getattr inside the sub-flux helpers. The test
exercises bit-identical parity by varying ``depth`` and the sibling-
cache values via direct attribute injection.
"""

from __future__ import annotations

from datetime import timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.alkalinity import Alkalinity


@pytest.fixture
def alk() -> Alkalinity:
    proc = Alkalinity(time_step=timedelta(minutes=5))
    proc.use_NH4 = True
    proc.use_NO3 = True
    proc.use_Algae = True
    proc.use_Balgae = True
    # No sibling-process wiring: the sub-flux helpers will all return 0
    # (use_nitrogen / use_floating_algae / use_benthic_algae default to
    # False).
    return proc


N_CELLS = 5


def _da(values: np.ndarray) -> xr.DataArray:
    return xr.DataArray(values, dims="cell")


def _depth(label: str) -> xr.DataArray:
    if label == "thin":
        return _da(np.full(N_CELLS, 0.05))
    if label == "deep":
        return _da(np.full(N_CELLS, 10.0))
    if label == "uniform":
        return _da(np.full(N_CELLS, 1.0))
    if label == "randomised":
        rng = np.random.default_rng(20260513)
        return _da(rng.uniform(0.3, 5.0, N_CELLS))
    raise ValueError(f"unknown depth label {label!r}")


DEPTH_SCENARIOS = ["uniform", "thin", "deep", "randomised"]


@pytest.mark.parametrize("depth_label", DEPTH_SCENARIOS)
def test_helper_matches_inline_no_siblings(
    alk: Alkalinity, depth_label: str
) -> None:
    """Helper and shadow produce bit-equal rates when no sibling
    processes are wired up (all sub-fluxes degrade to 0)."""
    depth = _depth(depth_label)

    rate_new, _components = alk._change_with_components(depth=depth)
    rate_old = alk._change_legacy_inline(depth=depth)

    np.testing.assert_array_equal(
        np.asarray(rate_new), np.asarray(rate_old),
        err_msg=f"Alkalinity rate differs for depth {depth_label!r}",
    )


@pytest.mark.parametrize("depth_label", DEPTH_SCENARIOS)
def test_components_dict_contains_all_registry_diagnostics(
    alk: Alkalinity, depth_label: str
) -> None:
    depth = _depth(depth_label)
    _, components = alk._change_with_components(depth=depth)
    assert set(components.keys()) == set(alk.REGISTRY_DIAGNOSTICS)


def test_helper_zero_state_produces_finite_components(alk: Alkalinity) -> None:
    depth = _da(np.full(N_CELLS, 1.0))
    _, components = alk._change_with_components(depth=depth)
    for name, value in components.items():
        arr = np.asarray(value)
        assert np.isfinite(arr).all()


def test_legacy_attribute_aliases_match_components(alk: Alkalinity) -> None:
    """The four legacy attribute names
    (``alk_nitrification_rate``, ``alk_denitrification_rate``,
    ``alk_benthic_algae_growth_rate``, ``alk_benthic_algae_respiration_rate``)
    must equal the corresponding Appendix A components after the helper
    runs. Pinned so a future refactor cannot silently drop the
    back-compat aliases that ``test_alkalinity_v1_parity_v3.py`` and
    ``test_alkalinity_tier1.py`` read."""
    depth = _da(np.full(N_CELLS, 1.0))
    _, components = alk._change_with_components(depth=depth)

    np.testing.assert_array_equal(
        np.asarray(alk.alk_nitrification_rate),
        np.asarray(components["alk_nitrification_sink_rate"]),
    )
    np.testing.assert_array_equal(
        np.asarray(alk.alk_denitrification_rate),
        np.asarray(components["alk_denitrification_source_rate"]),
    )
    np.testing.assert_array_equal(
        np.asarray(alk.alk_benthic_algae_growth_rate),
        np.asarray(components["alk_balgae_growth_rate"]),
    )
    np.testing.assert_array_equal(
        np.asarray(alk.alk_benthic_algae_respiration_rate),
        np.asarray(components["alk_balgae_respiration_rate"]),
    )


def test_helper_matches_inline_with_nitrogen_sibling_caches(
    alk: Alkalinity,
) -> None:
    """Inject mock nitrogen-process caches and verify bit-identical
    parity through the sub-flux helper path."""

    class _MockNitrogen:
        nitrification_flux_rate = _da(np.full(N_CELLS, 0.5))
        denitrification_flux_rate = _da(np.full(N_CELLS, 0.3))

    alk.use_nitrogen = True
    alk.nitrogen_process = _MockNitrogen()
    depth = _da(np.full(N_CELLS, 1.0))

    rate_new, _components = alk._change_with_components(depth=depth)
    rate_old = alk._change_legacy_inline(depth=depth)

    np.testing.assert_array_equal(
        np.asarray(rate_new), np.asarray(rate_old),
        err_msg="rate differs with Nitrogen sibling caches wired",
    )
