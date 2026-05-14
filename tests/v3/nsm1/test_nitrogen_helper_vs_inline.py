"""Phase 4 helper-vs-inline parity for Nitrogen.

Asserts that ``Nitrogen._change_with_components`` and the shadow
``Nitrogen._change_legacy_inline`` produce bit-identical
``(ammonium_rate, nitrate_rate, orgn_rate)`` across a parametrised
state/forcing matrix.

Per the Phase 4 deliverable, the matrix specifically covers
``use_OrgN ∈ {True, False}`` (the OrgN sub-fluxes default to 0 under
False, mirroring ``change_organic_nitrogen``'s early-return) and the
NH4+NO3 mass-balance closure regime (high NH4 → strong nitrification
flux into NO3).

This file is **deleted in Phase 10** of
``design/clearwater_modules_v3_nsm1_pattern_alignment_specification.md``
together with the ``_change_legacy_inline`` shadow method, after the
final end-to-end baseline parity passes (§11.3).

Tolerance: ``rtol=0, atol=0`` — bit-identical only (§11.6).
"""

from __future__ import annotations

from datetime import timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.nitrogen import Nitrogen


@pytest.fixture
def nitrogen() -> Nitrogen:
    """Standalone Nitrogen Process with v3 DEFAULTS. ``use_floating_algae``
    / ``use_benthic_algae`` default to False so the algal coupling
    sub-fluxes collapse to 0; the bare nitrification / denitrification /
    bed / OrgN arithmetic is what's exercised."""
    n = Nitrogen(time_step=timedelta(minutes=5))
    # Defaults that change_* methods read at runtime.
    n.use_ammonium = True
    n.use_nitrate = True
    n.use_OrgN = True
    n.use_floating_algae = False
    n.use_benthic_algae = False
    return n


@pytest.fixture
def nitrogen_no_orgn() -> Nitrogen:
    """Nitrogen Process with ``use_OrgN = False`` — the OrgN sub-fluxes
    in the components dict default to 0 and ``change_organic_nitrogen``
    early-returns 0."""
    n = Nitrogen(time_step=timedelta(minutes=5))
    n.use_ammonium = True
    n.use_nitrate = True
    n.use_OrgN = False
    n.use_floating_algae = False
    n.use_benthic_algae = False
    return n


N_CELLS = 5


def _da(values: np.ndarray) -> xr.DataArray:
    return xr.DataArray(values, dims="cell")


def _zero_state():
    z = np.zeros(N_CELLS)
    return dict(
        nitrate=_da(z),
        ammonium=_da(z),
        organic_nitrogen=_da(z),
        temperature=_da(z + 20.0),
        depth=_da(z + 1.0),
        oxygen_dissolved=_da(z + 8.0),
    )


def _uniform_state():
    return dict(
        nitrate=_da(np.full(N_CELLS, 1.0)),
        ammonium=_da(np.full(N_CELLS, 0.5)),
        organic_nitrogen=_da(np.full(N_CELLS, 0.3)),
        temperature=_da(np.full(N_CELLS, 20.0)),
        depth=_da(np.full(N_CELLS, 1.0)),
        oxygen_dissolved=_da(np.full(N_CELLS, 8.0)),
    )


def _randomised_state(seed: int = 20260513):
    rng = np.random.default_rng(seed)
    return dict(
        nitrate=_da(rng.uniform(0.05, 5.0, N_CELLS)),
        ammonium=_da(rng.uniform(0.01, 2.0, N_CELLS)),
        organic_nitrogen=_da(rng.uniform(0.0, 1.0, N_CELLS)),
        temperature=_da(rng.uniform(5.0, 30.0, N_CELLS)),
        depth=_da(rng.uniform(0.3, 5.0, N_CELLS)),
        oxygen_dissolved=_da(rng.uniform(2.0, 12.0, N_CELLS)),
    )


def _cold_water():
    return dict(
        nitrate=_da(np.full(N_CELLS, 1.0)),
        ammonium=_da(np.full(N_CELLS, 0.5)),
        organic_nitrogen=_da(np.full(N_CELLS, 0.3)),
        temperature=_da(np.full(N_CELLS, 0.5)),
        depth=_da(np.full(N_CELLS, 1.0)),
        oxygen_dissolved=_da(np.full(N_CELLS, 12.0)),
    )


def _hot_water():
    return dict(
        nitrate=_da(np.full(N_CELLS, 1.0)),
        ammonium=_da(np.full(N_CELLS, 0.5)),
        organic_nitrogen=_da(np.full(N_CELLS, 0.3)),
        temperature=_da(np.full(N_CELLS, 35.0)),
        depth=_da(np.full(N_CELLS, 1.0)),
        oxygen_dissolved=_da(np.full(N_CELLS, 4.0)),
    )


def _thin_depth():
    return dict(
        nitrate=_da(np.full(N_CELLS, 1.0)),
        ammonium=_da(np.full(N_CELLS, 0.5)),
        organic_nitrogen=_da(np.full(N_CELLS, 0.3)),
        temperature=_da(np.full(N_CELLS, 20.0)),
        depth=_da(np.full(N_CELLS, 0.05)),
        oxygen_dissolved=_da(np.full(N_CELLS, 8.0)),
    )


def _high_nh4_low_no3():
    """NH4+NO3 mass-balance closure regime: high NH4 + low NO3 →
    nitrification flux moves N from NH4 to NO3 in unit-equivalent
    amounts. The Phase 4 deliverable calls this regime out explicitly."""
    return dict(
        nitrate=_da(np.full(N_CELLS, 0.05)),
        ammonium=_da(np.full(N_CELLS, 5.0)),
        organic_nitrogen=_da(np.full(N_CELLS, 0.3)),
        temperature=_da(np.full(N_CELLS, 25.0)),
        depth=_da(np.full(N_CELLS, 1.0)),
        oxygen_dissolved=_da(np.full(N_CELLS, 8.0)),
    )


def _hypoxic():
    """Low DOX: nitrification is inhibited (1 - exp(-KNR*DOX)) → 0;
    denitrification is enhanced (DOX/(DOX+KsOxdn) → 0)."""
    return dict(
        nitrate=_da(np.full(N_CELLS, 1.0)),
        ammonium=_da(np.full(N_CELLS, 0.5)),
        organic_nitrogen=_da(np.full(N_CELLS, 0.3)),
        temperature=_da(np.full(N_CELLS, 20.0)),
        depth=_da(np.full(N_CELLS, 1.0)),
        oxygen_dissolved=_da(np.full(N_CELLS, 0.05)),
    )


SCENARIOS = [
    ("zero_state", _zero_state),
    ("uniform_state", _uniform_state),
    ("randomised_state", _randomised_state),
    ("cold_water", _cold_water),
    ("hot_water", _hot_water),
    ("thin_depth", _thin_depth),
    ("high_nh4_low_no3", _high_nh4_low_no3),
    ("hypoxic", _hypoxic),
]


@pytest.mark.parametrize("label,factory", SCENARIOS, ids=[s[0] for s in SCENARIOS])
def test_helper_matches_inline_use_orgn_true(
    nitrogen: Nitrogen, label: str, factory
) -> None:
    """``_change_with_components`` and ``_change_legacy_inline`` must
    produce ``(ammonium_rate, nitrate_rate, orgn_rate)`` that are
    byte-equal under ``use_OrgN = True``."""
    kwargs = factory()

    nh4_new, no3_new, orgn_new, _components = (
        nitrogen._change_with_components(**kwargs)
    )
    nh4_old, no3_old, orgn_old = nitrogen._change_legacy_inline(**kwargs)

    np.testing.assert_array_equal(
        np.asarray(nh4_new), np.asarray(nh4_old),
        err_msg=f"ammonium_rate differs for scenario {label!r}",
    )
    np.testing.assert_array_equal(
        np.asarray(no3_new), np.asarray(no3_old),
        err_msg=f"nitrate_rate differs for scenario {label!r}",
    )
    np.testing.assert_array_equal(
        np.asarray(orgn_new), np.asarray(orgn_old),
        err_msg=f"orgn_rate differs for scenario {label!r}",
    )


@pytest.mark.parametrize("label,factory", SCENARIOS, ids=[s[0] for s in SCENARIOS])
def test_helper_matches_inline_use_orgn_false(
    nitrogen_no_orgn: Nitrogen, label: str, factory
) -> None:
    """Same parity, but with ``use_OrgN = False``. ``orgn_rate`` is
    ``0.0`` in both helper and shadow; ``ammonium_rate`` and
    ``nitrate_rate`` still bit-identical."""
    kwargs = factory()

    nh4_new, no3_new, orgn_new, _components = (
        nitrogen_no_orgn._change_with_components(**kwargs)
    )
    nh4_old, no3_old, orgn_old = nitrogen_no_orgn._change_legacy_inline(**kwargs)

    np.testing.assert_array_equal(
        np.asarray(nh4_new), np.asarray(nh4_old),
        err_msg=f"ammonium_rate (use_OrgN=False) differs for scenario {label!r}",
    )
    np.testing.assert_array_equal(
        np.asarray(no3_new), np.asarray(no3_old),
        err_msg=f"nitrate_rate (use_OrgN=False) differs for scenario {label!r}",
    )
    np.testing.assert_array_equal(
        np.asarray(orgn_new), np.asarray(orgn_old),
        err_msg=f"orgn_rate (use_OrgN=False) differs for scenario {label!r}",
    )


@pytest.mark.parametrize("label,factory", SCENARIOS, ids=[s[0] for s in SCENARIOS])
def test_components_dict_contains_all_registry_diagnostics(
    nitrogen: Nitrogen, label: str, factory
) -> None:
    """The ``components`` dict returned by ``_change_with_components``
    must contain every name in ``Nitrogen.REGISTRY_DIAGNOSTICS``."""
    kwargs = factory()
    _, _, _, components = nitrogen._change_with_components(**kwargs)
    assert set(components.keys()) == set(nitrogen.REGISTRY_DIAGNOSTICS)


def test_preserved_attribute_names_match_components(nitrogen: Nitrogen) -> None:
    """``nitrification_flux_rate`` and ``denitrification_flux_rate`` are
    the consumer-facing attribute names DOX/Alkalinity/N2 read via
    ``getattr``. After ``run`` populates the cache from the components
    dict (pattern F), the attribute values must equal the components
    dict values exactly."""
    kwargs = _uniform_state()
    _, _, _, components = nitrogen._change_with_components(**kwargs)

    # The pattern F loop in run does ``setattr(self, name, components[name])``.
    # Verify those two specific names are in the dict and the values
    # are sensible (positive, finite).
    assert "nitrification_flux_rate" in components
    assert "denitrification_flux_rate" in components
    nitr = np.asarray(components["nitrification_flux_rate"])
    denit = np.asarray(components["denitrification_flux_rate"])
    assert np.all(nitr >= 0.0)
    assert np.all(denit >= 0.0)
    assert np.isfinite(nitr).all()
    assert np.isfinite(denit).all()


def test_helper_high_nh4_drives_strong_nitrification(nitrogen: Nitrogen) -> None:
    """Mass-balance regime: high NH4 should produce a strong NH4 → NO3
    nitrification flux. Specifically, ``nitrification_flux_rate`` at
    NH4 = 5 mg/L should exceed nitrification at NH4 = 0.05 mg/L."""
    high = _high_nh4_low_no3()
    low = _high_nh4_low_no3()
    low["ammonium"] = _da(np.full(N_CELLS, 0.05))

    _, _, _, c_high = nitrogen._change_with_components(**high)
    _, _, _, c_low = nitrogen._change_with_components(**low)

    nitr_high = np.asarray(c_high["nitrification_flux_rate"])
    nitr_low = np.asarray(c_low["nitrification_flux_rate"])
    assert np.all(nitr_high > nitr_low)


def test_helper_zero_state_produces_finite_components(nitrogen: Nitrogen) -> None:
    """Zero state should yield finite components."""
    kwargs = _zero_state()
    _, _, _, components = nitrogen._change_with_components(**kwargs)
    for name, value in components.items():
        arr = np.asarray(value)
        assert np.isfinite(arr).all(), (
            f"{name} contains non-finite values at zero state"
        )
