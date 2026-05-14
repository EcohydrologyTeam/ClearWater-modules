"""Phase 3 helper-vs-inline parity for DOX.

Asserts that ``DOX._change_with_components`` and the shadow
``DOX._change_legacy_inline`` produce bit-identical ``(delta_dox,
rate)`` across a parametrised state/forcing matrix.

This file is **deleted in Phase 10** of
``design/clearwater_modules_v3_nsm1_pattern_alignment_specification.md``
together with the ``_change_legacy_inline`` shadow method, after the
final end-to-end baseline parity passes (§11.3).

The DOX scenario matrix specifically covers (per Phase 3 deliverable):

- the hypoxic regime (``DOX → 0``) where SOD-Monod attenuation kicks in;
- the supersaturated regime (``DOX > DOX_sat``) where atmospheric
  reaeration reverses sign;
- DOX-saturated regime (no reaeration driver);
- standard zero / uniform / randomised / cold / hot / thin-depth.

Tolerance: ``rtol=0, atol=0`` — bit-identical only (§11.6).
"""

from __future__ import annotations

from datetime import timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.dox import DOX


@pytest.fixture(scope="module")
def dox() -> DOX:
    """Standalone DOX Process with v3 DEFAULTS. No sibling Processes
    wired in: the algae / nitrogen / carbon / cbod sub-fluxes collapse
    to 0 and only atmospheric reaeration + SOD + nitrification-via-NH4-
    forcing contribute. This isolates the bare integrator arithmetic."""
    return DOX(time_step=timedelta(minutes=5))


N_CELLS = 5


def _da(values: np.ndarray) -> xr.DataArray:
    return xr.DataArray(values, dims="cell")


def _zero_state():
    z = np.zeros(N_CELLS)
    return dict(
        dox=_da(z + 1e-12),  # avoid 0/0 in SOD Monod ratio
        t_water_c=_da(z + 20.0),
        depth=_da(z + 1.0),
        ammonium=_da(z),
        pressure_mb=1013.25,
    )


def _uniform_state():
    return dict(
        dox=_da(np.full(N_CELLS, 8.0)),
        t_water_c=_da(np.full(N_CELLS, 20.0)),
        depth=_da(np.full(N_CELLS, 1.0)),
        ammonium=_da(np.full(N_CELLS, 0.5)),
        pressure_mb=1013.25,
    )


def _randomised_state(seed: int = 20260513):
    rng = np.random.default_rng(seed)
    return dict(
        dox=_da(rng.uniform(0.5, 14.0, N_CELLS)),
        t_water_c=_da(rng.uniform(5.0, 30.0, N_CELLS)),
        depth=_da(rng.uniform(0.3, 5.0, N_CELLS)),
        ammonium=_da(rng.uniform(0.0, 2.0, N_CELLS)),
        pressure_mb=1013.25,
    )


def _cold_water():
    return dict(
        dox=_da(np.full(N_CELLS, 12.0)),
        t_water_c=_da(np.full(N_CELLS, 0.5)),
        depth=_da(np.full(N_CELLS, 1.0)),
        ammonium=_da(np.full(N_CELLS, 0.5)),
        pressure_mb=1013.25,
    )


def _hot_water():
    return dict(
        dox=_da(np.full(N_CELLS, 4.0)),
        t_water_c=_da(np.full(N_CELLS, 35.0)),
        depth=_da(np.full(N_CELLS, 1.0)),
        ammonium=_da(np.full(N_CELLS, 0.5)),
        pressure_mb=1013.25,
    )


def _thin_depth():
    return dict(
        dox=_da(np.full(N_CELLS, 8.0)),
        t_water_c=_da(np.full(N_CELLS, 20.0)),
        depth=_da(np.full(N_CELLS, 0.05)),
        ammonium=_da(np.full(N_CELLS, 0.5)),
        pressure_mb=1013.25,
    )


def _hypoxic():
    """DOX → 0: SOD Monod attenuation reduces sediment sink toward zero;
    nitrification ``(1 - exp(-KNR*DOX))`` term ≈ 0; reaeration drives
    DOX upward. Tests the hypoxic-regime parity explicitly."""
    return dict(
        dox=_da(np.full(N_CELLS, 0.05)),
        t_water_c=_da(np.full(N_CELLS, 20.0)),
        depth=_da(np.full(N_CELLS, 1.0)),
        ammonium=_da(np.full(N_CELLS, 1.0)),
        pressure_mb=1013.25,
    )


def _supersaturated():
    """DOX > DOX_sat: atmospheric reaeration term reverses sign
    (``ka_tc * (O2sat - DOX)`` becomes negative)."""
    return dict(
        dox=_da(np.full(N_CELLS, 15.0)),  # well above typical 9 mg/L at 20 C
        t_water_c=_da(np.full(N_CELLS, 20.0)),
        depth=_da(np.full(N_CELLS, 1.0)),
        ammonium=_da(np.full(N_CELLS, 0.5)),
        pressure_mb=1013.25,
    )


def _high_nh4():
    """High ammonium: nitrification sink dominates the budget. Cross-
    couples DOX with the NH4 input forcing."""
    return dict(
        dox=_da(np.full(N_CELLS, 8.0)),
        t_water_c=_da(np.full(N_CELLS, 20.0)),
        depth=_da(np.full(N_CELLS, 1.0)),
        ammonium=_da(np.full(N_CELLS, 5.0)),  # 10x typical
        pressure_mb=1013.25,
    )


SCENARIOS = [
    ("zero_state", _zero_state),
    ("uniform_state", _uniform_state),
    ("randomised_state", _randomised_state),
    ("cold_water", _cold_water),
    ("hot_water", _hot_water),
    ("thin_depth", _thin_depth),
    ("hypoxic", _hypoxic),
    ("supersaturated", _supersaturated),
    ("high_nh4", _high_nh4),
]


@pytest.mark.parametrize("label,factory", SCENARIOS, ids=[s[0] for s in SCENARIOS])
def test_helper_matches_inline_bit_identical(
    dox: DOX, label: str, factory
) -> None:
    """``_change_with_components`` and ``_change_legacy_inline`` must
    produce ``(delta_dox, rate)`` that are byte-equal."""
    kwargs = factory()

    delta_new, rate_new, _components = dox._change_with_components(**kwargs)
    delta_old, rate_old = dox._change_legacy_inline(**kwargs)

    np.testing.assert_array_equal(
        np.asarray(delta_new), np.asarray(delta_old),
        err_msg=f"delta_dox differs for scenario {label!r}",
    )
    np.testing.assert_array_equal(
        np.asarray(rate_new), np.asarray(rate_old),
        err_msg=f"rate differs for scenario {label!r}",
    )


@pytest.mark.parametrize("label,factory", SCENARIOS, ids=[s[0] for s in SCENARIOS])
def test_components_dict_contains_all_registry_diagnostics(
    dox: DOX, label: str, factory
) -> None:
    """The ``components`` dict returned by ``_change_with_components``
    must contain every name in ``DOX.REGISTRY_DIAGNOSTICS``."""
    kwargs = factory()
    _, _, components = dox._change_with_components(**kwargs)
    assert set(components.keys()) == set(dox.REGISTRY_DIAGNOSTICS)


def test_sod_and_dox_sod_rate_are_aliases(dox: DOX) -> None:
    """``sod_rate`` and ``dox_sod_rate`` are documented aliases for the
    same volumetric SOD sink (mg-O2/L/d). Phase 3 closeout documents
    this explicitly; this test pins it."""
    kwargs = _uniform_state()
    _, _, components = dox._change_with_components(**kwargs)
    np.testing.assert_array_equal(
        np.asarray(components["sod_rate"]),
        np.asarray(components["dox_sod_rate"]),
    )


def test_helper_hypoxic_sod_attenuates_toward_zero(dox: DOX) -> None:
    """As DOX → 0, the SOD sink (under Monod attenuation by DOX/(DOX+KsSOD))
    must approach zero. Verifies the hypoxic-regime physics in the
    helper output."""
    kwargs_normal = _uniform_state()
    _, _, c_normal = dox._change_with_components(**kwargs_normal)

    kwargs_low = _uniform_state()
    kwargs_low["dox"] = _da(np.full(N_CELLS, 0.01))
    _, _, c_low = dox._change_with_components(**kwargs_low)

    sod_normal = np.asarray(c_normal["dox_sod_rate"])
    sod_low = np.asarray(c_low["dox_sod_rate"])
    # Hypoxic SOD is strictly smaller than normoxic SOD at every cell.
    assert np.all(np.abs(sod_low) <= np.abs(sod_normal) + 1e-12)


def test_helper_zero_state_produces_finite_components(dox: DOX) -> None:
    """All sub-rate components must be finite at near-zero state (no
    NaN/inf from divisions or Monod-zero edges)."""
    kwargs = _zero_state()
    _, _, components = dox._change_with_components(**kwargs)
    for name, value in components.items():
        arr = np.asarray(value)
        assert np.isfinite(arr).all(), (
            f"{name} contains non-finite values at near-zero state"
        )
