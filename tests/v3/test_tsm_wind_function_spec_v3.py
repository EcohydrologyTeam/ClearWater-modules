"""v3 TSM wind-function spec implementation tests.

Implements the test plan from
``design/clearwater_modules_v3_tsm_wind_function_specification.md``:

1. Default-change regression: ``Temperature()`` constructs with
   ``wind_c = 2.0``.
2. Validator behavior at ``wind_c in {1.0, 2.0}`` (no warning),
   ``wind_c in {1.5, 3.0}`` (warns), and ``wind_c in {0.0, -1.0, 3.5}``
   (raises ``ValueError``).
3. Wind-input-height correction: at default 2.0 m the transform is a
   no-op; at ``wind_input_height = 10.0`` with ``surface_z0 = 0.001``
   the effective wind is ``0.825 * raw``, the log-law factor used by
   CE-QUAL-W2.
4. Wind shelter: scalar default 1.0 is a no-op; scalar 0.5 reduces
   the wind term ``b * W^c`` by ``0.5^c``; the registry forcing
   ``wind_shelter_coefficient`` overrides the constructor scalar.
5. Composition: ``wind_input_height = 10`` and ``wind_shelter = 0.65``
   compose multiplicatively as ``0.825 * 0.65`` applied to wind_speed
   (then squared at ``c = 2``).
6. Energy-conservation regression: handled by the existing test suite;
   no test added here.
"""
from __future__ import annotations

import warnings
from datetime import datetime, timedelta
from typing import Any

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.temperature import Temperature


class _FakeRegistry:
    """Minimal registry surrogate for ``__contains__`` + ``get_at_time``."""

    def __init__(self, values: dict[str, Any]) -> None:
        self._values = values

    def __contains__(self, name: str) -> bool:
        return name in self._values

    def get_at_time(self, name: str, time: Any) -> Any:
        return self._values[name]


# ---------------------------------------------------------------------------
# 1. Default-change regression
# ---------------------------------------------------------------------------


def test_default_constructor_wind_c_is_2():
    """``Temperature()`` with no kwargs constructs with ``wind_c = 2.0``."""
    t = Temperature()
    assert t.wind_c == 2.0


def test_default_constructor_wind_input_height_is_2_m():
    t = Temperature()
    assert t.wind_input_height == 2.0


def test_default_constructor_surface_z0_is_001_m():
    t = Temperature()
    assert t.surface_z0 == 0.001


def test_default_constructor_wind_shelter_is_unity():
    t = Temperature()
    assert t.wind_shelter == 1.0


# ---------------------------------------------------------------------------
# 2. Validator behavior
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("c", [1.0, 2.0])
def test_wind_c_supported_values_no_warning(c: float) -> None:
    """Supported values {1.0, 2.0} should construct silently."""
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        Temperature(wind_c=c)


@pytest.mark.parametrize("c", [0.5, 1.5, 2.5, 3.0])
def test_wind_c_in_range_but_nonstandard_warns(c: float) -> None:
    """Values inside (0, 3] but outside {1.0, 2.0} warn (with W2/QUAL2K
    consensus message)."""
    with pytest.warns(UserWarning, match="CE-QUAL-W2 explicitly defaults to CFW = 2.0"):
        Temperature(wind_c=c)


@pytest.mark.parametrize("c", [0.0, -1.0, 3.5, float("inf")])
def test_wind_c_out_of_range_raises(c: float) -> None:
    """Values outside (0.0, 3.0] raise ValueError."""
    with pytest.raises(ValueError, match="wind_c must be in"):
        Temperature(wind_c=c)


def test_wind_c_3_is_back_compat_opt_in_only_warns():
    """``c = 3.0`` is allowed at the upper bound for back-compat with
    explicit opt-ins; emits the validator warning but does not raise."""
    with pytest.warns(UserWarning, match="CE-QUAL-W2"):
        t = Temperature(wind_c=3.0)
    assert t.wind_c == 3.0


# ---------------------------------------------------------------------------
# 3. Wind-input-height correction
# ---------------------------------------------------------------------------


def test_wind_input_height_2m_is_no_op():
    """Default ``wind_input_height = 2.0`` produces effective_wind == raw
    for any wind_speed."""
    t = Temperature()
    for w in (0.0, 1.0, 3.7, 10.0):
        assert t._compute_effective_wind(w) == w


def test_wind_input_height_10m_to_2m_log_law_factor():
    """``wind_input_height = 10.0`` with ``surface_z0 = 0.001`` produces
    ``effective_wind = 0.825 * raw_wind``, the W2 log-law factor."""
    t = Temperature(wind_input_height=10.0, surface_z0=0.001)
    raw = 4.0
    effective = t._compute_effective_wind(raw)
    expected_factor = np.log(2.0 / 0.001) / np.log(10.0 / 0.001)
    assert effective == pytest.approx(raw * expected_factor)
    # Sanity check: factor is ~0.825 per the W2 reference.
    assert 0.82 < expected_factor < 0.83


def test_wind_input_height_quadratic_attenuation_at_c2():
    """At ``c = 2`` and ``wind_input_height = 10``, the wind term
    ``b * W^c`` attenuates by ``0.825^2 ~= 0.68``."""
    raw = 4.0  # m/s
    t_no_corr = Temperature(wind_a=0.3, wind_b=1.5, wind_c=2.0)
    t_at_10m = Temperature(wind_a=0.3, wind_b=1.5, wind_c=2.0,
                           wind_input_height=10.0, surface_z0=0.001)

    f_no_corr = t_no_corr.wind_function(raw, 1.0)
    f_at_10m = t_at_10m.wind_function(raw, 1.0)

    # The 'a' term is unaffected by W; only 'b * W^c' attenuates.
    a_term = 0.3 / 1e6
    b_term_no_corr = f_no_corr - a_term
    b_term_at_10m = f_at_10m - a_term
    factor = (np.log(2.0 / 0.001) / np.log(10.0 / 0.001)) ** 2
    assert b_term_at_10m == pytest.approx(b_term_no_corr * factor)


# ---------------------------------------------------------------------------
# 4. Wind shelter
# ---------------------------------------------------------------------------


def test_wind_shelter_unity_is_no_op():
    """Default ``wind_shelter = 1.0`` produces effective_wind == raw."""
    t = Temperature()
    assert t._compute_effective_wind(4.0) == 4.0


def test_wind_shelter_scalar_attenuation_at_c2():
    """At ``c = 2`` and ``wind_shelter = 0.5``, the wind term
    ``b * W^c`` attenuates by ``0.5^2 = 0.25``."""
    raw = 4.0
    t_open = Temperature(wind_a=0.3, wind_b=1.5, wind_c=2.0)
    t_sheltered = Temperature(wind_a=0.3, wind_b=1.5, wind_c=2.0,
                              wind_shelter=0.5)

    f_open = t_open.wind_function(raw, 1.0)
    f_sheltered = t_sheltered.wind_function(raw, 1.0)

    a_term = 0.3 / 1e6
    b_open = f_open - a_term
    b_sheltered = f_sheltered - a_term
    assert b_sheltered == pytest.approx(b_open * 0.25)


def test_registry_shelter_overrides_scalar():
    """When the cached registry shelter is set, it takes precedence
    over the constructor scalar."""
    t = Temperature(wind_shelter=0.5)
    # Simulate run() having found the per-cell forcing in the registry.
    t._cached_shelter = 0.7

    raw = 4.0
    effective = t._compute_effective_wind(raw)
    # 4.0 * 0.7 = 2.8 (per-cell wins, NOT 4.0 * 0.5 = 2.0)
    assert effective == pytest.approx(2.8)


def test_registry_shelter_per_cell_array():
    """Per-cell ``wind_shelter_coefficient`` as an xr.DataArray broadcasts
    against wind_speed."""
    t = Temperature(wind_shelter=1.0)
    t._cached_shelter = xr.DataArray([1.0, 0.7, 0.3])

    raw = xr.DataArray([4.0, 4.0, 4.0])
    effective = t._compute_effective_wind(raw)
    np.testing.assert_allclose(effective.values, [4.0, 2.8, 1.2])


def test_unset_cached_shelter_falls_back_to_scalar():
    """``self._cached_shelter is None`` (no registry forcing) falls back
    to the constructor scalar."""
    t = Temperature(wind_shelter=0.5)
    assert t._cached_shelter is None
    assert t._compute_effective_wind(4.0) == 2.0  # 4.0 * 0.5


# ---------------------------------------------------------------------------
# 5. Composition (shelter * height_factor)
# ---------------------------------------------------------------------------


def test_composition_shelter_and_height_correction():
    """Combined ``wind_input_height = 10``, ``wind_shelter = 0.65``,
    ``c = 2`` is mass-balance-equivalent to a single multiplicative
    factor ``0.825 * 0.65`` applied to wind_speed and squared in the
    wind term."""
    raw = 4.0
    t = Temperature(wind_a=0.3, wind_b=1.5, wind_c=2.0,
                    wind_input_height=10.0, surface_z0=0.001,
                    wind_shelter=0.65)

    effective = t._compute_effective_wind(raw)
    height_factor = np.log(2.0 / 0.001) / np.log(10.0 / 0.001)
    expected_effective = raw * 0.65 * height_factor
    assert effective == pytest.approx(expected_effective)

    # Wind term scales by (shelter * height_factor)^c.
    f_open = Temperature(wind_a=0.3, wind_b=1.5, wind_c=2.0).wind_function(raw, 1.0)
    f_combined = t.wind_function(raw, 1.0)
    a_term = 0.3 / 1e6
    composite_factor_squared = (0.65 * height_factor) ** 2
    assert (f_combined - a_term) == pytest.approx(
        (f_open - a_term) * composite_factor_squared
    )


# ---------------------------------------------------------------------------
# Validator coverage for the new params
# ---------------------------------------------------------------------------


def test_wind_input_height_must_be_positive():
    with pytest.raises(ValueError, match="wind_input_height must be > 0"):
        Temperature(wind_input_height=0.0)
    with pytest.raises(ValueError, match="wind_input_height must be > 0"):
        Temperature(wind_input_height=-1.0)


def test_surface_z0_must_be_positive():
    with pytest.raises(ValueError, match="surface_z0 must be > 0"):
        Temperature(surface_z0=0.0)
    with pytest.raises(ValueError, match="surface_z0 must be > 0"):
        Temperature(surface_z0=-0.001)


def test_surface_z0_must_be_less_than_input_height():
    with pytest.raises(ValueError, match="surface_z0 must be strictly less"):
        Temperature(wind_input_height=2.0, surface_z0=2.0)
    with pytest.raises(ValueError, match="surface_z0 must be strictly less"):
        Temperature(wind_input_height=2.0, surface_z0=3.0)


def test_wind_shelter_must_be_positive():
    with pytest.raises(ValueError, match="wind_shelter must be > 0"):
        Temperature(wind_shelter=0.0)
    with pytest.raises(ValueError, match="wind_shelter must be > 0"):
        Temperature(wind_shelter=-0.5)


def test_wind_shelter_above_one_warns():
    with pytest.warns(UserWarning, match="wind_shelter = 1.5 is greater than 1.0"):
        Temperature(wind_shelter=1.5)


def test_wind_shelter_unity_no_warning():
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        Temperature(wind_shelter=1.0)


# ---------------------------------------------------------------------------
# Composition with wind_function (integration with Richardson)
# ---------------------------------------------------------------------------


def test_wind_function_default_no_op():
    """At default constructor params, ``wind_function(W, Ri)`` reduces
    to ``Ri * (a + b * W^c) / 1e6`` exactly (no internal transforms)."""
    t = Temperature()
    raw = 4.0
    ri = 0.95
    expected = ri * (0.3 / 1e6 + 1.5 / 1e6 * raw ** 2.0)
    assert t.wind_function(raw, ri) == pytest.approx(expected)


def test_wind_function_applies_combined_transforms():
    """Non-default params: ``wind_function`` produces the formula
    evaluated at the effective wind, with Richardson factor."""
    t = Temperature(wind_a=0.3, wind_b=1.5, wind_c=2.0,
                    wind_input_height=10.0, surface_z0=0.001,
                    wind_shelter=0.65)
    raw = 4.0
    ri = 0.95
    height_factor = np.log(2.0 / 0.001) / np.log(10.0 / 0.001)
    eff = raw * 0.65 * height_factor
    expected = ri * (0.3 / 1e6 + 1.5 / 1e6 * eff ** 2.0)
    assert t.wind_function(raw, ri) == pytest.approx(expected)
