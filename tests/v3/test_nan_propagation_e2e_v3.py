"""Phase R-5 NaN propagation end-to-end test (review-findings 2026-05-04).

Goal: confirm a missing meteorology forcing (NaN ``wind_speed``,
``air_temperature``, or ``solar_radiation``) produces a visible NaN at
the end of a substep rather than a silent finite value. This guards
the M3 fix on ``richardson_number`` and ensures the surrounding chain
(Richardson stability function, latent + sensible flux, depth ramp,
rate cap) does not accidentally scrub the NaN.

Strategy: drive ``Temperature.run`` directly against a stub registry
pre-populated with one finite cell of forcings, then poison one
forcing variable with NaN per scenario and assert that the
post-substep ``water_temperature`` is NaN. A finite-baseline sanity
case confirms the harness does not produce NaN for the all-good path.

Stubs are inlined to avoid coupling to other test files in this
directory.
"""

from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np
import xarray as xr

from clearwater_modules_v3 import Temperature


# ---------------------------------------------------------------------------
# Stubs
# ---------------------------------------------------------------------------


class _StubRegistry:
    def __init__(self) -> None:
        self._data: dict[str, xr.DataArray] = {}

    def register(self, name: str, value) -> None:
        self._data[name] = value

    def get(self, name: str):
        return self._data[name]

    def get_at_time(self, name: str, time):
        if name not in self._data:
            raise KeyError(name)
        return self._data[name]

    def set_at_time(self, name: str, time, value) -> None:
        self._data[name] = value

    def get_variable(self, name: str):
        raise KeyError(name)

    def __contains__(self, name: str) -> bool:
        return name in self._data


def _make_temperature() -> Temperature:
    """Construct a Temperature whose ``__skip_first_time_step`` is False so
    the first ``run`` call actually exercises the kinetics path."""
    t = Temperature(wind_a=0.3, wind_b=1.5, wind_c=3.0)
    t.from_hotstart({})  # flips __skip_first_time_step to False
    return t


def _seed_single_cell_finite(registry: _StubRegistry) -> None:
    """Pre-populate a single-cell registry with all-finite forcings.
    Values pulled from the patterns used in
    ``tests/v3/test_wet_mask_scope_v3.py`` and
    ``tests/v3/test_5_tsm_calculations_v3.py``-style integration tests:
    realistic but unremarkable surface-water conditions."""
    registry.register(
        "water_temperature", xr.DataArray(np.array([20.0]))
    )
    registry.register(
        "wetted_surface_area", xr.DataArray(np.array([100.0]))
    )
    registry.register("volume", xr.DataArray(np.array([1000.0])))
    registry.register("cloudiness", xr.DataArray(np.array([0.0])))
    registry.register("air_temperature", xr.DataArray(np.array([25.0])))
    registry.register("solar_radiation", xr.DataArray(np.array([800.0])))
    registry.register("wind_speed", xr.DataArray(np.array([3.0])))
    registry.register(
        "atmospheric_pressure", xr.DataArray(np.array([1013.0]))
    )
    registry.register(
        "atmospheric_vapor_pressure", xr.DataArray(np.array([20.0]))
    )
    registry.register(
        "sediment_temperature", xr.DataArray(np.array([15.0]))
    )
    registry.register(
        "sediment_thickness", xr.DataArray(np.array([0.1]))
    )


def _run_one_substep(registry: _StubRegistry) -> float:
    t = _make_temperature()
    t.run(datetime(2026, 1, 1, 0, 0, 0), registry)
    return float(registry.get("water_temperature").values[0])


# ---------------------------------------------------------------------------
# Sanity baseline: all-finite forcings produce a finite output
# ---------------------------------------------------------------------------


def test_all_finite_forcings_produce_finite_output():
    """Sanity check for the harness: with every forcing finite, the
    post-substep ``water_temperature`` is finite. Pin only the
    finite-ness; the exact value is irrelevant to this test (kinetics
    coverage lives in test_5_tsm_calculations_v3.py)."""
    registry = _StubRegistry()
    _seed_single_cell_finite(registry)
    result = _run_one_substep(registry)
    assert np.isfinite(result), (
        f"Finite-baseline substep produced non-finite water_temperature "
        f"({result!r}); the harness or upstream kinetics is broken "
        f"independently of the NaN-propagation contract under test."
    )


# ---------------------------------------------------------------------------
# NaN forcing -> NaN output (one test per critical forcing)
# ---------------------------------------------------------------------------


def test_nan_wind_speed_propagates_to_water_temperature():
    """A NaN ``wind_speed`` poisons the Richardson number, the wind
    function, and through them both the latent and sensible heat fluxes.
    The M3 fix on ``richardson_number`` was designed to make this
    visible at the kinetics output. After one substep, the registry's
    ``water_temperature`` must be NaN."""
    registry = _StubRegistry()
    _seed_single_cell_finite(registry)
    registry.register("wind_speed", xr.DataArray(np.array([np.nan])))
    result = _run_one_substep(registry)
    assert np.isnan(result), (
        f"NaN wind_speed must propagate to NaN water_temperature "
        f"(M3 visible-defect contract); got {result!r}"
    )


def test_nan_air_temperature_propagates_to_water_temperature():
    """A NaN ``air_temperature`` enters via ``flux_atmospheric_longwave``
    (Swinbank polynomial in T_K), ``flux_sensible`` (T_air - T_water
    driving difference), and ``density_air``. After one substep, the
    registry's ``water_temperature`` must be NaN."""
    registry = _StubRegistry()
    _seed_single_cell_finite(registry)
    registry.register("air_temperature", xr.DataArray(np.array([np.nan])))
    result = _run_one_substep(registry)
    assert np.isnan(result), (
        f"NaN air_temperature must propagate to NaN water_temperature; "
        f"got {result!r}"
    )


def test_nan_solar_radiation_propagates_to_water_temperature():
    """A NaN ``solar_radiation`` adds directly to the net flux in
    ``flux_net``. After one substep, the registry's
    ``water_temperature`` must be NaN. Of the three forcings tested,
    this is the most direct path: solar_radiation is a summand in the
    energy balance and has no protective clamping."""
    registry = _StubRegistry()
    _seed_single_cell_finite(registry)
    registry.register("solar_radiation", xr.DataArray(np.array([np.nan])))
    result = _run_one_substep(registry)
    assert np.isnan(result), (
        f"NaN solar_radiation must propagate to NaN water_temperature; "
        f"got {result!r}"
    )
