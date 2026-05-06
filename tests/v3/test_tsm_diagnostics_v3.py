"""Regression tests for v3 TSM per-component flux diagnostics and the
equilibrium-temperature Newton-Raphson solver (audit 2026-05-05 open
questions 2 and 3).

Covered:

1. ``Temperature.flux_components`` returns a dict with the seven
   Fortran-A pathway-output keys and ``q_net`` matches the manual
   composition with the magnitudes-only sign convention.
2. After ``Temperature.run`` writes ``water_temperature``, the seven
   diagnostic values are cached on the process instance (``self.q_*``).
3. ``Temperature.run`` writes a diagnostic to the registry only when
   the user has pre-registered the variable (matching N2's
   ``total_dissolved_gas`` opt-in pattern).
4. ``Temperature.equilibrium_temperature`` converges to the
   physically-correct equilibrium for a uniform meteorological
   forcing, identified by ``q_net(T_eq) ~= 0``.

Refs:
    design/clearwater_modules_v3_tsm_audit_2026-05-05.md
    open questions 2 (TeqC) and 3 (per-component fluxes).
"""
from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3.processes.temperature import Temperature


# ---------------------------------------------------------------------------
# Stub registry (re-uses the pattern from test_tsm_sediment_v3.py)
# ---------------------------------------------------------------------------


class _StubRegistry:
    """Minimal in-memory registry compatible with ``Temperature.run``.

    Supports ``__contains__``, ``get_at_time``, ``set_at_time``,
    ``register``. Used to drive ``run`` without spinning up a full
    ``Model`` object.
    """

    def __init__(self) -> None:
        self._data: dict[str, xr.DataArray] = {}

    def __contains__(self, key: str) -> bool:
        return key in self._data

    def register(self, key: str, value: xr.DataArray) -> None:
        self._data[key] = value

    def get_at_time(self, key: str, time: datetime) -> xr.DataArray:
        return self._data[key]

    def set_at_time(self, key: str, time: datetime, value: xr.DataArray) -> None:
        self._data[key] = value


def _arr(value: float, n: int = 1) -> xr.DataArray:
    return xr.DataArray(np.full(n, value, dtype=float), dims=("nface",))


def _scalar(value) -> float:
    """Coerce a 0-D / 1-D xr.DataArray or numpy array to a Python float."""
    if isinstance(value, xr.DataArray):
        arr = value.values.reshape(-1)
    else:
        arr = np.asarray(value).reshape(-1)
    return float(arr[0])


def _populate_registry(reg: _StubRegistry, n: int = 1) -> None:
    """Populate the typical TSM forcing variables for ``run``."""
    reg.register("water_temperature", _arr(20.0, n))
    reg.register("wetted_surface_area", _arr(100.0, n))
    reg.register("volume", _arr(100.0, n))  # depth = 1 m -> ramp inactive
    reg.register("cloudiness", _arr(0.3, n))
    reg.register("air_temperature", _arr(22.0, n))
    reg.register("solar_radiation", _arr(400.0, n))
    reg.register("wind_speed", _arr(3.0, n))
    reg.register("atmospheric_pressure", _arr(1013.0, n))
    reg.register("atmospheric_vapor_pressure", _arr(15.0, n))
    reg.register("sediment_temperature", _arr(18.0, n))
    reg.register("sediment_thickness", _arr(0.1, n))


# ---------------------------------------------------------------------------
# (1) flux_components dict shape and net composition
# ---------------------------------------------------------------------------


def test_flux_components_returns_seven_keys_plus_net() -> None:
    """``flux_components`` returns the seven Fortran-A pathway outputs."""
    t = Temperature()
    components = t.flux_components(
        water_temperature=20.0,
        cloudiness=0.3,
        air_temperature=22.0,
        solar_flux=400.0,
        wind_speed=3.0,
        atmospheric_pressure=1013.0,
        atmospheric_vapor_pressure=15.0,
        sediment_temperature=18.0,
        sediment_thickness=0.1,
    )
    assert set(components.keys()) == {
        "q_sensible",
        "q_latent",
        "q_longwave_up",
        "q_longwave_down",
        "q_solar",
        "q_sediment",
        "q_net",
    }


def test_flux_components_q_net_matches_manual_composition() -> None:
    """``q_net = sensible + solar + sediment + LW_down - LW_up - latent``.

    Pin the magnitudes-only sign convention by reconstructing the net
    flux from the per-component values returned in the same dict.
    """
    t = Temperature()
    components = t.flux_components(
        water_temperature=20.0,
        cloudiness=0.3,
        air_temperature=22.0,
        solar_flux=400.0,
        wind_speed=3.0,
        atmospheric_pressure=1013.0,
        atmospheric_vapor_pressure=15.0,
        sediment_temperature=18.0,
        sediment_thickness=0.1,
    )
    expected_net = (
        float(components["q_sensible"])
        + float(components["q_solar"])
        + float(components["q_sediment"])
        + float(components["q_longwave_down"])
        - float(components["q_longwave_up"])
        - float(components["q_latent"])
    )
    np.testing.assert_allclose(
        float(components["q_net"]), expected_net, rtol=1e-12
    )


def test_flux_components_q_solar_passthrough() -> None:
    """``q_solar`` is the input ``solar_flux`` unchanged (not attenuated)."""
    t = Temperature()
    solar_in = 537.0
    components = t.flux_components(
        water_temperature=20.0,
        cloudiness=0.3,
        air_temperature=22.0,
        solar_flux=solar_in,
        wind_speed=3.0,
        atmospheric_pressure=1013.0,
        atmospheric_vapor_pressure=15.0,
        sediment_temperature=18.0,
        sediment_thickness=0.1,
    )
    assert float(components["q_solar"]) == solar_in


def test_flux_net_matches_flux_components_q_net() -> None:
    """The thin-wrapper ``flux_net`` matches ``flux_components['q_net']``."""
    t = Temperature()
    inputs = dict(
        water_temperature=20.0,
        cloudiness=0.3,
        air_temperature=22.0,
        solar_flux=400.0,
        wind_speed=3.0,
        atmospheric_pressure=1013.0,
        atmospheric_vapor_pressure=15.0,
        sediment_temperature=18.0,
        sediment_thickness=0.1,
    )
    np.testing.assert_allclose(
        float(t.flux_net(**inputs)),
        float(t.flux_components(**inputs)["q_net"]),
        rtol=1e-12,
    )


# ---------------------------------------------------------------------------
# (2) Per-component caching on the process instance
# ---------------------------------------------------------------------------


def test_run_caches_q_components_on_process() -> None:
    """``Temperature.run`` caches ``self.q_*`` after the substep."""
    t = Temperature(time_step=timedelta(minutes=5))
    # First substep is a no-op (v1-coupling skip); call run twice.
    reg = _StubRegistry()
    _populate_registry(reg)
    now = datetime(2026, 5, 5, 12, 0, 0)
    t.run(now, reg)
    t.run(now + timedelta(minutes=5), reg)

    # All seven attributes exist after a real substep.
    for name in (
        "q_sensible",
        "q_latent",
        "q_longwave_up",
        "q_longwave_down",
        "q_solar",
        "q_sediment",
        "q_net",
    ):
        assert hasattr(t, name), f"Temperature.run did not cache self.{name}"


# ---------------------------------------------------------------------------
# (3) Optional registry write — the N2 pattern
# ---------------------------------------------------------------------------


def test_run_writes_q_diagnostics_only_when_registered() -> None:
    """A registered diagnostic key receives the per-substep flux value;
    an unregistered key produces no registry change.
    """
    t = Temperature(time_step=timedelta(minutes=5))
    reg = _StubRegistry()
    _populate_registry(reg)

    # Pre-register two of the seven diagnostic keys.
    reg.register("q_net", _arr(np.nan))
    reg.register("q_latent", _arr(np.nan))

    now = datetime(2026, 5, 5, 12, 0, 0)
    t.run(now, reg)  # skip
    t.run(now + timedelta(minutes=5), reg)

    # Registered diagnostics are filled with finite values.
    assert np.isfinite(_scalar(reg.get_at_time("q_net", now))), (
        "q_net was registered but Temperature.run did not write a finite value"
    )
    assert np.isfinite(_scalar(reg.get_at_time("q_latent", now))), (
        "q_latent was registered but Temperature.run did not write a finite value"
    )

    # Unregistered diagnostics: not added to registry as a side effect.
    assert "q_sensible" not in reg, (
        "q_sensible must NOT be auto-registered when not pre-registered "
        "by the user"
    )
    assert "q_solar" not in reg
    assert "q_sediment" not in reg
    assert "q_longwave_up" not in reg
    assert "q_longwave_down" not in reg


# ---------------------------------------------------------------------------
# (4) equilibrium_temperature Newton-Raphson
# ---------------------------------------------------------------------------


def test_equilibrium_temperature_zero_q_net_at_solution() -> None:
    """``equilibrium_temperature`` returns ``T_eq`` such that
    ``q_net(T_eq) ~= 0`` under the same forcing.

    Direct check of the Newton-Raphson root: substitute the returned
    ``T_eq`` into ``flux_components`` with the same sediment
    temperature held during the solve, then assert ``|q_net|`` is
    small. Bound: the loop's per-iterate tolerance is 0.01 K, and
    ``|d_qnet/dT|`` for typical surface conditions is ~5-15 W/m^2/K,
    so the residual ``|q_net|`` should be below ~0.5 W/m^2 at
    convergence.
    """
    t = Temperature()
    forcing = dict(
        cloudiness=0.3,
        air_temperature=20.0,
        solar_flux=400.0,
        wind_speed=3.0,
        atmospheric_pressure=1013.0,
        atmospheric_vapor_pressure=15.0,
        sediment_temperature=20.0,
        sediment_thickness=0.1,
    )
    # Use a tight tolerance so the residual is well under 1 W/m^2.
    teq_c = float(
        t.equilibrium_temperature(
            max_iterations=20,
            tolerance_kelvin=1e-6,
            **forcing,
        )
    )

    # Re-evaluate q_net at the equilibrium with the SAME sediment
    # temperature used during the solve.
    components_at_teq = t.flux_components(
        water_temperature=teq_c,
        **forcing,
    )
    assert abs(float(components_at_teq["q_net"])) < 1.0, (
        f"|q_net(T_eq)| = {abs(float(components_at_teq['q_net'])):.4g} "
        f"W/m^2 (T_eq = {teq_c:.4f} C) -- equilibrium not found"
    )


def test_equilibrium_temperature_warmer_when_solar_higher() -> None:
    """Higher solar flux raises the equilibrium temperature.

    Sanity / monotonicity check on the Newton-Raphson convergence and
    the analytic derivatives that drive it.
    """
    t = Temperature()
    common = dict(
        cloudiness=0.3,
        air_temperature=20.0,
        wind_speed=3.0,
        atmospheric_pressure=1013.0,
        atmospheric_vapor_pressure=15.0,
        sediment_temperature=20.0,
        sediment_thickness=0.1,
    )
    teq_low = float(t.equilibrium_temperature(solar_flux=200.0, **common))
    teq_high = float(t.equilibrium_temperature(solar_flux=800.0, **common))
    assert teq_high > teq_low + 1.0, (
        f"Higher solar flux did not raise T_eq materially "
        f"(low={teq_low:.2f} C, high={teq_high:.2f} C)"
    )


def test_equilibrium_temperature_writes_to_registry_when_registered() -> None:
    """``Temperature.run`` writes ``equilibrium_temperature`` only when
    the user has pre-registered it. Otherwise the Newton-Raphson is
    skipped (off the hot path)."""
    t = Temperature(time_step=timedelta(minutes=5))
    reg = _StubRegistry()
    _populate_registry(reg)

    # Pre-register the diagnostic.
    reg.register("equilibrium_temperature", _arr(np.nan))

    now = datetime(2026, 5, 5, 12, 0, 0)
    t.run(now, reg)  # skip
    t.run(now + timedelta(minutes=5), reg)

    teq = _scalar(reg.get_at_time("equilibrium_temperature", now))
    assert np.isfinite(teq), (
        "equilibrium_temperature was registered but Temperature.run did "
        "not write a finite value"
    )
    # For 22 C air, 400 W/m^2 solar, light wind, T_eq sits in a
    # plausible band around the air temperature.
    assert 10.0 < teq < 50.0, f"T_eq = {teq:.2f} C is implausibly extreme"
