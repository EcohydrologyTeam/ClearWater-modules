"""v3 sediment heat-exchange tests.

Covers two related fixes (review findings doc, C3 and C10):

1. **C3** — sediment_diffusivity unit + default + docstring fix. v2/v3
   inherited a transcription error (default 0.0061, docstring "m^2/s")
   that does not match the Fortran TSM source
   (default 0.0432, units m^2/day). This test pins the corrected
   default and verifies the flux magnitude is consistent with the
   Fortran formula.

2. **C10** — dynamic sediment temperature evolution. The Fortran TSM
   updates ``T_sed`` each substep via
   ``dT_sed/dt = alphas / (0.5 * h2^2) * (T_water - T_sed)``. The
   Python ports (v1, v2) dropped this update, breaking energy
   conservation between the water and sediment heat reservoirs.
   v3 restores the Fortran behavior; tests below verify:

   - The relaxation time constant matches ``tau = 0.5 * h2^2 / alphas``.
   - The exchange is energy-conservative: heat lost by sediment per
     unit area equals heat gained by water per unit area, per substep.
   - ``evolve_sediment_temperature=False`` is a backward-compat opt-out.

References:
    Fortran source:
    HEC-RAS-WQ/RAS-1D-WQ/Kinetics Libraries/{TemperatureEnergyBudget,TemperatureEquilibrium}/Source files/modTemperature.f90
"""

from __future__ import annotations

from datetime import timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v3 import Temperature


# ---------------------------------------------------------------------------
# C3: sediment_diffusivity default + units
# ---------------------------------------------------------------------------


def test_sediment_diffusivity_default_matches_fortran():
    """The default sediment_diffusivity must match the Fortran TSM
    ``alphas = 0.0432`` (m^2/day). Pinning this catches future
    transcription errors."""
    t = Temperature(wind_a=0.3, wind_b=1.5, wind_c=3.0)
    assert t.sediment_diffusivity == pytest.approx(0.0432, abs=1e-12)


def test_sediment_density_and_specific_heat_match_fortran():
    """Fortran defaults: pb=1600, Cps=1673, h2=0.1 (h2 is a registry
    variable so not a constructor default). v3 should match."""
    t = Temperature(wind_a=0.3, wind_b=1.5, wind_c=3.0)
    assert t.sediment_density == 1600.0
    assert t.sediment_specific_heat == 1673.0


def test_flux_sediment_matches_fortran_formula_at_realistic_inputs():
    """For shallow-river-typical inputs, the v3 flux should match the
    Fortran formula
    ``q_sediment = pb * Cps * alphas / 0.5 / h2 * (T_sed - T_water) / 86400``.

    With pb=1600, Cps=1673, alphas=0.0432 m^2/day, h2=0.1 m,
    T_sed - T_water = 10 K::

        q_sediment = 1600 * 1673 * 0.0432 / 0.5 / 0.1 * 10 / 86400
                   ~= 267.7 W/m^2

    This is the canonical magnitude for a strong T_sed > T_water
    gradient in shallow water. v3's pre-fix default (0.0061) gave
    ~38 W/m^2, ~7x too low.
    """
    t = Temperature(wind_a=0.3, wind_b=1.5, wind_c=3.0)
    flux = t.flux_sediment(
        water_temperature=15.0,
        sediment_temperature=25.0,
        sediment_thickness=0.1,
    )
    expected = 1600.0 * 1673.0 * 0.0432 / 0.5 / 0.1 * 10.0 / 86400.0
    assert float(flux) == pytest.approx(expected, rel=1e-12)
    # Magnitude sanity: the corrected flux is ~267 W/m^2 at this gradient.
    # A 7x-too-low pre-fix flux would be ~38 W/m^2.
    assert float(flux) > 200.0


# ---------------------------------------------------------------------------
# C10: dynamic sediment temperature evolution
# ---------------------------------------------------------------------------


def test_sediment_temperature_change_zero_when_water_equals_sediment():
    """At thermal equilibrium (T_water == T_sed) the sediment
    temperature update is zero."""
    t = Temperature(wind_a=0.3, wind_b=1.5, wind_c=3.0)
    delta = t.sediment_temperature_change(
        water_temperature=20.0,
        sediment_temperature=20.0,
        sediment_thickness=0.1,
    )
    assert float(delta) == pytest.approx(0.0, abs=1e-15)


def test_sediment_temperature_change_relaxation_time_constant():
    """The Fortran formula ``dT_sed/dt = alphas / (0.5 * h2^2) * (T_w - T_s)``
    has time constant ``tau = 0.5 * h2^2 / alphas``. With
    alphas = 0.0432 m^2/day and h2 = 0.1 m,
    tau = 0.5 * 0.01 / 0.0432 ~ 0.1157 day ~ 9999.9 seconds ~ 2.78 hours.

    Per-substep: dT_sed = (T_w - T_s) * dt_seconds / tau_seconds.
    With T_w - T_s = 10 K and dt = 300 s::

        dT_sed = 10 * 300 / 9999.9 ~= 0.30 K

    Verifies both the formula and the unit chain.
    """
    t = Temperature(
        wind_a=0.3, wind_b=1.5, wind_c=3.0,
        time_step=timedelta(minutes=5),
    )
    delta = t.sediment_temperature_change(
        water_temperature=25.0,
        sediment_temperature=15.0,
        sediment_thickness=0.1,
    )
    tau_seconds = 0.5 * 0.1**2 / (0.0432 / 86400.0)
    expected = (25.0 - 15.0) * 300.0 / tau_seconds
    assert float(delta) == pytest.approx(expected, rel=1e-12)


def test_sediment_temperature_change_sign_drives_toward_water():
    """When T_water > T_sed, sediment warms (positive dT)."""
    t = Temperature(wind_a=0.3, wind_b=1.5, wind_c=3.0)
    delta_warming = t.sediment_temperature_change(
        water_temperature=25.0,
        sediment_temperature=15.0,
        sediment_thickness=0.1,
    )
    assert float(delta_warming) > 0.0

    delta_cooling = t.sediment_temperature_change(
        water_temperature=15.0,
        sediment_temperature=25.0,
        sediment_thickness=0.1,
    )
    assert float(delta_cooling) < 0.0


def test_water_sediment_energy_conservation_per_substep():
    """The water-sediment heat exchange must be energy-conservative.

    Per unit area per substep:
        dE_water    = q_sediment * dt
                    = pb * Cps * alphas / (0.5 * h2) * (T_s - T_w) * dt / 86400
        dE_sediment = pb * Cps * h2 * dT_sed
                    = pb * Cps * h2 * alphas / (0.5 * h2^2) * (T_w - T_s) * dt / 86400
                    = pb * Cps * alphas / (0.5 * h2) * (T_w - T_s) * dt / 86400
                    = -dE_water

    The two must sum to zero (heat lost by sediment = heat gained by water).
    """
    t = Temperature(
        wind_a=0.3, wind_b=1.5, wind_c=3.0,
        time_step=timedelta(minutes=5),
    )
    T_water, T_sed, h2 = 15.0, 25.0, 0.1
    dt = t.time_step_seconds

    flux_w_into_water_W_per_m2 = t.flux_sediment(
        water_temperature=T_water,
        sediment_temperature=T_sed,
        sediment_thickness=h2,
    )
    dE_water_J_per_m2 = float(flux_w_into_water_W_per_m2) * dt

    delta_T_sed = t.sediment_temperature_change(
        water_temperature=T_water,
        sediment_temperature=T_sed,
        sediment_thickness=h2,
    )
    dE_sediment_J_per_m2 = (
        t.sediment_density * t.sediment_specific_heat * h2 * float(delta_T_sed)
    )

    # Heat gained by water plus heat gained by sediment should sum to zero
    # (energy flows from sediment -> water in this case; T_sed > T_water).
    total = dE_water_J_per_m2 + dE_sediment_J_per_m2
    # Both terms are O(80,000 J/m^2) for these inputs.
    # Conservation should hold to floating-point precision (rel ~ 1e-12).
    reference_magnitude = max(abs(dE_water_J_per_m2), abs(dE_sediment_J_per_m2))
    assert abs(total) / reference_magnitude < 1e-12, (
        f"Energy non-conservation: dE_water={dE_water_J_per_m2:.6e}, "
        f"dE_sediment={dE_sediment_J_per_m2:.6e}, total={total:.6e}, "
        f"rel={abs(total)/reference_magnitude:.3e}"
    )


def test_evolve_sediment_temperature_default_is_true():
    t = Temperature(wind_a=0.3, wind_b=1.5, wind_c=3.0)
    assert t.evolve_sediment_temperature is True


def test_evolve_sediment_temperature_can_be_disabled():
    t = Temperature(
        wind_a=0.3, wind_b=1.5, wind_c=3.0,
        evolve_sediment_temperature=False,
    )
    assert t.evolve_sediment_temperature is False
    # The method itself still computes a non-zero value when called
    # directly — the gate is in `run()`, not in the method. This
    # makes the method usable for diagnostic / sensitivity-analysis
    # callers that want the rate without applying the update.
    delta = t.sediment_temperature_change(
        water_temperature=25.0,
        sediment_temperature=15.0,
        sediment_thickness=0.1,
    )
    assert float(delta) > 0.0


# ---------------------------------------------------------------------------
# Vectorized inputs
# ---------------------------------------------------------------------------


def test_sediment_temperature_change_works_on_multi_cell_dataarray():
    """v3 must support cell-shaped inputs since the kernel writes to
    a multi-cell registry."""
    t = Temperature(
        wind_a=0.3, wind_b=1.5, wind_c=3.0,
        time_step=timedelta(minutes=5),
    )
    # Three cells: equilibrium, +10 K gradient, -10 K gradient (mirror).
    T_w = xr.DataArray(np.array([20.0, 25.0, 15.0]))
    T_s = xr.DataArray(np.array([20.0, 15.0, 25.0]))
    h2 = xr.DataArray(np.array([0.1, 0.1, 0.1]))
    delta = t.sediment_temperature_change(
        water_temperature=T_w,
        sediment_temperature=T_s,
        sediment_thickness=h2,
    )
    arr = np.asarray(delta)
    assert arr.shape == (3,)
    assert arr[0] == pytest.approx(0.0, abs=1e-15)   # equilibrium
    assert arr[1] > 0.0                              # T_w > T_s, sediment warms
    assert arr[2] < 0.0                              # T_w < T_s, sediment cools
    # Mirror symmetry: cell 1 and cell 2 are exactly opposite (+/-10 K).
    assert arr[1] == pytest.approx(-arr[2], rel=1e-12)
