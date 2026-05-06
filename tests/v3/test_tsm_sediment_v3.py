"""v3 sediment heat-exchange tests.

Covers three related fixes (review findings doc, C3 and C10; audit
2026-05-05 finding F2):

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

3. **F2** (audit 2026-05-05) — when the depth ramp or rate cap is
   active on the water-side delta, the same scaling factors must be
   applied to the sediment-side delta to preserve the per-cell
   water-sediment energy pair-cancellation invariant. Without this,
   shallow cells (depth < q_net_depth_ramp_ref) leak energy
   monotonically because the water absorbs ``q_sed * ramp`` of energy
   per substep while the sediment loses an unattenuated ``q_sed * dt``.
   v3 fixes this by extracting the ramp factor and rate-cap clip ratio
   from ``_temperature_change_with_factors`` and applying them to the
   sediment delta in ``Temperature.run``.

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


# ---------------------------------------------------------------------------
# F2 (audit 2026-05-05): depth ramp and rate cap apply symmetrically to
# the sediment-side delta to preserve per-cell water-sediment energy
# conservation.
# ---------------------------------------------------------------------------


class _StubRegistry:
    """Minimal VariableRegistry stand-in for the F2 end-to-end tests.
    Mirrors the pattern in ``tests/v3/test_wet_mask_scope_v3.py``."""

    def __init__(self, initial: dict) -> None:
        self._data = dict(initial)

    def get_at_time(self, name, time):
        if name not in self._data:
            raise KeyError(name)
        return self._data[name]

    def set_at_time(self, name, time, value) -> None:
        self._data[name] = value

    def get(self, name):
        return self._data[name]

    def register(self, name, value) -> None:
        self._data[name] = value

    def get_variable(self, name):
        raise KeyError(name)

    def __contains__(self, name) -> bool:
        return name in self._data


def _arr(value, shape=(1,)):
    return xr.DataArray(np.full(shape, value, dtype=float))


def _scalar(da):
    """Extract a scalar from a (1,)-shape DataArray / ndarray result."""
    arr = np.asarray(da)
    return float(arr.reshape(-1)[0])


def _method_kwargs(
    *,
    T_water=15.0,
    T_sed=25.0,
    surface_area=100.0,
    volume=1000.0,           # depth = 10 m by default (deep -> ramp = 1.0)
    h2=0.1,
):
    """Return kwargs for ``temperature_change`` /
    ``_temperature_change_with_factors`` calls. The atmospheric forcings
    are chosen so the surface-flux contribution is small relative to
    the sediment-flux contribution, making the sediment-vs-non-sediment
    isolation in the conservation tests robust to surface-flux realism."""
    return {
        "water_temperature": _arr(T_water),
        "surface_area": _arr(surface_area),
        "volume": _arr(volume),
        "cloudiness": _arr(0.0),
        "air_temperature": _arr(T_water),    # T_air = T_water minimizes sensible
        "solar_flux": _arr(0.0),              # no solar
        "wind_speed": _arr(1.0),
        "atmospheric_pressure": _arr(1013.0),
        "atmospheric_vapor_pressure": _arr(15.0),  # near-saturation at 15 C
        "sediment_temperature": _arr(T_sed),
        "sediment_thickness": _arr(h2),
    }


def _registry_dict(
    *,
    T_water=15.0,
    T_sed=25.0,
    surface_area=100.0,
    volume=1000.0,
    h2=0.1,
):
    """Return registry-key dict for ``Temperature.run``-driven tests.
    Uses the registry variable names (e.g. ``wetted_surface_area``,
    ``solar_radiation``) rather than the method kwarg names."""
    return {
        "water_temperature": _arr(T_water),
        "wetted_surface_area": _arr(surface_area),
        "volume": _arr(volume),
        "cloudiness": _arr(0.0),
        "air_temperature": _arr(T_water),
        "solar_radiation": _arr(0.0),
        "wind_speed": _arr(1.0),
        "atmospheric_pressure": _arr(1013.0),
        "atmospheric_vapor_pressure": _arr(15.0),
        "sediment_temperature": _arr(T_sed),
        "sediment_thickness": _arr(h2),
    }


def test_F2_helper_returns_correct_ramp_when_active():
    """``_temperature_change_with_factors`` must return ``ramp`` equal to
    ``min(1, depth / q_net_depth_ramp_ref)`` for shallow cells."""
    t = Temperature(
        wind_a=0.3, wind_b=1.5, wind_c=3.0,
        time_step=timedelta(minutes=5),
        q_net_depth_ramp_ref=0.3,           # default; ramp active when depth < 0.3 m
        dTdt_max_per_hour=float("inf"),     # cap disabled to isolate ramp
    )
    # Shallow cell: depth = 0.05 m, ramp_ref = 0.3 m -> ramp = 1/6
    inputs = _method_kwargs(volume=5.0, surface_area=100.0)  # depth = 0.05 m
    delta, ramp, clip_ratio, _components = t._temperature_change_with_factors(**inputs)
    expected_ramp = 0.05 / 0.3
    assert _scalar(ramp) == pytest.approx(expected_ramp, rel=1e-12)
    # Cap disabled -> clip_ratio is 1.0 everywhere.
    assert _scalar(clip_ratio) == pytest.approx(1.0)


def test_F2_helper_returns_ramp_one_when_disabled_or_deep():
    """When the ramp is disabled (q_net_depth_ramp_ref = 0) or the cell
    is deep enough, the helper returns ``ramp = 1.0`` (no scaling)."""
    # Disabled
    t_disabled = Temperature(
        wind_a=0.3, wind_b=1.5, wind_c=3.0,
        time_step=timedelta(minutes=5),
        q_net_depth_ramp_ref=0.0,
        dTdt_max_per_hour=float("inf"),
    )
    inputs = _method_kwargs(volume=5.0, surface_area=100.0)  # would be shallow if active
    _, ramp, _, _components = t_disabled._temperature_change_with_factors(**inputs)
    assert _scalar(ramp) == pytest.approx(1.0)

    # Deep cell with ramp enabled
    t_active = Temperature(
        wind_a=0.3, wind_b=1.5, wind_c=3.0,
        time_step=timedelta(minutes=5),
        q_net_depth_ramp_ref=0.3,
        dTdt_max_per_hour=float("inf"),
    )
    deep = _method_kwargs(volume=1000.0, surface_area=100.0)  # depth = 10 m, ramp clamped to 1
    _, ramp, _, _components = t_active._temperature_change_with_factors(**deep)
    assert _scalar(ramp) == pytest.approx(1.0)


def test_F2_helper_returns_clip_ratio_when_cap_fires():
    """When the rate cap clips the water-side delta, the helper must
    return ``clip_ratio = cap / |delta_unclipped|`` so the sediment
    side gets the same proportional scaling."""
    # Tight cap: 0.01 K/hr -> per-substep cap = 0.01 * 5/60 = 0.000833 K
    # Deep cell so ramp = 1.0; T_water and T_sed differ a lot so the
    # uncapped delta exceeds 0.000833 K easily.
    t = Temperature(
        wind_a=0.3, wind_b=1.5, wind_c=3.0,
        time_step=timedelta(minutes=5),
        q_net_depth_ramp_ref=0.0,           # disable ramp to isolate cap
        dTdt_max_per_hour=0.01,             # very tight cap
    )
    inputs = _method_kwargs(T_water=10.0, T_sed=30.0, volume=100.0, surface_area=100.0)
    delta, ramp, clip_ratio, _components = t._temperature_change_with_factors(**inputs)
    cap_value = 0.01 * (5.0 / 60.0)
    delta_v = _scalar(delta)
    clip_v = _scalar(clip_ratio)
    # The cap fired on this cell:
    assert abs(delta_v) == pytest.approx(cap_value, rel=1e-9)
    # clip_ratio is in (0, 1) when the cap fires:
    assert 0.0 < clip_v < 1.0
    # ramp should be 1.0 (disabled).
    assert _scalar(ramp) == pytest.approx(1.0)


def _energy_change_water_J_per_substep(t: Temperature, T_water, volume, delta_T_water):
    """Heat content change of the water column over one substep (J)."""
    rho_w = _scalar(t.water_density(T_water))
    cp_w = _scalar(t.water_specific_heat(T_water))
    return rho_w * cp_w * volume * delta_T_water


def _energy_change_sediment_J_per_substep(
    t: Temperature, surface_area, h2, delta_T_sed,
):
    """Heat content change of the active sediment layer over one substep (J)."""
    return (
        t.sediment_density
        * t.sediment_specific_heat
        * h2
        * surface_area
        * delta_T_sed
    )


def _run_one_substep_capture_deltas(
    t: Temperature, *, T_water, T_sed, surface_area, volume, h2,
):
    """Run ``Temperature.run`` once with the given initial conditions,
    return the per-cell delta in water and sediment temperatures.
    Returns ``(delta_T_water, delta_T_sed)``."""
    from datetime import datetime
    inputs = _registry_dict(
        T_water=T_water, T_sed=T_sed,
        surface_area=surface_area, volume=volume, h2=h2,
    )
    registry = _StubRegistry(inputs)
    # Skip-first-step is on by default; bypass for the test by hand.
    t._Temperature__skip_first_time_step = False  # type: ignore[attr-defined]
    t.run(datetime(2026, 1, 1, 0, 0, 0), registry)  # type: ignore[arg-type]
    new_T_water = _scalar(registry.get("water_temperature"))
    new_T_sed = _scalar(registry.get("sediment_temperature"))
    return new_T_water - T_water, new_T_sed - T_sed


def test_F2_water_sediment_conservation_under_depth_ramp():
    """When the depth ramp is active (shallow cell), the per-cell
    water-sediment energy pair-cancellation invariant must still hold.

    Method: run two simulations through ``Temperature.run`` with
    identical setup; one with ``use_sediment_temperature=True``
    (sediment exchange active) and one with ``use_sediment_temperature=False``
    (no sediment exchange). The difference in water-side energy change
    between the two runs is the sediment-flux contribution to water.
    Pair this with the sediment-side energy change from the active run;
    they must sum to zero.

    Pre-F2-fix the sediment-side delta was unscaled while the
    water-side delta was scaled by ``ramp``, breaking the invariant
    by a factor of ``(1 - ramp)`` in shallow cells.
    """
    # Shallow cell: depth = 0.05 m, ramp_ref = 0.3 m -> ramp ~ 0.167
    surface_area = 100.0
    volume = 5.0           # depth = 0.05 m
    h2 = 0.1
    T_water = 15.0
    T_sed = 25.0

    common = dict(
        wind_a=0.3, wind_b=1.5, wind_c=3.0,
        time_step=timedelta(minutes=5),
        q_net_depth_ramp_ref=0.3,
        dTdt_max_per_hour=float("inf"),  # disable cap to isolate ramp effect
    )

    # Run with sediment exchange ON
    t_on = Temperature(**common, use_sediment_temperature=True)
    dT_water_with_sed, dT_sed = _run_one_substep_capture_deltas(
        t_on,
        T_water=T_water, T_sed=T_sed,
        surface_area=surface_area, volume=volume, h2=h2,
    )

    # Run with sediment exchange OFF
    t_off = Temperature(**common, use_sediment_temperature=False)
    dT_water_no_sed, _ = _run_one_substep_capture_deltas(
        t_off,
        T_water=T_water, T_sed=T_sed,
        surface_area=surface_area, volume=volume, h2=h2,
    )

    # Water-side energy contribution from sediment exchange = E(with) - E(without)
    dE_water_from_sed = _energy_change_water_J_per_substep(
        t_on, T_water, volume, dT_water_with_sed - dT_water_no_sed,
    )
    # Sediment-side energy change
    dE_sediment = _energy_change_sediment_J_per_substep(
        t_on, surface_area, h2, dT_sed,
    )

    # Conservation: dE_water_from_sed + dE_sediment = 0
    total = dE_water_from_sed + dE_sediment
    reference = max(abs(dE_water_from_sed), abs(dE_sediment))
    assert reference > 0.0, "Test setup produced no sediment exchange"
    rel_error = abs(total) / reference
    assert rel_error < 1e-9, (
        f"Energy non-conservation under depth ramp (F2): "
        f"dE_water_from_sed={dE_water_from_sed:.6e} J, "
        f"dE_sediment={dE_sediment:.6e} J, total={total:.6e} J, "
        f"rel_error={rel_error:.3e}. "
        f"Expected pair-cancellation; pre-F2-fix breaks by ~(1-ramp) where "
        f"ramp = depth/ramp_ref = 0.05/0.3 ~ 0.167 in this fixture."
    )


def test_F2_sediment_delta_scales_by_ramp_and_clip_ratio_in_run():
    """Direct verification of the F2 fix: the sediment-side delta
    written to the registry by ``Temperature.run`` equals
    ``sediment_temperature_change * ramp * clip_ratio`` from
    ``_temperature_change_with_factors``.

    A direct end-to-end conservation test on the rate-cap-fired path
    is unsuitable because the cap is non-linear in the total
    water-side delta: subtracting two runs (with vs without sediment
    exchange) would not isolate the sediment-attributable share when
    the cap fires on both runs. Instead, we verify the *implementation*
    of the F2 scaling at the run() boundary, which combined with the
    depth-ramp conservation test and the unguarded-path conservation
    test (``test_water_sediment_energy_conservation_per_substep``)
    pins the full F2 contract: scaling is correctly applied, scaled
    deltas conserve energy by the unguarded-pair derivation, and the
    fix holds for both ramp and cap.
    """
    from datetime import datetime
    surface_area = 100.0
    volume = 1000.0
    h2 = 0.1
    T_water = 10.0
    T_sed = 30.0           # large gradient drives a large unguarded delta

    # Tight cap so it fires; ramp disabled to isolate the cap effect.
    t = Temperature(
        wind_a=0.3, wind_b=1.5, wind_c=3.0,
        time_step=timedelta(minutes=5),
        q_net_depth_ramp_ref=0.0,
        dTdt_max_per_hour=0.001,
        use_sediment_temperature=True,
        evolve_sediment_temperature=True,
    )
    # Compute the unguarded sediment delta independently.
    delta_T_sed_unguarded = _scalar(t.sediment_temperature_change(
        water_temperature=_arr(T_water),
        sediment_temperature=_arr(T_sed),
        sediment_thickness=_arr(h2),
    ))
    # Compute the ramp + clip_ratio that the helper would return.
    method_inputs = _method_kwargs(
        T_water=T_water, T_sed=T_sed,
        surface_area=surface_area, volume=volume, h2=h2,
    )
    _, ramp, clip_ratio, _components = t._temperature_change_with_factors(**method_inputs)
    expected_scaled_delta_T_sed = (
        delta_T_sed_unguarded * _scalar(ramp) * _scalar(clip_ratio)
    )

    # Run end-to-end through Temperature.run; capture the actual
    # sediment delta written to the registry.
    registry = _StubRegistry(_registry_dict(
        T_water=T_water, T_sed=T_sed,
        surface_area=surface_area, volume=volume, h2=h2,
    ))
    t._Temperature__skip_first_time_step = False  # type: ignore[attr-defined]
    t.run(datetime(2026, 1, 1, 0, 0, 0), registry)  # type: ignore[arg-type]
    actual_delta_T_sed = _scalar(registry.get("sediment_temperature")) - T_sed

    # The actual delta must match the predicted scaled delta to
    # floating-point precision. Pre-F2-fix the actual would equal
    # `delta_T_sed_unguarded` (no scaling applied to sediment).
    rel_error = abs(actual_delta_T_sed - expected_scaled_delta_T_sed) / (
        abs(expected_scaled_delta_T_sed) + 1e-30
    )
    assert rel_error < 1e-12, (
        f"F2 scaling mismatch: actual={actual_delta_T_sed:.6e}, "
        f"expected={expected_scaled_delta_T_sed:.6e} "
        f"(=unguarded {delta_T_sed_unguarded:.6e} * ramp {_scalar(ramp):.4f} "
        f"* clip_ratio {_scalar(clip_ratio):.4f}), rel={rel_error:.3e}. "
        f"Pre-F2-fix the actual would match the unguarded value, ignoring "
        f"the cap; the scaled value above must hold for energy conservation."
    )
    # Sanity: the cap fired (scaling factor < 1 because clip_ratio < 1).
    assert _scalar(clip_ratio) < 1.0, (
        "Test setup did not actually trigger the cap; tighten the cap or "
        "increase the T_sed - T_water gradient."
    )
