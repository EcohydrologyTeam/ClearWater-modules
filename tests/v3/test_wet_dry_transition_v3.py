"""Phase R-5 wet/dry transition regression test (review-findings 2026-05-04).

Pins the documented behavior change between Phase 2 (per-process
``xr.where(volume > 0, delta, 0)`` guard inside ``Temperature.run``,
``T_water`` carries forward through dry intervals) and Phase 3 + the
m19 cleanup (registry-level wet-mask in ``Model.__apply_wet_mask`` is
the single point of dry-cell handling; ``T_water`` becomes NaN on dry
cells when a wet-mask is configured, and the per-process volume guard
has been removed).

Post-removal semantics that this file pins:

- **Without** ``wet_mask_variable`` configured: dry cells with
  ``volume == 0`` produce NaN/inf in ``temperature_change`` (division
  by zero); the kernel writes those NaN/inf back to the registry.
  Users are expected to either run on uniformly-wet meshes or
  configure ``wet_mask_variable`` for proper handling.
- **With** ``wet_mask_variable`` configured:
  ``Model.__apply_wet_mask`` writes NaN into ``output_variables`` on
  cells where the masking variable falls below threshold, regardless
  of what the kernel produced. Dry cells get clean NaN.

Stub patterns mirror those in
``tests/v3/test_model_orchestration_v3.py`` and
``tests/v3/test_wet_mask_scope_v3.py``; copies are inlined here so this
file does not depend on test-private state in those files.
"""

from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np
import xarray as xr

from clearwater_modules_v3 import Temperature
from clearwater_modules_v3.model import Model


# ---------------------------------------------------------------------------
# Stubs (inlined from test_model_orchestration_v3.py / test_wet_mask_scope_v3.py)
# ---------------------------------------------------------------------------


class _StubRegistry:
    """Minimal VariableRegistry stand-in: name -> DataArray dict."""

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


def _build_model(
    processes,
    *,
    wet_mask_variable: str | None = None,
    wet_mask_threshold: float = 0.0,
) -> Model:
    return Model(
        processes=tuple(processes),
        variable_registry=_StubRegistry(),
        variable_data_sources={},
        start_time=datetime(2026, 1, 1, 0, 0, 0),
        end_time=datetime(2026, 1, 1, 0, 5, 0),
        time_step=timedelta(minutes=5),
        output_variables=[],
        wet_mask_variable=wet_mask_variable,
        wet_mask_threshold=wet_mask_threshold,
    )


def _make_temperature() -> Temperature:
    """Construct a Temperature whose ``__skip_first_time_step`` is False so
    the first ``run`` call actually exercises the kinetics path. Matches
    the pattern used by the v3 hotstart tests."""
    t = Temperature(wind_a=0.3, wind_b=1.5, wind_c=3.0)
    t.from_hotstart({})  # flips __skip_first_time_step to False
    return t


def _seed_three_cell_registry(
    registry: _StubRegistry,
    *,
    water_temperature: np.ndarray,
    volume: np.ndarray,
) -> None:
    """Pre-populate a 3-cell registry with finite forcings; ``volume`` and
    ``water_temperature`` are caller-controlled so wet/dry mixes can be
    tested. Surface area follows volume (zero where dry) so the depth
    ramp behaves consistently with the wet/dry pattern."""
    surface_area = np.where(volume > 0.0, 100.0, 0.0)
    registry.register("water_temperature", xr.DataArray(water_temperature))
    registry.register("wetted_surface_area", xr.DataArray(surface_area))
    registry.register("volume", xr.DataArray(volume.astype(float)))
    registry.register("cloudiness", xr.DataArray(np.array([0.0, 0.0, 0.0])))
    registry.register("air_temperature", xr.DataArray(np.array([25.0, 25.0, 25.0])))
    registry.register("solar_radiation", xr.DataArray(np.array([800.0, 800.0, 800.0])))
    registry.register("wind_speed", xr.DataArray(np.array([3.0, 3.0, 3.0])))
    registry.register(
        "atmospheric_pressure", xr.DataArray(np.array([1013.0, 1013.0, 1013.0]))
    )
    registry.register(
        "atmospheric_vapor_pressure", xr.DataArray(np.array([20.0, 20.0, 20.0]))
    )
    registry.register(
        "sediment_temperature", xr.DataArray(np.array([15.0, 16.0, 17.0]))
    )
    registry.register(
        "sediment_thickness", xr.DataArray(np.array([0.1, 0.1, 0.1]))
    )


# ---------------------------------------------------------------------------
# Scenario 1: dry cell without wet-mask -> NaN/inf surfaces (defect visible)
# ---------------------------------------------------------------------------


def test_dry_cell_without_wet_mask_produces_nan_or_inf():
    """3-cell mesh: cell 1 dry (volume=0). One Temperature substep with
    NO ``wet_mask_variable`` configured. Post-m19 behavior: the per-process
    volume guard is removed, so ``temperature_change`` divides by zero on
    cell 1 and the resulting NaN/inf is written back to the registry. The
    wet cells must remain finite and unchanged in sign-of-magnitude."""
    t = _make_temperature()
    model = _build_model(processes=[t])  # wet_mask_variable=None
    registry = model._Model__registry

    initial_T = np.array([20.0, 18.0, 22.0])
    volume = np.array([1.0, 0.0, 1.0])
    _seed_three_cell_registry(
        registry, water_temperature=initial_T.copy(), volume=volume
    )

    t.run(datetime(2026, 1, 1, 0, 0, 0), registry)

    result = registry.get("water_temperature").values

    # Cell 1: division by volume=0 in temperature_change. The exact value
    # is NaN or inf depending on numerator sign; we assert "not finite"
    # rather than pinning a specific NaN-vs-inf outcome (the kernel does
    # not normalize between them, and the contract is "the defect is
    # visible at the output").
    assert not np.isfinite(result[1]), (
        f"Without wet-mask, dry-cell water_temperature should be NaN or "
        f"inf (the documented Phase-3-after-m19 defect-visibility "
        f"contract); got {result[1]!r}"
    )

    # Cells 0 and 2: wet, must remain finite. We do not pin the exact
    # post-substep temperature (kinetics depend on many constants); only
    # that the value did not become NaN/inf.
    assert np.isfinite(result[0])
    assert np.isfinite(result[2])


# ---------------------------------------------------------------------------
# Scenario 2: dry cell with wet-mask -> clean NaN at orchestration layer
# ---------------------------------------------------------------------------


def test_dry_cell_with_wet_mask_produces_clean_nan():
    """Same 3-cell mesh but with ``wet_mask_variable="volume"``. The
    Temperature kernel still produces NaN/inf on cell 1 (the per-process
    guard is gone), but ``Model.__apply_wet_mask`` overwrites
    ``output_variables`` on dry cells with NaN. The end-of-substep
    registry value on cell 1 must therefore be NaN, regardless of
    whether the kernel produced NaN or inf."""
    t = _make_temperature()
    model = _build_model(
        processes=[t], wet_mask_variable="volume", wet_mask_threshold=0.0
    )
    registry = model._Model__registry

    initial_T = np.array([20.0, 18.0, 22.0])
    volume = np.array([1.0, 0.0, 1.0])
    _seed_three_cell_registry(
        registry, water_temperature=initial_T.copy(), volume=volume
    )

    # Drive the orchestration path the way ``__process_loop_full`` does:
    # run the process, then apply the wet-mask.
    current_time = datetime(2026, 1, 1, 0, 0, 0)
    wet_mask = model._Model__compute_wet_mask(current_time)
    t.run(current_time, registry)
    model._Model__apply_wet_mask(t, current_time, wet_mask)

    result = registry.get("water_temperature").values

    # Wet cells: finite (kinetics ran).
    assert np.isfinite(result[0])
    assert np.isfinite(result[2])
    # Dry cell: clean NaN, not inf. The wet-mask overwrites whatever the
    # kernel produced (NaN or inf) with NaN.
    assert np.isnan(result[1]), (
        f"With wet-mask configured, dry-cell water_temperature must be "
        f"NaN (overwritten by Model.__apply_wet_mask); got {result[1]!r}"
    )


# ---------------------------------------------------------------------------
# Scenario 3: wet/dry/wet trajectory across three substeps
# ---------------------------------------------------------------------------


def test_wet_dry_wet_trajectory_with_wet_mask():
    """Single-cell trajectory volume = [1.0, 0.0, 1.0] across three
    substeps. With wet-mask configured:
      - substep 0 (volume=1.0): cell is wet -> finite output.
      - substep 1 (volume=0.0): cell is dry -> NaN output.
      - substep 2 (volume=1.0): cell is wet again. The kernel reads
        the previous substep's NaN ``water_temperature`` from the
        registry (NaN is not auto-restored by the wet-mask going
        wet->dry->wet); arithmetic on NaN produces NaN, so the result
        stays NaN. This is the documented behavior change versus
        Phase 2, where the per-process guard kept ``T_water``
        unchanged through dry intervals so the wet-again substep would
        resume from the pre-dry value.

    The test pins the post-m19 trajectory: NaN on the dry substep, and
    NaN persists into the wet-again substep until the caller restores
    a finite initial condition (which is the caller's responsibility
    under the documented contract)."""
    t = _make_temperature()
    model = _build_model(
        processes=[t], wet_mask_variable="volume", wet_mask_threshold=0.0
    )
    registry = model._Model__registry

    # Single-cell mesh.
    def _seed_one_cell(T_init: float, vol: float) -> None:
        sa = 100.0 if vol > 0.0 else 0.0
        registry.register("water_temperature", xr.DataArray(np.array([T_init])))
        registry.register("wetted_surface_area", xr.DataArray(np.array([sa])))
        registry.register("volume", xr.DataArray(np.array([float(vol)])))
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

    def _step_with_volume(time: datetime, vol: float) -> float:
        # Update only volume / surface_area; keep T_water as registry holds.
        sa = 100.0 if vol > 0.0 else 0.0
        registry.register(
            "wetted_surface_area", xr.DataArray(np.array([sa]))
        )
        registry.register("volume", xr.DataArray(np.array([float(vol)])))
        wet_mask = model._Model__compute_wet_mask(time)
        t.run(time, registry)
        model._Model__apply_wet_mask(t, time, wet_mask)
        return float(registry.get("water_temperature").values[0])

    _seed_one_cell(T_init=20.0, vol=1.0)

    # Substep 0: wet -> finite.
    val_0 = _step_with_volume(datetime(2026, 1, 1, 0, 0, 0), vol=1.0)
    assert np.isfinite(val_0)

    # Substep 1: dry -> NaN (wet-mask overwrites whatever kernel produced).
    val_1 = _step_with_volume(datetime(2026, 1, 1, 0, 5, 0), vol=0.0)
    assert np.isnan(val_1)

    # Substep 2: wet again, but the registry's ``water_temperature`` is
    # NaN. NaN + finite_delta == NaN. After the substep, the cell is
    # wet (volume=1.0), so the wet-mask leaves the kernel's output in
    # place. Result: NaN persists. This is the documented Phase-3
    # behavior; the caller is responsible for re-seeding the cell when
    # it transitions from dry back to wet.
    val_2 = _step_with_volume(datetime(2026, 1, 1, 0, 10, 0), vol=1.0)
    assert np.isnan(val_2), (
        f"Wet-again substep starting from a NaN T_water should remain "
        f"NaN under the documented post-m19 contract; got {val_2!r}. "
        f"Phase 2's per-process guard would have kept the pre-dry "
        f"value through the dry interval, but Phase 3 deliberately "
        f"removed that behavior."
    )


# ---------------------------------------------------------------------------
# Scenario 4: dry-throughout cell with wet-mask -> NaN every substep
# ---------------------------------------------------------------------------


def test_dry_throughout_cell_stays_nan_every_substep():
    """volume = [0, 0, 0] across three substeps with wet-mask configured.
    The cell is dry every substep; the wet-mask overwrites
    ``water_temperature`` with NaN on every substep, so the registry
    holds NaN throughout."""
    t = _make_temperature()
    model = _build_model(
        processes=[t], wet_mask_variable="volume", wet_mask_threshold=0.0
    )
    registry = model._Model__registry

    # Single-cell mesh, dry throughout.
    registry.register(
        "water_temperature", xr.DataArray(np.array([20.0]))
    )
    registry.register(
        "wetted_surface_area", xr.DataArray(np.array([0.0]))
    )
    registry.register("volume", xr.DataArray(np.array([0.0])))
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

    times = [
        datetime(2026, 1, 1, 0, 0, 0),
        datetime(2026, 1, 1, 0, 5, 0),
        datetime(2026, 1, 1, 0, 10, 0),
    ]
    for time in times:
        wet_mask = model._Model__compute_wet_mask(time)
        t.run(time, registry)
        model._Model__apply_wet_mask(t, time, wet_mask)
        result = registry.get("water_temperature").values
        assert np.isnan(result[0]), (
            f"Dry-throughout cell at {time}: water_temperature must be "
            f"NaN under the documented wet-mask contract; got {result[0]!r}"
        )
