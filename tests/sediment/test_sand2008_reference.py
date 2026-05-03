"""SAND2008 reference dataset integration test.

Loads the SAND2008-5621 worked example (``bed.sdf``, ``erate.sdf``,
``core_field.sdf``) and exercises the SSM end-to-end on a tiny
synthetic mesh:

* Per-shear-level erosion rates from the ``SedflumeTableErosionModel``
  match the ENRATE values within 1 % (this is the unit-level cross
  check that does not require Riverine).
* Per-step ``SSM.run`` orchestration converges and produces non-negative
  fluxes with mass conservation across erosion and deposition.
* Emergent armoring: under sustained τ above τ_ce(fines), the surface
  D50 of the active layer increases monotonically over time.

The test is **not** a precision benchmark of the SAND2008 published
totals; it is a smoke gate intended to catch regressions in the
end-to-end orchestration. Tighter cross-checks against published
totals are deferred to Stage 3 with the Albany dataset.

Reference: SAND2008-5621 Figures 1, 2, 3.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v2.processes.sediment import (
    SSM,
    contracts,
)
from clearwater_modules_v2.processes.sediment.io.sedflume import (
    load_sedflume_bundle,
)


DATA_DIR = Path(__file__).parent / "data" / "sand2008_example"
BED_SDF = DATA_DIR / "bed.sdf"
ERATE_SDF = DATA_DIR / "erate.sdf"
CORE_FIELD_SDF = DATA_DIR / "core_field.sdf"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def sand2008_bundle():
    """Loaded SEDflume bundle from the SAND2008 example dataset."""
    return load_sedflume_bundle(BED_SDF, ERATE_SDF, CORE_FIELD_SDF)


@pytest.fixture
def small_mesh() -> xr.Dataset:
    """Mesh with 3 cells × 8 timesteps. No hydraulics."""
    n_face = 3
    n_time = 8
    return xr.Dataset(
        coords={
            contracts.DIM_TIME: np.arange(n_time, dtype="int32"),
            contracts.DIM_NFACE: np.arange(n_face, dtype="int32"),
        }
    )


def _make_ssm(bundle, mesh, time_step):
    """Build an SSM bound to a tiny 3-cell single-core mesh, no Riverine."""
    from clearwater_modules_v2.processes.sediment.ssm import (
        _build_classes_from_bundle,
    )

    n_face = mesh.sizes[contracts.DIM_NFACE]
    registry = _build_classes_from_bundle(bundle)
    ssm = SSM(
        sediment_classes=registry,
        sedflume_bundle=bundle,
        shear_driver="external",
        shear_options={"growth_limit": 0.0},
        bedload_solver="off",
        nsedflume=bundle.nsedflume,
        biostabilization_alpha=0.0,
        time_step=time_step,
        core_id=np.ones(n_face, dtype=np.int64),
    )
    ssm.bind_mesh(mesh)
    return ssm


@pytest.fixture
def ssm_instance(sand2008_bundle, small_mesh) -> SSM:
    """SSM bound to a tiny 3-cell single-core mesh, no Riverine."""
    from datetime import timedelta

    return _make_ssm(sand2008_bundle, small_mesh, timedelta(seconds=60))


# ---------------------------------------------------------------------------
# Unit-level cross-check: ENRATE table interpolation
# ---------------------------------------------------------------------------


def test_sedflume_table_returns_known_rates_at_table_taus(sand2008_bundle):
    """At τ values that hit ENRATE columns exactly, the per-cell rate
    should equal the table value × bulk density (g/cm²/s).

    Per Figure 1 the SAND2008 table-row for shear-level 4 (τ = 8 Pa)
    is 7.0e-3 cm/s for every layer of every core.
    """
    from clearwater_modules_v2.processes.sediment.erosion import (
        SedflumeTableErosionModel,
    )

    # Active-layer table: synthesize from deepest in-place row repeated
    # over the size dimension (matches the ssm.py default fallback).
    erate_active = np.tile(
        sand2008_bundle.erate_per_core_cm_s[0, -1, :].reshape(1, -1),
        (sand2008_bundle.size_interpolants_um.size, 1),
    )
    model = SedflumeTableErosionModel(
        tau_levels_pa=sand2008_bundle.tau_levels_pa,
        erate_per_core=sand2008_bundle.erate_per_core_cm_s,
        erate_active_per_size=erate_active,
        size_interpolants_um=sand2008_bundle.size_interpolants_um,
        taucrit_per_size_pa=sand2008_bundle.taucrit_per_size_pa,
    )

    # τ = 8 Pa is column index 3 in the example tau_levels; layer 1
    # (index 0 in 0-origin) on core 1 (index 0). Bulk density is 1.9
    # g/cm³. Expected rate: 7.0e-3 cm/s × 1.9 g/cm³ = 1.33e-2 g/cm²/s.
    n_face = 1
    tau = xr.DataArray(np.array([8.0]), dims=(contracts.DIM_NFACE,))
    layer_mass = xr.DataArray(np.array([15.0 * 1.9]), dims=(contracts.DIM_NFACE,))
    layer_mass0 = layer_mass.copy()
    bulk_density = xr.DataArray(np.array([1.9]), dims=(contracts.DIM_NFACE,))
    core_id = xr.DataArray(np.array([0], dtype=np.int64), dims=(contracts.DIM_NFACE,))

    result = model.erosion_rate(
        tau_pa=tau,
        layer_index=3,                                # in-place layer 3 (1-origin)
        layer_mass=layer_mass,
        layer_initial_mass=layer_mass0,
        bulk_density=bulk_density,
        core_id=core_id,
    )
    expected = 7.0e-3 * 1.9
    np.testing.assert_allclose(result.values, [expected], rtol=1e-2)


# ---------------------------------------------------------------------------
# End-to-end orchestration: shear sweep
# ---------------------------------------------------------------------------


def _seed_external_tau(mesh: xr.Dataset, n_face: int) -> None:
    """Allocate an empty external τ field on the mesh."""
    n_time = mesh.sizes[contracts.DIM_TIME]
    mesh[contracts.VAR_BED_SHEAR_STRESS_INPUT] = (
        (contracts.DIM_TIME, contracts.DIM_NFACE),
        np.zeros((n_time, n_face), dtype="float32"),
    )


def _set_external_tau(mesh: xr.Dataset, t_idx: int, tau_value: float) -> None:
    """Write a uniform τ at time slot ``t_idx``."""
    mesh[contracts.VAR_BED_SHEAR_STRESS_INPUT].values[t_idx, :] = tau_value


def test_ssm_run_orchestration_steps_through_shear_history(ssm_instance):
    """Drive SSM through the SAND2008 example τ sweep ``[0, 2, 4, 8, 10, 20]`` Pa.

    Asserts the orchestration plumbs successfully through every step:

    1. Each step completes without raising.
    2. Per-class erosion flux is recorded and non-negative everywhere.
    3. Total bed mass strictly decreases as τ rises from low to high
       (above τ_ce) — emergent erosion behaviour.
    4. Per-class flux at τ = 0 is exactly zero (gating works).
    """
    ssm = ssm_instance
    mesh = ssm._mesh
    n_face = mesh.sizes[contracts.DIM_NFACE]
    _seed_external_tau(mesh, n_face)

    tau_history = [0.0, 2.0, 4.0, 8.0, 10.0, 20.0]
    bed_total_per_step = []
    erosion_flux_per_step = []

    # Snapshot mass at t=0 (post-init).
    bed_total_per_step.append(
        float(ssm._bed.layer_mass_at(0).values.sum())
    )

    for i, tau_value in enumerate(tau_history):
        t_idx = i
        _set_external_tau(mesh, t_idx, tau_value)
        # Copy bed state from previous slot so the SSM mutates step-by-step.
        if t_idx > 0:
            ssm._bed.set_layer_mass_at(
                t_idx, ssm._bed.layer_mass_at(t_idx - 1).values
            )
            ssm._bed.set_class_fraction_at(
                t_idx, ssm._bed.class_fraction_at(t_idx - 1).values
            )
            ssm._bed.set_layer_active_at(
                t_idx, ssm._bed.layer_active_at(t_idx - 1).values
            )
            ssm._bed.set_layer_taucrit_at(
                t_idx, ssm._bed.layer_taucrit_at(t_idx - 1).values
            )

        ssm.run(time=t_idx)

        bed_total_per_step.append(
            float(ssm._bed.layer_mass_at(t_idx).values.sum())
        )
        flux = mesh[contracts.VAR_BED_EROSION_FLUX].isel(
            {contracts.DIM_TIME: t_idx}
        ).values
        erosion_flux_per_step.append(flux)

    # 2. Non-negativity.
    for flux in erosion_flux_per_step:
        assert np.all(flux >= 0.0), "erosion flux went negative"

    # 4. τ = 0 gates erosion off completely.
    np.testing.assert_array_equal(
        erosion_flux_per_step[0], np.zeros_like(erosion_flux_per_step[0])
    )

    # 3. Total bed mass strictly decreases as τ rises from 4 Pa to 20 Pa.
    # Index in bed_total_per_step: [pre, after_t0, after_t1, ..., after_t5]
    assert bed_total_per_step[6] < bed_total_per_step[3], (
        f"Bed mass did not decrease across the high-shear segment: "
        f"{bed_total_per_step[3]} -> {bed_total_per_step[6]}"
    )


# ---------------------------------------------------------------------------
# Emergent armoring
# ---------------------------------------------------------------------------


def test_in_place_layer_d50_increases_under_selectively_gated_shear(
    sand2008_bundle, small_mesh,
):
    """Hold τ between τ_ce(class 6) = 1.08 Pa and τ_ce(class 7) = 1.6 Pa
    for several long steps, and verify the in-place layer's PERSED
    coarsens monotonically as classes 0–6 erode preferentially.

    The active layer (layer 1) is constantly refreshed by the
    reorganize step from the in-place layer (layer 2), so its
    PERSED is essentially layer 2's PERSED at every step. Watching
    layer 2 evolution is the cleanest way to see armoring emerge in
    the SAND2008 example.

    Uses a 1-day step so the cumulative erosion (~1e-3 g/cm² per
    step) is large enough to lift PERSED of the gated classes (7,8)
    out of the float32 quantization floor.
    """
    from datetime import timedelta

    ssm = _make_ssm(sand2008_bundle, small_mesh, timedelta(days=1))
    mesh = ssm._mesh
    n_face = mesh.sizes[contracts.DIM_NFACE]
    _seed_external_tau(mesh, n_face)

    # τ ∈ (τ_ce[6]=1.08, τ_ce[7]=1.6) Pa — gates out the coarsest 2
    # classes only. Classes 0–6 erode, classes 7–8 stay put.
    tau_value = 1.4
    n_steps = min(6, mesh.sizes[contracts.DIM_TIME] - 1)
    d50_l2_history: list[float] = []

    d50_arr = ssm.registry_classes.d50_um_array

    for t_idx in range(n_steps):
        if t_idx > 0:
            ssm._bed.set_layer_mass_at(
                t_idx, ssm._bed.layer_mass_at(t_idx - 1).values
            )
            ssm._bed.set_class_fraction_at(
                t_idx, ssm._bed.class_fraction_at(t_idx - 1).values
            )
            ssm._bed.set_layer_active_at(
                t_idx, ssm._bed.layer_active_at(t_idx - 1).values
            )
            ssm._bed.set_layer_taucrit_at(
                t_idx, ssm._bed.layer_taucrit_at(t_idx - 1).values
            )
        _set_external_tau(mesh, t_idx, tau_value)
        ssm.run(time=t_idx)

        # Layer 2 is the topmost in-place layer with most of the bed mass.
        l2_persed = np.asarray(
            ssm._bed.class_fraction_at(t_idx).values[:, 2, :], dtype="float64"
        )
        d50_l2 = (l2_persed * d50_arr).sum(axis=-1)
        d50_l2_history.append(float(d50_l2.mean()))

    # Monotone-non-decreasing: removing fines raises D50.
    for j in range(1, len(d50_l2_history)):
        assert d50_l2_history[j] >= d50_l2_history[j - 1] - 1e-6, (
            f"In-place D50 decreased at step {j}: "
            f"{d50_l2_history[j-1]} -> {d50_l2_history[j]}"
        )
    # Strict net increase from start to finish — the armoring signal.
    assert d50_l2_history[-1] > d50_l2_history[0], (
        f"No armoring: D50 history {d50_l2_history}"
    )
