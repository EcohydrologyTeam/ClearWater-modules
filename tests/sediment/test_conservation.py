"""Closed-domain mass-conservation integration test.

Builds a tiny synthetic 3-cell mesh with no flow boundaries and a
hand-crafted SedflumeBundle (1 cohesive class + 1 sand class, 5
layers), then drives SSM through three different shear regimes:

* τ = 0 Pa  — pure deposition (when seeded with non-zero suspended C).
* τ = high — strong erosion.
* τ = mid  — mixed erosion and deposition.

For each step, the total bed mass + suspended mass + bedload mass
must be invariant up to the per-step erosion / deposition fluxes
the SSM emits as Riverine source/sink terms. In a CLOSED domain
(no advection out, no boundary inflow), the algebraic identity is

    Δ(bed) + Δ(suspended_in_water_col)  =  -Δ(net injected source)

Since SSM is the only producer of those source terms here and we
account for them explicitly, the global mass invariant holds to
``_MASS_CONSERVATION_TOL`` (g/cm²).

Reference: bed.py ``_MASS_CONSERVATION_TOL``; design spec §5.8.
"""

from __future__ import annotations

from datetime import timedelta

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v2.processes.sediment import SSM, contracts
from clearwater_modules_v2.processes.sediment.classes import (
    SedimentClass,
    SedimentClassRegistry,
)
from clearwater_modules_v2.processes.sediment.io.sedflume import SedflumeBundle


_TOL_GCM2 = 1e-5  # matches bed.py's _MASS_CONSERVATION_TOL


# ---------------------------------------------------------------------------
# Hand-crafted bundle
# ---------------------------------------------------------------------------


def _make_synthetic_bundle() -> SedflumeBundle:
    """Build a tiny 1-core, 5-layer, 2-class bundle for unit testing.

    Mirrors the ``test_io_sedflume.test_load_yaml_config_round_trip``
    layout but is built directly from numpy arrays for full control.
    """
    n_layers = 5
    n_class = 2
    # layer indices [0,1] are active/deposition (start empty), [2,3,4] in-place
    layer_thickness = np.array([[0.0, 0.0, 5.0, 5.0, 5.0]])
    layer_taucrit_pa = np.array([[0.4, 0.4, 0.4, 0.5, 0.6]])
    bulk_density = np.array([[1.6, 1.6, 1.6, 1.7, 1.8]])
    # PSD per layer (0..4): 60 % silt / 40 % sand throughout.
    psd = np.array([[
        [60.0, 40.0],
        [60.0, 40.0],
        [60.0, 40.0],
        [60.0, 40.0],
        [60.0, 40.0],
    ]])

    # Power-law erosion (nsedflume = 2): same A, n, max for every layer.
    ea = np.full((1, n_layers), 5.0e-4)   # cm/s for τ in Pa
    en = np.full((1, n_layers), 2.0)
    max_rate = np.full((1, n_layers), 1.0)

    return SedflumeBundle(
        n_layers=n_layers,
        var_bed=1,
        icalc_bl=0,
        nsedflume=2,
        zb_skin_um=1500.0,
        tau_const_pa=0.0,
        bedload_cutoff_um=64.0,
        max_deposit_limit=1.0,
        d50_um=np.array([32.0, 250.0]),
        tau_ce_pa=np.array([0.15, 0.20]),
        tau_cs_pa=np.array([0.20, 0.30]),
        settling_cm_s=np.array([-1.0, -1.0]),
        size_interpolants_um=np.array([32.0, 250.0]),
        taucrit_per_size_pa=np.array([0.15, 0.20]),
        erate_active_table=None,
        actdep_a=np.array([5.0e-4, 5.0e-4]),
        actdep_n=np.array([2.0, 2.0]),
        actdep_max=np.array([1.0, 1.0]),
        n_cores=1,
        layer_thickness_cm=layer_thickness,
        layer_taucrit_pa=layer_taucrit_pa,
        bulk_density_g_cm3=bulk_density,
        water_density_g_cm3=1.0,
        solid_density_g_cm3=2.65,
        particle_size_distribution_pct=psd,
        tau_levels_pa=np.array([0.0, 1.0]),
        erate_per_core_cm_s=None,
        ea_per_core=ea,
        en_per_core=en,
        max_rate_per_core_cm_s=max_rate,
        core_field_ij=None,
    )


@pytest.fixture
def bundle() -> SedflumeBundle:
    return _make_synthetic_bundle()


@pytest.fixture
def mesh() -> xr.Dataset:
    """3 cells × 6 timesteps; no advection / boundary inflow."""
    n_face = 3
    n_time = 6
    return xr.Dataset(
        coords={
            contracts.DIM_TIME: np.arange(n_time, dtype="int32"),
            contracts.DIM_NFACE: np.arange(n_face, dtype="int32"),
        }
    )


@pytest.fixture
def ssm(bundle, mesh) -> SSM:
    """SSM bound to the synthetic mesh, no Riverine, no bedload."""
    n_face = mesh.sizes[contracts.DIM_NFACE]
    registry = SedimentClassRegistry.from_iterable(
        [
            SedimentClass(
                label="silt_fine",
                d50_um=32.0,
                tau_ce_pa=0.15,
                tau_cs_pa=0.20,
            ),
            SedimentClass(
                label="sand_med",
                d50_um=250.0,
                tau_ce_pa=0.20,
                tau_cs_pa=0.30,
            ),
        ]
    )
    instance = SSM(
        sediment_classes=registry,
        sedflume_bundle=bundle,
        shear_driver="external",
        shear_options={"growth_limit": 0.0},
        bedload_solver="off",
        nsedflume=2,
        biostabilization_alpha=0.0,
        time_step=timedelta(seconds=60),
        core_id=np.ones(n_face, dtype=np.int64),
    )
    instance.bind_mesh(mesh)
    return instance


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _seed_external_tau(mesh: xr.Dataset, n_face: int) -> None:
    n_time = mesh.sizes[contracts.DIM_TIME]
    if contracts.VAR_BED_SHEAR_STRESS_INPUT not in mesh.data_vars:
        mesh[contracts.VAR_BED_SHEAR_STRESS_INPUT] = (
            (contracts.DIM_TIME, contracts.DIM_NFACE),
            np.zeros((n_time, n_face), dtype="float32"),
        )


def _set_external_tau(mesh: xr.Dataset, t_idx: int, tau_value: float) -> None:
    mesh[contracts.VAR_BED_SHEAR_STRESS_INPUT].values[t_idx, :] = tau_value


def _bed_total_mass(ssm: SSM, t_idx: int) -> float:
    """Sum of per-cell, per-layer bed mass (g/cm²)."""
    return float(ssm._bed.layer_mass_at(t_idx).values.sum())


def _bedload_total_mass(mesh: xr.Dataset, t_idx: int) -> float:
    """Sum of per-cell per-class bedload mass; 0 when bedload is off."""
    if contracts.VAR_BEDLOAD_MASS not in mesh.data_vars:
        return 0.0
    da = mesh[contracts.VAR_BEDLOAD_MASS]
    if contracts.DIM_TIME in da.dims:
        return float(da.isel({contracts.DIM_TIME: t_idx}).values.sum())
    return float(da.values.sum())


def _net_source_total(
    mesh: xr.Dataset, registry: SedimentClassRegistry
) -> float:
    """Sum of net (E - D) source/sink staged on the mesh by the last run.

    SSM writes a per-class ``{cls.suspended_var}_source`` field on
    ``(nface,)`` per step. Returns the scalar total over cells × classes.
    """
    total = 0.0
    for cls in registry:
        name = f"{cls.suspended_var}_source"
        if name in mesh.data_vars:
            total += float(mesh[name].values.sum())
    return total


def _step_with_state_carried_forward(ssm: SSM, t_idx: int) -> None:
    """Advance one SSM step at slot ``t_idx``, copying state from t_idx-1."""
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


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_zero_shear_zero_concentration_is_a_no_op(ssm):
    """With τ = 0 and C = 0 (no suspended sediment), nothing happens.

    Verifies the no-flux baseline: total bed mass is exactly invariant.
    """
    mesh = ssm._mesh
    n_face = mesh.sizes[contracts.DIM_NFACE]
    _seed_external_tau(mesh, n_face)

    bed_mass_initial = _bed_total_mass(ssm, 0)
    for t_idx in range(3):
        _set_external_tau(mesh, t_idx, 0.0)
        _step_with_state_carried_forward(ssm, t_idx)
    bed_mass_final = _bed_total_mass(ssm, 2)

    np.testing.assert_allclose(
        bed_mass_final, bed_mass_initial, atol=_TOL_GCM2,
        err_msg="τ=0, C=0 should leave bed mass exactly invariant",
    )

    # Net source/sink should be exactly 0.
    net = _net_source_total(mesh, ssm.registry_classes)
    np.testing.assert_allclose(net, 0.0, atol=_TOL_GCM2)


def test_erosion_step_balances_with_injected_source(ssm):
    """Under high τ, erosion mass exactly equals the injected source.

    Initial suspended C = 0, so deposition is 0; net source = erosion.
    Mass balance: Δ(bed) = -Σ(source) per step.
    """
    mesh = ssm._mesh
    n_face = mesh.sizes[contracts.DIM_NFACE]
    _seed_external_tau(mesh, n_face)

    bed_mass_history: list[float] = [_bed_total_mass(ssm, 0)]
    src_history: list[float] = []

    tau_value = 1.0  # well above τ_ce(silt)=0.15 and τ_ce(sand)=0.20
    for t_idx in range(4):
        _set_external_tau(mesh, t_idx, tau_value)
        _step_with_state_carried_forward(ssm, t_idx)
        bed_mass_history.append(_bed_total_mass(ssm, t_idx))
        src_history.append(_net_source_total(mesh, ssm.registry_classes))

    # Bed mass strictly decreases under sustained erosion.
    for i in range(1, len(bed_mass_history)):
        assert bed_mass_history[i] <= bed_mass_history[i - 1] + _TOL_GCM2

    # Per-step balance: Δ(bed) ≈ -src (deposition is 0 with C=0).
    for i in range(1, len(bed_mass_history)):
        delta_bed = bed_mass_history[i] - bed_mass_history[i - 1]
        # The src field is overwritten each step (1-D, no time dim), so
        # src_history[i-1] is the most recent step's net source. Sign:
        # erosion → bed loses mass → src is positive.
        np.testing.assert_allclose(
            delta_bed, -src_history[i - 1], atol=_TOL_GCM2, rtol=1e-3,
            err_msg=f"step {i}: Δbed={delta_bed} ≠ -src={-src_history[i-1]}",
        )


def test_deposition_step_balances_with_injected_sink(ssm):
    """With τ < τ_cs and a non-zero suspended C, deposition removes mass
    from the water column and adds it to the bed.

    Mass balance: Δ(bed) = -Σ(source); the source array carries
    *negative* values (sink) for deposition, so bed gains and water loses.
    """
    mesh = ssm._mesh
    n_face = mesh.sizes[contracts.DIM_NFACE]
    _seed_external_tau(mesh, n_face)

    # Seed a nonzero suspended concentration on the mesh for both classes.
    # Override via run kwargs; the mesh itself does not need pre-allocated
    # constituent variables for this synthetic test.
    susp_c = xr.DataArray(
        np.full((n_face, len(ssm.registry_classes)), 100.0, dtype="float64"),
        dims=(contracts.DIM_NFACE, contracts.DIM_CLASS),
    )
    bottom_depth = xr.DataArray(
        np.full(n_face, 1.0, dtype="float64"),
        dims=(contracts.DIM_NFACE,),
    )

    bed_mass_history: list[float] = [_bed_total_mass(ssm, 0)]
    src_history: list[float] = []

    # τ = 0 ⇒ Krone P_d = 1 for cohesive, Gessler P_d → 1 for non-cohesive.
    for t_idx in range(3):
        _set_external_tau(mesh, t_idx, 0.0)
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
        ssm.run(
            time=t_idx,
            suspended_concentration=susp_c,
            bottom_water_layer_depth_m=bottom_depth,
        )
        bed_mass_history.append(_bed_total_mass(ssm, t_idx))
        src_history.append(_net_source_total(mesh, ssm.registry_classes))

    # Bed mass strictly increases under sustained deposition (no erosion).
    for i in range(1, len(bed_mass_history)):
        assert bed_mass_history[i] >= bed_mass_history[i - 1] - _TOL_GCM2

    # Per-step balance: Δ(bed) ≈ -src (erosion is 0 with τ=0).
    # Deposition → src is negative; bed gains mass.
    for i in range(1, len(bed_mass_history)):
        delta_bed = bed_mass_history[i] - bed_mass_history[i - 1]
        np.testing.assert_allclose(
            delta_bed, -src_history[i - 1], atol=_TOL_GCM2, rtol=1e-3,
            err_msg=f"step {i}: Δbed={delta_bed} ≠ -src={-src_history[i-1]}",
        )


def test_mixed_shear_history_preserves_global_mass_balance(ssm):
    """Drive the SSM through alternating erosion and deposition steps,
    and verify the cumulative balance: total bed mass change exactly
    cancels with the cumulative net source/sink.

    Sequence of (τ, C) per step:

        (0.0, 100)  — deposition
        (0.5, 50)   — mixed: erosion of fines, deposition of sand fraction
        (1.0, 0)    — pure erosion
        (0.0, 100)  — deposition again
    """
    mesh = ssm._mesh
    n_face = mesh.sizes[contracts.DIM_NFACE]
    _seed_external_tau(mesh, n_face)

    bottom_depth = xr.DataArray(
        np.full(n_face, 1.0, dtype="float64"),
        dims=(contracts.DIM_NFACE,),
    )

    bed_mass_initial = _bed_total_mass(ssm, 0)
    cumulative_src = 0.0

    schedule = [(0.0, 100.0), (0.5, 50.0), (1.0, 0.0), (0.0, 100.0)]
    for t_idx, (tau_value, c_value) in enumerate(schedule):
        _set_external_tau(mesh, t_idx, tau_value)
        susp_c = xr.DataArray(
            np.full(
                (n_face, len(ssm.registry_classes)), c_value, dtype="float64"
            ),
            dims=(contracts.DIM_NFACE, contracts.DIM_CLASS),
        )
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
        ssm.run(
            time=t_idx,
            suspended_concentration=susp_c,
            bottom_water_layer_depth_m=bottom_depth,
        )
        cumulative_src += _net_source_total(mesh, ssm.registry_classes)

    bed_mass_final = _bed_total_mass(ssm, len(schedule) - 1)
    delta_bed_total = bed_mass_final - bed_mass_initial

    # Global mass invariant: Δ(bed) + cumulative_src = 0
    # (everything that left the bed went into the source field as
    # "to-be-injected-into-water", and vice versa).
    np.testing.assert_allclose(
        delta_bed_total + cumulative_src, 0.0,
        atol=_TOL_GCM2, rtol=1e-3,
        err_msg=(
            f"Global mass balance violated: "
            f"Δbed_total={delta_bed_total}, cumulative_src={cumulative_src}, "
            f"sum={delta_bed_total + cumulative_src}"
        ),
    )


def test_per_layer_mass_conservation_invariant_holds(ssm):
    """The bed module's own mass-conservation invariant
    (``_MASS_CONSERVATION_TOL``) is checked inside
    :func:`reorganize_active_layer`; here we verify the global sum
    over all layers + steps stays bounded under repeated SSM steps.

    A regression that drops mass somewhere (e.g. miscoded PERSED
    renormalization) would manifest as a rapidly diverging total.
    """
    mesh = ssm._mesh
    n_face = mesh.sizes[contracts.DIM_NFACE]
    _seed_external_tau(mesh, n_face)

    bed_mass_history: list[float] = [_bed_total_mass(ssm, 0)]
    src_history: list[float] = []

    # Three steps under different shear levels, no suspended C.
    for t_idx, tau_value in enumerate([0.5, 1.0, 0.5]):
        _set_external_tau(mesh, t_idx, tau_value)
        _step_with_state_carried_forward(ssm, t_idx)
        bed_mass_history.append(_bed_total_mass(ssm, t_idx))
        src_history.append(_net_source_total(mesh, ssm.registry_classes))

    # Per-step Δ(bed) + src = 0 invariant.
    for i in range(1, len(bed_mass_history)):
        delta_bed = bed_mass_history[i] - bed_mass_history[i - 1]
        np.testing.assert_allclose(
            delta_bed + src_history[i - 1], 0.0,
            atol=_TOL_GCM2, rtol=1e-3,
            err_msg=f"step {i}: Δbed + src = {delta_bed + src_history[i-1]} ≠ 0",
        )
