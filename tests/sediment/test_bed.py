"""Unit tests for the multi-layer bed state and active-layer reorganization.

These exercise :mod:`clearwater_modules_v2.processes.sediment.bed`:

* allocation and IC population on a small mesh
* the three reorganization branches (deposition, borrow, promote)
* mass conservation across each branch
* bed-elevation update bookkeeping

Reference: design spec §5.6, §5.8, §7.4.
"""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v2.processes.sediment import (
    SedimentClass,
    SedimentClassRegistry,
)
from clearwater_modules_v2.processes.sediment import bed as bed_mod
from clearwater_modules_v2.processes.sediment import contracts


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def small_registry() -> SedimentClassRegistry:
    """2-class registry: one cohesive (silt), one non-cohesive (sand)."""
    return SedimentClassRegistry.from_iterable(
        [
            SedimentClass(label="silt_fine", d50_um=32.0),
            SedimentClass(label="sand_med", d50_um=250.0),
        ]
    )


@pytest.fixture
def small_mesh() -> xr.Dataset:
    """Empty mesh with 3 cells × 2 timesteps; no other variables."""
    return xr.Dataset(
        coords={
            contracts.DIM_TIME: np.arange(2, dtype="int32"),
            contracts.DIM_NFACE: np.arange(3, dtype="int32"),
        }
    )


@pytest.fixture
def initial_state_arrays() -> dict:
    """ICs for a 3-cell, 5-layer, 2-class bed.

    Layer 1 (active): light/empty in cell 0, deposited in cell 1, mid in cell 2.
    Layers 2..4: in-place core data.
    """
    n_face, n_layers, n_classes = 3, 5, 2
    layer_mass = np.array(
        [
            # cell 0
            [0.10, 1.0, 1.0, 1.0, 1.0],
            # cell 1: heavy layer-1 mass for branch (a) deposition test
            [0.50, 0.5, 1.0, 1.0, 1.0],
            # cell 2
            [0.20, 0.8, 1.0, 0.0, 0.0],
        ],
        dtype="float64",
    )
    # Equal split per layer per cell — but we want layer-2 fraction to differ
    # from layer-1 fraction so that the mass-weighted blend is non-trivial.
    class_fraction = np.zeros((n_face, n_layers, n_classes), dtype="float64")
    # Layer 1: silt-rich (0.8/0.2)
    class_fraction[:, 0, 0] = 0.8
    class_fraction[:, 0, 1] = 0.2
    # Layer 2: sand-rich (0.2/0.8)
    class_fraction[:, 1, 0] = 0.2
    class_fraction[:, 1, 1] = 0.8
    # Layers 3+: 50/50
    class_fraction[:, 2:, 0] = 0.5
    class_fraction[:, 2:, 1] = 0.5

    bulk_density = np.full((n_face, n_layers), 1.6, dtype="float64")

    layer_active = np.full((n_face, n_layers), bed_mod.LAYER_ACTIVE, dtype="int8")
    # Mark in-place layers as such (cell 2 has no layer 4/5 → mark absent).
    layer_active[:, 2:] = bed_mod.LAYER_IN_PLACE
    layer_active[2, 3:] = bed_mod.LAYER_ABSENT

    taucor = np.full((n_face, n_layers), 0.5, dtype="float64")  # Pa

    return dict(
        initial_layer_mass=layer_mass,
        initial_class_fraction=class_fraction,
        bulk_density=bulk_density,
        initial_layer_active=layer_active,
        taucor_initial=taucor,
    )


@pytest.fixture
def initialized_bed(small_mesh, small_registry, initial_state_arrays):
    """Bed state initialised on the small mesh."""
    return bed_mod.initialize_bed_state(
        mesh=small_mesh,
        registry=small_registry,
        n_layers=5,
        **initial_state_arrays,
    )


# ---------------------------------------------------------------------------
# initialize_bed_state
# ---------------------------------------------------------------------------


class TestInitializeBedState:
    def test_all_specs_present_with_correct_shape(
        self, small_mesh, small_registry, initial_state_arrays
    ):
        bed = bed_mod.initialize_bed_state(
            mesh=small_mesh,
            registry=small_registry,
            n_layers=5,
            **initial_state_arrays,
        )

        # Every BED_STATE_SPECS variable should be present with correct
        # dims, sizes, and dtype.
        n_face = small_mesh.sizes[contracts.DIM_NFACE]
        n_time = small_mesh.sizes[contracts.DIM_TIME]
        n_layers = 5
        n_classes = len(small_registry)
        size_lookup = {
            contracts.DIM_TIME: n_time,
            contracts.DIM_NFACE: n_face,
            contracts.DIM_LAYER: n_layers,
            contracts.DIM_CLASS: n_classes,
        }

        for spec in contracts.BED_STATE_SPECS:
            assert spec.name in small_mesh.data_vars, (
                f"missing bed-state variable {spec.name}"
            )
            da = small_mesh[spec.name]
            assert da.dims == spec.dims, (
                f"{spec.name} dims {da.dims} != expected {spec.dims}"
            )
            expected_shape = tuple(size_lookup[d] for d in spec.dims)
            assert da.shape == expected_shape, (
                f"{spec.name} shape {da.shape} != {expected_shape}"
            )
            assert da.dtype == np.dtype(spec.dtype), (
                f"{spec.name} dtype {da.dtype} != {spec.dtype}"
            )

        assert bed.n_layers == n_layers
        assert bed.n_classes == n_classes

    def test_initial_conditions_written_to_t0(self, initialized_bed, initial_state_arrays):
        np.testing.assert_allclose(
            initialized_bed.layer_mass_at(0).values,
            initial_state_arrays["initial_layer_mass"].astype("float32"),
        )
        np.testing.assert_allclose(
            initialized_bed.class_fraction_at(0).values,
            initial_state_arrays["initial_class_fraction"].astype("float32"),
        )
        np.testing.assert_array_equal(
            initialized_bed.layer_active_at(0).values,
            initial_state_arrays["initial_layer_active"],
        )
        np.testing.assert_allclose(
            initialized_bed.layer_bulk_density.values,
            initial_state_arrays["bulk_density"].astype("float32"),
        )

    def test_initial_thickness_consistent_with_mass_density(
        self, initialized_bed, initial_state_arrays
    ):
        expected = (
            0.01
            * initial_state_arrays["initial_layer_mass"]
            / initial_state_arrays["bulk_density"]
        ).astype("float32")
        np.testing.assert_allclose(
            initialized_bed.layer_thickness_at(0).values, expected, rtol=1e-6
        )

    def test_layer_and_class_dim_size_validation(
        self, small_mesh, small_registry, initial_state_arrays
    ):
        # Wrong layer count in IC → should raise.
        bad_mass = initial_state_arrays["initial_layer_mass"][:, :3]
        with pytest.raises(ValueError, match="initial_layer_mass shape"):
            bed_mod.initialize_bed_state(
                mesh=small_mesh.copy(),
                registry=small_registry,
                n_layers=5,
                initial_layer_mass=bad_mass,
                initial_class_fraction=initial_state_arrays["initial_class_fraction"],
                bulk_density=initial_state_arrays["bulk_density"],
                initial_layer_active=initial_state_arrays["initial_layer_active"],
                taucor_initial=initial_state_arrays["taucor_initial"],
            )


# ---------------------------------------------------------------------------
# reorganize_active_layer
# ---------------------------------------------------------------------------


def _make_single_cell_bed(
    layer_mass: np.ndarray,         # (n_layers,)
    class_fraction: np.ndarray,     # (n_layers, n_class)
    bulk_density: np.ndarray,       # (n_layers,)
    layer_active: np.ndarray,       # (n_layers,) int8
    layer_taucrit_pa: np.ndarray,   # (n_layers,) Pa
    registry: SedimentClassRegistry,
) -> bed_mod.BedState:
    """Build a 1-cell, n-layer BedState for branch-by-branch testing."""
    n_layers = layer_mass.size
    mesh = xr.Dataset(
        coords={
            contracts.DIM_TIME: np.arange(2, dtype="int32"),
            contracts.DIM_NFACE: np.arange(1, dtype="int32"),
        }
    )
    return bed_mod.initialize_bed_state(
        mesh=mesh,
        registry=registry,
        n_layers=n_layers,
        initial_layer_mass=layer_mass[None, :].copy(),
        initial_class_fraction=class_fraction[None, :, :].copy(),
        bulk_density=bulk_density[None, :].copy(),
        initial_layer_active=layer_active[None, :].copy(),
        taucor_initial=layer_taucrit_pa[None, :].copy(),
    )


class TestReorganizeActiveLayer:
    """Unit tests for the three branches of the SLLN reorganization."""

    def test_branch_a_deposition_pushes_excess_to_layer_2(self, small_registry):
        """m_1=0.5 g/cm², T_act → 0.3 g/cm² → 0.2 pushed to layer 2 (mass-weighted)."""
        # Setup: layer 1 has 0.5 g/cm², 80% silt / 20% sand.
        # Layer 2 has 1.0 g/cm², 20% silt / 80% sand.
        layer_mass = np.array([0.5, 1.0, 1.0, 0.0, 0.0])
        class_fraction = np.array(
            [
                [0.8, 0.2],   # layer 1
                [0.2, 0.8],   # layer 2
                [0.5, 0.5],
                [0.0, 0.0],
                [0.0, 0.0],
            ]
        )
        bulk_density = np.array([1.6, 1.6, 1.6, 1.6, 1.6])
        layer_active = np.array(
            [
                bed_mod.LAYER_ACTIVE,
                bed_mod.LAYER_ACTIVE,
                bed_mod.LAYER_IN_PLACE,
                bed_mod.LAYER_ABSENT,
                bed_mod.LAYER_ABSENT,
            ],
            dtype="int8",
        )
        taucrit = np.array([0.5, 0.5, 0.5, 1.0, 1.0])

        bed = _make_single_cell_bed(
            layer_mass, class_fraction, bulk_density, layer_active, taucrit, small_registry
        )

        # Choose tau, tau_crit, d50, bd1 so that T_act = 0.3 g/cm^2 exactly.
        # T_act = tactm * d50 * max(1, tau/tau_crit) * bd1 / 10000
        # For tau == tau_crit → factor = 1.
        # tactm=2, bd1=1.6 → T_act = 2 * d50 * 1 * 1.6 / 10000 = 3.2e-4 * d50
        # → d50 = 0.3 / 3.2e-4 = 937.5 μm
        tactm = 2.0
        d50 = 937.5
        tau = xr.DataArray(np.array([0.5]), dims=(contracts.DIM_NFACE,))
        tau_crit = xr.DataArray(np.array([0.5]), dims=(contracts.DIM_NFACE,))
        d50_arr = xr.DataArray(np.array([d50]), dims=(contracts.DIM_NFACE,))
        bd1 = xr.DataArray(np.array([1.6]), dims=(contracts.DIM_NFACE,))

        mass_before = bed.layer_mass_at(0).values.sum()

        bed_mod.reorganize_active_layer(
            bed=bed, t=0,
            tau_pa=tau, tau_crit_pa=tau_crit,
            d50_surface_um=d50_arr, bulk_density_layer1=bd1,
            tactm=tactm,
        )

        new_mass = bed.layer_mass_at(0).values[0, :]
        new_pers = bed.class_fraction_at(0).values[0, :, :]

        # Layer 1 capped at T_act = 0.3.
        np.testing.assert_allclose(new_mass[0], 0.3, atol=1e-6)
        # Layer 2 = original 1.0 + excess 0.2 = 1.2.
        np.testing.assert_allclose(new_mass[1], 1.2, atol=1e-6)
        # Layer-2 PERSED is mass-weighted blend of (1.0 of [0.2,0.8]) and
        # (0.2 of [0.8,0.2]):
        #   silt: (0.2*1.0 + 0.8*0.2) / 1.2 = (0.2 + 0.16)/1.2 = 0.3
        #   sand: (0.8*1.0 + 0.2*0.2) / 1.2 = (0.8 + 0.04)/1.2 = 0.7
        np.testing.assert_allclose(new_pers[1, 0], 0.3, atol=1e-6)
        np.testing.assert_allclose(new_pers[1, 1], 0.7, atol=1e-6)

        # Mass conservation: total layer mass unchanged.
        np.testing.assert_allclose(
            new_mass.sum(), mass_before, atol=1e-9
        )

    def test_branch_b_borrow_from_layer_2(self, small_registry):
        """m_1=0.1, T_act=0.5, m_2=1.0 with τ>τ_crit → borrow 0.4 from layer 2."""
        # Layer 1 = 0.1 g/cm², layer 2 = 1.0 g/cm² (silt 0.2 / sand 0.8),
        # layer 1 PERSED = silt 0.8 / sand 0.2.
        layer_mass = np.array([0.1, 1.0, 1.0, 0.0, 0.0])
        class_fraction = np.array(
            [
                [0.8, 0.2],   # layer 1
                [0.2, 0.8],   # layer 2 — donor
                [0.5, 0.5],
                [0.0, 0.0],
                [0.0, 0.0],
            ]
        )
        bulk_density = np.array([1.6, 1.6, 1.6, 1.6, 1.6])
        layer_active = np.array(
            [
                bed_mod.LAYER_ACTIVE,
                bed_mod.LAYER_ACTIVE,
                bed_mod.LAYER_IN_PLACE,
                bed_mod.LAYER_ABSENT,
                bed_mod.LAYER_ABSENT,
            ],
            dtype="int8",
        )
        # τ_crit at layer 2 = 0.3 Pa, so τ = 0.6 Pa → τ > τ_crit(SLLN).
        taucrit = np.array([0.5, 0.3, 0.5, 1.0, 1.0])

        bed = _make_single_cell_bed(
            layer_mass, class_fraction, bulk_density, layer_active, taucrit, small_registry
        )

        # Choose params so T_act = 0.5 g/cm^2 with the formula's τ/τ_crit
        # branch active. Use ratio = 1 for simplicity (set τ == τ_crit at
        # layer 1 = 0.5 Pa). Then T_act = 2 * d50 * 1 * 1.6 / 10000.
        # → d50 = 0.5 / (2 * 1.6 / 10000) = 0.5 / 3.2e-4 = 1562.5 μm.
        tactm = 2.0
        d50 = 1562.5
        tau = xr.DataArray(np.array([0.6]), dims=(contracts.DIM_NFACE,))
        # For T_act calc we want max(1, tau/tau_crit) = 1 → use surface
        # tau_crit = 0.6 (cancel with tau).
        tau_crit_surface = xr.DataArray(np.array([0.6]), dims=(contracts.DIM_NFACE,))
        d50_arr = xr.DataArray(np.array([d50]), dims=(contracts.DIM_NFACE,))
        bd1 = xr.DataArray(np.array([1.6]), dims=(contracts.DIM_NFACE,))

        mass_before = bed.layer_mass_at(0).values.sum()

        bed_mod.reorganize_active_layer(
            bed=bed, t=0,
            tau_pa=tau, tau_crit_pa=tau_crit_surface,
            d50_surface_um=d50_arr, bulk_density_layer1=bd1,
            tactm=tactm,
        )

        new_mass = bed.layer_mass_at(0).values[0, :]
        new_pers = bed.class_fraction_at(0).values[0, :, :]

        # Layer 1 should now hold T_act = 0.5; layer 2 should be 1.0 - 0.4 = 0.6.
        np.testing.assert_allclose(new_mass[0], 0.5, atol=1e-6)
        np.testing.assert_allclose(new_mass[1], 0.6, atol=1e-6)

        # Layer-1 PERSED is mass-weighted blend of (0.1 of [0.8,0.2]) +
        # (0.4 of [0.2,0.8]):
        #   silt: (0.8*0.1 + 0.2*0.4)/0.5 = (0.08 + 0.08)/0.5 = 0.32
        #   sand: (0.2*0.1 + 0.8*0.4)/0.5 = (0.02 + 0.32)/0.5 = 0.68
        np.testing.assert_allclose(new_pers[0, 0], 0.32, atol=1e-6)
        np.testing.assert_allclose(new_pers[0, 1], 0.68, atol=1e-6)

        # Mass conservation.
        np.testing.assert_allclose(new_mass.sum(), mass_before, atol=1e-9)

    def test_branch_c_promote_collapses_slln(self, small_registry):
        """m_1=0.1, T_act=0.5, m_2=0.1 (insufficient) → SLLN collapses into layer 1."""
        layer_mass = np.array([0.1, 0.1, 1.0, 0.0, 0.0])
        class_fraction = np.array(
            [
                [0.8, 0.2],   # layer 1
                [0.2, 0.8],   # layer 2 — donor, insufficient
                [0.5, 0.5],
                [0.0, 0.0],
                [0.0, 0.0],
            ]
        )
        bulk_density = np.array([1.6, 1.6, 1.6, 1.6, 1.6])
        layer_active = np.array(
            [
                bed_mod.LAYER_ACTIVE,
                bed_mod.LAYER_ACTIVE,
                bed_mod.LAYER_IN_PLACE,
                bed_mod.LAYER_ABSENT,
                bed_mod.LAYER_ABSENT,
            ],
            dtype="int8",
        )
        taucrit = np.array([0.5, 0.3, 0.5, 1.0, 1.0])

        bed = _make_single_cell_bed(
            layer_mass, class_fraction, bulk_density, layer_active, taucrit, small_registry
        )

        tactm = 2.0
        d50 = 1562.5  # gives T_act = 0.5 with bd1=1.6, factor=1
        tau = xr.DataArray(np.array([0.6]), dims=(contracts.DIM_NFACE,))
        tau_crit_surface = xr.DataArray(np.array([0.6]), dims=(contracts.DIM_NFACE,))
        d50_arr = xr.DataArray(np.array([d50]), dims=(contracts.DIM_NFACE,))
        bd1 = xr.DataArray(np.array([1.6]), dims=(contracts.DIM_NFACE,))

        mass_before = bed.layer_mass_at(0).values.sum()

        bed_mod.reorganize_active_layer(
            bed=bed, t=0,
            tau_pa=tau, tau_crit_pa=tau_crit_surface,
            d50_surface_um=d50_arr, bulk_density_layer1=bd1,
            tactm=tactm,
        )

        new_mass = bed.layer_mass_at(0).values[0, :]
        new_active = bed.layer_active_at(0).values[0, :]

        # Layer 1 = 0.1 + 0.1 = 0.2 (still less than T_act = 0.5, but no more
        # mass available immediately below).
        np.testing.assert_allclose(new_mass[0], 0.2, atol=1e-6)
        # Layer 2 zeroed out.
        np.testing.assert_allclose(new_mass[1], 0.0, atol=1e-12)
        assert new_active[1] == bed_mod.LAYER_ABSENT

        # Layer-1 PERSED is mass-weighted blend of (0.1 of [0.8,0.2]) +
        # (0.1 of [0.2,0.8]) = [0.5, 0.5].
        new_pers = bed.class_fraction_at(0).values[0, :, :]
        np.testing.assert_allclose(new_pers[0, 0], 0.5, atol=1e-6)
        np.testing.assert_allclose(new_pers[0, 1], 0.5, atol=1e-6)

        # Mass conservation.
        np.testing.assert_allclose(new_mass.sum(), mass_before, atol=1e-9)

    def test_mass_conservation_vectorized_across_cells(self, initialized_bed):
        """Run reorganization on the full 3-cell mesh and assert sum invariance."""
        n_face = initialized_bed.mesh.sizes[contracts.DIM_NFACE]

        tau = xr.DataArray(np.full(n_face, 0.6), dims=(contracts.DIM_NFACE,))
        tau_crit = xr.DataArray(np.full(n_face, 0.4), dims=(contracts.DIM_NFACE,))
        d50_arr = xr.DataArray(np.full(n_face, 500.0), dims=(contracts.DIM_NFACE,))
        bd1 = xr.DataArray(np.full(n_face, 1.6), dims=(contracts.DIM_NFACE,))

        mass_before_per_cell = (
            initialized_bed.layer_mass_at(0).values.sum(axis=-1).copy()
        )

        bed_mod.reorganize_active_layer(
            bed=initialized_bed, t=0,
            tau_pa=tau, tau_crit_pa=tau_crit,
            d50_surface_um=d50_arr, bulk_density_layer1=bd1,
            tactm=contracts.DEFAULT_TACTM,
        )

        mass_after_per_cell = initialized_bed.layer_mass_at(0).values.sum(axis=-1)
        np.testing.assert_allclose(
            mass_after_per_cell, mass_before_per_cell, atol=1e-9
        )

    def test_no_op_when_layer_1_matches_t_act(self, small_registry):
        """If m_1 == T_act and τ ≤ τ_crit, reorganization is a no-op."""
        layer_mass = np.array([0.5, 1.0, 1.0, 0.0, 0.0])
        class_fraction = np.array(
            [
                [0.8, 0.2],
                [0.2, 0.8],
                [0.5, 0.5],
                [0.0, 0.0],
                [0.0, 0.0],
            ]
        )
        bulk_density = np.array([1.6] * 5)
        layer_active = np.array(
            [
                bed_mod.LAYER_ACTIVE,
                bed_mod.LAYER_ACTIVE,
                bed_mod.LAYER_IN_PLACE,
                bed_mod.LAYER_ABSENT,
                bed_mod.LAYER_ABSENT,
            ],
            dtype="int8",
        )
        taucrit = np.array([0.5, 0.5, 0.5, 1.0, 1.0])

        bed = _make_single_cell_bed(
            layer_mass, class_fraction, bulk_density, layer_active, taucrit, small_registry
        )

        # Pick d50 → T_act = 0.5 exactly with factor=1.
        # T_act = 2 * d50 * 1 * 1.6 / 10000 → d50 = 1562.5
        tau = xr.DataArray(np.array([0.4]), dims=(contracts.DIM_NFACE,))   # τ < τ_crit
        tau_crit = xr.DataArray(np.array([0.5]), dims=(contracts.DIM_NFACE,))
        d50_arr = xr.DataArray(np.array([1562.5]), dims=(contracts.DIM_NFACE,))
        bd1 = xr.DataArray(np.array([1.6]), dims=(contracts.DIM_NFACE,))

        mass_before = bed.layer_mass_at(0).values.copy()
        pers_before = bed.class_fraction_at(0).values.copy()

        bed_mod.reorganize_active_layer(
            bed=bed, t=0,
            tau_pa=tau, tau_crit_pa=tau_crit,
            d50_surface_um=d50_arr, bulk_density_layer1=bd1,
            tactm=2.0,
        )

        np.testing.assert_allclose(bed.layer_mass_at(0).values, mass_before, atol=1e-7)
        np.testing.assert_allclose(bed.class_fraction_at(0).values, pers_before, atol=1e-7)


# ---------------------------------------------------------------------------
# update_bed_elevation
# ---------------------------------------------------------------------------


class TestUpdateBedElevation:
    def test_thickness_for_known_masses_and_density(
        self, initialized_bed, initial_state_arrays
    ):
        # Already verified by the IC test, but recompute explicitly here
        # using the function under test.
        bed_mod.update_bed_elevation(initialized_bed, t=0)

        layer_thickness = initialized_bed.layer_thickness_at(0).values
        expected = (
            0.01
            * initial_state_arrays["initial_layer_mass"]
            / initial_state_arrays["bulk_density"]
        ).astype("float32")
        np.testing.assert_allclose(layer_thickness, expected, rtol=1e-6)

        total_thickness = (
            initialized_bed.mesh[contracts.VAR_BED_TOTAL_THICKNESS]
            .isel({contracts.DIM_TIME: 0})
            .values
        )
        np.testing.assert_allclose(
            total_thickness, expected.sum(axis=-1), rtol=1e-6
        )

    def test_bed_change_zero_at_t0(self, initialized_bed):
        bed_mod.update_bed_elevation(initialized_bed, t=0)
        bed_change = (
            initialized_bed.mesh[contracts.VAR_BED_CHANGE]
            .isel({contracts.DIM_TIME: 0})
            .values
        )
        np.testing.assert_allclose(bed_change, 0.0, atol=1e-12)

        cum = (
            initialized_bed.mesh[contracts.VAR_BED_CUMULATIVE_CHANGE]
            .isel({contracts.DIM_TIME: 0})
            .values
        )
        np.testing.assert_allclose(cum, 0.0, atol=1e-12)

    def test_bed_change_at_t1_matches_delta(self, initialized_bed, initial_state_arrays):
        # Initialize t=0 elevation, then nudge mass at t=1 and recompute.
        bed_mod.update_bed_elevation(initialized_bed, t=0)

        # Copy t=0 mass to t=1 with cell-0 layer-1 mass increased by 0.1 g/cm².
        m1 = initial_state_arrays["initial_layer_mass"].copy()
        m1[0, 0] += 0.1
        initialized_bed.set_layer_mass_at(1, m1)

        bed_mod.update_bed_elevation(initialized_bed, t=1)

        bed_change = (
            initialized_bed.mesh[contracts.VAR_BED_CHANGE]
            .isel({contracts.DIM_TIME: 1})
            .values
        )

        # Per-cell delta: 0.01 * 0.1 / 1.6 = 6.25e-4 m on cell 0 only.
        expected = np.zeros(3, dtype="float32")
        expected[0] = 0.01 * 0.1 / 1.6
        np.testing.assert_allclose(bed_change, expected, rtol=1e-5, atol=1e-9)

        # Cumulative at t=1 equals the delta (since t=0 was zero).
        cum = (
            initialized_bed.mesh[contracts.VAR_BED_CUMULATIVE_CHANGE]
            .isel({contracts.DIM_TIME: 1})
            .values
        )
        np.testing.assert_allclose(cum, expected, rtol=1e-5, atol=1e-9)
