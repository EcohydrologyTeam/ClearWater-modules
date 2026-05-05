"""Unit tests for the cohesive-bed consolidation module.

Covers :mod:`clearwater_modules_v2.processes.sediment.consolidation`:

* The Sanford-Maa (2001) τ_ce(age) envelope returns the freshly-deposited
  value at age 0, the fully-consolidated value at large age, and the
  expected (1-1/e) intermediate at age ``T_c``.
* Per-class application correctly leaves non-cohesive classes alone.
* Layer-age advances by ``dt`` each call to :func:`update_bed_elevation`.
* Age-dilution on deposition gives the mass-weighted result.
* Borrow / promote / collapse correctly propagate ages.
* Integration: a single-cell single-class run with sustained low shear
  plus periodic deposition shows τ_ce rising over time.

References
----------
* Sanford, L. P., and Maa, J. P.-Y. (2001). *Marine Geology* 179, 9-23.
* Mehta, A. J., and Partheniades, E. (1975). *J. Hydraul. Res.* 13, 361-381.
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
from clearwater_modules_v2.processes.sediment import consolidation as consol_mod
from clearwater_modules_v2.processes.sediment import erosion as erosion_mod


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def sm_model() -> consol_mod.SanfordMaaConsolidation:
    """Sanford-Maa model with the package-default calibration."""
    return consol_mod.SanfordMaaConsolidation(
        tau_ce_zero_pa=contracts.DEFAULT_CONSOLIDATION_TAU_CE_ZERO_PA,
        tau_ce_inf_pa=contracts.DEFAULT_CONSOLIDATION_TAU_CE_INF_PA,
        consolidation_time_s=contracts.DEFAULT_CONSOLIDATION_TIME_S,
    )


@pytest.fixture
def cohesive_only_registry() -> SedimentClassRegistry:
    """Single-class cohesive registry (silt)."""
    return SedimentClassRegistry.from_iterable(
        [SedimentClass(label="silt_fine", d50_um=16.0, tau_ce_pa=0.10)]
    )


@pytest.fixture
def two_class_registry() -> SedimentClassRegistry:
    """One cohesive class + one non-cohesive (sand)."""
    return SedimentClassRegistry.from_iterable(
        [
            SedimentClass(label="silt_fine", d50_um=16.0, tau_ce_pa=0.10),
            SedimentClass(label="sand_med", d50_um=250.0, tau_ce_pa=0.30),
        ]
    )


# ---------------------------------------------------------------------------
# SanfordMaaConsolidation — unit tests
# ---------------------------------------------------------------------------


class TestSanfordMaaTauCe:
    def test_age_zero_returns_tau_ce_zero(self, sm_model):
        ages = xr.DataArray(np.array([[0.0, 0.0, 0.0]]), dims=("nface", "ssm_layer"))
        eff = sm_model.effective_tau_ce(ages)
        np.testing.assert_allclose(
            eff.values, sm_model.tau_ce_zero_pa, atol=1e-12
        )

    def test_age_large_approaches_tau_ce_inf(self, sm_model):
        # 50 × T_c → exp(-50) ≈ 1.9e-22; effectively the asymptote.
        ages = xr.DataArray(
            np.array([[50.0 * sm_model.consolidation_time_s]]),
            dims=("nface", "ssm_layer"),
        )
        eff = sm_model.effective_tau_ce(ages)
        np.testing.assert_allclose(
            eff.values, sm_model.tau_ce_inf_pa, atol=1e-10
        )

    def test_age_at_tc_matches_one_minus_one_over_e(self, sm_model):
        """At t = T_c the recovery is (1 - 1/e) of the (∞ - 0) gap."""
        ages = xr.DataArray(
            np.array([[sm_model.consolidation_time_s]]),
            dims=("nface", "ssm_layer"),
        )
        eff = sm_model.effective_tau_ce(ages)
        delta = sm_model.tau_ce_inf_pa - sm_model.tau_ce_zero_pa
        expected = sm_model.tau_ce_zero_pa + delta * (1.0 - 1.0 / np.e)
        np.testing.assert_allclose(eff.values, expected, atol=1e-10)

    def test_negative_ages_clamped_to_zero(self, sm_model):
        ages = xr.DataArray(
            np.array([[-1.0e6, -10.0, 0.0]]), dims=("nface", "ssm_layer")
        )
        eff = sm_model.effective_tau_ce(ages)
        np.testing.assert_allclose(
            eff.values, sm_model.tau_ce_zero_pa, atol=1e-12
        )

    def test_monotonic_increasing_in_age(self, sm_model):
        ages = xr.DataArray(
            np.linspace(0.0, 30.0 * sm_model.consolidation_time_s, 50)[None, :],
            dims=("nface", "ssm_layer"),
        )
        eff = sm_model.effective_tau_ce(ages).values[0]
        assert np.all(np.diff(eff) >= -1e-12), "τ_ce(age) must be non-decreasing"

    def test_validation_rejects_nonpositive_tc(self):
        with pytest.raises(ValueError, match="consolidation_time_s"):
            consol_mod.SanfordMaaConsolidation(
                tau_ce_zero_pa=0.1, tau_ce_inf_pa=0.5, consolidation_time_s=0.0
            )

    def test_validation_rejects_inf_below_zero(self):
        with pytest.raises(ValueError, match="tau_ce_inf_pa"):
            consol_mod.SanfordMaaConsolidation(
                tau_ce_zero_pa=0.5, tau_ce_inf_pa=0.1, consolidation_time_s=86400.0
            )


# ---------------------------------------------------------------------------
# apply_consolidation_per_class — unit tests
# ---------------------------------------------------------------------------


class TestApplyConsolidationPerClass:
    def test_cohesive_class_uses_aged_value(self, sm_model):
        # 1 face, 2 layers, 2 classes; class 0 cohesive, class 1 not.
        n_face, n_layer, n_class = 1, 2, 2
        baseline = xr.DataArray(
            np.full((n_face, n_layer, n_class), 0.30, dtype="float64"),
            dims=("nface", "ssm_layer", "ssm_class"),
        )
        ages = xr.DataArray(
            np.array([[0.0, sm_model.consolidation_time_s]]),
            dims=("nface", "ssm_layer"),
        )
        is_cohesive = np.array([True, False])
        out = consol_mod.apply_consolidation_per_class(
            baseline, ages, is_cohesive, sm_model
        )

        # Cohesive class (idx 0): layer 0 → tau_ce_zero, layer 1 → 1-1/e form.
        np.testing.assert_allclose(
            out.values[0, 0, 0], sm_model.tau_ce_zero_pa, atol=1e-12
        )
        delta = sm_model.tau_ce_inf_pa - sm_model.tau_ce_zero_pa
        expected_layer1 = sm_model.tau_ce_zero_pa + delta * (1.0 - 1.0 / np.e)
        np.testing.assert_allclose(out.values[0, 1, 0], expected_layer1, atol=1e-10)

        # Non-cohesive class (idx 1): unchanged at both layers.
        np.testing.assert_allclose(out.values[0, :, 1], 0.30, atol=1e-12)

    def test_apply_consolidation_helper_passthrough_when_model_none(self):
        baseline = xr.DataArray(
            np.full((2, 3, 2), 0.5),
            dims=("nface", "ssm_layer", "ssm_class"),
        )
        ages = xr.DataArray(np.zeros((2, 3)), dims=("nface", "ssm_layer"))
        out = erosion_mod.apply_consolidation(
            baseline, ages, np.array([True, False]), consolidation_model=None
        )
        # Same object returned (no copy needed).
        assert out is baseline


# ---------------------------------------------------------------------------
# Age-dilution on deposition — unit test
# ---------------------------------------------------------------------------


def _make_single_cell_two_layer_bed(
    layer_mass: np.ndarray,
    age0: np.ndarray,
    registry: SedimentClassRegistry,
) -> bed_mod.BedState:
    n_layers = layer_mass.size
    n_class = len(registry)
    mesh = xr.Dataset(
        coords={
            contracts.DIM_TIME: np.arange(2, dtype="int32"),
            contracts.DIM_NFACE: np.arange(1, dtype="int32"),
        }
    )
    pers = np.zeros((1, n_layers, n_class), dtype="float64")
    pers[..., 0] = 1.0  # all silt
    bd = np.full((1, n_layers), 1.6, dtype="float64")
    active = np.where(
        layer_mass > 0.0, bed_mod.LAYER_ACTIVE, bed_mod.LAYER_ABSENT
    ).astype("int8")[None, :]
    taucor = np.full((1, n_layers), 0.5, dtype="float64")

    bed = bed_mod.initialize_bed_state(
        mesh=mesh,
        registry=registry,
        n_layers=n_layers,
        initial_layer_mass=layer_mass[None, :].astype("float64"),
        initial_class_fraction=pers,
        bulk_density=bd,
        initial_layer_active=active,
        taucor_initial=taucor,
    )
    bed.set_layer_age_at(0, age0[None, :].astype("float64"))
    return bed


class TestAgeDilutionOnDeposition:
    def test_simple_two_to_one_dilution(self, cohesive_only_registry):
        # m1 = 1.0 g/cm² at age 100 s; deposit Δm = 1.0 g/cm² with age 0.
        # New age = 100 * (1 / (1 + 1)) = 50 s.
        bed = _make_single_cell_two_layer_bed(
            np.array([1.0, 0.0]),
            np.array([100.0, 0.0]),
            cohesive_only_registry,
        )
        bed_mod.dilute_layer1_age_on_deposition(
            bed, t=0,
            layer1_mass_before=np.array([1.0]),
            deposited_mass=np.array([1.0]),
        )
        new_age = bed.layer_age_at(0).values[0, 0]
        assert new_age == pytest.approx(50.0, rel=1e-6)

    def test_zero_deposit_leaves_age_unchanged(self, cohesive_only_registry):
        bed = _make_single_cell_two_layer_bed(
            np.array([2.0, 0.0]),
            np.array([300.0, 0.0]),
            cohesive_only_registry,
        )
        bed_mod.dilute_layer1_age_on_deposition(
            bed, t=0,
            layer1_mass_before=np.array([2.0]),
            deposited_mass=np.array([0.0]),
        )
        assert bed.layer_age_at(0).values[0, 0] == pytest.approx(300.0)

    def test_deposit_into_empty_layer_yields_zero_age(self, cohesive_only_registry):
        bed = _make_single_cell_two_layer_bed(
            np.array([0.0, 1.0]),
            np.array([0.0, 100.0]),
            cohesive_only_registry,
        )
        bed_mod.dilute_layer1_age_on_deposition(
            bed, t=0,
            layer1_mass_before=np.array([0.0]),
            deposited_mass=np.array([0.5]),
        )
        # Old age 0 × (0/0.5) = 0.
        assert bed.layer_age_at(0).values[0, 0] == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# Age advancement in update_bed_elevation
# ---------------------------------------------------------------------------


class TestAgeAdvancement:
    def test_age_advances_by_dt_for_layers_with_mass(self, cohesive_only_registry):
        bed = _make_single_cell_two_layer_bed(
            np.array([1.0, 0.5]),
            np.array([10.0, 20.0]),
            cohesive_only_registry,
        )
        dt = 60.0
        bed_mod.update_bed_elevation(bed, t=0, dt_seconds=dt)
        new_age = bed.layer_age_at(0).values[0, :]
        np.testing.assert_allclose(new_age, [10.0 + dt, 20.0 + dt], atol=1e-9)

    def test_empty_layer_age_pinned_to_zero(self, cohesive_only_registry):
        bed = _make_single_cell_two_layer_bed(
            np.array([1.0, 0.0]),
            np.array([10.0, 999.0]),  # ghost age on empty layer
            cohesive_only_registry,
        )
        bed_mod.update_bed_elevation(bed, t=0, dt_seconds=60.0)
        new_age = bed.layer_age_at(0).values[0, :]
        # Layer 0 advanced; layer 1 reset to zero (pin).
        assert new_age[0] == pytest.approx(70.0)
        assert new_age[1] == pytest.approx(0.0)

    def test_dt_zero_is_a_no_op(self, cohesive_only_registry):
        bed = _make_single_cell_two_layer_bed(
            np.array([1.0, 0.5]),
            np.array([10.0, 20.0]),
            cohesive_only_registry,
        )
        bed_mod.update_bed_elevation(bed, t=0, dt_seconds=0.0)
        new_age = bed.layer_age_at(0).values[0, :]
        np.testing.assert_allclose(new_age, [10.0, 20.0], atol=1e-12)


# ---------------------------------------------------------------------------
# Borrow / promote / collapse age inheritance
# ---------------------------------------------------------------------------


def _make_small_bed_for_branch(
    registry: SedimentClassRegistry,
    layer_mass: np.ndarray,
    layer_age: np.ndarray,
    layer_active: np.ndarray,
    layer_taucrit: np.ndarray,
) -> bed_mod.BedState:
    n_layers = layer_mass.size
    n_class = len(registry)
    mesh = xr.Dataset(
        coords={
            contracts.DIM_TIME: np.arange(2, dtype="int32"),
            contracts.DIM_NFACE: np.arange(1, dtype="int32"),
        }
    )
    pers = np.zeros((1, n_layers, n_class), dtype="float64")
    pers[..., 0] = 1.0  # all silt to keep the per-class blend trivial
    bd = np.full((1, n_layers), 1.6, dtype="float64")

    bed = bed_mod.initialize_bed_state(
        mesh=mesh,
        registry=registry,
        n_layers=n_layers,
        initial_layer_mass=layer_mass[None, :].astype("float64"),
        initial_class_fraction=pers,
        bulk_density=bd,
        initial_layer_active=layer_active[None, :].astype("int8"),
        taucor_initial=layer_taucrit[None, :].astype("float64"),
    )
    bed.set_layer_age_at(0, layer_age[None, :].astype("float64"))
    return bed


class TestAgeInheritanceInReorganization:
    def test_branch_a_pushes_layer1_age_into_layer2(self, cohesive_only_registry):
        # m1 = 0.5 (age 100), m2 = 1.0 (age 1000); excess 0.2 → layer 2.
        # Layer 2 new age = (1000*1.0 + 100*0.2) / 1.2 = 1020/1.2 = 850.
        layer_mass = np.array([0.5, 1.0, 1.0, 0.0, 0.0])
        layer_age = np.array([100.0, 1000.0, 0.0, 0.0, 0.0])
        layer_active = np.array(
            [
                bed_mod.LAYER_ACTIVE,
                bed_mod.LAYER_ACTIVE,
                bed_mod.LAYER_IN_PLACE,
                bed_mod.LAYER_ABSENT,
                bed_mod.LAYER_ABSENT,
            ]
        )
        taucrit = np.array([0.5, 0.5, 0.5, 1.0, 1.0])
        bed = _make_small_bed_for_branch(
            cohesive_only_registry, layer_mass, layer_age, layer_active, taucrit
        )

        # T_act = 0.3 → excess 0.2 (matches test_bed branch (a) calibration).
        tau = xr.DataArray(np.array([0.5]), dims=(contracts.DIM_NFACE,))
        tau_crit = xr.DataArray(np.array([0.5]), dims=(contracts.DIM_NFACE,))
        d50 = xr.DataArray(np.array([937.5]), dims=(contracts.DIM_NFACE,))
        bd1 = xr.DataArray(np.array([1.6]), dims=(contracts.DIM_NFACE,))

        bed_mod.reorganize_active_layer(
            bed=bed, t=0,
            tau_pa=tau, tau_crit_pa=tau_crit,
            d50_surface_um=d50, bulk_density_layer1=bd1,
            tactm=2.0,
        )

        new_age = bed.layer_age_at(0).values[0, :]
        # Layer 1 unchanged (uniform-aged mass removed off the top).
        assert new_age[0] == pytest.approx(100.0, rel=1e-6)
        # Layer 2 = (1000 × 1 + 100 × 0.2) / 1.2 = 850.
        assert new_age[1] == pytest.approx(850.0, rel=1e-5)

    def test_branch_b_borrow_blends_slln_age_into_layer1(self, cohesive_only_registry):
        # m1 = 0.1 (age 50), m2 = 1.0 (age 5000); deficit 0.4 borrowed.
        # New layer-1 age = (50 × 0.1 + 5000 × 0.4) / 0.5 = (5 + 2000) / 0.5 = 4010.
        layer_mass = np.array([0.1, 1.0, 1.0, 0.0, 0.0])
        layer_age = np.array([50.0, 5000.0, 0.0, 0.0, 0.0])
        layer_active = np.array(
            [
                bed_mod.LAYER_ACTIVE,
                bed_mod.LAYER_ACTIVE,
                bed_mod.LAYER_IN_PLACE,
                bed_mod.LAYER_ABSENT,
                bed_mod.LAYER_ABSENT,
            ]
        )
        taucrit = np.array([0.5, 0.3, 0.5, 1.0, 1.0])
        bed = _make_small_bed_for_branch(
            cohesive_only_registry, layer_mass, layer_age, layer_active, taucrit
        )

        # Same calibration as test_bed branch (b).
        tau = xr.DataArray(np.array([0.6]), dims=(contracts.DIM_NFACE,))
        tau_crit_surface = xr.DataArray(
            np.array([0.6]), dims=(contracts.DIM_NFACE,)
        )
        d50 = xr.DataArray(np.array([1562.5]), dims=(contracts.DIM_NFACE,))
        bd1 = xr.DataArray(np.array([1.6]), dims=(contracts.DIM_NFACE,))

        bed_mod.reorganize_active_layer(
            bed=bed, t=0,
            tau_pa=tau, tau_crit_pa=tau_crit_surface,
            d50_surface_um=d50, bulk_density_layer1=bd1,
            tactm=2.0,
        )

        new_age = bed.layer_age_at(0).values[0, :]
        # Layer 1: (50*0.1 + 5000*0.4)/0.5 = 4010.
        assert new_age[0] == pytest.approx(4010.0, rel=1e-6)
        # Layer 2 unchanged at 5000 (uniform-aged donor).
        assert new_age[1] == pytest.approx(5000.0, rel=1e-6)

    def test_branch_c_collapse_blends_full_slln_into_layer1(self, cohesive_only_registry):
        # m1 = 0.1 (age 50), m2 = 0.1 (age 5000); collapse → m1 = 0.2.
        # New age = (50*0.1 + 5000*0.1)/0.2 = (5 + 500)/0.2 = 2525.
        layer_mass = np.array([0.1, 0.1, 1.0, 0.0, 0.0])
        layer_age = np.array([50.0, 5000.0, 0.0, 0.0, 0.0])
        layer_active = np.array(
            [
                bed_mod.LAYER_ACTIVE,
                bed_mod.LAYER_ACTIVE,
                bed_mod.LAYER_IN_PLACE,
                bed_mod.LAYER_ABSENT,
                bed_mod.LAYER_ABSENT,
            ]
        )
        taucrit = np.array([0.5, 0.3, 0.5, 1.0, 1.0])
        bed = _make_small_bed_for_branch(
            cohesive_only_registry, layer_mass, layer_age, layer_active, taucrit
        )

        tau = xr.DataArray(np.array([0.6]), dims=(contracts.DIM_NFACE,))
        tau_crit_surface = xr.DataArray(
            np.array([0.6]), dims=(contracts.DIM_NFACE,)
        )
        d50 = xr.DataArray(np.array([1562.5]), dims=(contracts.DIM_NFACE,))
        bd1 = xr.DataArray(np.array([1.6]), dims=(contracts.DIM_NFACE,))

        bed_mod.reorganize_active_layer(
            bed=bed, t=0,
            tau_pa=tau, tau_crit_pa=tau_crit_surface,
            d50_surface_um=d50, bulk_density_layer1=bd1,
            tactm=2.0,
        )

        new_age = bed.layer_age_at(0).values[0, :]
        assert new_age[0] == pytest.approx(2525.0, rel=1e-6)
        # Layer 2 zeroed (absent now).
        assert new_age[1] == pytest.approx(0.0, abs=1e-12)


# ---------------------------------------------------------------------------
# Integration: τ_ce rises with sustained low-shear deposition
# ---------------------------------------------------------------------------


class TestConsolidationIntegration:
    def test_tau_ce_rises_over_many_steps(self, sm_model, cohesive_only_registry):
        """Single-cell, single-class run: layer 1 ages over many no-erosion steps;
        consolidation-aged τ_ce rises monotonically toward the asymptote."""
        n_steps = 30  # 30 days at dt = 1 day
        n_layers = 2
        dt = 86400.0

        mesh = xr.Dataset(
            coords={
                contracts.DIM_TIME: np.arange(n_steps + 1, dtype="int32"),
                contracts.DIM_NFACE: np.arange(1, dtype="int32"),
            }
        )
        layer_mass = np.array([[1.0, 0.0]])
        pers = np.zeros((1, n_layers, 1))
        pers[..., 0] = 1.0
        bd = np.full((1, n_layers), 1.6)
        active = np.array([[bed_mod.LAYER_ACTIVE, bed_mod.LAYER_ABSENT]], dtype="int8")
        taucor = np.full((1, n_layers), 0.10)

        bed = bed_mod.initialize_bed_state(
            mesh=mesh,
            registry=cohesive_only_registry,
            n_layers=n_layers,
            initial_layer_mass=layer_mass,
            initial_class_fraction=pers,
            bulk_density=bd,
            initial_layer_active=active,
            taucor_initial=taucor,
        )
        # Start age at zero everywhere (default).

        tau_ce_history = []
        for t_idx in range(n_steps + 1):
            if t_idx > 0:
                # Carry layer state forward (mimic ssm orchestration).
                bed.set_layer_mass_at(t_idx, bed.layer_mass_at(t_idx - 1).values)
                bed.set_class_fraction_at(
                    t_idx, bed.class_fraction_at(t_idx - 1).values
                )
                bed.set_layer_active_at(
                    t_idx, bed.layer_active_at(t_idx - 1).values
                )
                bed.set_layer_age_at(t_idx, bed.layer_age_at(t_idx - 1).values)
            # Advance age by dt.
            bed_mod.update_bed_elevation(bed, t=t_idx, dt_seconds=dt)
            # Read age at layer 0 and compute effective τ_ce.
            age = bed.layer_age_at(t_idx)
            eff = sm_model.effective_tau_ce(age)
            tau_ce_history.append(float(eff.values[0, 0]))

        tau_arr = np.asarray(tau_ce_history)
        # Monotonic non-decreasing.
        assert np.all(np.diff(tau_arr) >= -1e-12)
        # First step (just after dt added) is > tau_ce_zero.
        assert tau_arr[1] > sm_model.tau_ce_zero_pa
        # Final step (30 days, with T_c = 7 days → 30/7 ≈ 4.3) should be
        # well above τ_ce_zero and approaching τ_ce_inf.
        assert tau_arr[-1] > sm_model.tau_ce_zero_pa + 0.9 * (
            sm_model.tau_ce_inf_pa - sm_model.tau_ce_zero_pa
        )
        assert tau_arr[-1] < sm_model.tau_ce_inf_pa + 1e-9
