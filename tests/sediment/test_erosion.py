"""Unit tests for SSM erosion-rate models.

Covers:
* :class:`PowerLawErosionModel` — single-class, single-cell sanity checks
  against the closed-form ``E = A·τ^n·ρ_b`` for a fully-fresh layer, plus
  vectorization across multiple cells.
* :class:`SedflumeTableErosionModel` — bilinear (log-depth, linear-τ)
  interpolation against a synthetic 2-shear × 3-layer × 1-core table with
  values chosen so that exp/log midpoint identities give analytic answers.
* :func:`apply_vegetation_cohesion` — broadcasting of biostabilization and
  root-cohesion across the class dimension.

The expected values are derived from the SAND2008-5621 §S_SEDZLJ.f90
formulae, not from the EFDC+ Fortran source.
"""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v2.processes.sediment.erosion import (
    PowerLawErosionModel,
    SedflumeTableErosionModel,
    apply_vegetation_cohesion,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _da(values, dims=("nface",)):
    """Wrap a Python iterable / scalar as an xarray DataArray on ``dims``."""
    arr = np.asarray(values, dtype=np.float64)
    return xr.DataArray(arr, dims=dims)


# ---------------------------------------------------------------------------
# PowerLawErosionModel
# ---------------------------------------------------------------------------

class TestPowerLawErosionModel:
    """Tests for ``nsedflume=2`` per-core per-layer power law."""

    @pytest.fixture
    def synthetic_model(self) -> PowerLawErosionModel:
        """1 core × 3 layers, A=1e-4, n=2.0 everywhere, cap=1.0 g/cm²/s."""
        n_cores, k_b = 1, 3
        ea = np.full((n_cores, k_b), 1.0e-4)
        en = np.full((n_cores, k_b), 2.0)
        max_rate = np.full((n_cores, k_b), 1.0)
        # Active/deposited tables — not exercised here but required by ctor.
        nsicm = 2
        actdep_a = np.full((nsicm,), 5.0e-5)
        actdep_n = np.full((nsicm,), 2.0)
        actdep_max = np.full((nsicm,), 0.5)
        return PowerLawErosionModel(
            ea_per_core=ea,
            en_per_core=en,
            max_rate_per_core=max_rate,
            actdep_a=actdep_a,
            actdep_n=actdep_n,
            actdep_max=actdep_max,
        )

    def test_full_layer_dominates_e_top(self, synthetic_model):
        """Spec example: τ=2 Pa, A=1e-4, n=2.0, ρ_b=1.6, full layer.

        ``E_top = A · τ^n · ρ_b = 1e-4 × 4 × 1.6 = 6.4e-4 g/cm²/s``.
        With m_K = m_K0 the SN11 weighting collapses to E_top exactly.
        """
        tau = _da([2.0])
        layer_mass = _da([15.0])           # m_K
        layer_initial_mass = _da([15.0])   # m_K0 = m_K → SN11 = 0
        bulk_density = _da([1.6])
        core_id = xr.DataArray(np.array([0], dtype=np.int64), dims=("nface",))

        result = synthetic_model.erosion_rate(
            tau_pa=tau,
            layer_index=1,
            layer_mass=layer_mass,
            layer_initial_mass=layer_initial_mass,
            bulk_density=bulk_density,
            core_id=core_id,
        )

        expected = 1.0e-4 * (2.0 ** 2.0) * 1.6  # = 6.4e-4 g/cm²/s
        np.testing.assert_allclose(result.values, [expected], rtol=1e-12)
        assert result.dims == ("nface",)

    def test_deepest_layer_no_through_bottom_erosion(self, synthetic_model):
        """For the deepest layer, E_bottom = 0 → an empty layer erodes at 0."""
        tau = _da([2.0])
        # m_K = 0 → SN11 = 1 → E = E_bottom = 0 (deepest layer).
        layer_mass = _da([0.0])
        layer_initial_mass = _da([15.0])
        bulk_density = _da([1.6])
        core_id = xr.DataArray(np.array([0], dtype=np.int64), dims=("nface",))

        result = synthetic_model.erosion_rate(
            tau_pa=tau,
            layer_index=3,  # K_B
            layer_mass=layer_mass,
            layer_initial_mass=layer_initial_mass,
            bulk_density=bulk_density,
            core_id=core_id,
        )
        # Empty mass also zeroes the active mask, so result must be 0.
        np.testing.assert_allclose(result.values, [0.0])

    def test_zero_shear_zero_rate(self, synthetic_model):
        """τ ≤ 0 gates erosion off."""
        tau = _da([0.0, -0.1])
        layer_mass = _da([10.0, 10.0])
        layer_initial_mass = _da([10.0, 10.0])
        bulk_density = _da([1.6, 1.6])
        core_id = xr.DataArray(np.array([0, 0], dtype=np.int64), dims=("nface",))

        result = synthetic_model.erosion_rate(
            tau_pa=tau,
            layer_index=1,
            layer_mass=layer_mass,
            layer_initial_mass=layer_initial_mass,
            bulk_density=bulk_density,
            core_id=core_id,
        )
        np.testing.assert_allclose(result.values, [0.0, 0.0])

    def test_max_rate_cap_enforced(self, synthetic_model):
        """Cap clamps the rate to ``max_rate_per_core``."""
        # Override the cap so the spec τ=2 Pa case would otherwise exceed it.
        synthetic_model.max_rate_per_core[:] = 1.0e-4  # well below 6.4e-4

        tau = _da([2.0])
        layer_mass = _da([15.0])
        layer_initial_mass = _da([15.0])
        bulk_density = _da([1.6])
        core_id = xr.DataArray(np.array([0], dtype=np.int64), dims=("nface",))

        result = synthetic_model.erosion_rate(
            tau_pa=tau,
            layer_index=1,
            layer_mass=layer_mass,
            layer_initial_mass=layer_initial_mass,
            bulk_density=bulk_density,
            core_id=core_id,
        )
        np.testing.assert_allclose(result.values, [1.0e-4])

    def test_vectorized_across_cells(self, synthetic_model):
        """Multiple cells with different τ produce the analytic per-cell rates."""
        tau = _da([1.0, 2.0, 3.0])
        layer_mass = _da([10.0, 10.0, 10.0])
        layer_initial_mass = _da([10.0, 10.0, 10.0])
        bulk_density = _da([1.6, 1.6, 1.6])
        core_id = xr.DataArray(np.array([0, 0, 0], dtype=np.int64), dims=("nface",))

        result = synthetic_model.erosion_rate(
            tau_pa=tau,
            layer_index=1,
            layer_mass=layer_mass,
            layer_initial_mass=layer_initial_mass,
            bulk_density=bulk_density,
            core_id=core_id,
        )
        expected = 1.0e-4 * (np.array([1.0, 2.0, 3.0]) ** 2.0) * 1.6
        np.testing.assert_allclose(result.values, expected, rtol=1e-12)

    def test_partial_layer_blends_top_and_bottom(self, synthetic_model):
        """Half-eroded layer blends E_top and E_bottom equally.

        Both layers have the same (A, n) here, so the blend is degenerate
        and equals E_top regardless. Set distinct A_K+1 to exercise the blend.
        """
        # K=1 -> A=1e-4; K=2 -> A=2e-4. Both n=2.0.
        synthetic_model.ea_per_core[0, 1] = 2.0e-4

        tau = _da([2.0])
        # m_K = m_K0/2 → SN11 = 0.5 → E = 0.5*(E_bot + E_top)
        layer_mass = _da([5.0])
        layer_initial_mass = _da([10.0])
        bulk_density = _da([1.0])  # easier to read
        core_id = xr.DataArray(np.array([0], dtype=np.int64), dims=("nface",))

        result = synthetic_model.erosion_rate(
            tau_pa=tau,
            layer_index=1,  # not deepest → uses K=2 below
            layer_mass=layer_mass,
            layer_initial_mass=layer_initial_mass,
            bulk_density=bulk_density,
            core_id=core_id,
        )

        e_top = 1.0e-4 * 4.0   # 4e-4
        e_bot = 2.0e-4 * 4.0   # 8e-4
        expected = 0.5 * (e_bot - e_top) + e_top  # = 0.5*(e_top + e_bot) = 6e-4
        np.testing.assert_allclose(result.values, [expected], rtol=1e-12)


# ---------------------------------------------------------------------------
# SedflumeTableErosionModel
# ---------------------------------------------------------------------------

class TestSedflumeTableErosionModel:
    """Tests for ``nsedflume=1`` bilinear (log-depth, linear-τ) interpolation."""

    @pytest.fixture
    def synthetic_model(self) -> SedflumeTableErosionModel:
        """One core, K_B = 3 layers, ITBM = 2 shear levels.

        ERATE values (cm/s) are chosen so the log-depth blend at the
        midpoint of layer fractional remaining mass produces a clean
        geometric mean that's easy to verify analytically.

        Table layout: erate_per_core[core, layer-1, tau_level].
        """
        tau_levels_pa = np.array([0.5, 2.0])  # ITBM = 2

        # Per layer K, per τ level: pick values so log-blends are easy.
        # K=1 (top): [low_τ=1e-4, high_τ=1e-2]
        # K=2:        [low_τ=1e-3, high_τ=1e-1]
        # K=3 (deepest): [low_τ=1e-5, high_τ=1e-3]
        erate_per_core = np.array([[
            [1.0e-4, 1.0e-2],   # K=1
            [1.0e-3, 1.0e-1],   # K=2
            [1.0e-5, 1.0e-3],   # K=3
        ]])  # shape (1, 3, 2)

        size_interpolants_um = np.array([10.0, 100.0])
        erate_active_per_size = np.array([
            [1.0e-4, 1.0e-2],
            [1.0e-3, 1.0e-1],
        ])
        taucrit_per_size_pa = np.array([0.1, 0.5])

        return SedflumeTableErosionModel(
            tau_levels_pa=tau_levels_pa,
            erate_per_core=erate_per_core,
            erate_active_per_size=erate_active_per_size,
            size_interpolants_um=size_interpolants_um,
            taucrit_per_size_pa=taucrit_per_size_pa,
        )

    def test_full_fresh_layer_at_tau_low_endpoint(self, synthetic_model):
        """At τ = tau_low, full layer (m_K = m_K0): result = ERATE[K, low] · ρ_b.

        SN00 = 1, SN10 = 0, SN01 = 1, SN11 = 0:
            E = 1 · exp(0·ln(E_K+1) + 1·ln(E_K)) · ρ_b = E_K(low) · ρ_b
        """
        tau = _da([0.5])  # = tau_low
        layer_mass = _da([10.0])
        layer_initial_mass = _da([10.0])
        bulk_density = _da([1.5])
        core_id = xr.DataArray(np.array([0], dtype=np.int64), dims=("nface",))

        result = synthetic_model.erosion_rate(
            tau_pa=tau,
            layer_index=1,  # K=1 → ERATE[0,0,0] = 1e-4
            layer_mass=layer_mass,
            layer_initial_mass=layer_initial_mass,
            bulk_density=bulk_density,
            core_id=core_id,
        )
        expected = 1.0e-4 * 1.5
        np.testing.assert_allclose(result.values, [expected], rtol=1e-12)

    def test_full_fresh_layer_at_tau_high_endpoint(self, synthetic_model):
        """At τ = tau_high, full layer: E = ERATE[K, high] · ρ_b."""
        tau = _da([2.0])
        layer_mass = _da([10.0])
        layer_initial_mass = _da([10.0])
        bulk_density = _da([1.0])
        core_id = xr.DataArray(np.array([0], dtype=np.int64), dims=("nface",))

        result = synthetic_model.erosion_rate(
            tau_pa=tau,
            layer_index=1,
            layer_mass=layer_mass,
            layer_initial_mass=layer_initial_mass,
            bulk_density=bulk_density,
            core_id=core_id,
        )
        expected = 1.0e-2 * 1.0  # ERATE[K=1, high] · ρ_b
        np.testing.assert_allclose(result.values, [expected], rtol=1e-12)

    def test_tau_midpoint_full_layer(self, synthetic_model):
        """At τ midpoint (1.25 Pa), full layer:

        SN00 = (2.0 - 1.25)/(2.0 - 0.5) = 0.5, SN10 = 0.5
        SN01 = 1, SN11 = 0
        E = 0.5·E_K(low)·ρ_b + 0.5·E_K(high)·ρ_b
        """
        tau = _da([1.25])
        layer_mass = _da([10.0])
        layer_initial_mass = _da([10.0])
        bulk_density = _da([1.0])
        core_id = xr.DataArray(np.array([0], dtype=np.int64), dims=("nface",))

        result = synthetic_model.erosion_rate(
            tau_pa=tau,
            layer_index=1,
            layer_mass=layer_mass,
            layer_initial_mass=layer_initial_mass,
            bulk_density=bulk_density,
            core_id=core_id,
        )
        expected = 0.5 * 1.0e-4 + 0.5 * 1.0e-2
        np.testing.assert_allclose(result.values, [expected], rtol=1e-12)

    def test_depth_midpoint_at_tau_low(self, synthetic_model):
        """At τ = tau_low, half-eroded layer:

        SN00 = 1, SN10 = 0, SN01 = 0.5, SN11 = 0.5
        E = 1·exp(0.5·ln(E_K+1,low) + 0.5·ln(E_K,low))·ρ_b
          = sqrt(E_K+1,low · E_K,low) · ρ_b   (geometric mean)

        For K=1: sqrt(1e-3 · 1e-4) = sqrt(1e-7) ≈ 3.162277e-4 cm/s.
        """
        tau = _da([0.5])
        layer_mass = _da([5.0])
        layer_initial_mass = _da([10.0])
        bulk_density = _da([1.0])
        core_id = xr.DataArray(np.array([0], dtype=np.int64), dims=("nface",))

        result = synthetic_model.erosion_rate(
            tau_pa=tau,
            layer_index=1,
            layer_mass=layer_mass,
            layer_initial_mass=layer_initial_mass,
            bulk_density=bulk_density,
            core_id=core_id,
        )
        expected = np.sqrt(1.0e-3 * 1.0e-4)  # geometric mean of E_K and E_K+1
        np.testing.assert_allclose(result.values, [expected], rtol=1e-12)

    def test_bilinear_midpoint(self, synthetic_model):
        """Both midpoints: τ midpoint AND depth midpoint, K=1.

        SN00 = SN10 = 0.5; SN01 = SN11 = 0.5.
        E = 0.5·sqrt(E_K+1,low · E_K,low) + 0.5·sqrt(E_K+1,high · E_K,high)
          = 0.5·sqrt(1e-3·1e-4) + 0.5·sqrt(1e-1·1e-2)
        """
        tau = _da([1.25])
        layer_mass = _da([5.0])
        layer_initial_mass = _da([10.0])
        bulk_density = _da([1.0])
        core_id = xr.DataArray(np.array([0], dtype=np.int64), dims=("nface",))

        result = synthetic_model.erosion_rate(
            tau_pa=tau,
            layer_index=1,
            layer_mass=layer_mass,
            layer_initial_mass=layer_initial_mass,
            bulk_density=bulk_density,
            core_id=core_id,
        )
        expected = 0.5 * np.sqrt(1.0e-3 * 1.0e-4) + 0.5 * np.sqrt(1.0e-1 * 1.0e-2)
        np.testing.assert_allclose(result.values, [expected], rtol=1e-12)

    def test_deepest_layer_uses_floor_for_kp1(self, synthetic_model):
        """Deepest layer K=K_B replaces E_K+1 with the 1e-9 floor.

        Half-eroded → E ≈ sqrt(floor · E_K,low) · ρ_b.
        For K=3 at τ=tau_low: sqrt(1e-9 · 1e-5) = sqrt(1e-14) = 1e-7 cm/s.
        """
        tau = _da([0.5])
        layer_mass = _da([5.0])
        layer_initial_mass = _da([10.0])
        bulk_density = _da([1.0])
        core_id = xr.DataArray(np.array([0], dtype=np.int64), dims=("nface",))

        result = synthetic_model.erosion_rate(
            tau_pa=tau,
            layer_index=3,  # K_B
            layer_mass=layer_mass,
            layer_initial_mass=layer_initial_mass,
            bulk_density=bulk_density,
            core_id=core_id,
        )
        expected = np.sqrt(1.0e-9 * 1.0e-5)  # = 1e-7
        np.testing.assert_allclose(result.values, [expected], rtol=1e-10)

    def test_out_of_range_tau_clamps_with_warning(self, synthetic_model):
        """τ above the highest table level clamps and emits a RuntimeWarning."""
        tau = _da([10.0])  # well above tau_high = 2.0
        layer_mass = _da([10.0])
        layer_initial_mass = _da([10.0])
        bulk_density = _da([1.0])
        core_id = xr.DataArray(np.array([0], dtype=np.int64), dims=("nface",))

        with pytest.warns(RuntimeWarning):
            result = synthetic_model.erosion_rate(
                tau_pa=tau,
                layer_index=1,
                layer_mass=layer_mass,
                layer_initial_mass=layer_initial_mass,
                bulk_density=bulk_density,
                core_id=core_id,
            )
        # Should equal the rate at tau_high.
        np.testing.assert_allclose(result.values, [1.0e-2], rtol=1e-12)

    def test_vectorized_across_cells(self, synthetic_model):
        """Three cells: tau_low endpoint, midpoint, tau_high endpoint."""
        tau = _da([0.5, 1.25, 2.0])
        layer_mass = _da([10.0, 10.0, 10.0])
        layer_initial_mass = _da([10.0, 10.0, 10.0])
        bulk_density = _da([1.0, 1.0, 1.0])
        core_id = xr.DataArray(np.zeros(3, dtype=np.int64), dims=("nface",))

        result = synthetic_model.erosion_rate(
            tau_pa=tau,
            layer_index=1,
            layer_mass=layer_mass,
            layer_initial_mass=layer_initial_mass,
            bulk_density=bulk_density,
            core_id=core_id,
        )
        expected = np.array([
            1.0e-4,
            0.5 * 1.0e-4 + 0.5 * 1.0e-2,
            1.0e-2,
        ])
        np.testing.assert_allclose(result.values, expected, rtol=1e-12)

    def test_validation_rejects_nonincreasing_tau(self):
        with pytest.raises(ValueError, match="strictly increasing"):
            SedflumeTableErosionModel(
                tau_levels_pa=np.array([1.0, 1.0]),
                erate_per_core=np.ones((1, 1, 2)),
                erate_active_per_size=np.ones((1, 2)),
                size_interpolants_um=np.array([10.0]),
                taucrit_per_size_pa=np.array([0.1]),
            )

    def test_validation_rejects_shape_mismatch(self):
        with pytest.raises(ValueError, match="last dim"):
            SedflumeTableErosionModel(
                tau_levels_pa=np.array([0.5, 2.0]),
                erate_per_core=np.ones((1, 2, 3)),  # last dim ≠ 2
                erate_active_per_size=np.ones((1, 2)),
                size_interpolants_um=np.array([10.0]),
                taucrit_per_size_pa=np.array([0.1]),
            )


# ---------------------------------------------------------------------------
# apply_vegetation_cohesion
# ---------------------------------------------------------------------------

class TestApplyVegetationCohesion:
    """Tests for the vegetation-cohesion feedback on τ_ce."""

    def test_no_inputs_returns_unchanged(self):
        tau_ce = xr.DataArray(
            np.array([[0.2, 0.4], [0.3, 0.5]]),
            dims=("nface", "ssm_class"),
        )
        out = apply_vegetation_cohesion(tau_ce, None, None)
        # Should be the same object (identity), not a copy.
        assert out is tau_ce

    def test_biostabilization_only_spec_example(self):
        """τ_ce_base = 0.2 Pa, B = 1.0, α = 0.5 → τ_ce_eff = 0.3 Pa."""
        tau_ce = xr.DataArray(np.array([[0.2]]), dims=("nface", "ssm_class"))
        biostab = xr.DataArray(np.array([1.0]), dims=("nface",))

        out = apply_vegetation_cohesion(
            tau_ce_pa=tau_ce,
            biostabilization=biostab,
            biostabilization_alpha=0.5,
        )
        np.testing.assert_allclose(out.values, [[0.3]], rtol=1e-12)

    def test_biostabilization_and_root_cohesion_spec_example(self):
        """τ_ce_base = 0.2, B = 1.0, α = 0.5, root = 0.1 → τ_ce_eff = 0.4 Pa."""
        tau_ce = xr.DataArray(np.array([[0.2]]), dims=("nface", "ssm_class"))
        biostab = xr.DataArray(np.array([1.0]), dims=("nface",))
        root = xr.DataArray(np.array([0.1]), dims=("nface",))

        out = apply_vegetation_cohesion(
            tau_ce_pa=tau_ce,
            biostabilization=biostab,
            root_cohesion_pa=root,
            biostabilization_alpha=0.5,
        )
        np.testing.assert_allclose(out.values, [[0.4]], rtol=1e-12)

    def test_root_cohesion_only(self):
        tau_ce = xr.DataArray(np.array([[0.2, 0.5]]), dims=("nface", "ssm_class"))
        root = xr.DataArray(np.array([0.1]), dims=("nface",))

        out = apply_vegetation_cohesion(
            tau_ce_pa=tau_ce,
            root_cohesion_pa=root,
        )
        np.testing.assert_allclose(out.values, [[0.3, 0.6]], rtol=1e-12)

    def test_broadcasts_across_class_dim(self):
        """Per-cell biostab and root_cohesion broadcast across (nface, ssm_class)."""
        tau_ce = xr.DataArray(
            np.array([[0.2, 0.4],   # cell 0, two classes
                      [0.3, 0.6]]),  # cell 1, two classes
            dims=("nface", "ssm_class"),
        )
        biostab = xr.DataArray(np.array([1.0, 0.0]), dims=("nface",))  # cell-0 vegetated, cell-1 bare
        root = xr.DataArray(np.array([0.1, 0.0]), dims=("nface",))

        out = apply_vegetation_cohesion(
            tau_ce_pa=tau_ce,
            biostabilization=biostab,
            root_cohesion_pa=root,
            biostabilization_alpha=0.5,
        )
        # Cell 0: τ_ce × (1 + 0.5·1.0) + 0.1 = 1.5·τ_ce + 0.1
        #   → [0.3 + 0.1, 0.6 + 0.1] = [0.4, 0.7]
        # Cell 1: τ_ce × (1 + 0) + 0 = τ_ce → [0.3, 0.6]
        expected = np.array([[0.4, 0.7], [0.3, 0.6]])
        np.testing.assert_allclose(out.values, expected, rtol=1e-12)

    def test_custom_alpha(self):
        tau_ce = xr.DataArray(np.array([[1.0]]), dims=("nface", "ssm_class"))
        biostab = xr.DataArray(np.array([0.5]), dims=("nface",))

        out = apply_vegetation_cohesion(
            tau_ce_pa=tau_ce,
            biostabilization=biostab,
            biostabilization_alpha=2.0,
        )
        # 1.0 × (1 + 2.0 × 0.5) = 1.0 × 2.0 = 2.0
        np.testing.assert_allclose(out.values, [[2.0]], rtol=1e-12)
