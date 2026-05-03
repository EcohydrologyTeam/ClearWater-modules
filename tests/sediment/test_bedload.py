"""Unit tests for SSM bedload transport (van Rijn 1984).

Covers:

* Closed-form van Rijn helpers (velocity, height, equilibrium concentration).
* The standalone explicit-upwind solver on a synthetic 3-cell mesh, with
  a closed-domain mass-conservation check.
* The Riverine-constituent solver's advection-coefficient field write.
* A scaffolded parity test deferred until Stage 2 wire-up makes both
  solvers runnable on the same 1-D channel.

The tolerances follow the issue spec: ``5%`` for the closed-form helpers
(targets given to two figures), tighter for mass conservation (which is
algebraic up to round-off in the explicit upwind step on a closed domain).
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest
import xarray as xr

from clearwater_modules_v2.processes.sediment import bedload, contracts
from clearwater_modules_v2.processes.sediment.classes import (
    SedimentClass,
    SedimentClassRegistry,
)


# ---------------------------------------------------------------------------
# Reference values for the van Rijn closed-form helpers
# ---------------------------------------------------------------------------
#
# All three closed-form tests share the same physical setup:
#
#   D50      = 250 μm  (sand_medium)
#   τ_ce     = 0.2 Pa
#   τ        = 0.4 Pa  ⇒  T_R = 1.0
#   s_s      = 2.65 (quartz)
#   ν        = 0.01 cm²/s, g = 980 cm/s²
#
# Hand-derived expected values (see issue spec):
#   d_*      = 0.025 cm × ((1.65·980)/(0.01²))^(1/3) ≈ 6.32
#   u_BL     = 1.5 × 1^0.6 × √(1.65·980·0.025)       ≈ 9.54 cm/s
#   δ_BL    = 0.3 × 0.025 × 6.32^0.7 × √1            ≈ 0.0273 cm
#   C_eq     = 0.117 × 2.65 × 1.0 / 6.32             ≈ 0.0490 g/cm³
# ---------------------------------------------------------------------------

D50_UM = 250.0
TAU_CE = 0.2
TAU = 0.4
EXPECTED_TR = 1.0
EXPECTED_D_STAR = 6.32
EXPECTED_UBL_CM_S = 9.54
EXPECTED_DELTA_BL_CM = 0.0273
EXPECTED_C_EQ_G_CM3 = 0.0490


def _tau_array(value: float, nface: int = 3) -> xr.DataArray:
    """Helper: build a constant-tau xarray DataArray of length ``nface``."""
    return xr.DataArray(
        np.full(nface, value, dtype="float64"),
        dims=("nface",),
    )


# ---------------------------------------------------------------------------
# Closed-form van Rijn helpers
# ---------------------------------------------------------------------------

class TestVanRijnVelocity:
    def test_matches_hand_derivation_within_5pct(self):
        tau = _tau_array(TAU)
        u_bl = bedload.van_rijn_bedload_velocity_cm_s(
            tau, tau_ce_pa=TAU_CE, d50_um=D50_UM
        )
        # All cells equal under uniform forcing.
        assert u_bl.dims == ("nface",)
        np.testing.assert_allclose(u_bl.values, EXPECTED_UBL_CM_S, rtol=0.05)

    def test_subcritical_returns_zero(self):
        # τ < τ_ce ⇒ T_R clipped to 0 ⇒ u_BL = 0.
        tau = _tau_array(0.1)  # below τ_ce = 0.2
        u_bl = bedload.van_rijn_bedload_velocity_cm_s(
            tau, tau_ce_pa=TAU_CE, d50_um=D50_UM
        )
        np.testing.assert_array_equal(u_bl.values, np.zeros(3))

    def test_invalid_inputs_raise(self):
        tau = _tau_array(TAU)
        with pytest.raises(ValueError):
            bedload.van_rijn_bedload_velocity_cm_s(
                tau, tau_ce_pa=0.0, d50_um=D50_UM
            )
        with pytest.raises(ValueError):
            bedload.van_rijn_bedload_velocity_cm_s(
                tau, tau_ce_pa=TAU_CE, d50_um=-1.0
            )


class TestVanRijnHeight:
    def test_matches_hand_derivation_within_5pct(self):
        tau = _tau_array(TAU)
        delta_bl = bedload.van_rijn_bedload_height_cm(
            tau, tau_ce_pa=TAU_CE, d50_um=D50_UM, d_star=EXPECTED_D_STAR
        )
        np.testing.assert_allclose(delta_bl.values, EXPECTED_DELTA_BL_CM, rtol=0.05)

    def test_subcritical_returns_zero(self):
        tau = _tau_array(0.1)
        delta_bl = bedload.van_rijn_bedload_height_cm(
            tau, tau_ce_pa=TAU_CE, d50_um=D50_UM, d_star=EXPECTED_D_STAR
        )
        np.testing.assert_array_equal(delta_bl.values, np.zeros(3))


class TestVanRijnEquilibriumConcentration:
    def test_matches_hand_derivation_within_5pct(self):
        tau = _tau_array(TAU)
        c_eq = bedload.van_rijn_equilibrium_concentration(
            tau, tau_ce_pa=TAU_CE, d_star=EXPECTED_D_STAR
        )
        np.testing.assert_allclose(c_eq.values, EXPECTED_C_EQ_G_CM3, rtol=0.05)

    def test_invalid_d_star_raises(self):
        tau = _tau_array(TAU)
        with pytest.raises(ValueError):
            bedload.van_rijn_equilibrium_concentration(
                tau, tau_ce_pa=TAU_CE, d_star=0.0
            )


# ---------------------------------------------------------------------------
# Synthetic 3-cell mesh fixtures
# ---------------------------------------------------------------------------
#
# Topology: a 1-D chain of 3 cells joined by 2 internal edges. No boundary
# edges (closed domain) — this is the key property that lets us check
# bedload mass conservation exactly.
#
#     cell 0 ── edge 0 ── cell 1 ── edge 1 ── cell 2
#
# We deliberately *omit* boundary edges so every face exchange is internal
# and the per-step mass flux sums to zero across all cells.
# ---------------------------------------------------------------------------

def _make_closed_3cell_mesh() -> xr.Dataset:
    nface = 3
    nedge = 2
    edges_face1 = np.array([0, 1], dtype=np.int64)
    edges_face2 = np.array([1, 2], dtype=np.int64)
    edge_length_m = np.array([10.0, 10.0], dtype="float64")
    face_area_m2 = np.array([100.0, 100.0, 100.0], dtype="float64")

    return xr.Dataset(
        data_vars={
            "edges_face1": (("nedge",), edges_face1),
            "edges_face2": (("nedge",), edges_face2),
            "edge_length": (("nedge",), edge_length_m),
            "faces_surface_area": (("nface",), face_area_m2),
        },
        coords={
            "nface": np.arange(nface),
            "nedge": np.arange(nedge),
        },
    )


def _make_registry() -> SedimentClassRegistry:
    """One bedload-eligible (sand_medium) and one cohesive (silt) class.

    Forces the standalone solver to filter the silt class out by D50.
    """
    classes = [
        SedimentClass(label="silt_fine", d50_um=20.0, tau_ce_pa=0.05),
        SedimentClass(label="sand_medium", d50_um=D50_UM, tau_ce_pa=TAU_CE),
    ]
    return SedimentClassRegistry.from_iterable(classes)


# ---------------------------------------------------------------------------
# Standalone explicit-upwind solver
# ---------------------------------------------------------------------------

class TestBedloadStandaloneExplicit:
    def test_filters_to_bedload_eligible_classes(self):
        registry = _make_registry()
        solver = bedload.BedloadStandaloneExplicit(registry)
        # silt_fine (20 μm) is below the 64 μm cutoff; sand_medium (250 μm) above.
        assert solver.eligible_class_indices == [1]
        assert pytest.approx(solver._d_star[1], rel=0.05) == EXPECTED_D_STAR

    def test_one_step_conserves_mass_in_closed_domain(self):
        registry = _make_registry()
        solver = bedload.BedloadStandaloneExplicit(registry)
        mesh = _make_closed_3cell_mesh()

        # Seed CBL: only sand_medium (class index 1) carries bedload mass.
        # Initial concentration: cell 0 has 1.0 g/cm², others zero.
        # Under positive u_BL this should propagate downstream while
        # leaving the closed-domain total invariant.
        nface = mesh.sizes["nface"]
        n_classes = len(registry)
        cbl0 = np.zeros((nface, n_classes), dtype="float32")
        cbl0[0, 1] = 1.0
        mesh[contracts.VAR_BEDLOAD_MASS] = (
            (contracts.DIM_NFACE, contracts.DIM_CLASS),
            cbl0.copy(),
        )

        # Uniform shear well above τ_ce so u_BL > 0 in every cell.
        tau_pa = _tau_array(TAU, nface=nface)

        total_before = float(
            mesh[contracts.VAR_BEDLOAD_MASS].sum().values
        )

        solver.step(
            mesh=mesh,
            time=None,                  # no time dim on this synthetic mesh
            tau_pa=tau_pa,
            dt_seconds=1.0,
        )

        total_after = float(mesh[contracts.VAR_BEDLOAD_MASS].sum().values)
        # Closed-domain explicit-upwind step is mass-conservative up to
        # float32 round-off (the bedload array is stored as float32).
        np.testing.assert_allclose(total_after, total_before, rtol=1e-5, atol=1e-7)

        # Sanity: mass moved from cell 0 toward cells 1 and 2 (downstream).
        cbl_after = mesh[contracts.VAR_BEDLOAD_MASS].values[:, 1]
        assert cbl_after[0] < cbl0[0, 1]   # cell 0 lost mass
        assert cbl_after[1] > 0.0          # cell 1 gained mass
        # cell 2 may still be at 0 after one step (no upstream mass yet)
        # because the upwind from edge 1 is cell 1 which started at 0.
        assert cbl_after[2] == pytest.approx(0.0)

    def test_zero_dt_is_a_noop(self):
        registry = _make_registry()
        solver = bedload.BedloadStandaloneExplicit(registry)
        mesh = _make_closed_3cell_mesh()
        nface = mesh.sizes["nface"]
        n_classes = len(registry)

        cbl0 = np.zeros((nface, n_classes), dtype="float32")
        cbl0[0, 1] = 1.0
        mesh[contracts.VAR_BEDLOAD_MASS] = (
            (contracts.DIM_NFACE, contracts.DIM_CLASS),
            cbl0.copy(),
        )
        tau_pa = _tau_array(TAU, nface=nface)

        solver.step(mesh=mesh, time=None, tau_pa=tau_pa, dt_seconds=0.0)
        np.testing.assert_array_equal(
            mesh[contracts.VAR_BEDLOAD_MASS].values, cbl0
        )


# ---------------------------------------------------------------------------
# Riverine-constituent solver
# ---------------------------------------------------------------------------

class _RiverineStub:
    """Minimal stand-in for the Riverine instance used by Bedload Mode B.

    Exposes the ``constituent_dict`` attribute the solver writes to.
    """

    def __init__(self):
        self.constituent_dict: dict = {}


class TestBedloadRiverineConstituent:
    def test_registers_bedload_companion_constituents(self):
        registry = _make_registry()
        riverine = _RiverineStub()
        solver = bedload.BedloadRiverineConstituent(registry, riverine)

        # Only sand_medium is bedload-eligible.
        expected_name = (
            contracts.suspended_var_name("sand_medium") + "_bedload"
        )
        assert expected_name in riverine.constituent_dict
        cfg = riverine.constituent_dict[expected_name]
        assert cfg["decay_rate"] == 0.0
        assert (
            cfg["advection_coefficient_var"]
            == contracts.advection_coef_var_name("sand_medium")
        )
        # silt_fine is below cutoff and must NOT be registered.
        silt_name = contracts.suspended_var_name("silt_fine") + "_bedload"
        assert silt_name not in riverine.constituent_dict

    def test_step_writes_advection_coefficient_field(self):
        registry = _make_registry()
        riverine = _RiverineStub()
        solver = bedload.BedloadRiverineConstituent(registry, riverine)
        mesh = _make_closed_3cell_mesh()

        var_name = contracts.advection_coef_var_name("sand_medium")
        assert var_name not in mesh.variables  # not yet allocated

        tau_pa = _tau_array(TAU, nface=mesh.sizes["nface"])
        solver.step(mesh=mesh, time=None, tau_pa=tau_pa, dt_seconds=60.0)

        # Field is written, lives on the edge dim, and carries m/s units.
        assert var_name in mesh.variables
        field = mesh[var_name]
        assert field.dims == ("nedge",)
        assert field.shape == (mesh.sizes["nedge"],)
        assert field.attrs.get("units") == "m s-1"

        # Magnitude check: u_BL ≈ 9.54 cm/s ⇒ 0.0954 m/s on every edge
        # (uniform shear, uniform D50). Allow 5% to match the helper test.
        np.testing.assert_allclose(
            field.values, 0.01 * EXPECTED_UBL_CM_S, rtol=0.05
        )


# ---------------------------------------------------------------------------
# Parity test (placeholder)
# ---------------------------------------------------------------------------

@pytest.mark.skip(
    reason=(
        "requires Stage 2 wire-up — both solvers run a 1D channel and "
        "CBL agrees within tol"
    )
)
def test_standalone_vs_riverine_parity_on_1d_channel():
    """Run BedloadStandaloneExplicit and BedloadRiverineConstituent on
    the same 1-D channel and verify the per-class CBL fields agree to
    within numerical tolerance.

    Deferred: depends on Batch C, which extends Riverine's implicit
    advection solver to read per-constituent advection-coefficient
    fields. Once Batch C lands, populate this test with a sand-only
    1-D channel, run both solvers for N steps with the same ``dt`` and
    boundary conditions, and assert ``np.allclose(cbl_a, cbl_b, ...)``.
    """
    raise AssertionError("placeholder — should be skipped by the marker above")
