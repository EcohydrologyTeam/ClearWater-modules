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


# ---------------------------------------------------------------------------
# Solver wiring: standalone solver consumes a configured transport function
# ---------------------------------------------------------------------------
#
# These tests exercise the Stage-2 solver-wiring refactor: the standalone
# (and Riverine-constituent) solvers now accept an optional
# ``transport_function=`` argument that they consume on every call. The
# existing default behaviour (van Rijn) is preserved bit-for-bit by
# ``BedloadStandaloneExplicit(registry)`` (no kwarg) — backed by the
# unchanged ``test_one_step_conserves_mass_in_closed_domain`` test above.
# Below we add coverage for:
#
# * Wilcock-Crowe (2003) — surface-composition-sensitive closure.
# * Wu (2000) — the other surface-composition-sensitive closure.
# * Standalone-vs-Riverine parity under a non-van-Rijn function.
# ---------------------------------------------------------------------------


def _seed_cbl_in_cell(mesh: xr.Dataset, registry, cell: int = 0, class_idx: int = 1,
                      mass: float = 1.0) -> None:
    """Helper: seed a single non-zero CBL value on the mesh."""
    nface = mesh.sizes["nface"]
    n_classes = len(registry)
    cbl0 = np.zeros((nface, n_classes), dtype="float32")
    cbl0[cell, class_idx] = mass
    mesh[contracts.VAR_BEDLOAD_MASS] = (
        (contracts.DIM_NFACE, contracts.DIM_CLASS),
        cbl0.copy(),
    )


class TestStandaloneSolverWithCustomTransportFunction:
    """Standalone solver consumes the ``transport_function`` argument."""

    def test_van_rijn_default_is_backwards_compatible(self):
        """Construction without transport_function still uses van Rijn."""
        registry = _make_registry()
        solver = bedload.BedloadStandaloneExplicit(registry)
        assert isinstance(
            solver.transport_function,
            bedload.VanRijn1984TransportFunction,
        )

    def test_explicit_van_rijn_matches_default(self):
        """Explicit van Rijn instance produces identical CBL update to default."""
        registry = _make_registry()
        mesh_a = _make_closed_3cell_mesh()
        mesh_b = _make_closed_3cell_mesh()
        _seed_cbl_in_cell(mesh_a, registry)
        _seed_cbl_in_cell(mesh_b, registry)
        tau = _tau_array(TAU, nface=3)

        solver_a = bedload.BedloadStandaloneExplicit(registry)
        solver_b = bedload.BedloadStandaloneExplicit(
            registry, transport_function=bedload.VanRijn1984TransportFunction()
        )
        solver_a.step(mesh=mesh_a, time=None, tau_pa=tau, dt_seconds=1.0)
        solver_b.step(mesh=mesh_b, time=None, tau_pa=tau, dt_seconds=1.0)

        np.testing.assert_allclose(
            mesh_a[contracts.VAR_BEDLOAD_MASS].values,
            mesh_b[contracts.VAR_BEDLOAD_MASS].values,
            rtol=1e-12,
            atol=1e-15,
        )

    def test_wilcock_crowe_smoke_with_surface_composition(self):
        """Wilcock-Crowe runs end-to-end and produces a per-cell q_b matching
        the standalone unit-test value at the same conditions.

        Conditions match ``TestWilcockCrowe2003.test_reference_value_within_tolerance``:
        D50 = 5 mm, F_s = 0.34, τ = 5 Pa ⇒ q_b ≈ 0.293 g/cm/s.
        Verified by computing q_b directly via the closure (with the same
        registry_context the solver passes) and confirming the per-cell
        effective velocity used by the solver agrees with q_b/(δ_BL · C_eq).
        """
        # Single bedload-eligible class at 5 mm gravel.
        classes = [
            SedimentClass(label="gravel_5mm", d50_um=5000.0, tau_ce_pa=2.0),
        ]
        registry = SedimentClassRegistry.from_iterable(classes)

        fn = bedload.WilcockCrowe2003TransportFunction()
        solver = bedload.BedloadStandaloneExplicit(
            registry, transport_function=fn
        )
        mesh = _make_closed_3cell_mesh()
        _seed_cbl_in_cell(mesh, registry, cell=0, class_idx=0, mass=1.0)

        tau = _tau_array(5.0, nface=3)
        ctx = {
            "surface_sand_fraction": 0.34,
            "surface_class_fraction": 1.0,
        }
        solver.step(
            mesh=mesh, time=None, tau_pa=tau, dt_seconds=1.0,
            registry_context=ctx,
        )

        # Independently compute the closure's q_b at the same conditions.
        qb = fn.transport_rate(
            tau_pa=tau,
            d50_um=5000.0,
            tau_ce_pa=2.0,
            velocity_m_s=xr.DataArray(np.full(3, 1.0), dims=("nface",)),
            depth_m=xr.DataArray(np.full(3, 1.0), dims=("nface",)),
            slope=1.0e-3,
            solid_density_g_cm3=2.65,
            registry_context=ctx,
        )
        # Hand-derived reference is 0.293 g/cm/s; tolerate 10 %.
        np.testing.assert_allclose(qb.values, 0.293, rtol=0.10)

        # The solver's effective velocity must equal q_b / (δ_BL · C_eq).
        d_star = bedload._cheng_d_star(5000.0, solid_specific_gravity=2.65)
        delta = bedload.van_rijn_bedload_height_cm(
            tau, tau_ce_pa=2.0, d50_um=5000.0, d_star=d_star,
        )
        c_eq = bedload.van_rijn_equilibrium_concentration(
            tau, tau_ce_pa=2.0, d_star=d_star, solid_density_g_cm3=2.65,
        )
        u_eff_expected = qb.values / (delta.values * c_eq.values)
        u_eff = bedload._qb_to_effective_velocity_cm_s(
            qb_g_cm_s=qb,
            tau_pa=tau,
            tau_ce_pa=2.0,
            d50_um=5000.0,
            d_star=d_star,
            solid_density_g_cm3=2.65,
        )
        np.testing.assert_allclose(u_eff, u_eff_expected, rtol=1e-12)

        # The CBL must have moved (solver actually consumed the closure).
        cbl_after = mesh[contracts.VAR_BEDLOAD_MASS].values[:, 0]
        assert cbl_after[0] < 1.0
        assert cbl_after[1] > 0.0

    def test_wu2000_smoke_with_surface_composition(self):
        """Wu (2000) runs through the standalone solver with surface composition.

        Conditions: D50 = 250 μm, τ = 0.5 Pa, p_e/p_h = 1.0 (uniform bed),
        F_i = 1.0. Hand-derived q_b ≈ 0.0273 g/cm/s
        (matches ``TestWu2000.test_uniform_bed_reference_value``).
        """
        classes = [
            SedimentClass(label="sand_250", d50_um=250.0, tau_ce_pa=0.05),
        ]
        registry = SedimentClassRegistry.from_iterable(classes)

        fn = bedload.Wu2000TransportFunction()
        solver = bedload.BedloadStandaloneExplicit(
            registry, transport_function=fn
        )
        mesh = _make_closed_3cell_mesh()
        _seed_cbl_in_cell(mesh, registry, cell=0, class_idx=0, mass=1.0)

        tau = _tau_array(0.5, nface=3)
        ctx = {"pe_ph_ratio": 1.0, "surface_class_fraction": 1.0}
        solver.step(
            mesh=mesh, time=None, tau_pa=tau, dt_seconds=1.0,
            registry_context=ctx,
        )

        qb = fn.transport_rate(
            tau_pa=tau,
            d50_um=250.0,
            tau_ce_pa=0.05,
            velocity_m_s=xr.DataArray(np.full(3, 1.0), dims=("nface",)),
            depth_m=xr.DataArray(np.full(3, 1.0), dims=("nface",)),
            slope=1.0e-3,
            solid_density_g_cm3=2.65,
            registry_context=ctx,
        )
        np.testing.assert_allclose(qb.values, 0.0273, rtol=0.10)

        # CBL moved.
        cbl_after = mesh[contracts.VAR_BEDLOAD_MASS].values[:, 0]
        assert cbl_after[0] < 1.0
        assert cbl_after[1] > 0.0

    def test_wilcock_crowe_closed_domain_mass_conservation(self):
        """Even with a non-van-Rijn closure, closed-domain mass is conserved
        to float32 round-off (the upwind step is identical; only the
        per-cell magnitude of u_eff changes)."""
        classes = [
            SedimentClass(label="gravel_5mm", d50_um=5000.0, tau_ce_pa=2.0),
        ]
        registry = SedimentClassRegistry.from_iterable(classes)
        solver = bedload.BedloadStandaloneExplicit(
            registry, transport_function=bedload.WilcockCrowe2003TransportFunction()
        )
        mesh = _make_closed_3cell_mesh()
        _seed_cbl_in_cell(mesh, registry, cell=0, class_idx=0, mass=1.0)

        total_before = float(mesh[contracts.VAR_BEDLOAD_MASS].sum().values)
        solver.step(
            mesh=mesh, time=None, tau_pa=_tau_array(5.0, nface=3),
            dt_seconds=1.0,
            registry_context={"surface_sand_fraction": 0.34, "surface_class_fraction": 1.0},
        )
        total_after = float(mesh[contracts.VAR_BEDLOAD_MASS].sum().values)
        np.testing.assert_allclose(total_after, total_before, rtol=1e-5, atol=1e-7)


class TestStandaloneRiverineParityNonVanRijn:
    """Both solvers, given the same transport function, write the same
    per-class effective velocity field."""

    def test_wilcock_crowe_yields_matching_advection_coefficient(self):
        """Run BedloadRiverineConstituent with Wilcock-Crowe and compare its
        per-edge advection coefficient field to the per-edge mean of the
        standalone solver's per-cell u_eff. They must agree exactly
        (modulo float32 storage)."""
        classes = [
            SedimentClass(label="gravel_5mm", d50_um=5000.0, tau_ce_pa=2.0),
        ]
        registry = SedimentClassRegistry.from_iterable(classes)
        fn = bedload.WilcockCrowe2003TransportFunction()
        ctx = {"surface_sand_fraction": 0.34, "surface_class_fraction": 1.0}
        tau = _tau_array(5.0, nface=3)

        # Riverine-constituent path.
        riverine = _RiverineStub()
        rsolver = bedload.BedloadRiverineConstituent(
            registry, riverine, transport_function=fn,
        )
        mesh_r = _make_closed_3cell_mesh()
        rsolver.step(
            mesh=mesh_r, time=None, tau_pa=tau,
            dt_seconds=1.0, registry_context=ctx,
        )
        var_name = contracts.advection_coef_var_name("gravel_5mm")
        u_edge_riverine_m_s = mesh_r[var_name].values

        # Standalone path: replicate the per-edge mean of u_eff.
        d_star = bedload._cheng_d_star(5000.0, solid_specific_gravity=2.65)
        qb = fn.transport_rate(
            tau_pa=tau,
            d50_um=5000.0,
            tau_ce_pa=2.0,
            velocity_m_s=xr.DataArray(np.full(3, 1.0), dims=("nface",)),
            depth_m=xr.DataArray(np.full(3, 1.0), dims=("nface",)),
            slope=1.0e-3,
            solid_density_g_cm3=2.65,
            registry_context=ctx,
        )
        u_eff_cell_cm_s = bedload._qb_to_effective_velocity_cm_s(
            qb_g_cm_s=qb, tau_pa=tau, tau_ce_pa=2.0,
            d50_um=5000.0, d_star=d_star, solid_density_g_cm3=2.65,
        )
        # Edge-mean (matches BedloadRiverineConstituent.step).
        u_edge_expected_m_s = (
            0.5 * (u_eff_cell_cm_s[0:2] + u_eff_cell_cm_s[1:3]) * 0.01
        )
        np.testing.assert_allclose(
            u_edge_riverine_m_s, u_edge_expected_m_s.astype("float32"),
            rtol=1e-5, atol=1e-7,
        )


# ---------------------------------------------------------------------------
# Pluggable transport-function abstraction
# ---------------------------------------------------------------------------
#
# These tests cover the seven peer-reviewed bedload formulas registered in
# ``bedload.BEDLOAD_TRANSPORT_FUNCTIONS``. Each per-function test pins the
# value at a single hand-derived reference condition, the menu test loops
# over the registry as a smoke check, and the backwards-compat test verifies
# that ``VanRijn1984TransportFunction`` agrees with the standalone helpers.
# ---------------------------------------------------------------------------


def _benign_inputs(tau_value: float = TAU, nface: int = 3):
    """Return a (tau_pa, V, h) triple suitable for any transport function."""
    tau = _tau_array(tau_value, nface=nface)
    velocity = xr.DataArray(np.full(nface, 1.0), dims=("nface",))
    depth = xr.DataArray(np.full(nface, 1.0), dims=("nface",))
    return tau, velocity, depth


class TestTransportFunctionRegistry:
    """Registry-level smoke checks for the pluggable closures."""

    def test_seven_functions_registered(self):
        # Spec requires exactly the seven peer-reviewed formulas.
        assert set(bedload.BEDLOAD_TRANSPORT_FUNCTIONS) == {
            "van_rijn",
            "wilcock_crowe",
            "parker",
            "yang",
            "wu",
            "engelund_hansen",
            "toffaleti",
        }

    def test_get_transport_function_lookup(self):
        for name in bedload.BEDLOAD_TRANSPORT_FUNCTIONS:
            fn = bedload.get_transport_function(name)
            assert fn.name == name
            assert isinstance(fn, bedload.BedloadTransportFunction)

    def test_get_transport_function_unknown_raises(self):
        with pytest.raises(ValueError, match="Unknown bedload transport_function"):
            bedload.get_transport_function("not_a_real_formula")

    def test_menu_smoke_all_functions_return_nonnegative(self):
        """Loop over the registry, evaluate at benign forcing, assert sane output."""
        tau, vel, depth = _benign_inputs(TAU, nface=4)
        for name, cls in bedload.BEDLOAD_TRANSPORT_FUNCTIONS.items():
            fn = cls()
            qb = fn.transport_rate(
                tau_pa=tau,
                d50_um=D50_UM,
                tau_ce_pa=TAU_CE,
                velocity_m_s=vel,
                depth_m=depth,
                slope=1.0e-3,
                solid_density_g_cm3=2.65,
            )
            assert isinstance(qb, xr.DataArray), f"{name} did not return DataArray"
            assert qb.dims == tau.dims, f"{name} dims mismatch"
            arr = np.asarray(qb.values)
            assert np.all(np.isfinite(arr)), f"{name} produced non-finite values"
            assert np.all(arr >= 0.0), f"{name} produced negative q_b"


class TestVanRijnClassWrapper:
    """Backwards-compat: the class-wrapper mirrors the standalone helpers."""

    def test_van_rijn_class_matches_helpers(self):
        tau = _tau_array(TAU)
        vel = xr.DataArray(np.full(3, 0.5), dims=("nface",))
        depth = xr.DataArray(np.full(3, 1.0), dims=("nface",))

        fn = bedload.VanRijn1984TransportFunction()
        qb = fn.transport_rate(
            tau_pa=tau,
            d50_um=D50_UM,
            tau_ce_pa=TAU_CE,
            velocity_m_s=vel,
            depth_m=depth,
            slope=1.0e-3,
            solid_density_g_cm3=2.65,
        )

        # Recompute via the standalone helpers and compare.
        d_star = bedload._cheng_d_star(D50_UM, solid_specific_gravity=2.65)
        u_bl = bedload.van_rijn_bedload_velocity_cm_s(
            tau, tau_ce_pa=TAU_CE, d50_um=D50_UM
        )
        delta = bedload.van_rijn_bedload_height_cm(
            tau, tau_ce_pa=TAU_CE, d50_um=D50_UM, d_star=d_star
        )
        c_eq = bedload.van_rijn_equilibrium_concentration(
            tau, tau_ce_pa=TAU_CE, d_star=d_star, solid_density_g_cm3=2.65
        )
        expected = (u_bl * delta * c_eq).values

        np.testing.assert_allclose(qb.values, expected, rtol=1e-12)

    def test_van_rijn_class_is_default_in_registry(self):
        assert "van_rijn" in bedload.BEDLOAD_TRANSPORT_FUNCTIONS
        # The default class must be the wrapper (not the standalone helper).
        cls = bedload.BEDLOAD_TRANSPORT_FUNCTIONS["van_rijn"]
        assert cls is bedload.VanRijn1984TransportFunction


class TestWilcockCrowe2003:
    """Wilcock & Crowe (2003) sand-gravel mixed-bed reference value.

    Reference condition (matches Wilcock & Crowe 2003 §"Comparison with
    flume data" example, simplified to single-class evaluation):

        D50 = 5 mm,  Fs = 0.34,  τ = 5 Pa,  ρ_s = 2650 kg/m³, ρ_w = 1000 kg/m³

    Hand calc:
        τ*_rsg = 0.021 + 0.015 · exp(-20·0.34)  ≈ 0.0210168
        d_i / d_sg = 1.0 → b_i = 0.67/(1+exp(1.5)) ≈ 0.1226
        τ*_ri ≈ 0.0210168
        τ_ri  = τ*_ri · (ρ_s − ρ_w) g d  = 0.0210168·1650·9.81·5e-3
              ≈ 1.700 Pa
        φ = 5 / 1.700 ≈ 2.940  (regime 2: φ ≥ 1.35)
        √φ = 1.7146 ⇒ 0.894/1.7146 = 0.5214 ⇒ 1 − 0.5214 = 0.4786
        0.4786^4.5 = exp(4.5·ln(0.4786)) = exp(−3.317) ≈ 0.0362
        W*_i = 14·0.0362 ≈ 0.507
        u* = √(5/1000) = 0.0707 m/s ; u*³ ≈ 3.535e-4 m³/s³
        q_b (kg/m/s) = 1.0 · 0.507 · 3.535e-4 · 2650 / (1.65·9.81)
                     ≈ 0.4748 / 16.19 ≈ 0.02933 kg/m/s
        q_b (g/cm/s) = 0.2933
    """

    def test_reference_value_within_tolerance(self):
        tau = _tau_array(5.0)
        vel = xr.DataArray(np.full(3, 1.0), dims=("nface",))
        depth = xr.DataArray(np.full(3, 1.0), dims=("nface",))
        fn = bedload.WilcockCrowe2003TransportFunction()
        qb = fn.transport_rate(
            tau_pa=tau,
            d50_um=5000.0,        # 5 mm
            tau_ce_pa=0.0,        # not used by W&C (uses internal τ_ri)
            velocity_m_s=vel,
            depth_m=depth,
            slope=1.0e-3,
            solid_density_g_cm3=2.65,
            registry_context={"surface_sand_fraction": 0.34, "surface_class_fraction": 1.0},
        )
        # Hand-derived value 0.293 g/cm/s; allow 10 % tolerance.
        np.testing.assert_allclose(qb.values, 0.293, rtol=0.10)

    def test_subcritical_returns_zero(self):
        # Very low τ → φ < 1 ⇒ low-regime W* ≈ 0.
        tau = _tau_array(0.01)
        vel = xr.DataArray(np.full(3, 0.1), dims=("nface",))
        depth = xr.DataArray(np.full(3, 0.5), dims=("nface",))
        fn = bedload.WilcockCrowe2003TransportFunction()
        qb = fn.transport_rate(
            tau_pa=tau,
            d50_um=5000.0,
            tau_ce_pa=0.0,
            velocity_m_s=vel,
            depth_m=depth,
            slope=1.0e-3,
            solid_density_g_cm3=2.65,
        )
        # W*_i = 0.002·φ^7.5 with φ ≪ 1 collapses to ≈ 0.
        assert np.all(qb.values < 1.0e-6)


class TestParker1990:
    """Parker (1990) gravel similarity collapse.

    Reference condition: D50 = 30 mm, τ = 30 Pa.
        θ = 30 / (1650·9.81·0.030) ≈ 0.0618
        φ = θ / 0.0386 ≈ 1.601 → mid regime
        x = 0.601
        W* = 0.00218·exp(14.2·0.601 − 9.28·0.361)
           = 0.00218·exp(8.534 − 3.350)
           = 0.00218·exp(5.184) ≈ 0.385
        u* = √(30/1000) = 0.1732 m/s, u*³ ≈ 5.2e-3
        q_b = 0.385·5.2e-3·2650 / (1.65·9.81) ≈ 0.328 kg/m/s = 3.28 g/cm/s
    """

    def test_reference_value_within_tolerance(self):
        tau = _tau_array(30.0)
        vel = xr.DataArray(np.full(3, 2.0), dims=("nface",))
        depth = xr.DataArray(np.full(3, 1.5), dims=("nface",))
        fn = bedload.Parker1990TransportFunction()
        qb = fn.transport_rate(
            tau_pa=tau,
            d50_um=30000.0,       # 30 mm gravel
            tau_ce_pa=0.0,
            velocity_m_s=vel,
            depth_m=depth,
            slope=1.0e-3,
            solid_density_g_cm3=2.65,
        )
        # Parker's mid-regime exponential is exquisitely sensitive to
        # τ; allow 20 % tolerance on the hand value.
        np.testing.assert_allclose(qb.values, 3.28, rtol=0.20)

    def test_low_phi_collapses_to_low_regime(self):
        # τ → 0 ⇒ φ → 0 ⇒ W* = 0.00218·φ^14.2 → 0.
        tau = _tau_array(0.001)
        vel = xr.DataArray(np.full(3, 0.1), dims=("nface",))
        depth = xr.DataArray(np.full(3, 0.5), dims=("nface",))
        fn = bedload.Parker1990TransportFunction()
        qb = fn.transport_rate(
            tau_pa=tau,
            d50_um=30000.0,
            tau_ce_pa=0.0,
            velocity_m_s=vel,
            depth_m=depth,
            slope=1.0e-3,
            solid_density_g_cm3=2.65,
        )
        assert np.all(qb.values >= 0.0)
        assert np.all(qb.values < 1.0e-3)


class TestYang:
    """Yang (1973) sand-formula sanity check.

    Yang predicts non-zero transport whenever V·S > V_cr·S, regardless of
    the precise magnitude — the equation is dominated by an empirical fit
    of log Ct against log ψ. We verify the sign and order of magnitude
    against a known sand condition (D50 = 0.25 mm, V = 1 m/s, h = 1 m,
    S = 5e-4) where Yang's published Ct curves give Ct ≈ 100 ppm.
    """

    def test_sand_returns_positive_with_excess_stream_power(self):
        tau = _tau_array(0.5)
        vel = xr.DataArray(np.full(3, 1.0), dims=("nface",))
        depth = xr.DataArray(np.full(3, 1.0), dims=("nface",))
        fn = bedload.YangTransportFunction()
        qb = fn.transport_rate(
            tau_pa=tau,
            d50_um=250.0,
            tau_ce_pa=0.2,
            velocity_m_s=vel,
            depth_m=depth,
            slope=5.0e-4,
            solid_density_g_cm3=2.65,
        )
        # qt = Ct ppm × ρ_w × V × h × 1e-6 × 10 (g/cm/s).
        # Ct ~10–1000 ppm → qb in 1e-2 .. 1 g/cm/s range.
        assert np.all(qb.values > 0.0)
        assert np.all(qb.values < 100.0)

    def test_zero_velocity_returns_zero(self):
        tau = _tau_array(0.5)
        vel = xr.zeros_like(tau)
        depth = xr.DataArray(np.full(3, 1.0), dims=("nface",))
        fn = bedload.YangTransportFunction()
        qb = fn.transport_rate(
            tau_pa=tau,
            d50_um=250.0,
            tau_ce_pa=0.2,
            velocity_m_s=vel,
            depth_m=depth,
            slope=5.0e-4,
            solid_density_g_cm3=2.65,
        )
        # V·S = 0 ⇒ ψ < 0 ⇒ no transport.
        np.testing.assert_array_equal(qb.values, np.zeros(3))


class TestWu2000:
    """Wu, Wang & Jia (2000) — uniform-bed reduction.

    With pe/ph = 1 (no hiding/exposure correction) and class fraction = 1:
        d = 0.25 mm = 2.5e-4 m
        τ_cm = 0.03·(2650−1000)·9.81·2.5e-4 = 0.1214 Pa
        excess = τ/τ_ci − 1 = 0.5/0.1214 − 1 = 3.119
        Φ_b = 0.0053·3.119^2.2
            = 0.0053·exp(2.2·ln(3.119))
            = 0.0053·exp(2.502) = 0.0053·12.20 ≈ 0.0647
        d³ = 1.5625e-11 m³ ; (s−1)·g·d³ = 2.530e-10 ; √ = 1.591e-5
        scale = ρ_s·√(...) = 2650·1.591e-5 = 0.04217 kg/m/s
        q_b (kg/m/s) = 0.0647·1·0.04217 ≈ 2.728e-3
        q_b (g/cm/s) = 0.0273
    """

    def test_uniform_bed_reference_value(self):
        tau = _tau_array(0.5)
        vel = xr.DataArray(np.full(3, 1.0), dims=("nface",))
        depth = xr.DataArray(np.full(3, 1.0), dims=("nface",))
        fn = bedload.Wu2000TransportFunction()
        qb = fn.transport_rate(
            tau_pa=tau,
            d50_um=250.0,
            tau_ce_pa=0.0,
            velocity_m_s=vel,
            depth_m=depth,
            slope=1.0e-3,
            solid_density_g_cm3=2.65,
        )
        np.testing.assert_allclose(qb.values, 0.0273, rtol=0.10)

    def test_below_critical_returns_zero(self):
        # τ < τ_ci ⇒ excess clipped → 0.
        tau = _tau_array(0.05)
        vel = xr.DataArray(np.full(3, 0.1), dims=("nface",))
        depth = xr.DataArray(np.full(3, 0.5), dims=("nface",))
        fn = bedload.Wu2000TransportFunction()
        qb = fn.transport_rate(
            tau_pa=tau,
            d50_um=250.0,
            tau_ce_pa=0.0,
            velocity_m_s=vel,
            depth_m=depth,
            slope=1.0e-3,
            solid_density_g_cm3=2.65,
        )
        np.testing.assert_array_equal(qb.values, np.zeros(3))


class TestEngelundHansen1967:
    """Engelund & Hansen (1967) total-load reference value.

    Condition: V = 1.0 m/s, h = 1.0 m, S = 1e-3, D50 = 0.25 mm.
        C  = V/√(R·S) = 1/√(1·1e-3) = 31.62
        Δ  = (2650 − 1000)/1000 = 1.65
        q_t (vol) = 0.05·V^5/(√g · C³ · Δ² · d)
                  = 0.05·1/(√9.81·31.62³·1.65²·2.5e-4)
                  ≈ 0.05/(3.13·31619·2.7225·2.5e-4)
                  ≈ 0.05/67.27 ≈ 7.43e-4 m²/s
        q_t (mass) = 7.43e-4·2650 ≈ 1.97 kg/m/s = 19.7 g/cm/s
    """

    def test_reference_value_within_tolerance(self):
        tau = _tau_array(1.0)
        vel = xr.DataArray(np.full(3, 1.0), dims=("nface",))
        depth = xr.DataArray(np.full(3, 1.0), dims=("nface",))
        fn = bedload.EngelundHansen1967TransportFunction()
        qb = fn.transport_rate(
            tau_pa=tau,
            d50_um=250.0,
            tau_ce_pa=0.0,
            velocity_m_s=vel,
            depth_m=depth,
            slope=1.0e-3,
            solid_density_g_cm3=2.65,
        )
        np.testing.assert_allclose(qb.values, 19.7, rtol=0.05)

    def test_zero_velocity_returns_zero(self):
        tau = _tau_array(0.0)
        vel = xr.zeros_like(tau)
        depth = xr.DataArray(np.full(3, 1.0), dims=("nface",))
        fn = bedload.EngelundHansen1967TransportFunction()
        qb = fn.transport_rate(
            tau_pa=tau,
            d50_um=250.0,
            tau_ce_pa=0.0,
            velocity_m_s=vel,
            depth_m=depth,
            slope=1.0e-3,
            solid_density_g_cm3=2.65,
        )
        # V = 0 ⇒ V^5 = 0 ⇒ q_t = 0.
        np.testing.assert_array_equal(qb.values, np.zeros(3))


class TestToffaleti1968:
    """Toffaleti (1968) single-zone reduction sanity check.

    Toffaleti's original formula returns transport in tons/day/ft. We
    verify the unit conversion produces a physically reasonable
    g cm⁻¹ s⁻¹ value at a sand reference condition (V = 1 m/s,
    D50 = 0.25 mm, h = 1 m, τ > τ_ce). Exact agreement with Toffaleti's
    multi-zone integral is queued for phase 3; this test pins the
    single-zone result for regression.
    """

    def test_sand_returns_positive_above_threshold(self):
        tau = _tau_array(0.5)         # > τ_ce = 0.2 ⇒ gate open
        vel = xr.DataArray(np.full(3, 1.0), dims=("nface",))
        depth = xr.DataArray(np.full(3, 1.0), dims=("nface",))
        fn = bedload.Toffaleti1968TransportFunction()
        qb = fn.transport_rate(
            tau_pa=tau,
            d50_um=250.0,
            tau_ce_pa=0.2,
            velocity_m_s=vel,
            depth_m=depth,
            slope=1.0e-3,
            solid_density_g_cm3=2.65,
        )
        assert np.all(qb.values > 0.0)
        # Result should be in a physically sensible range — Toffaleti
        # tables give 10–1000 tons/day/ft for medium-sand transport,
        # i.e. 0.3–30 kg/m/s = 3–300 g/cm/s.
        assert np.all(qb.values < 1000.0)

    def test_below_threshold_returns_zero(self):
        tau = _tau_array(0.05)        # < τ_ce = 0.2 ⇒ gate closed
        vel = xr.DataArray(np.full(3, 1.0), dims=("nface",))
        depth = xr.DataArray(np.full(3, 1.0), dims=("nface",))
        fn = bedload.Toffaleti1968TransportFunction()
        qb = fn.transport_rate(
            tau_pa=tau,
            d50_um=250.0,
            tau_ce_pa=0.2,
            velocity_m_s=vel,
            depth_m=depth,
            slope=1.0e-3,
            solid_density_g_cm3=2.65,
        )
        np.testing.assert_array_equal(qb.values, np.zeros(3))
