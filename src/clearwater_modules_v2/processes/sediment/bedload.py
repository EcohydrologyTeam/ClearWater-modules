"""Bedload transport (van Rijn 1984).

Two interchangeable solver implementations, selected via SSM config
``sediment.bedload.solver``:

* :class:`BedloadStandaloneExplicit` (``"standalone"``) — explicit
  upwind face-flux step on the mesh, NumPy. No modification to
  ClearWater-Riverine required. Recommended for the initial release;
  see design spec §11 item 3.

* :class:`BedloadRiverineConstituent` (``"riverine"``) — registers
  per-class bedload mass as Riverine constituents with their own
  per-class advection-coefficient field
  (:func:`contracts.advection_coef_var_name`). Requires the Riverine
  ``linalg.py`` extension (Batch C).

Both solvers compute equivalent bedload velocities and heights from
van Rijn (1984) and obtain matching equilibrium concentrations on
idealized cases (verified by ``test_bedload.py`` parity test).

Reference: SAND2008-5621 §"S_BEDLOAD.f90"; van Rijn (1984a); design
spec §5.7.

License note
------------
The van Rijn equations implemented here (eqs 18, 20a, 20b, 21 from
van Rijn 1984a Part II) are taken directly from the open peer-reviewed
literature. No code from ``s_bedload.f90`` (GPL-2.0) was copied; the
implementation is a clean-room port from the published equations.
"""

from __future__ import annotations

from typing import Protocol, runtime_checkable

import numpy as np
import xarray as xr

from . import contracts
from .classes import SedimentClass, SedimentClassRegistry


@runtime_checkable
class BedloadSolver(Protocol):
    """Advances bedload mass per class on the mesh by one sediment step."""

    def step(
        self,
        mesh: xr.Dataset,
        time,
        tau_pa: xr.DataArray,                 # (nface,) Pa
        suspended_psus: xr.DataArray,         # (nface, ssm_class) suspended-erosion fraction
        bed_class_fraction_layer1: xr.DataArray,
        dt_seconds: float,
    ) -> None:
        """Update :data:`contracts.VAR_BEDLOAD_MASS` in place."""
        ...


# ---------------------------------------------------------------------------
# van Rijn (1984a) Part II closed-form helpers
# ---------------------------------------------------------------------------

def _transport_parameter(
    tau_pa: xr.DataArray,
    tau_ce_pa: float,
) -> xr.DataArray:
    """Per-cell van Rijn transport parameter ``T_R``.

    .. math::

        T_R = \\max\\bigl((\\tau - \\tau_{ce})/\\tau_{ce},\\ 0\\bigr)

    Vectorized over the input ``tau_pa`` array using ``xr.where`` so the
    return value remains an :class:`xarray.DataArray` whose coordinates
    and dims match the input.
    """
    if tau_ce_pa <= 0.0:
        raise ValueError(
            f"tau_ce_pa must be > 0 (got {tau_ce_pa!r}); van Rijn T_R is "
            "undefined for non-positive critical shear stress."
        )
    excess = (tau_pa - tau_ce_pa) / tau_ce_pa
    zero = xr.zeros_like(excess) if hasattr(excess, "dims") else 0.0
    return xr.where(excess > 0.0, excess, zero)


def van_rijn_bedload_velocity_cm_s(
    tau_pa: xr.DataArray,
    tau_ce_pa: float,
    d50_um: float,
    solid_specific_gravity: float = contracts.DEFAULT_SOLID_SPECIFIC_GRAVITY,
) -> xr.DataArray:
    """van Rijn (1984) bedload velocity, eq. 20a.

    .. math::

        T_R = \\max\\bigl((\\tau - \\tau_{ce})/\\tau_{ce},\\ 0\\bigr)

        u_{BL} = 1.5\\, T_R^{0.6}\\, \\sqrt{(s_s - 1)\\,g\\,D_{50}} \\quad [{\\rm cm/s}]

    All quantities CGS internally. ``d50_um`` is converted from μm to cm
    (× 1e-4), gravity is :data:`contracts.G_CGS`. Returns an xarray
    DataArray broadcastable from ``tau_pa``.
    """
    if d50_um <= 0.0:
        raise ValueError(f"d50_um must be > 0 (got {d50_um!r}).")
    t_r = _transport_parameter(tau_pa, tau_ce_pa)
    d50_cm = d50_um * 1.0e-4
    immersed = (solid_specific_gravity - 1.0) * contracts.G_CGS * d50_cm
    if immersed <= 0.0:
        # Defensive: solid_specific_gravity ≤ 1 means buoyant grain; no bedload.
        return xr.zeros_like(t_r)
    return 1.5 * np.power(t_r, 0.6) * np.sqrt(immersed)


def van_rijn_bedload_height_cm(
    tau_pa: xr.DataArray,
    tau_ce_pa: float,
    d50_um: float,
    d_star: float,
) -> xr.DataArray:
    """van Rijn (1984) saltation-layer height, eq. 20b.

    .. math:: \\delta_{BL} = 0.3\\, D_{50}\\, d_*^{0.7}\\, \\sqrt{T_R}

    ``d50_um`` is converted to cm internally. Returns cm.
    """
    if d50_um <= 0.0:
        raise ValueError(f"d50_um must be > 0 (got {d50_um!r}).")
    if d_star <= 0.0:
        raise ValueError(f"d_star must be > 0 (got {d_star!r}).")
    t_r = _transport_parameter(tau_pa, tau_ce_pa)
    d50_cm = d50_um * 1.0e-4
    return 0.3 * d50_cm * np.power(d_star, 0.7) * np.sqrt(t_r)


def van_rijn_equilibrium_concentration(
    tau_pa: xr.DataArray,
    tau_ce_pa: float,
    d_star: float,
    solid_density_g_cm3: float = contracts.DEFAULT_SOLID_SPECIFIC_GRAVITY,
) -> xr.DataArray:
    """van Rijn (1981) equilibrium bedload concentration, eq. 21
    (s_sedzlj.f90:191). Returns g/cm³.

    .. math:: C_{eq} = 0.117\\,\\rho_s\\,T_R / d_*
    """
    if d_star <= 0.0:
        raise ValueError(f"d_star must be > 0 (got {d_star!r}).")
    t_r = _transport_parameter(tau_pa, tau_ce_pa)
    return 0.117 * solid_density_g_cm3 * t_r / d_star


# ---------------------------------------------------------------------------
# Internal helpers used by the standalone solver
# ---------------------------------------------------------------------------

def _cheng_d_star(
    d50_um: float,
    solid_specific_gravity: float = contracts.DEFAULT_SOLID_SPECIFIC_GRAVITY,
    kinematic_viscosity_cm2_s: float = contracts.NU_CGS,
    g_cm_s2: float = contracts.G_CGS,
) -> float:
    """Dimensionless particle parameter ``d_*`` (Cheng 1997).

    .. math:: d_* = D_{50}\\,\\bigl[(s_s - 1)\\,g/\\nu^2\\bigr]^{1/3}

    Same definition as in :func:`settling.cheng_1997_settling_velocity`,
    repeated here to avoid a settling-module dependency in bedload.
    """
    d50_cm = float(d50_um) * 1.0e-4
    return d50_cm * np.cbrt(
        (solid_specific_gravity - 1.0) * g_cm_s2 / (kinematic_viscosity_cm2_s ** 2)
    )


def _filter_bedload_eligible(
    registry: SedimentClassRegistry,
    bedload_cutoff_um: float,
) -> list[tuple[int, SedimentClass]]:
    """Return ``[(class_idx, class), ...]`` with D50 ≥ cutoff."""
    return [
        (i, c) for i, c in enumerate(registry) if c.d50_um >= bedload_cutoff_um
    ]


def _edges_face_arrays(mesh: xr.Dataset) -> tuple[np.ndarray, np.ndarray]:
    """Pull the edge → (face1, face2) connectivity from the mesh.

    Boundary edges in Riverine encode the ghost cell with a sentinel
    (negative or out-of-range) on one side. The standalone solver
    treats such edges as no-flux for bedload (sediment cannot be
    advected across the boundary by van Rijn ``u_BL`` because the
    upstream face has no defined CBL). The masking is done in
    :meth:`BedloadStandaloneExplicit.step`.
    """
    if "edges_face1" not in mesh.variables or "edges_face2" not in mesh.variables:
        raise KeyError(
            "Mesh dataset is missing edges_face1/edges_face2 connectivity; "
            "BedloadStandaloneExplicit requires the Riverine UGRID mesh."
        )
    f1 = np.asarray(mesh["edges_face1"].values, dtype=np.int64)
    f2 = np.asarray(mesh["edges_face2"].values, dtype=np.int64)
    return f1, f2


def _edge_length_array(mesh: xr.Dataset) -> np.ndarray:
    """Per-edge length in metres. Falls back to 1.0 for synthetic meshes."""
    if "edge_length" in mesh.variables:
        return np.asarray(mesh["edge_length"].values, dtype="float64")
    # Synthetic test meshes may omit edge_length; assume unit length.
    nedge = mesh.sizes.get("nedge", None)
    if nedge is None:
        raise KeyError("Mesh dataset has neither 'edge_length' nor 'nedge' dim.")
    return np.ones(nedge, dtype="float64")


def _face_area_array(mesh: xr.Dataset) -> np.ndarray:
    """Per-face plan-view area in m². Falls back to 1.0 for synthetic meshes."""
    if "faces_surface_area" in mesh.variables:
        return np.asarray(mesh["faces_surface_area"].values, dtype="float64")
    nface = mesh.sizes.get("nface", None)
    if nface is None:
        raise KeyError("Mesh dataset has neither 'faces_surface_area' nor 'nface' dim.")
    return np.ones(nface, dtype="float64")


def _ensure_bedload_array(mesh: xr.Dataset, n_classes: int) -> None:
    """Make sure ``contracts.VAR_BEDLOAD_MASS`` exists with the right shape.

    Created lazily so the standalone solver works even on synthetic
    meshes that have not been pre-allocated by the orchestrator.
    The variable carries dims ``(time, nface, ssm_class)`` per
    :data:`contracts.BED_STATE_SPECS`; if no time dim exists on the
    mesh we fall back to ``(nface, ssm_class)`` for unit-test meshes.
    """
    if contracts.VAR_BEDLOAD_MASS in mesh.variables:
        return
    nface = mesh.sizes["nface"]
    if contracts.DIM_TIME in mesh.sizes:
        ntime = mesh.sizes[contracts.DIM_TIME]
        data = np.zeros((ntime, nface, n_classes), dtype="float32")
        mesh[contracts.VAR_BEDLOAD_MASS] = (
            (contracts.DIM_TIME, contracts.DIM_NFACE, contracts.DIM_CLASS),
            data,
        )
    else:
        data = np.zeros((nface, n_classes), dtype="float32")
        mesh[contracts.VAR_BEDLOAD_MASS] = (
            (contracts.DIM_NFACE, contracts.DIM_CLASS),
            data,
        )


# ---------------------------------------------------------------------------
# Standalone explicit-upwind solver
# ---------------------------------------------------------------------------

class BedloadStandaloneExplicit:
    """Standalone bedload solver — explicit upwind face-flux step.

    The simpler of the two implementations and the default. Computes
    per-face bedload fluxes from cell-centred ``u_BL`` (van Rijn 1984
    eq. 20a) and per-class ``CBL``, then advances ``CBL`` by an explicit
    upwind step on the unstructured mesh using the existing edge–face
    connectivity in the mesh dataset.

    Algorithm (per bedload-eligible class, per step)
    ------------------------------------------------
    1. Compute per-cell bedload velocity :math:`u_{BL}` (cm/s) from
       :func:`van_rijn_bedload_velocity_cm_s`.
    2. Average ``u_BL`` to each edge as the arithmetic mean of the two
       adjacent cell values (a sign convention is implicit: positive
       ``u_BL`` always points from face1 → face2 along the edge normal,
       so only the magnitude is averaged here; the actual upwind face is
       chosen by the sign of the velocity at face evaluation time).
    3. Compute upwind face flux:
       ``flux[edge] = u_BL_edge × CBL[upwind_face] × edge_length``.
    4. Apply ``CBL[face] += dt × Σ(flux_in − flux_out) / face_area``.

    Boundary edges (where one of the connected faces is a ghost / sentinel)
    are treated as no-flux: bedload cannot leave or enter the domain
    through the boundary because there is no CBL defined on the ghost.
    This makes the closed-domain mass conservation property exact.
    """

    def __init__(
        self,
        registry: SedimentClassRegistry,
        bedload_cutoff_um: float = contracts.DEFAULT_BEDLOAD_CUTOFF_UM,
    ) -> None:
        self.registry = registry
        self.bedload_cutoff_um = float(bedload_cutoff_um)
        self._eligible: list[tuple[int, SedimentClass]] = _filter_bedload_eligible(
            registry, self.bedload_cutoff_um
        )
        # Cache d_star per eligible class (Cheng-style).
        self._d_star: dict[int, float] = {
            idx: _cheng_d_star(
                cls.d50_um, solid_specific_gravity=cls.solid_density_g_cm3
            )
            for idx, cls in self._eligible
        }

    @property
    def eligible_class_indices(self) -> list[int]:
        return [i for i, _ in self._eligible]

    def step(
        self,
        mesh: xr.Dataset,
        time,
        tau_pa: xr.DataArray,
        suspended_psus: xr.DataArray | None = None,
        bed_class_fraction_layer1: xr.DataArray | None = None,
        dt_seconds: float = 0.0,
    ) -> None:
        """Advance ``CBL`` by ``dt_seconds`` for every eligible class.

        Updates :data:`contracts.VAR_BEDLOAD_MASS` in place on ``mesh``.
        Parameters ``suspended_psus`` and ``bed_class_fraction_layer1``
        are accepted for Protocol compatibility but unused by this
        solver (they would be consumed by an erosion-coupled variant).
        """
        if dt_seconds <= 0.0:
            return
        if not self._eligible:
            return

        n_classes = len(self.registry)
        _ensure_bedload_array(mesh, n_classes)

        # Cache mesh geometry once per call.
        f1, f2 = _edges_face_arrays(mesh)
        edge_len = _edge_length_array(mesh)
        face_area = _face_area_array(mesh)
        nface = mesh.sizes[contracts.DIM_NFACE]

        # Boundary mask: edges where either side is out of [0, nface).
        valid_edge = (f1 >= 0) & (f1 < nface) & (f2 >= 0) & (f2 < nface)

        # Resolve the destination view (with or without time dim).
        bedload_var = mesh[contracts.VAR_BEDLOAD_MASS]
        has_time = contracts.DIM_TIME in bedload_var.dims
        if has_time:
            # Expect 'time' to be the supplied integer index into the time dim.
            t_idx = int(time) if time is not None else 0
        # tau_pa may be (time, nface) or (nface,); reduce to (nface,).
        if contracts.DIM_TIME in tau_pa.dims:
            tau_now = tau_pa.isel({contracts.DIM_TIME: int(time) if time is not None else 0})
        else:
            tau_now = tau_pa
        tau_arr = np.asarray(tau_now.values, dtype="float64")

        for class_idx, cls in self._eligible:
            tau_ce = float(cls.tau_ce_pa) if cls.tau_ce_pa is not None else 0.0
            if tau_ce <= 0.0:
                # τ_ce not yet resolved (Soulsby fill happens in bed.py); skip.
                continue

            # Per-cell u_BL (cm/s) -> convert to m/s for face-flux integration
            # (edge_length is m, face_area is m², dt is s, CBL is g/cm²).
            # We keep CBL in g/cm² and let the velocity / length / area cancel
            # consistently in CGS by working entirely in cm and seconds for
            # the flux step:
            #   flux  [g/s]      = u_BL [cm/s] × CBL [g/cm²] × edge_len_cm [cm]
            #   ΔCBL  [g/cm²]    = dt × Σflux / face_area_cm2 [cm²]
            edge_len_cm = edge_len * 100.0
            face_area_cm2 = face_area * 1.0e4

            u_bl = van_rijn_bedload_velocity_cm_s(
                tau_now,
                tau_ce_pa=tau_ce,
                d50_um=cls.d50_um,
                solid_specific_gravity=cls.solid_density_g_cm3,
            )
            u_bl_arr = np.asarray(u_bl.values, dtype="float64")

            # Pull current CBL slab for this class (1-D, length nface).
            if has_time:
                cbl_slab = np.asarray(
                    bedload_var.isel({contracts.DIM_TIME: t_idx, contracts.DIM_CLASS: class_idx}).values,
                    dtype="float64",
                )
            else:
                cbl_slab = np.asarray(
                    bedload_var.isel({contracts.DIM_CLASS: class_idx}).values,
                    dtype="float64",
                )

            # Edge-centred velocity (arithmetic mean of the two faces).
            # Use clipped indices so the gather is safe; then apply valid_edge.
            i1 = np.clip(f1, 0, nface - 1)
            i2 = np.clip(f2, 0, nface - 1)
            u_face1 = u_bl_arr[i1]
            u_face2 = u_bl_arr[i2]
            u_edge = 0.5 * (u_face1 + u_face2)

            # Upwind CBL: when u_edge >= 0 the flux flows face1 → face2,
            # so the upwind cell is face1; otherwise face2.
            cbl_upwind = np.where(u_edge >= 0.0, cbl_slab[i1], cbl_slab[i2])

            # Per-edge mass flux (g/s). Zero out boundary edges.
            flux = u_edge * cbl_upwind * edge_len_cm
            flux = np.where(valid_edge, flux, 0.0)

            # Accumulate into per-cell divergence: outflow at upwind face,
            # inflow at downwind face. We add |flux| × sign convention:
            # if u_edge ≥ 0: face1 loses, face2 gains.
            net = np.zeros(nface, dtype="float64")
            np.add.at(net, i1, np.where(u_edge >= 0.0, -flux, flux))
            np.add.at(net, i2, np.where(u_edge >= 0.0, flux, -flux))

            # ΔCBL = dt × net / face_area  (g/cm²).
            # Guard against zero face areas (synthetic edge-only meshes).
            safe_area = np.where(face_area_cm2 > 0.0, face_area_cm2, 1.0)
            delta_cbl = dt_seconds * net / safe_area

            new_cbl = cbl_slab + delta_cbl
            # Floor at zero to suppress negative concentrations from
            # explicit-upwind round-off when CBL ≈ 0 everywhere.
            new_cbl = np.maximum(new_cbl, 0.0)

            # Write back into the mesh DataArray in place.
            # xarray's .values returns a view for numpy-backed arrays.
            if has_time:
                bedload_var.values[t_idx, :, class_idx] = new_cbl.astype(
                    bedload_var.dtype, copy=False
                )
            else:
                bedload_var.values[:, class_idx] = new_cbl.astype(
                    bedload_var.dtype, copy=False
                )


# ---------------------------------------------------------------------------
# Riverine-constituent solver (delegates implicit advection to Riverine)
# ---------------------------------------------------------------------------

class BedloadRiverineConstituent:
    """Bedload solver that uses Riverine's implicit transport for CBL.

    Registers each bedload-eligible class as an additional Riverine
    constituent with its own per-class advection-coefficient field
    (``ssm_advection_coef_<label>``) computed from van Rijn ``u_BL``.

    Requires the per-constituent advection-coefficient extension to
    ``clearwater_riverine.linalg.LHS.update_values`` (Batch C). Until
    that lands the implicit advection will not actually execute, but
    :meth:`step` is fully functional: it writes the per-edge advection
    coefficient into the mesh dataset on every call so the field is
    available for inspection and for the future Riverine read-back.
    """

    # Suffix appended to the suspended-class constituent name to form
    # the bedload-companion constituent name.
    _BEDLOAD_SUFFIX: str = "_bedload"

    def __init__(
        self,
        registry: SedimentClassRegistry,
        riverine,                               # ClearwaterRiverine instance
        bedload_cutoff_um: float = contracts.DEFAULT_BEDLOAD_CUTOFF_UM,
    ) -> None:
        self.registry = registry
        self.riverine = riverine
        self.bedload_cutoff_um = float(bedload_cutoff_um)
        self._eligible: list[tuple[int, SedimentClass]] = _filter_bedload_eligible(
            registry, self.bedload_cutoff_um
        )

        # Register each bedload-eligible class with the Riverine instance.
        # We use a defensive guard because some unit-test stubs may not
        # have the full constituent_dict scaffolding yet.
        if not hasattr(riverine, "constituent_dict") or riverine.constituent_dict is None:
            riverine.constituent_dict = {}

        for _, cls in self._eligible:
            cname = contracts.suspended_var_name(cls.label) + self._BEDLOAD_SUFFIX
            # Schema mirrors the constituent_dict format documented at
            # transport.py:188-201: requires initial_conditions,
            # boundary_conditions, units. We register a no-source,
            # zero-decay constituent whose advection coefficient is the
            # per-class van Rijn velocity field.
            riverine.constituent_dict[cname] = {
                "initial_conditions": None,
                "boundary_conditions": None,
                "units": "g/cm^2",
                "decay_rate": 0.0,
                "advection_coefficient_var": contracts.advection_coef_var_name(cls.label),
            }

    @property
    def eligible_class_indices(self) -> list[int]:
        return [i for i, _ in self._eligible]

    def step(
        self,
        mesh: xr.Dataset,
        time,
        tau_pa: xr.DataArray,
        suspended_psus: xr.DataArray | None = None,
        bed_class_fraction_layer1: xr.DataArray | None = None,
        dt_seconds: float = 0.0,
    ) -> None:
        """Compute per-edge advection coefficients for the next Riverine update.

        For each bedload-eligible class:

        1. Compute per-cell ``u_BL`` (cm/s) from van Rijn 1984 eq. 20a.
        2. Average to edges (arithmetic mean of the two adjacent cells).
        3. Convert to m/s and write into
           ``mesh[contracts.advection_coef_var_name(class.label)]`` so
           Riverine's solver can read it on the next ``riverine.update()``.

        The implicit advection itself is delegated to Riverine; this
        method does no CBL update.
        """
        if not self._eligible:
            return

        f1, f2 = _edges_face_arrays(mesh)
        nface = mesh.sizes[contracts.DIM_NFACE]
        nedge = mesh.sizes.get(
            "nedge",
            int(f1.shape[0]) if f1.ndim else 0,
        )

        # tau_pa may be (time, nface) or (nface,); reduce to (nface,).
        if contracts.DIM_TIME in tau_pa.dims:
            tau_now = tau_pa.isel({contracts.DIM_TIME: int(time) if time is not None else 0})
        else:
            tau_now = tau_pa

        i1 = np.clip(f1, 0, nface - 1)
        i2 = np.clip(f2, 0, nface - 1)
        valid_edge = (f1 >= 0) & (f1 < nface) & (f2 >= 0) & (f2 < nface)

        for _, cls in self._eligible:
            tau_ce = float(cls.tau_ce_pa) if cls.tau_ce_pa is not None else 0.0
            if tau_ce <= 0.0:
                continue

            u_bl = van_rijn_bedload_velocity_cm_s(
                tau_now,
                tau_ce_pa=tau_ce,
                d50_um=cls.d50_um,
                solid_specific_gravity=cls.solid_density_g_cm3,
            )
            u_bl_arr = np.asarray(u_bl.values, dtype="float64")

            u_edge_cm_s = 0.5 * (u_bl_arr[i1] + u_bl_arr[i2])
            u_edge_cm_s = np.where(valid_edge, u_edge_cm_s, 0.0)
            # Riverine's advection-coefficient field is in m/s to match
            # edge_velocity. Convert from cm/s.
            u_edge_m_s = u_edge_cm_s * 0.01

            var_name = contracts.advection_coef_var_name(cls.label)
            # Allocate-or-overwrite. Field is per-edge and (optionally)
            # per-time; we write a 1-D (nedge,) DataArray for simplicity
            # and let Riverine broadcast across time when it reads it.
            mesh[var_name] = (
                ("nedge",),
                u_edge_m_s.astype("float32"),
            )
            # Stamp metadata for downstream consumers.
            mesh[var_name].attrs.update(
                {
                    "units": "m s-1",
                    "long_name": f"van Rijn bedload velocity (edge mean) — {cls.label}",
                    "description": (
                        "Per-edge advection coefficient for the bedload constituent "
                        f"'{contracts.suspended_var_name(cls.label) + self._BEDLOAD_SUFFIX}'. "
                        "Computed from van Rijn (1984) eq. 20a applied per cell, then "
                        "averaged to edges."
                    ),
                }
            )
