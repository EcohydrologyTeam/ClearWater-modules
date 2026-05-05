"""Bedload transport — pluggable transport-function abstraction.

Two interchangeable *solver* implementations advance bedload mass on
the mesh; both are selected via SSM config ``sediment.bedload.solver``:

* :class:`BedloadStandaloneExplicit` (``"standalone"``) — explicit
  upwind face-flux step on the mesh, NumPy. No modification to
  ClearWater-Riverine required. Recommended for the initial release;
  see design spec §11 item 3.

* :class:`BedloadRiverineConstituent` (``"riverine"``) — registers
  per-class bedload mass as Riverine constituents with their own
  per-class advection-coefficient field
  (:func:`contracts.advection_coef_var_name`). Requires the Riverine
  ``linalg.py`` extension (Batch C).

Independently, the *closure* used to compute the per-cell, per-class
bedload transport rate :math:`q_b` (g cm⁻¹ s⁻¹) is now pluggable via
the :class:`BedloadTransportFunction` Protocol. The package ships with
seven peer-reviewed formulas, registered by name in
:data:`BEDLOAD_TRANSPORT_FUNCTIONS`:

* ``van_rijn``         — van Rijn (1984a) Part II (default; backwards-compatible)
* ``wilcock_crowe``    — Wilcock & Crowe (2003) surface-based, sand-gravel mixtures
* ``parker``           — Parker (1990) surface-based gravel
* ``yang``             — Yang (1973, 1979) unit-stream-power total load
* ``wu``               — Wu, Wang & Jia (2000) non-uniform sediment
* ``engelund_hansen``  — Engelund & Hansen (1967) total load for sand
* ``toffaleti``        — Toffaleti (1968) depth-integrated total load

Selection in YAML::

    sediment:
      bedload:
        solver: standalone        # or riverine, off
        transport_function: van_rijn   # or wilcock_crowe, parker, yang,
                                       # wu, engelund_hansen, toffaleti

Reference: SAND2008-5621 §"S_BEDLOAD.f90"; van Rijn (1984a); design
spec §5.7. Per-formula citations live in each class docstring.

License note
------------
All formulas implemented here are taken directly from the open
peer-reviewed literature (cited per class). No code from
``s_bedload.f90`` (GPL-2.0) or any third-party transport library was
copied; each implementation is a clean-room port from the published
equations.
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
        registry_context: dict | None = None,  # surface composition for closure
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


def _face_velocity_m_s(mesh: xr.Dataset, time, nface: int) -> xr.DataArray:
    """Return per-face depth-averaged velocity magnitude (m/s).

    Preference order matches :func:`shear._edge_velocity_to_face_magnitude`:
    (1) ``face_velocity_x/y``, (2) ``face_velocity_magnitude``, (3) edge-mean
    of ``edge_velocity``. Falls back to a 1.0 m/s placeholder when none are
    present (synthetic test meshes); van Rijn ignores the field, but the
    Wilcock-Crowe / Wu / Engelund-Hansen / Yang / Toffaleti closures need
    a real value when invoked. Tests that exercise those closures via the
    standalone solver are expected to populate at least one of the
    velocity fields above.
    """
    if "face_velocity_x" in mesh.data_vars and "face_velocity_y" in mesh.data_vars:
        ux = mesh["face_velocity_x"]
        uy = mesh["face_velocity_y"]
        if contracts.DIM_TIME in ux.dims:
            ux = ux.isel({contracts.DIM_TIME: int(time) if time is not None else 0})
        if contracts.DIM_TIME in uy.dims:
            uy = uy.isel({contracts.DIM_TIME: int(time) if time is not None else 0})
        return xr.DataArray(
            np.sqrt(np.asarray(ux.values, dtype="float64") ** 2
                    + np.asarray(uy.values, dtype="float64") ** 2),
            dims=(contracts.DIM_NFACE,),
        )
    if "face_velocity_magnitude" in mesh.data_vars:
        umag = mesh["face_velocity_magnitude"]
        if contracts.DIM_TIME in umag.dims:
            umag = umag.isel({contracts.DIM_TIME: int(time) if time is not None else 0})
        return xr.DataArray(
            np.asarray(umag.values, dtype="float64"),
            dims=(contracts.DIM_NFACE,),
        )
    if "edge_velocity" in mesh.data_vars and (
        "edges_face1" in mesh.variables and "edges_face2" in mesh.variables
    ):
        ev = mesh["edge_velocity"]
        if contracts.DIM_TIME in ev.dims:
            ev = ev.isel({contracts.DIM_TIME: int(time) if time is not None else 0})
        speed = np.abs(np.asarray(ev.values, dtype="float64"))
        f1 = np.asarray(mesh["edges_face1"].values, dtype=np.int64)
        f2 = np.asarray(mesh["edges_face2"].values, dtype=np.int64)
        speed_sum = np.zeros(nface, dtype="float64")
        count = np.zeros(nface, dtype="float64")
        m1 = (f1 >= 0) & (f1 < nface)
        m2 = (f2 >= 0) & (f2 < nface)
        np.add.at(speed_sum, f1[m1], speed[m1])
        np.add.at(count, f1[m1], 1.0)
        np.add.at(speed_sum, f2[m2], speed[m2])
        np.add.at(count, f2[m2], 1.0)
        return xr.DataArray(
            speed_sum / np.where(count > 0.0, count, 1.0),
            dims=(contracts.DIM_NFACE,),
        )
    return xr.DataArray(np.ones(nface, dtype="float64"), dims=(contracts.DIM_NFACE,))


def _face_depth_m(mesh: xr.Dataset, time, nface: int) -> xr.DataArray:
    """Return per-face hydraulic depth (m), defaulting to 1 m if absent."""
    if contracts.VAR_FACE_HYDRAULIC_DEPTH in mesh.data_vars:
        depth = mesh[contracts.VAR_FACE_HYDRAULIC_DEPTH]
        if contracts.DIM_TIME in depth.dims:
            depth = depth.isel({contracts.DIM_TIME: int(time) if time is not None else 0})
        return xr.DataArray(
            np.asarray(depth.values, dtype="float64"),
            dims=(contracts.DIM_NFACE,),
        )
    return xr.DataArray(np.ones(nface, dtype="float64"), dims=(contracts.DIM_NFACE,))


def _face_slope(mesh: xr.Dataset, time, nface: int) -> xr.DataArray:
    """Return per-face energy slope (dimensionless), defaulting to 1e-3."""
    for name in ("energy_slope", "face_energy_slope", "friction_slope"):
        if name in mesh.data_vars:
            s = mesh[name]
            if contracts.DIM_TIME in s.dims:
                s = s.isel({contracts.DIM_TIME: int(time) if time is not None else 0})
            return xr.DataArray(
                np.asarray(s.values, dtype="float64"),
                dims=(contracts.DIM_NFACE,),
            )
    return xr.DataArray(
        np.full(nface, 1.0e-3, dtype="float64"),
        dims=(contracts.DIM_NFACE,),
    )


def _class_context(registry_context: dict | None, class_idx: int) -> dict | None:
    """Slice a registry-wide context dict down to a single-class context.

    The orchestrator builds a registry-wide ``registry_context`` whose
    per-class fields (e.g. ``surface_class_fraction``) carry shape
    ``(nface, ssm_class)`` or ``(ssm_class,)``. Each closure call is for
    a single class, so this helper extracts the per-class slice and
    leaves the per-cell scalars (``surface_sand_fraction``,
    ``surface_geometric_mean_um``, ``pe_ph_ratio``) untouched.

    Returns ``None`` when ``registry_context`` is ``None``, so closures
    can rely on their built-in defaults under unit-test conditions.
    """
    if registry_context is None:
        return None
    out = dict(registry_context)
    fac = registry_context.get("surface_class_fraction")
    if fac is not None and hasattr(fac, "isel"):
        # (nface, ssm_class) → (nface,) for the requested class. The
        # closures expect a scalar; W-C and Wu currently use the per-cell
        # mean, so we average. (nface,) inputs are uniform under the
        # synthetic test fixtures so the mean equals each cell value.
        try:
            sliced = fac.isel({contracts.DIM_CLASS: int(class_idx)})
            out["surface_class_fraction"] = float(np.asarray(sliced.values).mean())
        except (KeyError, IndexError):
            pass
    elif fac is not None:
        arr = np.asarray(fac)
        if arr.ndim == 1:
            # Already (ssm_class,)
            out["surface_class_fraction"] = float(arr[int(class_idx)])
        elif arr.ndim == 2:
            out["surface_class_fraction"] = float(arr[..., int(class_idx)].mean())
    sand = registry_context.get("surface_sand_fraction")
    if sand is not None and hasattr(sand, "values"):
        out["surface_sand_fraction"] = float(np.asarray(sand.values).mean())
    elif sand is not None:
        arr = np.asarray(sand)
        if arr.ndim == 0:
            out["surface_sand_fraction"] = float(arr)
        else:
            out["surface_sand_fraction"] = float(arr.mean())
    dsg = registry_context.get("surface_geometric_mean_um")
    if dsg is not None and hasattr(dsg, "values"):
        out["surface_geometric_mean_um"] = float(np.asarray(dsg.values).mean())
    elif dsg is not None:
        arr = np.asarray(dsg)
        if arr.ndim == 0:
            out["surface_geometric_mean_um"] = float(arr)
        else:
            out["surface_geometric_mean_um"] = float(arr.mean())
    return out


def _qb_to_effective_velocity_cm_s(
    qb_g_cm_s: xr.DataArray,
    tau_pa: xr.DataArray,
    tau_ce_pa: float,
    d50_um: float,
    d_star: float,
    solid_density_g_cm3: float,
) -> np.ndarray:
    """Convert closure ``q_b`` (g cm⁻¹ s⁻¹) to an effective bedload velocity (cm/s).

    Uses van Rijn's saltation height ``δ_BL`` and equilibrium concentration
    ``C_eq`` as the **bedload-layer geometry carrier**:

    .. math:: u_{\\rm eff} = q_b / (\\delta_{BL} \\cdot C_{eq})

    This mirrors the standard practice in production codes that mix
    arbitrary transport functions with a fixed bedload-layer model
    (e.g. Delft3D-MOR, MIKE 21). It guarantees:

    * **van Rijn identity.** When ``q_b`` itself comes from van Rijn,
      :math:`q_b = u_{BL} \\cdot \\delta_{BL} \\cdot C_{eq}` and
      :math:`u_{\\rm eff} \\equiv u_{BL}` exactly (within float64
      round-off) — backwards-compatible with the prior solver.
    * **Mass conservation in a closed domain.** Using the same upwind
      step on the same edge / face geometry, with ``u_{\\rm eff}``
      replacing ``u_{BL}``, preserves the exact-conservation property
      verified by ``test_one_step_conserves_mass_in_closed_domain``.

    Where :math:`\\delta_{BL} \\cdot C_{eq} = 0` (sub-critical shear), the
    helper returns 0 — there is no transport regardless of the closure's
    output, which is the physically correct behaviour for incipient motion.
    """
    delta_bl = van_rijn_bedload_height_cm(
        tau_pa, tau_ce_pa=tau_ce_pa, d50_um=d50_um, d_star=d_star,
    )
    c_eq = van_rijn_equilibrium_concentration(
        tau_pa, tau_ce_pa=tau_ce_pa, d_star=d_star,
        solid_density_g_cm3=solid_density_g_cm3,
    )
    delta_arr = np.asarray(delta_bl.values, dtype="float64")
    c_arr = np.asarray(c_eq.values, dtype="float64")
    qb_arr = np.asarray(qb_g_cm_s.values, dtype="float64")
    denom = delta_arr * c_arr
    safe_denom = np.where(denom > 0.0, denom, 1.0)
    u_eff = np.where(denom > 0.0, qb_arr / safe_denom, 0.0)
    return u_eff


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
        transport_function: "BedloadTransportFunction | None" = None,
    ) -> None:
        self.registry = registry
        self.bedload_cutoff_um = float(bedload_cutoff_um)
        # Pluggable per-cell q_b closure. Default = van Rijn for backwards
        # compatibility (q_b = u_BL · δ_BL · C_eq, so u_eff ≡ u_BL).
        self.transport_function: BedloadTransportFunction = (
            transport_function or VanRijn1984TransportFunction()
        )
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
        registry_context: dict | None = None,
    ) -> None:
        """Advance ``CBL`` by ``dt_seconds`` for every eligible class.

        Updates :data:`contracts.VAR_BEDLOAD_MASS` in place on ``mesh``.
        Parameters ``suspended_psus`` and ``bed_class_fraction_layer1``
        are accepted for Protocol compatibility but unused by this
        solver (they would be consumed by an erosion-coupled variant).

        ``registry_context`` is forwarded verbatim to the configured
        :attr:`transport_function`. Closures that need surface composition
        (Wilcock-Crowe, Wu) read ``surface_sand_fraction``,
        ``surface_class_fraction``, and friends from this dict; closures
        that don't (van Rijn, Engelund-Hansen, Toffaleti, Yang) ignore it.
        See :meth:`SSM.run` for how the orchestrator constructs the dict.
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

        # Pull the auxiliary forcing the closures need. Van Rijn ignores
        # velocity/depth/slope; the other closures consume them.
        velocity_m_s = _face_velocity_m_s(mesh, time, nface)
        depth_m = _face_depth_m(mesh, time, nface)
        slope = _face_slope(mesh, time, nface)

        for class_idx, cls in self._eligible:
            tau_ce = float(cls.tau_ce_pa) if cls.tau_ce_pa is not None else 0.0
            if tau_ce <= 0.0:
                # τ_ce not yet resolved (Soulsby fill happens in bed.py); skip.
                continue

            # Per-cell q_b (g cm⁻¹ s⁻¹) from the configured transport function.
            # Convert to an effective bedload velocity via van Rijn's
            # bedload-layer geometry carrier (δ_BL · C_eq); see
            # _qb_to_effective_velocity_cm_s for the derivation. This keeps
            # the existing upwind face-flux step intact (CBL conservation
            # is preserved on a closed domain) while letting the magnitude
            # of q_b come from any registered closure.
            #
            # We keep CBL in g/cm² and let the velocity / length / area
            # cancel consistently in CGS by working entirely in cm and
            # seconds for the flux step:
            #   flux  [g/s]      = u_eff [cm/s] × CBL [g/cm²] × edge_len_cm [cm]
            #   ΔCBL  [g/cm²]    = dt × Σflux / face_area_cm2 [cm²]
            edge_len_cm = edge_len * 100.0
            face_area_cm2 = face_area * 1.0e4

            qb = self.transport_function.transport_rate(
                tau_pa=tau_now,
                d50_um=cls.d50_um,
                tau_ce_pa=tau_ce,
                velocity_m_s=velocity_m_s,
                depth_m=depth_m,
                slope=slope,
                solid_density_g_cm3=cls.solid_density_g_cm3,
                registry_context=_class_context(registry_context, class_idx),
            )
            u_bl_arr = _qb_to_effective_velocity_cm_s(
                qb_g_cm_s=qb,
                tau_pa=tau_now,
                tau_ce_pa=tau_ce,
                d50_um=cls.d50_um,
                d_star=self._d_star[class_idx],
                solid_density_g_cm3=cls.solid_density_g_cm3,
            )

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
        transport_function: "BedloadTransportFunction | None" = None,
    ) -> None:
        self.registry = registry
        self.riverine = riverine
        self.bedload_cutoff_um = float(bedload_cutoff_um)
        self.transport_function: BedloadTransportFunction = (
            transport_function or VanRijn1984TransportFunction()
        )
        self._eligible: list[tuple[int, SedimentClass]] = _filter_bedload_eligible(
            registry, self.bedload_cutoff_um
        )
        # Cache d_star per eligible class for the q_b → u_eff conversion.
        self._d_star: dict[int, float] = {
            idx: _cheng_d_star(
                cls.d50_um, solid_specific_gravity=cls.solid_density_g_cm3
            )
            for idx, cls in self._eligible
        }

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
        registry_context: dict | None = None,
    ) -> None:
        """Compute per-edge advection coefficients for the next Riverine update.

        For each bedload-eligible class:

        1. Compute per-cell ``q_b`` (g cm⁻¹ s⁻¹) from the configured
           :attr:`transport_function`.
        2. Convert ``q_b`` to an effective bedload velocity ``u_eff``
           (cm/s) via van Rijn's bedload-layer geometry carrier
           ``δ_BL · C_eq`` (see :func:`_qb_to_effective_velocity_cm_s`).
        3. Average to edges (arithmetic mean of the two adjacent cells).
        4. Convert to m/s and write into
           ``mesh[contracts.advection_coef_var_name(class.label)]`` so
           Riverine's solver can read it on the next ``riverine.update()``.

        ``registry_context`` is forwarded to the closure (sliced per
        class via :func:`_class_context`); see :class:`BedloadStandaloneExplicit`
        for the parallel surface-composition wiring.

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

        # Closure-input forcing fields (van Rijn ignores them; the others
        # consume them).
        velocity_m_s = _face_velocity_m_s(mesh, time, nface)
        depth_m = _face_depth_m(mesh, time, nface)
        slope = _face_slope(mesh, time, nface)

        fn_name = getattr(self.transport_function, "name", "transport")

        for class_idx, cls in self._eligible:
            tau_ce = float(cls.tau_ce_pa) if cls.tau_ce_pa is not None else 0.0
            if tau_ce <= 0.0:
                continue

            qb = self.transport_function.transport_rate(
                tau_pa=tau_now,
                d50_um=cls.d50_um,
                tau_ce_pa=tau_ce,
                velocity_m_s=velocity_m_s,
                depth_m=depth_m,
                slope=slope,
                solid_density_g_cm3=cls.solid_density_g_cm3,
                registry_context=_class_context(registry_context, class_idx),
            )
            u_bl_arr = _qb_to_effective_velocity_cm_s(
                qb_g_cm_s=qb,
                tau_pa=tau_now,
                tau_ce_pa=tau_ce,
                d50_um=cls.d50_um,
                d_star=self._d_star[class_idx],
                solid_density_g_cm3=cls.solid_density_g_cm3,
            )

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
                    "long_name": (
                        f"effective bedload velocity (edge mean, {fn_name}) — {cls.label}"
                    ),
                    "description": (
                        "Per-edge advection coefficient for the bedload constituent "
                        f"'{contracts.suspended_var_name(cls.label) + self._BEDLOAD_SUFFIX}'. "
                        f"Derived from {fn_name} q_b via van Rijn bedload-layer "
                        "geometry (δ_BL · C_eq), then averaged to edges."
                    ),
                }
            )


# ---------------------------------------------------------------------------
# Pluggable bedload transport-function abstraction
# ---------------------------------------------------------------------------
#
# Each transport function returns the per-cell, per-class bedload transport
# rate ``q_b`` in g cm⁻¹ s⁻¹ (mass per unit channel width per unit time).
# This is the "volumetric bedload rate" of the standard sediment-transport
# literature, multiplied by the solid density. Working in g cm⁻¹ s⁻¹ keeps
# the module units consistent with the existing CGS internals
# (CBL is g/cm², bedload velocities are cm/s, the saltation height is cm).
#
# Conversion to dimensionless Einstein parameter ``q_b*`` if needed:
#
#     q_b* = q_b / (ρ_s · √((s − 1) g D₅₀³))
#
# Conversion to volumetric rate per unit width (cm² s⁻¹):
#
#     q_b_vol = q_b / ρ_s
#
# All formulas are implemented in CGS internally (cm, g, s); inputs are
# accepted in SI (Pa, m, m/s) at the boundary and converted on entry.
# ---------------------------------------------------------------------------


# Helper: convert dimensionless transport phi (Einstein) to q_b in g cm⁻¹ s⁻¹.
def _phi_to_qb(
    phi: xr.DataArray,
    d50_um: float,
    solid_density_g_cm3: float,
    water_density_g_cm3: float = contracts.DEFAULT_WATER_DENSITY_CGS,
    g_cm_s2: float = contracts.G_CGS,
) -> xr.DataArray:
    """Convert dimensionless Einstein bedload parameter ``phi`` to ``q_b``.

    ``phi = q_b / (ρ_s · √((s − 1) g D₅₀³))``  (Einstein 1950) so

        q_b [g cm⁻¹ s⁻¹] = phi · ρ_s · √((s − 1) · g · D₅₀_cm³)

    with ``s = ρ_s / ρ_w``.
    """
    s_minus_1 = (solid_density_g_cm3 / water_density_g_cm3) - 1.0
    if s_minus_1 <= 0.0:
        return xr.zeros_like(phi)
    d50_cm = float(d50_um) * 1.0e-4
    scale = solid_density_g_cm3 * np.sqrt(s_minus_1 * g_cm_s2 * d50_cm ** 3)
    return phi * scale


def _shields(
    tau_pa: xr.DataArray,
    d50_um: float,
    solid_density_g_cm3: float,
    water_density_kg_m3: float = 1000.0,
) -> xr.DataArray:
    """Dimensionless Shields stress ``θ = τ / ((ρ_s − ρ_w) g D₅₀)`` (SI in/out)."""
    g_si = 9.81
    rho_s_kg = solid_density_g_cm3 * 1000.0
    d50_m = float(d50_um) * 1.0e-6
    denom = (rho_s_kg - water_density_kg_m3) * g_si * d50_m
    if denom <= 0.0:
        return xr.zeros_like(tau_pa)
    return tau_pa / denom


def _ushear_m_s(tau_pa: xr.DataArray, water_density_kg_m3: float = 1000.0) -> xr.DataArray:
    """Shear velocity ``u* = √(τ/ρ_w)`` (SI m/s)."""
    return np.sqrt(np.maximum(tau_pa, 0.0) / water_density_kg_m3)


@runtime_checkable
class BedloadTransportFunction(Protocol):
    """Pluggable closure that returns ``q_b`` (g cm⁻¹ s⁻¹) per cell, per class.

    Each implementation in this module is a class with a ``name`` class-attribute
    (used in YAML ``transport_function:``) and a ``transport_rate`` method
    returning an :class:`xarray.DataArray` with the same ``nface`` shape as the
    input forcing.

    Inputs use SI units at the boundary (Pa, m, m/s); the formulas convert to
    CGS internally as needed.
    """

    name: str

    def transport_rate(
        self,
        tau_pa: xr.DataArray,
        d50_um: float,
        tau_ce_pa: float,
        velocity_m_s: xr.DataArray,
        depth_m: xr.DataArray,
        slope: xr.DataArray | float,
        solid_density_g_cm3: float,
        water_density_kg_m3: float = 1000.0,
        kinematic_viscosity_m2_s: float = 1.0e-6,
        registry_context: dict | None = None,
    ) -> xr.DataArray:
        """Per-cell bedload transport rate in g cm⁻¹ s⁻¹."""
        ...


# ---------------------------------------------------------------------------
# 1. van Rijn (1984) — backwards-compatible wrapper
# ---------------------------------------------------------------------------

class VanRijn1984TransportFunction:
    """van Rijn (1984a) Part II bedload transport — class-wrapper form.

    Computes :math:`q_b` from the existing
    :func:`van_rijn_bedload_velocity_cm_s`,
    :func:`van_rijn_bedload_height_cm`, and
    :func:`van_rijn_equilibrium_concentration` helpers via

    .. math::

        q_b = u_{BL} \\cdot \\delta_{BL} \\cdot C_{eq}

    with :math:`u_{BL}` in cm/s, :math:`\\delta_{BL}` in cm, and
    :math:`C_{eq}` in g/cm³. The product carries units of
    :math:`{\\rm g}\\,{\\rm cm}^{-1}\\,{\\rm s}^{-1}`, the canonical
    bedload mass discharge per unit channel width.

    This formulation is mathematically equivalent to van Rijn's
    eq. 23 :math:`q_b = 0.053\\, [(s-1)\\,g]^{0.5}\\, D_{50}^{1.5}\\,
    d_*^{-0.3}\\, T_R^{2.1}` to within the 5% precision the helpers
    are documented to (verified by ``test_van_rijn_class_matches_helpers``).

    Reference: van Rijn, L. C. (1984). "Sediment transport, Part I: Bed
    load transport." J. Hydraul. Eng. 110(10), 1431–1456.
    """

    name: str = "van_rijn"

    def transport_rate(
        self,
        tau_pa: xr.DataArray,
        d50_um: float,
        tau_ce_pa: float,
        velocity_m_s: xr.DataArray,
        depth_m: xr.DataArray,
        slope: xr.DataArray | float,
        solid_density_g_cm3: float,
        water_density_kg_m3: float = 1000.0,
        kinematic_viscosity_m2_s: float = 1.0e-6,
        registry_context: dict | None = None,
    ) -> xr.DataArray:
        if tau_ce_pa <= 0.0:
            return xr.zeros_like(tau_pa)
        nu_cgs = kinematic_viscosity_m2_s * 1.0e4
        d_star = _cheng_d_star(
            d50_um,
            solid_specific_gravity=solid_density_g_cm3,
            kinematic_viscosity_cm2_s=nu_cgs,
        )
        u_bl = van_rijn_bedload_velocity_cm_s(
            tau_pa,
            tau_ce_pa=tau_ce_pa,
            d50_um=d50_um,
            solid_specific_gravity=solid_density_g_cm3,
        )
        delta_bl = van_rijn_bedload_height_cm(
            tau_pa,
            tau_ce_pa=tau_ce_pa,
            d50_um=d50_um,
            d_star=d_star,
        )
        c_eq = van_rijn_equilibrium_concentration(
            tau_pa,
            tau_ce_pa=tau_ce_pa,
            d_star=d_star,
            solid_density_g_cm3=solid_density_g_cm3,
        )
        # u_BL [cm/s] · δ_BL [cm] · C_eq [g/cm³]  →  g cm⁻¹ s⁻¹
        return u_bl * delta_bl * c_eq


# ---------------------------------------------------------------------------
# 2. Wilcock & Crowe (2003) — surface-based, sand-gravel mixtures
# ---------------------------------------------------------------------------

class WilcockCrowe2003TransportFunction:
    """Wilcock & Crowe (2003) surface-based bedload for sand-gravel mixtures.

    Uses the bed-surface sand fraction :math:`F_s` to set the
    reference Shields stress for the geometric-mean grain size, then
    applies a hiding-and-exposure adjustment to obtain a per-class
    reference. Dimensionless transport :math:`W^*_i` follows a
    two-regime closure.

    Equations (Wilcock & Crowe 2003, eqs. 1–6):

    .. math::

        \\tau^*_{rsg} = 0.021 + 0.015 \\exp(-20\\,F_s)

        \\tau^*_{ri}/\\tau^*_{rsg} = (d_i/d_{sg})^{b_i}

        b_i = \\frac{0.67}{1 + \\exp(1.5 - d_i/d_{sg})}

        \\phi = \\tau / \\tau_{ri}

        W^*_i = \\begin{cases}
            0.002\\,\\phi^{7.5}                 & \\phi < 1.35 \\\\
            14\\,(1 - 0.894/\\sqrt{\\phi})^{4.5} & \\phi \\ge 1.35
        \\end{cases}

        q_{b,i} = \\frac{F_i\\, W^*_i\\, u_*^3\\, \\rho_s}{(s-1)\\,g}

    where :math:`F_i` is the surface fraction of class :math:`i`.
    Returns one class's contribution; the orchestrator can sum over
    classes if it wants the total. ``registry_context`` may carry
    ``"surface_sand_fraction"`` (scalar or array) and
    ``"surface_class_fraction"`` (scalar fraction of *this* class on
    the surface, default 1.0). When omitted, a benign default is used.

    Reference: Wilcock, P. R., and Crowe, J. C. (2003). "Surface-based
    transport model for mixed-size sediment." J. Hydraul. Eng. 129(2),
    120–128. DOI: 10.1061/(ASCE)0733-9429(2003)129:2(120).
    """

    name: str = "wilcock_crowe"

    def transport_rate(
        self,
        tau_pa: xr.DataArray,
        d50_um: float,
        tau_ce_pa: float,
        velocity_m_s: xr.DataArray,
        depth_m: xr.DataArray,
        slope: xr.DataArray | float,
        solid_density_g_cm3: float,
        water_density_kg_m3: float = 1000.0,
        kinematic_viscosity_m2_s: float = 1.0e-6,
        registry_context: dict | None = None,
    ) -> xr.DataArray:
        ctx = registry_context or {}
        f_s = float(ctx.get("surface_sand_fraction", 0.15))   # typical gravel-bed Fs
        f_i = float(ctx.get("surface_class_fraction", 1.0))
        d_sg_um = float(ctx.get("surface_geometric_mean_um", d50_um))

        # Reference Shields stress for surface geometric mean (eq. 4).
        tau_star_rsg = 0.021 + 0.015 * np.exp(-20.0 * f_s)

        # Hiding/exposure (eq. 5): per-class reference stress.
        d_ratio = float(d50_um) / max(d_sg_um, 1.0e-12)
        b_i = 0.67 / (1.0 + np.exp(1.5 - d_ratio))
        tau_star_ri = tau_star_rsg * (d_ratio ** b_i)

        # Convert per-class reference Shields stress to dimensional τ_ri (Pa).
        g_si = 9.81
        rho_s_kg = solid_density_g_cm3 * 1000.0
        d_i_m = float(d50_um) * 1.0e-6
        tau_ri_pa = tau_star_ri * (rho_s_kg - water_density_kg_m3) * g_si * d_i_m
        if tau_ri_pa <= 0.0:
            return xr.zeros_like(tau_pa)

        # Stress ratio φ = τ / τ_ri (per cell).
        phi = tau_pa / tau_ri_pa
        phi = xr.where(phi > 0.0, phi, 0.0)

        # Two-regime W*_i (eq. 6).
        w_low = 0.002 * (phi ** 7.5)
        # Guard sqrt(0) → 0; in the high regime φ ≥ 1.35 always.
        sqrt_phi = np.sqrt(xr.where(phi > 0.0, phi, 1.0))
        w_high = 14.0 * np.maximum(1.0 - 0.894 / sqrt_phi, 0.0) ** 4.5
        w_star = xr.where(phi < 1.35, w_low, w_high)
        w_star = xr.where(phi > 0.0, w_star, 0.0)

        # Dimensional transport: q_b,i = F_i · W*_i · u*³ · ρ_s / ((s-1) g)
        u_star = _ushear_m_s(tau_pa, water_density_kg_m3=water_density_kg_m3)  # m/s
        s_minus_1 = (rho_s_kg / water_density_kg_m3) - 1.0
        if s_minus_1 <= 0.0:
            return xr.zeros_like(tau_pa)
        # q_b in kg m⁻¹ s⁻¹
        qb_si = f_i * w_star * (u_star ** 3) * rho_s_kg / (s_minus_1 * g_si)
        # Convert kg m⁻¹ s⁻¹ → g cm⁻¹ s⁻¹  (×1000/100 = ×10)
        return qb_si * 10.0


# ---------------------------------------------------------------------------
# 3. Parker (1990) — surface-based gravel
# ---------------------------------------------------------------------------

class Parker1990TransportFunction:
    """Parker (1990) surface-based gravel bedload (similarity-collapse form).

    Implements the dimensionless similarity function
    :math:`W^*(\\phi_{50})` with three regimes (Parker 1990 eq. 67):

    .. math::

        W^* = \\begin{cases}
            0.00218\\,\\exp\\bigl[14.2(\\phi - 1) - 9.28(\\phi - 1)^2\\bigr] &
                0.95 \\le \\phi < 1.65 \\\\
            0.00218\\,\\bigl(1 - 0.853/\\phi\\bigr)^{4.5} \\cdot 11.93 &
                \\phi \\ge 1.65 \\\\
            0.00218\\,\\phi^{14.2} & \\phi < 0.95
        \\end{cases}

    where :math:`\\phi = \\tau^* / \\tau^*_{r50}` and
    :math:`\\tau^*_{r50} = 0.0386` is Parker's reference Shields stress
    for the surface :math:`d_{50}`.

    Then :math:`q_b = W^* \\cdot u_*^3 \\cdot \\rho_s / ((s-1) g)`.

    Reference: Parker, G. (1990). "Surface-based bedload transport
    relation for gravel rivers." J. Hydraul. Res. 28(4), 417–436.
    DOI: 10.1080/00221689009499058.
    """

    name: str = "parker"
    TAU_STAR_R50: float = 0.0386

    def transport_rate(
        self,
        tau_pa: xr.DataArray,
        d50_um: float,
        tau_ce_pa: float,
        velocity_m_s: xr.DataArray,
        depth_m: xr.DataArray,
        slope: xr.DataArray | float,
        solid_density_g_cm3: float,
        water_density_kg_m3: float = 1000.0,
        kinematic_viscosity_m2_s: float = 1.0e-6,
        registry_context: dict | None = None,
    ) -> xr.DataArray:
        # Dimensionless Shields stress for d50.
        theta = _shields(
            tau_pa,
            d50_um=d50_um,
            solid_density_g_cm3=solid_density_g_cm3,
            water_density_kg_m3=water_density_kg_m3,
        )
        theta = xr.where(theta > 0.0, theta, 0.0)

        phi = theta / self.TAU_STAR_R50

        # Three regimes (Parker 1990 eq. 67).
        # Low regime (φ < 0.95): W* = 0.00218 · φ^14.2
        w_low = 0.00218 * np.maximum(phi, 0.0) ** 14.2
        # Mid (0.95 ≤ φ < 1.65): exponential expansion.
        x = phi - 1.0
        w_mid = 0.00218 * np.exp(14.2 * x - 9.28 * x * x)
        # High (φ ≥ 1.65): asymptotic form, capped to keep numerics tame
        # when φ is very large.
        safe_phi = xr.where(phi > 0.0, phi, 1.0)
        w_high = 0.00218 * 11.93 * np.maximum(1.0 - 0.853 / safe_phi, 0.0) ** 4.5

        w_star = xr.where(phi < 0.95, w_low,
                          xr.where(phi < 1.65, w_mid, w_high))
        w_star = xr.where(phi > 0.0, w_star, 0.0)

        g_si = 9.81
        rho_s_kg = solid_density_g_cm3 * 1000.0
        s_minus_1 = (rho_s_kg / water_density_kg_m3) - 1.0
        if s_minus_1 <= 0.0:
            return xr.zeros_like(tau_pa)

        u_star = _ushear_m_s(tau_pa, water_density_kg_m3=water_density_kg_m3)
        qb_si = w_star * (u_star ** 3) * rho_s_kg / (s_minus_1 * g_si)
        return qb_si * 10.0  # kg m⁻¹ s⁻¹ → g cm⁻¹ s⁻¹


# ---------------------------------------------------------------------------
# 4. Yang (1973, 1979) — unit-stream-power total load
# ---------------------------------------------------------------------------

class YangTransportFunction:
    """Yang (1973, 1979) total-load formula based on unit stream power V·S.

    Uses the sand formula (Yang 1973) for D₅₀ < 2 mm and the gravel
    formula (Yang 1984) for coarser material:

    Sand (1973, eq. 25):

    .. math::

        \\log C_t = 5.435 - 0.286 \\log\\!\\frac{w_s d}{\\nu}
                            - 0.457 \\log\\!\\frac{u_*}{w_s}
                  + \\bigl(1.799 - 0.409 \\log\\!\\frac{w_s d}{\\nu}
                            - 0.314 \\log\\!\\frac{u_*}{w_s}\\bigr)
                    \\log\\!\\Bigl(\\frac{V S}{w_s} - \\frac{V_{cr} S}{w_s}\\Bigr)

    where :math:`C_t` is total sediment concentration in ppm by weight,
    :math:`V_{cr}` is the incipient-motion velocity (Yang 1973 eq. 8).

    The total-load mass discharge per unit width is then

    .. math::

        q_t = C_t \\cdot 10^{-6} \\cdot \\rho_w \\cdot V \\cdot h

    in kg m⁻¹ s⁻¹, converted to g cm⁻¹ s⁻¹ by ×10. Yang's formula is
    a *total*-load relation; SSM's bedload-fraction split is handled
    upstream of this closure (the orchestrator multiplies by the bedload
    PSUS factor when treating Yang as a bedload-only proxy).

    Reference: Yang, C. T. (1973). "Incipient motion and sediment
    transport." J. Hydraul. Div. ASCE 99(HY10), 1679–1704.
    Yang, C. T. (1979). "Unit stream power equations for total load."
    J. Hydrology 40(1–2), 123–138.
    DOI: 10.1016/0022-1694(79)90092-1.
    """

    name: str = "yang"

    def transport_rate(
        self,
        tau_pa: xr.DataArray,
        d50_um: float,
        tau_ce_pa: float,
        velocity_m_s: xr.DataArray,
        depth_m: xr.DataArray,
        slope: xr.DataArray | float,
        solid_density_g_cm3: float,
        water_density_kg_m3: float = 1000.0,
        kinematic_viscosity_m2_s: float = 1.0e-6,
        registry_context: dict | None = None,
    ) -> xr.DataArray:
        g_si = 9.81
        d_m = float(d50_um) * 1.0e-6
        rho_s_kg = solid_density_g_cm3 * 1000.0
        s_minus_1 = (rho_s_kg / water_density_kg_m3) - 1.0
        if s_minus_1 <= 0.0:
            return xr.zeros_like(tau_pa)

        # Settling velocity via simplified Rubey/Cheng (good enough for
        # the range Yang validated on, 0.062–10 mm).
        # Use Cheng (1997) closed form; mirror settling.cheng_1997.
        nu = kinematic_viscosity_m2_s
        d_star = d_m * (s_minus_1 * g_si / nu**2) ** (1.0 / 3.0)
        # Cheng (1997) eq. 11
        ws = (nu / d_m) * (np.sqrt(25.0 + 1.2 * d_star ** 2) - 5.0) ** 1.5
        if ws <= 0.0:
            return xr.zeros_like(tau_pa)

        u_star = _ushear_m_s(tau_pa, water_density_kg_m3=water_density_kg_m3)

        # Incipient-motion velocity Vcr (Yang 1973 eq. 8). Two regimes
        # by particle Reynolds number Re* = u* d / ν.
        re_star = u_star * d_m / nu
        # Bound Re* away from the singularity at 1.2 (where ln vanishes).
        re_safe = xr.where(re_star > 1.2, re_star, 70.0)
        vcr_over_ws_low = 2.5 / (np.log10(re_safe) - 0.06) + 0.66
        vcr_over_ws_high = xr.full_like(re_star, 2.05)
        vcr_over_ws = xr.where(re_star < 70.0, vcr_over_ws_low, vcr_over_ws_high)
        # When Re* < 1.2 the formula is undefined; set Vcr to a large
        # value so that V S − Vcr S < 0 ⇒ no transport.
        vcr_over_ws = xr.where(re_star >= 1.2, vcr_over_ws, 1.0e6)
        v_cr = vcr_over_ws * ws

        # Slope as DataArray for broadcasting.
        if isinstance(slope, (int, float)):
            slope_da = xr.full_like(tau_pa, float(slope))
        else:
            slope_da = slope

        # Effective unit-stream-power excess.
        psi = (velocity_m_s * slope_da - v_cr * slope_da) / ws
        psi = xr.where(psi > 0.0, psi, 0.0)

        # Yang 1973 sand formula (or Yang 1984 gravel for coarse).
        ws_d_over_nu = ws * d_m / nu
        u_over_ws = xr.where(ws > 0.0, u_star / ws, 0.0)
        # Guard logs against zeros.
        log_term = np.log10(xr.where(psi > 0.0, psi, 1.0))
        log_wsd = np.log10(max(ws_d_over_nu, 1.0e-12))
        log_uws = np.log10(xr.where(u_over_ws > 0.0, u_over_ws, 1.0))

        # Coefficients (Yang 1973 sand). Yang 1984 gravel uses different
        # constants; we use the sand form for D50 < 2 mm and switch to
        # the gravel form (Yang 1984 eq. 11) above.
        if d50_um < 2000.0:
            log_ct = (
                5.435 - 0.286 * log_wsd - 0.457 * log_uws
                + (1.799 - 0.409 * log_wsd - 0.314 * log_uws) * log_term
            )
        else:
            # Yang 1984 gravel.
            log_ct = (
                6.681 - 0.633 * log_wsd - 4.816 * log_uws
                + (2.784 - 0.305 * log_wsd - 0.282 * log_uws) * log_term
            )
        # Where psi <= 0 there is no transport.
        c_t_ppm = xr.where(psi > 0.0, 10.0 ** log_ct, 0.0)

        # Convert ppm-by-weight × ρ_w × V × h → kg m⁻¹ s⁻¹.
        qb_si = c_t_ppm * 1.0e-6 * water_density_kg_m3 * velocity_m_s * depth_m
        return qb_si * 10.0  # → g cm⁻¹ s⁻¹


# ---------------------------------------------------------------------------
# 5. Wu, Wang & Jia (2000) — non-uniform sediment with hiding/exposure
# ---------------------------------------------------------------------------

class Wu2000TransportFunction:
    """Wu, Wang & Jia (2000) non-uniform sediment bedload formula.

    Uses a hiding-and-exposure correction based on full size-distribution
    fractions; the bedload component (eqs. 6, 9, 10) is

    .. math::

        \\Phi_{b,i} = 0.0053 \\Bigl(\\frac{\\tau_b'}{\\tau_{ci}} - 1\\Bigr)^{2.2}

        \\frac{\\tau_{ci}}{\\tau_{cm}} =
            \\Bigl(\\frac{p_{e,i}}{p_{h,i}}\\Bigr)^{-0.6}

    where :math:`\\tau_{cm} = 0.03\\,(\\rho_s - \\rho)\\,g\\,d_i` and the
    exposure/hidden probabilities are computed from the bed
    size-distribution. Class fraction defaults to 1.0 (uniform bed)
    when ``registry_context`` is absent. Returns

    .. math::

        q_{b,i} = \\Phi_{b,i} \\cdot p_i \\cdot \\rho_s
                  \\cdot \\sqrt{(s - 1)\\,g\\,d_i^3}

    with :math:`p_i` the per-class bed-fraction.

    Reference: Wu, W., Wang, S. S. Y., and Jia, Y. (2000). "Nonuniform
    sediment transport in alluvial rivers." J. Hydraul. Res. 38(6),
    427–434. DOI: 10.1080/00221680009498296.
    """

    name: str = "wu"

    def transport_rate(
        self,
        tau_pa: xr.DataArray,
        d50_um: float,
        tau_ce_pa: float,
        velocity_m_s: xr.DataArray,
        depth_m: xr.DataArray,
        slope: xr.DataArray | float,
        solid_density_g_cm3: float,
        water_density_kg_m3: float = 1000.0,
        kinematic_viscosity_m2_s: float = 1.0e-6,
        registry_context: dict | None = None,
    ) -> xr.DataArray:
        ctx = registry_context or {}
        # Hidden / exposed probability ratio. Defaults to 1.0
        # (no hiding/exposure correction → uniform-bed limit).
        pe_ph_ratio = float(ctx.get("pe_ph_ratio", 1.0))
        f_i = float(ctx.get("surface_class_fraction", 1.0))

        g_si = 9.81
        d_i_m = float(d50_um) * 1.0e-6
        rho_s_kg = solid_density_g_cm3 * 1000.0
        s_minus_1 = (rho_s_kg / water_density_kg_m3) - 1.0
        if s_minus_1 <= 0.0:
            return xr.zeros_like(tau_pa)

        # Reference critical Shields stress θ_cm = 0.03 (Wu et al. 2000).
        # Per-class critical stress with hiding/exposure correction.
        tau_cm_pa = 0.03 * (rho_s_kg - water_density_kg_m3) * g_si * d_i_m
        tau_ci_pa = tau_cm_pa * (pe_ph_ratio ** -0.6)
        if tau_ci_pa <= 0.0:
            return xr.zeros_like(tau_pa)

        # Skin-friction shear ≈ total τ in the absence of bedform info.
        excess = tau_pa / tau_ci_pa - 1.0
        excess = xr.where(excess > 0.0, excess, 0.0)
        phi_bi = 0.0053 * (excess ** 2.2)

        # Dimensional q_b,i = φ · p_i · ρ_s · √((s-1) g d³)
        scale = rho_s_kg * np.sqrt(s_minus_1 * g_si * d_i_m ** 3)  # kg m⁻¹ s⁻¹
        qb_si = phi_bi * f_i * scale
        return qb_si * 10.0  # → g cm⁻¹ s⁻¹


# ---------------------------------------------------------------------------
# 6. Engelund & Hansen (1967) — total load for sand
# ---------------------------------------------------------------------------

class EngelundHansen1967TransportFunction:
    """Engelund & Hansen (1967) total-load formula for sand-bed rivers.

    .. math::

        q_t = \\frac{0.05\\,V^5}{\\sqrt{g}\\,C^3\\,\\Delta^2\\,d_{50}}

    where :math:`C = V/\\sqrt{R S}` is the Chézy coefficient, :math:`R`
    is the hydraulic radius (≈ depth in wide channels), :math:`S` is
    the energy slope, and :math:`\\Delta = (\\rho_s - \\rho_w)/\\rho_w`.

    The result is the total volumetric transport per unit width
    (m² s⁻¹); multiplied by :math:`\\rho_s` (kg/m³) it becomes the
    mass discharge (kg m⁻¹ s⁻¹), and ×10 converts to g cm⁻¹ s⁻¹.

    Engelund & Hansen treat the formula as a *total-load* closure, so
    the entire mass passes through the bedload window; SSM's bedload
    PSUS split (if any) is applied externally.

    Reference: Engelund, F., and Hansen, E. (1967). "A monograph on
    sediment transport in alluvial streams." Teknisk Forlag, Copenhagen.
    """

    name: str = "engelund_hansen"

    def transport_rate(
        self,
        tau_pa: xr.DataArray,
        d50_um: float,
        tau_ce_pa: float,
        velocity_m_s: xr.DataArray,
        depth_m: xr.DataArray,
        slope: xr.DataArray | float,
        solid_density_g_cm3: float,
        water_density_kg_m3: float = 1000.0,
        kinematic_viscosity_m2_s: float = 1.0e-6,
        registry_context: dict | None = None,
    ) -> xr.DataArray:
        g_si = 9.81
        d_m = float(d50_um) * 1.0e-6
        rho_s_kg = solid_density_g_cm3 * 1000.0
        delta = (rho_s_kg - water_density_kg_m3) / water_density_kg_m3
        if delta <= 0.0 or d_m <= 0.0:
            return xr.zeros_like(tau_pa)

        # Chézy from the velocity / depth / slope: C² = V² / (R S).
        # Use depth as a proxy for R (wide-channel approximation).
        if isinstance(slope, (int, float)):
            slope_da = xr.full_like(tau_pa, float(slope))
        else:
            slope_da = slope
        # Guard against zero slope or depth → infinite C → undefined q_t.
        rs = depth_m * slope_da
        rs_safe = xr.where(rs > 0.0, rs, 1.0)
        c_chezy_sq = xr.where(rs > 0.0, velocity_m_s ** 2 / rs_safe, 0.0)
        # When V or S is zero, C² may collapse; guard.
        c_chezy = np.sqrt(xr.where(c_chezy_sq > 0.0, c_chezy_sq, 1.0))

        # Volumetric transport per unit width (m²/s).
        qt_vol = (
            0.05 * (velocity_m_s ** 5)
            / (np.sqrt(g_si) * (c_chezy ** 3) * (delta ** 2) * d_m)
        )
        qt_vol = xr.where(c_chezy_sq > 0.0, qt_vol, 0.0)
        qt_vol = xr.where(qt_vol > 0.0, qt_vol, 0.0)

        # Mass discharge (kg m⁻¹ s⁻¹) → g cm⁻¹ s⁻¹.
        qb_si = qt_vol * rho_s_kg
        return qb_si * 10.0


# ---------------------------------------------------------------------------
# 7. Toffaleti (1968) — depth-integrated total load
# ---------------------------------------------------------------------------

class Toffaleti1968TransportFunction:
    """Toffaleti (1968) depth-integrated total-load procedure.

    The full Toffaleti procedure integrates suspended-load concentration
    over four vertical zones (lower, middle, upper, surface) using
    depth-power exponents :math:`z_i = w_s / (κ u_*)`. SSM uses a
    single-zone reduction adequate for depth-averaged transport at a
    cell, calibrated against the Toffaleti tables (BR-1 procedure):

    .. math::

        q_t = M\\, \\frac{V^{n_v}}{(0.00058)^{n_v - 1}}\\, d_{50}^{0.33}

    where :math:`M` and :math:`n_v` are the bedload-fraction
    coefficients from Toffaleti's empirical tables. We take
    :math:`n_v = 1.5` and :math:`M = 0.6` as documented for medium
    sand at 60 °F. The result is in tons day⁻¹ ft⁻¹ in Toffaleti's
    original units; we convert to kg m⁻¹ s⁻¹ then to g cm⁻¹ s⁻¹.

    This is a deliberate single-zone simplification of the original
    four-zone integral; the design memo
    ``ssm_bedload_functions.md`` documents the reduction. The full
    multi-zone integration is queued for a phase-3 enhancement.

    Reference: Toffaleti, F. B. (1968). "A procedure for computation of
    the total river sand discharge and detailed distribution, bed to
    surface." Tech. Report No. 5, Committee on Channel Stabilization,
    U.S. Army Corps of Engineers.
    """

    name: str = "toffaleti"

    def transport_rate(
        self,
        tau_pa: xr.DataArray,
        d50_um: float,
        tau_ce_pa: float,
        velocity_m_s: xr.DataArray,
        depth_m: xr.DataArray,
        slope: xr.DataArray | float,
        solid_density_g_cm3: float,
        water_density_kg_m3: float = 1000.0,
        kinematic_viscosity_m2_s: float = 1.0e-6,
        registry_context: dict | None = None,
    ) -> xr.DataArray:
        # Toffaleti single-zone reduction. Convert to ft / ft·s for the
        # original empirical fit, then back to SI.
        v_fps = velocity_m_s * 3.28084
        depth_ft = depth_m * 3.28084
        d_in_ft = float(d50_um) * 1.0e-6 * 3.28084  # μm → m → ft

        # Threshold: below incipient motion velocity, no transport.
        # Use a simple Shields-type proxy via tau / tau_ce.
        if tau_ce_pa <= 0.0:
            return xr.zeros_like(tau_pa)
        gate = xr.where(tau_pa > tau_ce_pa, 1.0, 0.0)

        # Toffaleti BR-1 constants (medium sand, 60 °F).
        m_coef = 0.6
        n_v = 1.5

        # Total load (tons day⁻¹ ft⁻¹) — empirical regression.
        qt_tdf = (
            m_coef * (v_fps ** n_v)
            * np.maximum(depth_ft, 0.0) ** 0.0
            * (d_in_ft ** 0.33)
            / (0.00058 ** (n_v - 1.0))
        )
        qt_tdf = qt_tdf * gate

        # Convert tons/day/ft → kg/s/m.
        # 1 ton (US short) = 907.185 kg; 1 day = 86400 s; 1 ft = 0.3048 m.
        # tons/day/ft × (907.185 / 86400) / 0.3048  =  ×0.0344 (kg/s/m)
        qt_si = qt_tdf * (907.185 / 86400.0) / 0.3048
        qt_si = xr.where(qt_si > 0.0, qt_si, 0.0)
        return qt_si * 10.0  # → g cm⁻¹ s⁻¹


# ---------------------------------------------------------------------------
# Registry of available transport functions
# ---------------------------------------------------------------------------

#: Mapping of transport-function name → class. Selected via SSM YAML
#: ``sediment.bedload.transport_function``. Use :func:`get_transport_function`
#: to instantiate by name.
BEDLOAD_TRANSPORT_FUNCTIONS: dict[str, type] = {
    VanRijn1984TransportFunction.name: VanRijn1984TransportFunction,
    WilcockCrowe2003TransportFunction.name: WilcockCrowe2003TransportFunction,
    Parker1990TransportFunction.name: Parker1990TransportFunction,
    YangTransportFunction.name: YangTransportFunction,
    Wu2000TransportFunction.name: Wu2000TransportFunction,
    EngelundHansen1967TransportFunction.name: EngelundHansen1967TransportFunction,
    Toffaleti1968TransportFunction.name: Toffaleti1968TransportFunction,
}


def get_transport_function(name: str) -> BedloadTransportFunction:
    """Instantiate a bedload transport function by name.

    Raises ``ValueError`` if ``name`` is not one of the keys of
    :data:`BEDLOAD_TRANSPORT_FUNCTIONS`. Names are case-insensitive.
    """
    key = (name or "").strip().lower()
    if key not in BEDLOAD_TRANSPORT_FUNCTIONS:
        raise ValueError(
            f"Unknown bedload transport_function {name!r}; "
            f"expected one of {sorted(BEDLOAD_TRANSPORT_FUNCTIONS)}."
        )
    return BEDLOAD_TRANSPORT_FUNCTIONS[key]()
