"""Bed shear-stress drivers.

Three implementations behind a common Protocol:

* :class:`ExternalShearDriver` — read τ_b from a mesh variable supplied
  by the user / wave model / measurement.
* :class:`CurrentOnlyShearDriver` — compute from RAS face velocities via
  Parker (2004) log-law (default), or alternatively from Manning's n.
* :class:`WaveCurrentShearDriver` — Christoffersen & Jonsson (1985)
  combined wave-current iteration. **Stub in v1; full impl deferred.**

Reference: SAND2008-5621 §"S_SHEAR.f90"; design spec §5.2 and §8;
Parker (2004); Christoffersen & Jonsson (1985); EFDC s_shear.f90.

Implementation note (license/provenance): the Parker log-law and
Christoffersen-Jonsson formulations are reproduced directly from the
peer-reviewed literature. The EFDC ``s_shear.f90`` source is treated as
a verification reference only; no Fortran code has been ported.
"""

from __future__ import annotations

from datetime import datetime
from typing import Protocol, runtime_checkable

import numpy as np
import xarray as xr

from . import contracts, coupling


# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

#: Water density used for τ = ρ f_c U² in SI Pa.
_RHO_W_SI: float = 1000.0  # kg/m^3

#: Gravitational acceleration in SI for the Manning's-n formulation.
_G_SI: float = 9.81  # m/s^2

#: Floor on hydraulic depth (m) used to guard the log-law denominator and
#: any 1/h^(1/3) factor against vanishing dry cells.
_DEPTH_FLOOR_M: float = 1.0e-6

#: Floor on roughness k_n (m) used to keep ln(11h/(2k_n)) well-defined.
_KN_FLOOR_M: float = 1.0e-6


# ---------------------------------------------------------------------------
# Protocol
# ---------------------------------------------------------------------------


@runtime_checkable
class ShearStressDriver(Protocol):
    """Returns bed shear stress in Pa on (nface,) for a single time step."""

    def compute(
        self,
        mesh: xr.Dataset,
        time: datetime,
        d50_surface_um: xr.DataArray,
        previous_tau_pa: xr.DataArray,
    ) -> xr.DataArray:
        """Compute τ_b in Pa.

        Parameters
        ----------
        mesh : xr.Dataset
            Shared mesh dataset with hydraulics variables already populated.
        time : datetime
            Time stamp of the step being computed.
        d50_surface_um : xr.DataArray  (nface,)
            Surface-layer mean D50 in μm. Used as the skin-friction roughness
            fallback when ``zb_skin`` is not configured.
        previous_tau_pa : xr.DataArray  (nface,)
            τ_b from the previous step. Used by the growth limiter.

        Returns
        -------
        xr.DataArray
            (nface,) bed shear stress in Pa.
        """
        ...


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def apply_growth_limiter(
    tau_new: xr.DataArray,
    previous_tau_pa: xr.DataArray,
    growth_limit: float,
) -> xr.DataArray:
    """Limit the per-step growth of τ_b to ``growth_limit`` (fraction).

    SEDZLJ-style stability device (mirrors ``s_shear.f90:315``). When
    ``τ_new`` exceeds ``τ_prev × (1 + growth_limit)``, the rise is clamped
    to ``growth_limit × (τ_new − τ_prev)`` above ``τ_prev``. Drops in τ
    are passed through unchanged. Pure xarray; broadcasts naturally.

    Parameters
    ----------
    tau_new : xr.DataArray
        Newly computed τ_b (Pa).
    previous_tau_pa : xr.DataArray
        τ_b from the previous step (Pa). Same dims/shape as ``tau_new``.
    growth_limit : float
        Maximum allowed fractional increase per step (e.g. ``0.10`` for 10 %).

    Returns
    -------
    xr.DataArray
        Limited τ_b in Pa.
    """
    if growth_limit <= 0.0:
        return tau_new
    threshold = previous_tau_pa * (1.0 + growth_limit)
    clamped = previous_tau_pa + growth_limit * (tau_new - previous_tau_pa)
    return xr.where(tau_new > threshold, clamped, tau_new)


def _edge_velocity_to_face_magnitude(mesh: xr.Dataset, time) -> xr.DataArray:
    """Reconstruct cell-centroid velocity magnitude from edge velocities.

    Strategy (preference order):

    1. If RAS exposed cell-velocity components on faces
       (``face_velocity_x``/``face_velocity_y``), use them directly —
       they already live on the cell centroid and are by far the most
       accurate source.
    2. If a precomputed cell-velocity-magnitude field is present
       (``face_velocity_magnitude``), use it.
    3. Otherwise fall back to averaging the *absolute* values of the
       edge-normal velocities incident on each face. This is a coarse
       reconstruction (it does not recover the vector velocity), but it
       is the same approach used by EFDC SEDZLJ when face vector velocity
       is not available, and it preserves shear scaling
       :math:`\\tau \\propto |U|^2`.

    Returns an :class:`xr.DataArray` on ``(nface,)`` with values in m/s.
    """

    # Path 1 — face vector components (best)
    if "face_velocity_x" in mesh.data_vars and "face_velocity_y" in mesh.data_vars:
        ux = mesh["face_velocity_x"]
        uy = mesh["face_velocity_y"]
        if contracts.DIM_TIME in ux.dims:
            ux = ux.sel({contracts.DIM_TIME: time})
        if contracts.DIM_TIME in uy.dims:
            uy = uy.sel({contracts.DIM_TIME: time})
        return np.sqrt(ux * ux + uy * uy)

    # Path 2 — face velocity magnitude already on the mesh
    if "face_velocity_magnitude" in mesh.data_vars:
        umag = mesh["face_velocity_magnitude"]
        if contracts.DIM_TIME in umag.dims:
            umag = umag.sel({contracts.DIM_TIME: time})
        return umag

    # Path 3 — average |edge_velocity| onto each face using the
    # edges_face1/edges_face2 connectivity. We accumulate Σ|U_e| and a
    # count per face, then divide. This is the same edge→face averaging
    # pattern used by clearwater_riverine.utilities for diffusion.
    if contracts.VAR_EDGE_VELOCITY not in mesh.data_vars:
        raise ValueError(
            "CurrentOnlyShearDriver requires either "
            "(face_velocity_x, face_velocity_y), face_velocity_magnitude, or "
            f"{contracts.VAR_EDGE_VELOCITY!r} on the mesh dataset."
        )
    if "edges_face1" not in mesh.variables or "edges_face2" not in mesh.variables:
        raise ValueError(
            "Edge-to-face velocity averaging requires 'edges_face1' and "
            "'edges_face2' connectivity arrays on the mesh dataset."
        )

    edge_v = mesh[contracts.VAR_EDGE_VELOCITY]
    if contracts.DIM_TIME in edge_v.dims:
        edge_v = edge_v.sel({contracts.DIM_TIME: time})
    edge_speed = np.abs(edge_v.values).astype("float64")  # (nedge,)

    face1 = np.asarray(mesh["edges_face1"].values, dtype=np.int64)
    face2 = np.asarray(mesh["edges_face2"].values, dtype=np.int64)

    nface = mesh.sizes[contracts.DIM_NFACE]
    speed_sum = np.zeros(nface, dtype="float64")
    count = np.zeros(nface, dtype="float64")

    # face1 / face2 may include ghost-cell indices >= nface; mask those.
    f1_valid = (face1 >= 0) & (face1 < nface)
    f2_valid = (face2 >= 0) & (face2 < nface)
    np.add.at(speed_sum, face1[f1_valid], edge_speed[f1_valid])
    np.add.at(count, face1[f1_valid], 1.0)
    np.add.at(speed_sum, face2[f2_valid], edge_speed[f2_valid])
    np.add.at(count, face2[f2_valid], 1.0)

    safe_count = np.where(count > 0, count, 1.0)
    face_speed = speed_sum / safe_count
    return xr.DataArray(
        face_speed,
        dims=(contracts.DIM_NFACE,),
        name="face_velocity_magnitude",
    )


def _face_depth(mesh: xr.Dataset, time) -> xr.DataArray:
    """Return hydraulic depth on (nface,) at the requested time, floored
    at :data:`_DEPTH_FLOOR_M` to keep log-law and 1/h^(1/3) finite for
    dry cells."""
    if contracts.VAR_FACE_HYDRAULIC_DEPTH not in mesh.data_vars:
        raise ValueError(
            "CurrentOnlyShearDriver requires "
            f"{contracts.VAR_FACE_HYDRAULIC_DEPTH!r} on the mesh dataset."
        )
    depth = mesh[contracts.VAR_FACE_HYDRAULIC_DEPTH]
    if contracts.DIM_TIME in depth.dims:
        depth = depth.sel({contracts.DIM_TIME: time})
    return xr.where(depth > _DEPTH_FLOOR_M, depth, _DEPTH_FLOOR_M)


def _roughness_kn(d50_surface_um: xr.DataArray, zb_skin_m: float) -> xr.DataArray:
    """Skin-friction roughness ``k_n = max(D50_surface (m), zb_skin_m)`` (m).

    ``d50_surface_um`` arrives in micrometres (μm) per the contracts
    schema; convert to metres before the max. A small floor is applied
    so the log-law denominator can never blow up.
    """
    d50_m = d50_surface_um * 1.0e-6
    floored = xr.where(d50_m > zb_skin_m, d50_m, zb_skin_m)
    return xr.where(floored > _KN_FLOOR_M, floored, _KN_FLOOR_M)


# ---------------------------------------------------------------------------
# Driver implementations
# ---------------------------------------------------------------------------


class ExternalShearDriver:
    """Mode A: read τ_b directly from :data:`contracts.VAR_BED_SHEAR_STRESS_INPUT`.

    Useful when τ_b has been precomputed by an external wave or
    coupled-model tool and stored on the mesh.

    Parameters
    ----------
    growth_limit : float
        Per-step fractional cap on τ growth (default 10 %; design spec §5.2).
        Set to ``0.0`` to disable.
    """

    def __init__(
        self,
        growth_limit: float = contracts.DEFAULT_SHEAR_GROWTH_LIMIT,
    ) -> None:
        self.growth_limit = float(growth_limit)

    def compute(
        self,
        mesh: xr.Dataset,
        time: datetime,
        d50_surface_um: xr.DataArray,
        previous_tau_pa: xr.DataArray,
    ) -> xr.DataArray:
        if contracts.VAR_BED_SHEAR_STRESS_INPUT not in mesh.data_vars:
            raise ValueError(
                "ExternalShearDriver requires "
                f"{contracts.VAR_BED_SHEAR_STRESS_INPUT!r} on the mesh dataset. "
                "Either populate this field from your external τ source, or "
                "switch the configured shear driver to 'current_only'."
            )
        tau = mesh[contracts.VAR_BED_SHEAR_STRESS_INPUT]
        if contracts.DIM_TIME in tau.dims:
            tau = tau.sel({contracts.DIM_TIME: time})
        tau = tau.astype("float64")
        return apply_growth_limiter(tau, previous_tau_pa, self.growth_limit)


class CurrentOnlyShearDriver:
    """Mode B (default): compute τ_b from face velocities.

    Two sub-formulations selectable via ``formulation``:

    * ``"log_law"`` (default) — Parker (2004) log-law:
      :math:`f_c = (0.42 / \\ln(11 h / (2 k_n)))^2`, then
      :math:`\\tau = \\rho_w f_c |U|^2`.
    * ``"manning"`` — :math:`f_c = g n^2 / h^{1/3}`. Useful when ESM
      supplies a calibrated composite Manning's n.

    Applies a per-step growth limiter (default 10 %) for stability.

    Parameters
    ----------
    formulation : {"log_law", "manning"}
        Shear formulation. ``"log_law"`` ignores Manning's n; ``"manning"``
        ignores the skin-friction k_n.
    zb_skin_m : float
        Skin-friction roughness fallback (m); used in ``log_law`` whenever
        the surface-layer D50 is smaller than this value.
    growth_limit : float
        Per-step fractional cap on τ growth (default 10 %).
    use_composite_manning : bool
        If True (and ``formulation="manning"``), prefer ESM's
        :data:`contracts.VAR_COMPOSITE_MANNINGS_N` over the static RAS
        Manning's n when both are present.
    """

    _ALLOWED_FORMULATIONS = ("log_law", "manning")

    def __init__(
        self,
        formulation: str = "log_law",
        zb_skin_m: float = contracts.DEFAULT_ZB_SKIN_M,
        growth_limit: float = contracts.DEFAULT_SHEAR_GROWTH_LIMIT,
        use_composite_manning: bool = True,
    ) -> None:
        if formulation not in self._ALLOWED_FORMULATIONS:
            raise ValueError(
                f"formulation must be one of {self._ALLOWED_FORMULATIONS}, "
                f"got {formulation!r}"
            )
        self.formulation = formulation
        self.zb_skin_m = float(zb_skin_m)
        self.growth_limit = float(growth_limit)
        self.use_composite_manning = bool(use_composite_manning)

    # -- helpers --------------------------------------------------------

    def _f_c_log_law(
        self,
        depth_m: xr.DataArray,
        d50_surface_um: xr.DataArray,
    ) -> xr.DataArray:
        """Parker (2004) log-law friction factor (dimensionless).

        :math:`f_c = (0.42 / \\ln(11 h / (2 k_n)))^2`, with
        :math:`k_n = \\max(D_{50}\\,[\\mathrm{m}], z_{b,\\mathrm{skin}})`.
        """
        kn = _roughness_kn(d50_surface_um, self.zb_skin_m)
        # Guard the argument of the logarithm to be > 1 so ln > 0 and
        # the friction factor remains well-defined for very shallow
        # cells (h ~ k_n). The 1.0e-3 floor keeps the divisor bounded.
        ln_arg = (11.0 * depth_m) / (2.0 * kn)
        ln_arg = xr.where(ln_arg > 1.0 + 1.0e-9, ln_arg, 1.0 + 1.0e-9)
        denom = np.log(ln_arg)
        return (0.42 / denom) ** 2

    def _f_c_manning(
        self,
        depth_m: xr.DataArray,
        mannings_n: xr.DataArray,
    ) -> xr.DataArray:
        """Manning friction factor :math:`f_c = g n^2 / h^{1/3}`."""
        return _G_SI * (mannings_n ** 2) / (depth_m ** (1.0 / 3.0))

    def _resolve_mannings_n(
        self,
        mesh: xr.Dataset,
        time,
    ) -> xr.DataArray:
        """Return Manning's n on (nface,), preferring ESM composite if enabled."""
        if contracts.VAR_MANNINGS_N not in mesh.data_vars:
            raise ValueError(
                "CurrentOnlyShearDriver(formulation='manning') requires "
                f"{contracts.VAR_MANNINGS_N!r} on the mesh dataset."
            )
        static_n = mesh[contracts.VAR_MANNINGS_N]
        if self.use_composite_manning:
            return coupling.read_composite_manning_n(mesh, time, static_n)
        return static_n

    # -- main entry point ----------------------------------------------

    def compute(
        self,
        mesh: xr.Dataset,
        time: datetime,
        d50_surface_um: xr.DataArray,
        previous_tau_pa: xr.DataArray,
    ) -> xr.DataArray:
        depth_m = _face_depth(mesh, time)
        u_mag = _edge_velocity_to_face_magnitude(mesh, time)

        if self.formulation == "log_law":
            f_c = self._f_c_log_law(depth_m, d50_surface_um)
        else:  # "manning" (validated in __init__)
            mannings_n = self._resolve_mannings_n(mesh, time)
            f_c = self._f_c_manning(depth_m, mannings_n)

        tau_new = _RHO_W_SI * f_c * (u_mag ** 2)
        # Strip auxiliary names/coords so downstream code can broadcast freely.
        tau_new = tau_new.rename(None) if tau_new.name is not None else tau_new
        return apply_growth_limiter(tau_new, previous_tau_pa, self.growth_limit)


class WaveCurrentShearDriver:
    """Mode C (stub): combined wave-current shear via Christoffersen & Jonsson (1985).

    Full implementation deferred — see design spec §15.
    """

    def __init__(self, **kwargs) -> None:
        # Constructor is permitted (so YAML config can still resolve the
        # class), but any attempt to actually compute τ raises.
        self._kwargs = kwargs

    def compute(
        self,
        mesh: xr.Dataset,
        time: datetime,
        d50_surface_um: xr.DataArray,
        previous_tau_pa: xr.DataArray,
    ) -> xr.DataArray:
        raise NotImplementedError(
            "Mode C wave-current shear driver deferred to phase 5; "
            "see design spec §15"
        )
