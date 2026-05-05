"""Erosion-rate models.

Two formulations behind a common Protocol, selected by
``nsedflume`` config flag:

* :class:`SedflumeTableErosionModel` (``nsedflume=1``) — bilinear
  interpolation in (τ, fractional layer-mass remaining) of the
  per-core ERATE table. Log-linear in depth, linear in τ. Mirrors
  s_sedzlj.f90:498 and :535.

* :class:`PowerLawErosionModel` (``nsedflume=2``) — per-core, per-layer
  power law :math:`E = A \\tau^n` with cap :math:`E_{\\max}`. Mirrors
  s_sedzlj.f90:508 and :539.

Vegetation-cohesion feedback (`tau_ce_eff`) is applied at the class
gate: per-class erosion is suppressed when τ < τ_ce × (1 + α B), where
B is the vegetation biostabilization factor.

Reference: SAND2008-5621 §"S_SEDZLJ.f90"; design spec §5.4.

Implementation note: this module is a clean-room re-implementation from
the published equations in SAND2008-5621 (Thanh, Grace & James 2008).
The EFDC+ Fortran source (GPL-2.0) was consulted only as a behavioural
reference; no source lines were copied. The variable names
(SN00/SN10/SN01/SN11, ERATEMOD, etc.) follow the SAND manual's
notation, which is not copyrightable.
"""

from __future__ import annotations

import warnings
from typing import Optional, Protocol, runtime_checkable

import numpy as np
import xarray as xr

from . import contracts  # noqa: F401  (re-exported for downstream submodules)
from .consolidation import ConsolidationModel


# Floor value used when the layer below the deepest in-place layer would be
# referenced; matches SEDZLJ behaviour to prevent through-bottom erosion.
_DEEPEST_LAYER_ERATE_FLOOR_CM_S: float = 1.0e-9


@runtime_checkable
class ErosionRateModel(Protocol):
    """Returns per-cell erosion rate in g/cm²/s for a single layer."""

    def erosion_rate(
        self,
        tau_pa: xr.DataArray,                  # (nface,)
        layer_index: int,                      # K index, 1-origin top-down
        layer_mass: xr.DataArray,              # (nface,) g/cm^2  TSED(K)
        layer_initial_mass: xr.DataArray,      # (nface,) g/cm^2  TSED0(K)
        bulk_density: xr.DataArray,            # (nface,) g/cm^3  BULKDENS(K)
        core_id: xr.DataArray,                 # (nface,) int
    ) -> xr.DataArray:
        """g/cm²/s erosion rate, before per-class fractionation."""
        ...


def _to_dataarray_like(values: np.ndarray, template: xr.DataArray) -> xr.DataArray:
    """Wrap a numpy array as a DataArray with the same dims/coords as template."""
    return xr.DataArray(values, dims=template.dims, coords=template.coords)


class SedflumeTableErosionModel:
    """``nsedflume=1`` — interpolate the SEDflume ERATE table.

    Constructor takes the loaded SEDflume tables; per-core, per-layer,
    per-shear-level erosion rates (cm/s).

    The interpolation is bilinear in (τ, layer fractional mass remaining):

    * Linear in τ between adjacent ``tau_levels_pa`` brackets.
    * Log-linear in fractional remaining mass ``m_K / m_K0`` between the
      current layer (K) and the next deeper layer (K+1). For the deepest
      layer (K = K_B) the next-layer rate is replaced by a small floor
      value (``1e-9 cm/s``) so erosion smoothly tapers to zero rather
      than punching through the bottom.

    Out-of-range τ values are clamped to the table endpoints (with a
    one-time warning), mirroring SEDZLJ's tabular extrapolation guard.
    """

    def __init__(
        self,
        tau_levels_pa: np.ndarray,                  # (ITBM,) shear-stress interpolants
        erate_per_core: np.ndarray,                 # (n_cores, K_B, ITBM) cm/s
        erate_active_per_size: np.ndarray,          # (NSICM, ITBM) cm/s — for active/deposited layers
        size_interpolants_um: np.ndarray,           # (NSICM,) μm — SCND
        taucrit_per_size_pa: np.ndarray,            # (NSICM,) Pa — TAUCRITE
        consolidation_model: Optional[ConsolidationModel] = None,
    ) -> None:
        tau_levels_pa = np.asarray(tau_levels_pa, dtype=np.float64)
        erate_per_core = np.asarray(erate_per_core, dtype=np.float64)
        erate_active_per_size = np.asarray(erate_active_per_size, dtype=np.float64)
        size_interpolants_um = np.asarray(size_interpolants_um, dtype=np.float64)
        taucrit_per_size_pa = np.asarray(taucrit_per_size_pa, dtype=np.float64)

        if tau_levels_pa.ndim != 1:
            raise ValueError(
                f"tau_levels_pa must be 1-D (ITBM,); got shape {tau_levels_pa.shape}"
            )
        itbm = tau_levels_pa.shape[0]
        if itbm < 2:
            raise ValueError(
                f"tau_levels_pa must have at least 2 levels for interpolation; got {itbm}"
            )
        if not np.all(np.diff(tau_levels_pa) > 0):
            raise ValueError("tau_levels_pa must be strictly increasing")

        if erate_per_core.ndim != 3:
            raise ValueError(
                f"erate_per_core must be 3-D (n_cores, K_B, ITBM); got shape {erate_per_core.shape}"
            )
        if erate_per_core.shape[2] != itbm:
            raise ValueError(
                f"erate_per_core last dim ({erate_per_core.shape[2]}) "
                f"must equal len(tau_levels_pa)={itbm}"
            )
        if np.any(erate_per_core <= 0):
            # Log interpolation requires strictly positive rates; clip very-low
            # values to the floor instead of erroring (SEDflume tables sometimes
            # contain zeros for the lowest shear level).
            erate_per_core = np.where(
                erate_per_core <= 0, _DEEPEST_LAYER_ERATE_FLOOR_CM_S, erate_per_core
            )

        if size_interpolants_um.ndim != 1:
            raise ValueError(
                "size_interpolants_um must be 1-D (NSICM,); "
                f"got shape {size_interpolants_um.shape}"
            )
        nsicm = size_interpolants_um.shape[0]
        if erate_active_per_size.shape != (nsicm, itbm):
            raise ValueError(
                f"erate_active_per_size shape {erate_active_per_size.shape} "
                f"must equal (NSICM={nsicm}, ITBM={itbm})"
            )
        if taucrit_per_size_pa.shape != (nsicm,):
            raise ValueError(
                f"taucrit_per_size_pa shape {taucrit_per_size_pa.shape} "
                f"must equal (NSICM={nsicm},)"
            )

        self.tau_levels_pa = tau_levels_pa
        self.erate_per_core = erate_per_core
        self.erate_active_per_size = erate_active_per_size
        self.size_interpolants_um = size_interpolants_um
        self.taucrit_per_size_pa = taucrit_per_size_pa

        self.n_cores = erate_per_core.shape[0]
        self.k_b = erate_per_core.shape[1]
        self.itbm = itbm
        self.nsicm = nsicm
        # Optional consolidation model used by the SSM driver to age the
        # per-(layer, class) τ_ce gate; this class itself does not apply
        # the gate. Stored as an attribute so the driver can introspect.
        self.consolidation_model: Optional[ConsolidationModel] = consolidation_model

    def erosion_rate(
        self,
        tau_pa: xr.DataArray,
        layer_index: int,
        layer_mass: xr.DataArray,
        layer_initial_mass: xr.DataArray,
        bulk_density: xr.DataArray,
        core_id: xr.DataArray,
    ) -> xr.DataArray:
        """Bilinear (log-depth, linear-τ) interpolation in the SEDflume ERATE table.

        Returns
        -------
        xr.DataArray
            Per-cell erosion rate in g/cm²/s (= cm/s × g/cm³).
        """
        if not (1 <= layer_index <= self.k_b):
            raise IndexError(
                f"layer_index {layer_index} out of range [1, {self.k_b}]"
            )

        # Convert inputs to numpy for vectorized work; keep template for the
        # final wrap-back.
        template = tau_pa
        tau = np.asarray(tau_pa.values, dtype=np.float64)
        m_k = np.asarray(layer_mass.values, dtype=np.float64)
        m_k0 = np.asarray(layer_initial_mass.values, dtype=np.float64)
        rho_b = np.asarray(bulk_density.values, dtype=np.float64)
        cores = np.asarray(core_id.values, dtype=np.int64)

        # Validate core indices.
        if cores.size and (cores.min() < 0 or cores.max() >= self.n_cores):
            raise IndexError(
                f"core_id values must be in [0, {self.n_cores - 1}]; "
                f"got range [{cores.min()}, {cores.max()}]"
            )

        # τ-bracket: find indices i such that tau_levels[i-1] <= tau < tau_levels[i].
        # Clamp out-of-range values to the endpoints with a one-time warning.
        tau_min = self.tau_levels_pa[0]
        tau_max = self.tau_levels_pa[-1]
        if np.any(tau < tau_min) or np.any(tau > tau_max):
            warnings.warn(
                f"tau_pa contains values outside SEDflume table range "
                f"[{tau_min:.4g}, {tau_max:.4g}] Pa; clamping for interpolation.",
                RuntimeWarning,
                stacklevel=2,
            )
        tau_clamped = np.clip(tau, tau_min, tau_max)

        # searchsorted with side='right' on a strictly increasing array gives
        # i = number of levels <= tau. We want bracket [i_lo, i_hi] with
        # i_lo = max(0, i-1), i_hi = min(itbm-1, i).
        idx_hi = np.searchsorted(self.tau_levels_pa, tau_clamped, side="right")
        idx_hi = np.clip(idx_hi, 1, self.itbm - 1)
        idx_lo = idx_hi - 1

        tau_low = self.tau_levels_pa[idx_lo]
        tau_high = self.tau_levels_pa[idx_hi]

        # τ weighting (matches SEDZLJ SN00/SN10 with tau_high > tau_low):
        #   SN00 = (tau_high - tau) / (tau_high - tau_low)  -> weight on tau_low
        #   SN10 = 1 - SN00                                  -> weight on tau_high
        denom = tau_high - tau_low
        # Shouldn't divide by zero given strictly-increasing levels.
        sn00 = (tau_high - tau_clamped) / denom
        sn10 = 1.0 - sn00

        # Depth weighting from layer fractional remaining mass.
        # Guard against zero or negative initial mass (empty layer): if m_K0 <= 0,
        # there is no sediment to erode -> rate = 0.
        valid_mass = m_k0 > 0.0
        # Compute ratio safely; for invalid cells we'll zero the result later.
        with np.errstate(divide="ignore", invalid="ignore"):
            sn01 = np.where(valid_mass, m_k / m_k0, 0.0)
        sn01 = np.clip(sn01, 0.0, 1.0)
        sn11 = 1.0 - sn01

        # Per-cell rate at the shear bracket endpoints. Index ERATE[core, K-1, level].
        k0 = layer_index - 1  # convert to 0-origin
        # Gather per-cell vectors of shape (nface,) for the four corners of the
        # (τ-bracket, depth) bilinear cell.
        e_k_low = self.erate_per_core[cores, k0, idx_lo]
        e_k_high = self.erate_per_core[cores, k0, idx_hi]

        if layer_index < self.k_b:
            k1 = layer_index  # K+1 in 0-origin
            e_kp1_low = self.erate_per_core[cores, k1, idx_lo]
            e_kp1_high = self.erate_per_core[cores, k1, idx_hi]
        else:
            # Deepest layer: replace K+1 rate with the floor to taper erosion
            # smoothly to zero as the layer empties (s_sedzlj.f90:501 behaviour).
            floor = np.full_like(e_k_low, _DEEPEST_LAYER_ERATE_FLOOR_CM_S)
            e_kp1_low = floor
            e_kp1_high = floor

        # Log-linear depth interpolation at each shear bracket endpoint.
        log_e_low = sn11 * np.log(e_kp1_low) + sn01 * np.log(e_k_low)
        log_e_high = sn11 * np.log(e_kp1_high) + sn01 * np.log(e_k_high)

        # NB: SEDZLJ writes the formula as
        #   ERATEMOD = (SN00 * exp(SN11*ln(E_{K+1,low})  + SN01*ln(E_{K,low}))
        #             + SN10 * exp(SN11*ln(E_{K+1,high}) + SN01*ln(E_{K,high}))) * BULKDENS
        # i.e. SN00 weights the low-τ corner (which interpolates depth at tau_low)
        # and SN10 weights the high-τ corner. We follow that convention exactly.
        e_rate_cm_s = sn00 * np.exp(log_e_low) + sn10 * np.exp(log_e_high)
        e_rate_g_cm2_s = e_rate_cm_s * rho_b

        # Zero out cells with no remaining (or never-existed) layer mass.
        e_rate_g_cm2_s = np.where(valid_mass & (m_k > 0.0), e_rate_g_cm2_s, 0.0)

        return _to_dataarray_like(e_rate_g_cm2_s, template)


class PowerLawErosionModel:
    """``nsedflume=2`` — per-core per-layer power law E = A τ^n, capped at E_max.

    Per layer K, per core, the surface erosion rate (cm/s) at the top of
    the layer is

    .. math::

        E_{\\rm top}(\\tau) = A_K\\, \\tau^{n_K}

    with τ in **Pascals** (the EFDC code converts ``0.1 * TAU`` because it
    carries τ internally in dynes/cm²; we receive τ in Pa directly so no
    conversion is needed). The effective rate within the layer is a
    linear blend with the next-deeper layer's rate, weighted by how much
    of the layer has been removed:

    .. math::

        E = (E_{\\rm bottom} - E_{\\rm top}) \\, \\frac{m_{K,0} - m_K}{m_{K,0}}
            + E_{\\rm top}

    so a fully-fresh layer (``m_K = m_{K,0}``) erodes at ``E_top`` and a
    nearly-emptied layer transitions toward ``E_bottom``. For the deepest
    layer ``E_bottom = 0`` to prevent through-bottom erosion. The final
    rate is multiplied by the dry bulk density to convert cm/s → g/cm²/s
    and capped at ``max_rate_per_core[core, K]``.
    """

    def __init__(
        self,
        ea_per_core: np.ndarray,        # (n_cores, K_B) coefficient A
        en_per_core: np.ndarray,        # (n_cores, K_B) exponent n
        max_rate_per_core: np.ndarray,  # (n_cores, K_B) cap g/cm^2/s
        actdep_a: np.ndarray,           # (NSICM,) for active/deposited layers
        actdep_n: np.ndarray,
        actdep_max: np.ndarray,
        consolidation_model: Optional[ConsolidationModel] = None,
    ) -> None:
        ea_per_core = np.asarray(ea_per_core, dtype=np.float64)
        en_per_core = np.asarray(en_per_core, dtype=np.float64)
        max_rate_per_core = np.asarray(max_rate_per_core, dtype=np.float64)
        actdep_a = np.asarray(actdep_a, dtype=np.float64)
        actdep_n = np.asarray(actdep_n, dtype=np.float64)
        actdep_max = np.asarray(actdep_max, dtype=np.float64)

        if ea_per_core.ndim != 2:
            raise ValueError(
                f"ea_per_core must be 2-D (n_cores, K_B); got shape {ea_per_core.shape}"
            )
        if en_per_core.shape != ea_per_core.shape:
            raise ValueError(
                f"en_per_core shape {en_per_core.shape} "
                f"must match ea_per_core shape {ea_per_core.shape}"
            )
        if max_rate_per_core.shape != ea_per_core.shape:
            raise ValueError(
                f"max_rate_per_core shape {max_rate_per_core.shape} "
                f"must match ea_per_core shape {ea_per_core.shape}"
            )

        if not (actdep_a.shape == actdep_n.shape == actdep_max.shape):
            raise ValueError(
                "actdep_a, actdep_n, and actdep_max must all have the same shape; "
                f"got {actdep_a.shape}, {actdep_n.shape}, {actdep_max.shape}"
            )
        if actdep_a.ndim != 1:
            raise ValueError(
                f"actdep_* arrays must be 1-D (NSICM,); got shape {actdep_a.shape}"
            )

        self.ea_per_core = ea_per_core
        self.en_per_core = en_per_core
        self.max_rate_per_core = max_rate_per_core
        self.actdep_a = actdep_a
        self.actdep_n = actdep_n
        self.actdep_max = actdep_max

        self.n_cores = ea_per_core.shape[0]
        self.k_b = ea_per_core.shape[1]
        self.nsicm = actdep_a.shape[0]
        # Optional consolidation model — see SedflumeTableErosionModel.
        self.consolidation_model: Optional[ConsolidationModel] = consolidation_model

    def erosion_rate(
        self,
        tau_pa: xr.DataArray,
        layer_index: int,
        layer_mass: xr.DataArray,
        layer_initial_mass: xr.DataArray,
        bulk_density: xr.DataArray,
        core_id: xr.DataArray,
    ) -> xr.DataArray:
        """Power-law in-place-layer erosion rate, in g/cm²/s.

        See class docstring for the equation. Vectorized across cells.
        """
        if not (1 <= layer_index <= self.k_b):
            raise IndexError(
                f"layer_index {layer_index} out of range [1, {self.k_b}]"
            )

        template = tau_pa
        tau = np.asarray(tau_pa.values, dtype=np.float64)
        m_k = np.asarray(layer_mass.values, dtype=np.float64)
        m_k0 = np.asarray(layer_initial_mass.values, dtype=np.float64)
        rho_b = np.asarray(bulk_density.values, dtype=np.float64)
        cores = np.asarray(core_id.values, dtype=np.int64)

        if cores.size and (cores.min() < 0 or cores.max() >= self.n_cores):
            raise IndexError(
                f"core_id values must be in [0, {self.n_cores - 1}]; "
                f"got range [{cores.min()}, {cores.max()}]"
            )

        k0 = layer_index - 1
        a_top = self.ea_per_core[cores, k0]
        n_top = self.en_per_core[cores, k0]
        cap = self.max_rate_per_core[cores, k0]

        # τ in Pa raised to n_top, with a guard for τ <= 0 (no erosion).
        # Use a where to avoid 0**negative -> inf.
        positive_tau = tau > 0.0
        # Use safe substitute (1.0) where tau<=0 to avoid runtime warnings;
        # we'll mask back to zero at the end.
        tau_safe = np.where(positive_tau, tau, 1.0)

        e_top = a_top * np.power(tau_safe, n_top)

        if layer_index < self.k_b:
            k1 = layer_index
            a_bot = self.ea_per_core[cores, k1]
            n_bot = self.en_per_core[cores, k1]
            e_bot = a_bot * np.power(tau_safe, n_bot)
        else:
            # Deepest layer: no erosion through the bottom (s_sedzlj.f90:513).
            e_bot = np.zeros_like(e_top)

        # Mass weighting: full layer (m_k = m_k0) → SN11 = 0 → E = E_top.
        # Empty layer (m_k = 0)            → SN11 = 1 → E = E_bottom.
        valid_mass = m_k0 > 0.0
        with np.errstate(divide="ignore", invalid="ignore"):
            sn11 = np.where(valid_mass, (m_k0 - m_k) / m_k0, 0.0)
        sn11 = np.clip(sn11, 0.0, 1.0)

        e_rate_cm_s = (e_bot - e_top) * sn11 + e_top
        e_rate_g_cm2_s = e_rate_cm_s * rho_b

        # Cap at per-(core, layer) maximum.
        e_rate_g_cm2_s = np.minimum(e_rate_g_cm2_s, cap)

        # Zero out: no shear, no remaining mass, or never-existed layer.
        active = positive_tau & valid_mass & (m_k > 0.0)
        e_rate_g_cm2_s = np.where(active, e_rate_g_cm2_s, 0.0)
        # Final non-negativity guard (cap could be negative if mis-specified;
        # erosion never runs backwards).
        e_rate_g_cm2_s = np.maximum(e_rate_g_cm2_s, 0.0)

        return _to_dataarray_like(e_rate_g_cm2_s, template)


def apply_consolidation(
    tau_ce_layer_class_pa: xr.DataArray,                      # (nface, ssm_layer, ssm_class)
    layer_age_s: xr.DataArray,                                # (nface, ssm_layer)
    is_cohesive: np.ndarray,                                  # (n_class,) bool
    consolidation_model: Optional[ConsolidationModel] = None,
) -> xr.DataArray:
    """Apply Sanford-Maa-style consolidation to per-(layer, class) τ_ce.

    For cohesive classes only (``is_cohesive[c] == True``), the per-layer
    effective τ_ce is replaced by ``consolidation_model.effective_tau_ce(
    layer_age_s)``. Non-cohesive classes retain the static value.

    If ``consolidation_model`` is ``None`` the input is returned
    unchanged — opt-in semantics matching SEDZLJ defaults.
    """
    if consolidation_model is None:
        return tau_ce_layer_class_pa
    # Delegate to consolidation module (handles broadcast + masking).
    from .consolidation import apply_consolidation_per_class
    return apply_consolidation_per_class(
        tau_ce_layer_class_pa=tau_ce_layer_class_pa,
        layer_age_s=layer_age_s,
        is_cohesive=is_cohesive,
        model=consolidation_model,
    )


def apply_vegetation_cohesion(
    tau_ce_pa: xr.DataArray,                                # (nface, ssm_class) base τ_ce
    biostabilization: xr.DataArray | None = None,           # (nface,) [0,1]
    root_cohesion_pa: xr.DataArray | None = None,           # (nface,) Pa
    biostabilization_alpha: float = 0.5,
) -> xr.DataArray:
    """Effective τ_ce given vegetation feedback.

    .. math::

        \\tau_{ce}^{\\rm eff} = \\tau_{ce}\\,(1 + \\alpha B) + \\tau_{\\rm root}

    Both ``biostabilization`` and ``root_cohesion_pa`` are broadcast across
    the class dimension if it is present in ``tau_ce_pa`` (they are typically
    per-cell quantities, while ``tau_ce_pa`` is per-cell per-class).

    If both inputs are ``None``, returns ``tau_ce_pa`` unchanged. If only
    one is supplied, the missing term is treated as zero.
    """
    if biostabilization is None and root_cohesion_pa is None:
        return tau_ce_pa

    result = tau_ce_pa
    if biostabilization is not None:
        # xarray broadcasts on shared dim names: biostabilization on (nface,)
        # combined with tau_ce_pa on (nface, ssm_class) yields (nface, ssm_class).
        result = result * (1.0 + biostabilization_alpha * biostabilization)
    if root_cohesion_pa is not None:
        result = result + root_cohesion_pa

    return result
