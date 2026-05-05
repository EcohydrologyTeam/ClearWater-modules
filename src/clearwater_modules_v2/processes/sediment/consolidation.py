"""Cohesive-bed consolidation: age-dependent critical shear stress.

Distinguishes SSM from the SEDZLJ baseline (which holds τ_ce constant
once initialised) by letting the effective τ_ce grow with the age of a
cohesive bed layer. This addresses a known SEDZLJ limitation flagged in
the source itself (``s_sedzlj.f90:707`` "SEDZLJ DOES NOT HAVE
CONSOLIDATION") and brings SSM closer to the cohesive-bed treatment in
MIKE 21, Delft3D, and TELEMAC.

Formulation
-----------
For each cohesive class (D50 < bedload cutoff), the layer-effective
critical shear stress at age :math:`t_{\\rm age}` is:

.. math::

    \\tau_{ce}^{\\rm eff}(t_{\\rm age}) =
        \\tau_{ce,\\infty}
        - (\\tau_{ce,\\infty} - \\tau_{ce,0})
          \\, \\exp(-t_{\\rm age} / T_c)

where:

* :math:`\\tau_{ce,0}` is the freshly-deposited critical shear stress
  (lower limit; equal to or close to the unconsolidated τ_ce from the
  bed-input file).
* :math:`\\tau_{ce,\\infty}` is the fully-consolidated critical shear
  stress (upper limit; typically 3×–5× the unconsolidated value, per
  Sanford & Maa 2001 calibration).
* :math:`T_c` is the consolidation time scale (typically days to weeks
  for cohesive sediments).
* :math:`t_{\\rm age}` is the layer's mean age at the current step.

References
----------
* Sanford, L. P., and Maa, J. P.-Y. (2001). "A unified erosion
  formulation for fine sediments." *Marine Geology* 179(1–2), 9–23.
  DOI: 10.1016/S0025-3227(01)00201-8.
* Mehta, A. J., and Partheniades, E. (1975). "An investigation of the
  depositional properties of flocculated fine sediments."
  *J. Hydraul. Res.* 13(4), 361–381.
  DOI: 10.1080/00221687509499694.

Limitations
-----------
This first release scopes consolidation to the τ_ce(age) formulation
**only**. A more complete model would additionally evolve:

* time-varying porosity / bulk density (Gibson, England & Hussey 1967),
* gel-point dynamics (Toorman 1999),
* finite-strain (large deformation) self-weight consolidation.

Those are deferred to a follow-on release. See
``ClearWater-Riverine-streaming/design/ssm_consolidation.md`` for the
full design memo.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np
import xarray as xr


@runtime_checkable
class ConsolidationModel(Protocol):
    """Protocol for layer-age → effective τ_ce models.

    Implementations must provide :meth:`effective_tau_ce` returning a
    DataArray broadcastable across ``(nface, ssm_layer)``.
    """

    def effective_tau_ce(self, layer_age_s: xr.DataArray) -> xr.DataArray:
        """Return effective τ_ce(age) in Pa, same dims as ``layer_age_s``."""
        ...


@dataclass
class SanfordMaaConsolidation:
    """Sanford & Maa (2001) single-mode age-dependent τ_ce.

    Parameters
    ----------
    tau_ce_zero_pa : float
        Freshly-deposited critical shear stress (Pa). Lower limit of the
        consolidation envelope; recovered as :math:`t_{\\rm age} \\to 0`.
    tau_ce_inf_pa : float
        Fully-consolidated critical shear stress (Pa). Upper limit;
        approached as :math:`t_{\\rm age} \\to \\infty`. Must be
        ``>= tau_ce_zero_pa``.
    consolidation_time_s : float
        Consolidation time scale :math:`T_c` (s). The age at which
        :math:`\\tau_{ce}^{\\rm eff}` reaches roughly
        :math:`\\tau_{ce,0} + (\\tau_{ce,\\infty}-\\tau_{ce,0})(1-1/e)
        \\approx \\tau_{ce,0} + 0.632(\\tau_{ce,\\infty}-\\tau_{ce,0})`.
        Must be strictly positive.
    """

    tau_ce_zero_pa: float
    tau_ce_inf_pa: float
    consolidation_time_s: float

    def __post_init__(self) -> None:
        if self.consolidation_time_s <= 0.0:
            raise ValueError(
                "consolidation_time_s must be > 0; got "
                f"{self.consolidation_time_s!r}"
            )
        if self.tau_ce_zero_pa < 0.0:
            raise ValueError(
                f"tau_ce_zero_pa must be >= 0; got {self.tau_ce_zero_pa!r}"
            )
        if self.tau_ce_inf_pa < self.tau_ce_zero_pa:
            raise ValueError(
                "tau_ce_inf_pa must be >= tau_ce_zero_pa; got "
                f"tau_ce_inf={self.tau_ce_inf_pa!r}, "
                f"tau_ce_zero={self.tau_ce_zero_pa!r}"
            )

    def effective_tau_ce(self, layer_age_s: xr.DataArray) -> xr.DataArray:
        """Effective τ_ce(age) in Pa.

        .. math::

            \\tau_{ce}^{\\rm eff}(t) =
                \\tau_{ce,\\infty}
                - (\\tau_{ce,\\infty} - \\tau_{ce,0})
                  \\, \\exp(-t / T_c)

        Parameters
        ----------
        layer_age_s : xr.DataArray
            Per-layer age (s); typically dims ``(nface, ssm_layer)``.
            Negative ages (numerically possible after a borrow that
            re-blends) are clamped to zero.

        Returns
        -------
        xr.DataArray
            Effective τ_ce in Pa with the same dims/coords as ``layer_age_s``.
        """
        # Clamp any small negative values from float drift.
        age = xr.where(layer_age_s < 0.0, 0.0, layer_age_s)
        tau_inf = float(self.tau_ce_inf_pa)
        tau_zero = float(self.tau_ce_zero_pa)
        tc = float(self.consolidation_time_s)
        delta = tau_inf - tau_zero
        return tau_inf - delta * np.exp(-age / tc)


def apply_consolidation_per_class(
    tau_ce_layer_class_pa: xr.DataArray,
    layer_age_s: xr.DataArray,
    is_cohesive: np.ndarray,
    model: ConsolidationModel,
) -> xr.DataArray:
    """Apply consolidation to per-(layer, class) τ_ce, cohesive classes only.

    For cohesive classes, the per-layer effective τ_ce becomes the
    consolidation-aged value (broadcast across the class dimension);
    non-cohesive classes retain their static τ_ce.

    Parameters
    ----------
    tau_ce_layer_class_pa : xr.DataArray
        Baseline per-layer per-class τ_ce, dims ``(nface, ssm_layer,
        ssm_class)``.
    layer_age_s : xr.DataArray
        Per-layer age, dims ``(nface, ssm_layer)``.
    is_cohesive : np.ndarray
        Per-class boolean flag; ``True`` for cohesive classes (D50
        below bedload cutoff). Shape ``(n_class,)``.
    model : ConsolidationModel
        The age → τ_ce mapping. Typically
        :class:`SanfordMaaConsolidation`.

    Returns
    -------
    xr.DataArray
        Effective per-layer per-class τ_ce with the same dims as the
        input, holding ``model.effective_tau_ce(age)`` on cohesive
        classes and the input on non-cohesive ones.
    """
    eff_layer = model.effective_tau_ce(layer_age_s)  # (nface, ssm_layer)
    # Broadcast (nface, ssm_layer) across ssm_class by manual tile so we
    # don't rely on xarray's broadcast rules (which won't add a new
    # dimension by themselves).
    n_class = tau_ce_layer_class_pa.sizes["ssm_class"]
    eff_layer_class_arr = np.broadcast_to(
        eff_layer.values[..., None], eff_layer.values.shape + (n_class,)
    ).copy()
    eff_layer_class = xr.DataArray(
        eff_layer_class_arr,
        dims=tuple(eff_layer.dims) + ("ssm_class",),
    )
    # Per-class boolean mask aligned to ssm_class.
    cohesive_mask = xr.DataArray(
        np.asarray(is_cohesive, dtype=bool),
        dims=("ssm_class",),
    )
    # xr.where re-orders dims to put broadcast-only dims first; transpose
    # back so the result has the same dim order as the input baseline.
    blended = xr.where(cohesive_mask, eff_layer_class, tau_ce_layer_class_pa)
    return blended.transpose(*tau_ce_layer_class_pa.dims)
