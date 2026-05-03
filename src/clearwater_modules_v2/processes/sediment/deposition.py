"""Deposition probabilities and rates.

* Cohesive (D50 < bedload cutoff) — Krone (1962) linear probability.
* Non-cohesive (D50 >= bedload cutoff) — Gessler (1965) erfc-based
  probability. The implementation uses :func:`scipy.special.erfc`
  directly (numerically equivalent to the Abramowitz-Stegun rational
  approximation used in EFDC ``s_sedzlj.f90:141-148``, but cleaner).

References
----------
- Gessler, J. (1965). The beginning of bedload movement of mixtures
  investigated as natural armoring in channels.
- Krone, R. B. (1962). Flume Studies of the Transport of Sediment in
  Estuarial Shoaling Processes.
- SAND2008-5621 (Thanh, Grace & James 2008); SSM design spec sec 5.5.
"""

from __future__ import annotations

import math

import numpy as np
import xarray as xr
from scipy.special import erfc

from . import contracts
from .classes import SedimentClassRegistry


# Numerical floor for tau when computing 1/tau ratios. Below this we treat
# tau as effectively zero (no shear) and force P_d -> 1.
_TAU_FLOOR_PA: float = 1.0e-12


def gessler_probability(
    tau_pa: xr.DataArray,
    tau_cs_pa: float,
) -> xr.DataArray:
    r"""Gessler (1965) deposition probability for non-cohesive sediment.

    .. math::

        P_y = \frac{1}{0.57}\left(\frac{\tau_{cs}}{\tau} - 1\right)

        P_d = \tfrac{1}{2}\,\mathrm{erfc}(-P_y / \sqrt{2})

    The second equation is the standard form of the cumulative normal
    distribution (mean 0, variance 1) evaluated at ``P_y``, which yields a
    probability in [0, 1].

    Limits:
        * ``tau -> 0``  =>  ``P_y -> +inf``  =>  ``P_d -> 1`` (deposition certain).
        * ``tau == tau_cs``  =>  ``P_y = 0``  =>  ``P_d = 0.5``.
        * ``tau -> +inf``  =>  ``P_y -> -inf``  =>  ``P_d -> 0``.

    Parameters
    ----------
    tau_pa : xr.DataArray  (nface,)
        Bed shear stress (Pa).
    tau_cs_pa : float
        Class-level critical shear stress for suspension (Pa).

    Returns
    -------
    xr.DataArray  (nface,)
        Deposition probability in [0, 1].
    """
    if tau_cs_pa <= 0.0:
        raise ValueError(
            f"tau_cs_pa must be > 0; got {tau_cs_pa}"
        )

    # Branch on tau ~ 0 (no shear) to avoid 1/0 and force P_d = 1.
    safe_tau = xr.where(tau_pa > _TAU_FLOOR_PA, tau_pa, _TAU_FLOOR_PA)
    p_y = (1.0 / 0.57) * (tau_cs_pa / safe_tau - 1.0)

    # P_d = 0.5 * erfc(-P_y / sqrt(2))  -- standard normal CDF at P_y.
    sqrt2 = math.sqrt(2.0)
    p_d = 0.5 * xr.apply_ufunc(erfc, -p_y / sqrt2, dask="parallelized")

    # Force P_d = 1 where tau is effectively zero (deposition certain).
    p_d = xr.where(tau_pa > _TAU_FLOOR_PA, p_d, 1.0)

    # Guard numerical noise: clip to [0, 1].
    return p_d.clip(0.0, 1.0)


def krone_probability(
    tau_pa: xr.DataArray,
    tau_cs_pa: float,
) -> xr.DataArray:
    r"""Krone (1962) deposition probability for cohesive sediment.

    .. math::

        P_d = \max\bigl(1 - \tau / \tau_{cs},\, 0\bigr)

    Limits:
        * ``tau == 0``         =>  ``P_d = 1``.
        * ``tau == tau_cs``    =>  ``P_d = 0``.
        * ``tau == 0.5 tau_cs``=>  ``P_d = 0.5``.

    Parameters
    ----------
    tau_pa : xr.DataArray  (nface,)
        Bed shear stress (Pa).
    tau_cs_pa : float
        Class-level critical shear stress for suspension (Pa).

    Returns
    -------
    xr.DataArray  (nface,)
        Deposition probability in [0, 1].
    """
    if tau_cs_pa <= 0.0:
        raise ValueError(
            f"tau_cs_pa must be > 0; got {tau_cs_pa}"
        )
    return xr.where(tau_pa < tau_cs_pa, 1.0 - tau_pa / tau_cs_pa, 0.0)


def compute_deposition_flux(
    registry: SedimentClassRegistry,
    suspended_concentration: xr.DataArray,        # (nface, ssm_class) mg/L
    tau_pa: xr.DataArray,                         # (nface,) Pa
    settling_velocity_cm_s: np.ndarray,           # (ssm_class,) cm/s
    bottom_water_layer_depth_m: xr.DataArray,     # (nface,) m
    dt_seconds: float,
    max_deposit_fraction: float = 1.0,
) -> xr.DataArray:
    r"""Per-class deposition mass flux (g/cm^2) over one sediment time step.

    For each class ``s`` in the registry:

    .. math::

        D_s = P_{d,s}\, C_s\, w_{s,s}\, \Delta t

    where :math:`C_s` is the suspended concentration converted from mg/L to
    g/cm^3 (factor ``1e-6``), :math:`w_{s,s}` is the settling velocity in
    cm/s, and :math:`\Delta t` is in seconds. ``P_{d,s}`` is supplied by
    :func:`gessler_probability` for non-cohesive classes and
    :func:`krone_probability` for cohesive classes.

    Each cell-class deposition is capped at the available suspended mass in
    the bottom water layer (mirrors SEDZLJ ``MAXDEPLIMIT`` logic at
    ``s_sedzlj.f90:157``):

    .. math::

        D_s \le \mathrm{max\_deposit\_fraction}\, \cdot\, C_s\,
                \cdot\, h_{\mathrm{bot}}\,\cdot\, 100

    where ``h_bot`` is the bottom-water-layer depth in metres (factor 100
    converts to cm so that ``C_s [g/cm^3] x h [cm] = mass/area [g/cm^2]``).

    Parameters
    ----------
    registry : SedimentClassRegistry
        Ordered sediment classes. Selection of Gessler vs. Krone is keyed on
        :attr:`SedimentClass.is_cohesive`.
    suspended_concentration : xr.DataArray  (nface, ssm_class)
        Bottom-cell suspended concentration in mg/L.
    tau_pa : xr.DataArray  (nface,)
        Bed shear stress (Pa).
    settling_velocity_cm_s : np.ndarray  (ssm_class,)
        Per-class settling velocity in cm/s.
    bottom_water_layer_depth_m : xr.DataArray  (nface,)
        Depth of the bottom water layer (m).
    dt_seconds : float
        Sediment time step (s).
    max_deposit_fraction : float, default 1.0
        Fraction of available bottom-cell mass that may deposit in one step.

    Returns
    -------
    xr.DataArray  (nface, ssm_class)
        Per-class deposition mass per unit bed area (g/cm^2) for this step.
    """
    if dt_seconds <= 0.0:
        raise ValueError(f"dt_seconds must be > 0; got {dt_seconds}")
    if max_deposit_fraction <= 0.0:
        raise ValueError(
            f"max_deposit_fraction must be > 0; got {max_deposit_fraction}"
        )
    if len(settling_velocity_cm_s) != len(registry):
        raise ValueError(
            "settling_velocity_cm_s length "
            f"({len(settling_velocity_cm_s)}) does not match registry size "
            f"({len(registry)})"
        )
    if contracts.DIM_CLASS not in suspended_concentration.dims:
        raise ValueError(
            f"suspended_concentration missing dim '{contracts.DIM_CLASS}'; "
            f"got dims {suspended_concentration.dims}"
        )

    # mg/L  ->  g/cm^3 conversion factor: 1 mg/L = 1e-6 g/cm^3.
    MG_PER_L_TO_G_PER_CM3: float = 1.0e-6
    # m  ->  cm conversion factor (for the bottom-layer depth).
    M_TO_CM: float = 100.0

    # Build a per-class probability stack on (nface, ssm_class) by computing
    # each class's P_d as a function of tau and concatenating along the
    # ssm_class dim. Pure xarray, no python for-loops over cells.
    p_d_per_class: list[xr.DataArray] = []
    for s, sed_class in enumerate(registry):
        if sed_class.tau_cs_pa is None:
            raise ValueError(
                f"Class '{sed_class.label}' has tau_cs_pa = None; "
                "deposition probability requires a defined value."
            )
        if sed_class.is_cohesive:
            p_d_s = krone_probability(tau_pa, sed_class.tau_cs_pa)
        else:
            p_d_s = gessler_probability(tau_pa, sed_class.tau_cs_pa)
        # Tag with class-dim coord so concat aligns cleanly.
        p_d_s = p_d_s.expand_dims({contracts.DIM_CLASS: [s]})
        p_d_per_class.append(p_d_s)
    p_d = xr.concat(p_d_per_class, dim=contracts.DIM_CLASS)

    # Settling-velocity DataArray on (ssm_class,) so broadcasting picks up
    # the right axis automatically.
    w_s = xr.DataArray(
        np.asarray(settling_velocity_cm_s, dtype="float64"),
        dims=(contracts.DIM_CLASS,),
        coords={contracts.DIM_CLASS: np.arange(len(registry))},
    )

    # Suspended C in g/cm^3.
    c_g_per_cm3 = suspended_concentration * MG_PER_L_TO_G_PER_CM3

    # D_s (g/cm^2/step) = P_d * C [g/cm^3] * w_s [cm/s] * dt [s].
    deposition = p_d * c_g_per_cm3 * w_s * dt_seconds

    # Mass cap: D_s <= max_deposit_fraction * C [g/cm^3] * h [cm].
    # h_cm = bottom_water_layer_depth_m * 100.
    h_cm = bottom_water_layer_depth_m * M_TO_CM
    cap = max_deposit_fraction * c_g_per_cm3 * h_cm

    deposition = xr.where(deposition > cap, cap, deposition)

    # Floor at 0 to suppress any tiny negative numerical artefacts (none
    # are possible analytically since all factors are non-negative).
    deposition = xr.where(deposition > 0.0, deposition, 0.0)

    # Make sure the class dim ends up last for downstream consumers.
    return deposition.transpose(..., contracts.DIM_CLASS)
