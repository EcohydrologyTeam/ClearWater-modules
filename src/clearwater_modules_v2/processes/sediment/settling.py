"""Settling velocities (Cheng 1997).

Reference: SAND2008-5621 eq. 3; Cheng N.-S. (1997) "Simplified settling
velocity formula for sediment particle." J. Hydraul. Eng. 123(2), 149–152.
Design spec §5.1.

Public API:
    cheng_1997_settling_velocity(d50_um, ...) -> float | np.ndarray
    resolve_settling_velocities(registry) -> np.ndarray  (cm/s, indexed by class)
"""

from __future__ import annotations

from typing import Union

import numpy as np

from . import contracts
from .classes import SedimentClassRegistry

ArrayLike = Union[float, np.ndarray]


def cheng_1997_settling_velocity(
    d50_um: ArrayLike,
    solid_specific_gravity: float = contracts.DEFAULT_SOLID_SPECIFIC_GRAVITY,
    water_density_g_cm3: float = contracts.DEFAULT_WATER_DENSITY_CGS,
    kinematic_viscosity_cm2_s: float = contracts.NU_CGS,
    g_cm_s2: float = contracts.G_CGS,
) -> ArrayLike:
    """Settling velocity in cm/s via Cheng (1997).

    .. math::

        d_* = D_{50}\\bigl[(s_s - 1)\\,g / \\nu^2\\bigr]^{1/3}

        w_s = (\\nu / D_{50})\\bigl(\\sqrt{25 + 1.2\\,d_*^2} - 5\\bigr)^{1.5}

    Parameters
    ----------
    d50_um : float or array
        Median grain diameter in micrometres.
    solid_specific_gravity : float
        Sediment specific gravity (default 2.65, quartz).
    water_density_g_cm3 : float
        Water density (default 1.0).
    kinematic_viscosity_cm2_s : float
        Kinematic viscosity (default 0.01, ~20 °C).
    g_cm_s2 : float
        Gravitational acceleration in CGS (default 980).

    Returns
    -------
    float or array
        Settling velocity in cm/s, same shape as ``d50_um``.
    """
    # Convert μm -> cm for the algebraic form below.
    d_cm = np.asarray(d50_um, dtype="float64") * 1.0e-4

    # Avoid division-by-zero downstream; D50 must be strictly positive.
    if np.any(d_cm <= 0.0):
        raise ValueError(
            "cheng_1997_settling_velocity requires D50 > 0 (got non-positive value)."
        )

    s_s = solid_specific_gravity / water_density_g_cm3
    nu = kinematic_viscosity_cm2_s

    # Dimensionless particle parameter d* = D * [(s_s - 1) g / nu^2]^(1/3)
    d_star = d_cm * np.cbrt((s_s - 1.0) * g_cm_s2 / (nu * nu))

    # Cheng (1997) closed-form settling velocity.
    inner = np.sqrt(25.0 + 1.2 * d_star * d_star) - 5.0
    # Mathematically inner > 0 for any d_star > 0; clip defensively to avoid
    # NaNs from floating-point noise when d_star is extremely small.
    inner = np.maximum(inner, 0.0)
    w_s = (nu / d_cm) * np.power(inner, 1.5)

    # Preserve scalar input -> scalar output.
    if np.ndim(d50_um) == 0:
        return float(w_s)
    return w_s


def resolve_settling_velocities(registry: SedimentClassRegistry) -> np.ndarray:
    """Build the per-class settling velocity vector (cm/s).

    For each class in ``registry``, returns the user-supplied
    ``settling_cm_s`` if set and positive, otherwise computes it via
    :func:`cheng_1997_settling_velocity` from ``d50_um``.
    """
    out = np.empty(len(registry), dtype="float64")
    for i, cls in enumerate(registry):
        ws = cls.settling_cm_s
        if ws is not None and ws > 0.0:
            out[i] = float(ws)
        else:
            out[i] = float(
                cheng_1997_settling_velocity(
                    cls.d50_um,
                    solid_specific_gravity=cls.solid_density_g_cm3,
                )
            )
    return out
