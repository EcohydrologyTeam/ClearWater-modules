"""
File contains process to calculate new algae biomass concentration and associated dependent variables
"""

import math
from clearwater_modules.shared.processes import arrhenius_correction
import numba
import xarray as xr

@numba.njit
def kdx_tc(
    TwaterC : xr.DataArray,
    kdx_20: xr.DataArray
) -> xr.DataArray :

    """Calculate kdx_tc: pathogen death rate (1/d).

    Args:
        TwaterC: Water temperature (C)
        kdx_20: Pathogen death rate at 20 degree (1/d)
    """

    return arrhenius_correction(TwaterC, kdx_20, 1.07)

@numba.njit
def PathogenDeath(
    kdx_tc : xr.DataArray,
    PX: xr.DataArray
) -> xr.DataArray :

    """Calculate PathogenDeath: pathogen natural death (cfu/100mL/d).

    Args:
      kdx_tc: pathogen death rate with temperature correction (1/d),
      PX: pathogen concentration (cfu/100mL)

    """
    return kdx_tc * PX

@numba.njit
def PathogenDecay(
  apx: xr.DataArray,
  q_solar: xr.DataArray,
  L: xr.DataArray,
  depth: xr.DataArray,
  PX: xr.DataArray
) -> xr.DataArray :

    """Calculate PathogenDecay: pathogen death due to light (cfu/100mL/d).

    Args:
      apx: light efficiency factor for pathogen decay,
      q_solar: solar radiation (1/d),
      L: lambda (1/m),
      depth: water depth (m),
      PX: Pathogen concentration (cfu/100mL)

    """
    return apx * q_solar / (L*depth) * (1-math.exp(-L*depth)) * PX

@numba.njit
def PathogenSettling(
  vx: xr.DataArray,
  depth: xr.DataArray,
  PX: xr.DataArray
) -> xr.DataArray :

    """Calculate PathogenSettling: pathogen settling (cfu/100mL/d).

    Args:
      vx: pathogen net settling velocity (m)
      depth: water depth (m),
      PX: Pathogen concentration (cfu/100mL)
    """
    return vx/depth*PX

@numba.njit
def dPXdt(
  PathogenDeath: xr.DataArray,
  PathogenDecay: xr.DataArray,
  PathogenSettling: xr.DataArray,

) -> xr.DataArray :

    """Calculate dPXdt: change in pathogen concentration (cfu/100mL/d).

    Args:
      PathogenSettling: pathogen settling (cfu/100mL/d)
      PathogenDecay: pathogen death due to light (cfu/100mL/d)
      PathogenDeath: pathogen natural death (cfu/100mL/d)

    """
    return -PathogenDeath - PathogenDecay - PathogenSettling

@numba.njit
def PX_new(
  PX:xr.DataArray,
  dPXdt: xr.DataArray

) -> xr.DataArray :

    """Calculate PX_new: New pathogen concentration (cfu/100mL).

    Args:
      dPXdt: change in pathogen concentration (cfu/100mL/d)
      PX: Pathogen concentration (cfu/100mL)
    """
    return PX+dPXdt