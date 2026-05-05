"""v3 NSM1 light extinction utilities.

Stateless light primitives ported from v1 NSM1's shared library
(``clearwater_modules/shared/processes.py``). Provides the Beer-Lambert
extinction coefficient with contributions from inorganic suspended solids,
particulate organic matter (via POC and the C/OM ratio), and chlorophyll-a
self-shading, plus the photosynthetically active radiation conversion.
"""

import xarray as xr


def L(
    lambda0: xr.DataArray,
    lambda1: xr.DataArray,
    lambda2: xr.DataArray,
    lambdas: xr.DataArray,
    lambdam: xr.DataArray,
    Solid: xr.DataArray,
    POC: xr.DataArray,
    fcom: xr.DataArray,
    Ap: xr.DataArray,
    use_Algae: xr.DataArray,
    use_POC: xr.DataArray,
) -> xr.DataArray:
    """Beer-Lambert light extinction coefficient.

    Args:
        lambda0 | 1/m | background extinction coefficient.
        lambda1 | 1/m / (ug-Chla/L) | linear chlorophyll self-shading
            coefficient.
        lambda2 | 1/m / (ug-Chla/L)^(2/3) | non-linear chlorophyll
            self-shading coefficient.
        lambdas | L/mg/m | inorganic suspended solids extinction coefficient.
        lambdam | L/mg/m | particulate organic matter extinction coefficient.
        Solid | mg/L | inorganic suspended solids concentration.
        POC | mg-C/L | particulate organic carbon concentration.
        fcom | mg-C/mg-D | carbon-to-organic-matter mass ratio.
        Ap | ug-Chla/L | floating algae chlorophyll-a concentration.
        use_Algae | bool | floating algae module switch.
        use_POC | bool | particulate organic carbon module switch.

    Returns:
        1/m | total light extinction coefficient.
    """
    extinction = lambda0 + lambdas * Solid
    extinction = xr.where(use_POC, extinction + lambdam * POC / fcom, extinction)
    extinction = xr.where(
        use_Algae,
        extinction + lambda1 * Ap + lambda2 * Ap ** 0.66667,
        extinction,
    )
    return extinction


def PAR(
    q_solar: xr.DataArray,
    Fr_PAR: xr.DataArray,
) -> xr.DataArray:
    """Photosynthetically active radiation at the water surface.

    Args:
        q_solar | W/m^2 | total incident solar radiation at the water surface.
        Fr_PAR | dimensionless | fraction of incident solar radiation in the
            PAR band (typical 0.45-0.50).

    Returns:
        W/m^2 | photosynthetically active radiation at the water surface.
    """
    return q_solar * Fr_PAR
