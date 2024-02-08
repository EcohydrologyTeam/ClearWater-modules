"""Process functions used by one or more modules"""
import numba
import xarray as xr


@numba.njit
def celsius_to_kelvin(tempc: xr.DataArray) -> xr.DataArray:
    return tempc + 273.16


@numba.njit
def kelvin_to_celsius(tempk: xr.DataArray) -> xr.DataArray:
    return tempk - 273.16

@numba.njit
def arrhenius_correction(
    water_temp_c: xr.DataArray,
    rc20: xr.DataArray,
    theta: xr.DataArray,
) -> xr.DataArray:
    """
    Computes an adjusted kinetics reaction rate coefficient for the specified water 
    temperature using the van't Hoff form of the Arrhenius equation

    Parameters
    ----------
    water_temp_c : xr.DataArray
        Water temperature in degrees Celsius
    rc20 : xr.DataArray
        Kinetics reaction (decay) coefficient at 20 degrees Celsius
    theta : xr.DataArray
        Temperature correction factor

    Returns
    ----------
    float
        Adjusted kinetics rate for the specified water temperature
    """
    return rc20 * theta**(water_temp_c - 20.0)

@numba.njit
def compute_depth(
    surface_area: xr.DataArray,
    volume: xr.DataArray
) -> xr.DataArray:
    """Compute depth of a computation cell

    Args:
        surface_area: state variable for surface area of computational cell provided by CWR engine
        volume: state variable for volume of computational cell provided by CWR engine
    """
    return volume / surface_area

@numba.njit
def L(
    lambda0: xr.DataArray,
    lambda1: xr.DataArray,
    lambda2: xr.DataArray,
    lambdas: xr.DataArray,
    lambdam: xr.DataArray,
    Solid: xr.DataArray,
    POC: xr.DataArray,
    focm: xr.DataArray,
    use_Algae: xr.DataArray,
    use_POC: xr.DataArray,
    Ap: xr.DataArray,

) -> xr.DataArray:
    """Compute L: lambda

    Args:
        lambda0: xr.DataArray,
        lambda1: xr.DataArray,
        lambda2: xr.DataArray,
        lambdas: xr.DataArray,
        lambdam: xr.DataArray,
    """
    L=lambdas * Solid
    L=xr.where (use_POC, L=L+lambdam*POC/focm,
                xr.where (use_Algae, L+ lambda1*Ap + lambda2*Ap**0.66667))

    return L

@numba.njit
def TwaterK(
    TwaterC : xr.DataArray,
) -> xr.DataArray :
    """Calculate temperature in kelvin (K)
    Args:
        TwaterC: water temperature celcius (C)
    """
    return celsius_to_kelvin(TwaterK)