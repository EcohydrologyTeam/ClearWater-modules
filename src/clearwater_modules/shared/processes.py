"""Process functions used by one or more modules"""
import numba


@numba.njit
def celsius_to_kelvin(tempc: float) -> float:
    return tempc + 273.16


@numba.njit
def kelvin_to_celsius(tempk: float) -> float:
    return tempk - 273.16

@numba.njit
def arrhenius_correction(
    water_temp_c: float,
    rc20: float,
    theta: float,
) -> float:
    """
    Computes an adjusted kinetics reaction rate coefficient for the specified water 
    temperature using the van't Hoff form of the Arrhenius equation

    Parameters
    ----------
    water_temp_c : float
        Water temperature in degrees Celsius
    rc20 : float
        Kinetics reaction (decay) coefficient at 20 degrees Celsius
    theta : float
        Temperature correction factor

    Returns
    ----------
    float
        Adjusted kinetics rate for the specified water temperature
    """
    return rc20 * theta**(water_temp_c - 20.0)

@numba.njit
def compute_depth(
    surface_area: float,
    volume: float
) -> float:
    """Compute depth of a computation cell

    Args:
        surface_area: state variable for surface area of computational cell provided by CWR engine
        volume: state variable for volume of computational cell provided by CWR engine
    """
    return volume / surface_area