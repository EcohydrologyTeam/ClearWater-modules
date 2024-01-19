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


def kah_20(
    kah_20_user: xr.DataArray,
    hydraulic_reaeration_option: xr.DataArray,
    velocity: xr.DataArray,
    depth: xr.DataArray,
    flow: xr.DataArray,
    topwidth: xr.DataArray,
    slope: xr.DataArray,
    shear_velocity: xr.DataArray
) -> xr.DataArray:
    """Calculate hydraulic oxygen reaeration rate based on flow parameters in different cells

    Args:
        kah_20_user: User defined O2 reaeration rate at 20 degrees (1/d)
        hydraulic_reaeration_option: Integer value which selects method for computing O2 reaeration rate 
        velocity: Average water velocity in cell (m/s)
        depth: Average water depth in cell (m)
        flow: Average flow rate in cell (m3/s)
        topwidth: Average topwidth of cell (m)
        slope: Average slope of bottom surface 
        shear_velocity: Average shear velocity on bottom surface (m/s)
    """

    da: xr.DataArray = xr.where(hydraulic_reaeration_option == 1, kah_20_user,
                        xr.where(hydraulic_reaeration_option == 2, (3.93 * velocity**0.5) / (depth**1.5),
                        xr.where(hydraulic_reaeration_option == 3, (5.32 * velocity**0.67) / (depth**1.85),
                        xr.where(hydraulic_reaeration_option == 4, (5.026 * velocity) / (depth**1.67),
                        xr.where(hydraulic_reaeration_option == 5, xr.where(depth < 0.61, (5.32 * velocity**0.67) / (depth**1.85), xr.where(depth > 0.61, (3.93 * velocity**0.5) / (depth**1.5), (5.026 * velocity) / (depth**1.67))),
                        xr.where(hydraulic_reaeration_option == 6, xr.where(flow < 0.556, 517 * (velocity * slope)**0.524 * flow**-0.242, 596 * (velocity * slope)**0.528 * flow**-0.136),
                        xr.where(hydraulic_reaeration_option == 7, xr.where(flow < 0.556, 88 * (velocity * slope)**0.313 * depth**-0.353, 142 * (velocity * slope)**0.333 * depth**-0.66 * topwidth**-0.243),
                        xr.where(hydraulic_reaeration_option == 8, xr.where(flow < 0.425, 31183 * velocity * slope, 15308 * velocity * slope),
                        xr.where(hydraulic_reaeration_option == 9, 2.16 * (1 + 9 * (velocity / (9.81 * depth)**0.5)**0.25) * shear_velocity / depth, -9999
                                 )))))))))
    return da


@numba.njit
def kah_T(
    water_temp_c: xr.DataArray,
    kah_20: xr.DataArray,
    theta: xr.DataArray
) -> xr.DataArray:
    """Calculate the temperature adjusted hydraulic oxygen reaeration rate (/d)

    Args:
        water_temp_c: Water temperature in Celsius
        kah_20: Hydraulic oxygen reaeration rate at 20 degrees Celsius
        theta: Arrhenius coefficient
    """
    return arrhenius_correction(water_temp_c, kah_20, theta)


def kaw_20(
    kaw_20_user: xr.DataArray,
    wind_speed: xr.DataArray,
    wind_reaeration_option: xr.DataArray
) -> xr.DataArray:
    """Calculate the wind oxygen reaeration velocity (m/d) based on wind speed, r stands for regional

    Args:
        kaw_20_user: User defined wind oxygen reaeration velocity at 20 degrees C (m/d)
        wind_speed: Wind speed at 10 meters above the water surface (m/s)
        wind_reaeration_option: Integer value which selects method for computing wind oxygen reaeration velocity
    """
    Uw10 = wind_speed * (10 / 2)**0.143

    da: xr.DataArray = xr.where(wind_reaeration_option == 1, kaw_20_user,
                        xr.where(wind_reaeration_option == 2, 0.864 * Uw10,
                        xr.where(wind_reaeration_option == 3, xr.where(Uw10 <= 3.5, 0.2 * Uw10, 0.057 * Uw10**2),
                        xr.where(wind_reaeration_option == 4, 0.728 * Uw10**0.5 - 0.317 * Uw10 + 0.0372 * Uw10**2,
                        xr.where(wind_reaeration_option == 5, 0.0986 * Uw10**1.64,
                        xr.where(wind_reaeration_option == 6, 0.5 + 0.05 * Uw10**2,
                        xr.where(wind_reaeration_option == 7, xr.where(Uw10 <= 5.5, 0.362 * Uw10**0.5, 0.0277 * Uw10**2),
                        xr.where(wind_reaeration_option == 8, 0.64 + 0.128 * Uw10**2,
                        xr.where(wind_reaeration_option == 9, xr.where(Uw10 <= 4.1, 0.156 * Uw10**0.63, 0.0269 * Uw10**1.9),
                        xr.where(wind_reaeration_option == 10, 0.0276 * Uw10**2,
                        xr.where(wind_reaeration_option == 11, 0.0432 * Uw10**2,
                        xr.where(wind_reaeration_option == 12, 0.319 * Uw10,
                        xr.where(wind_reaeration_option == 13, xr.where(Uw10 < 1.6, 0.398, 0.155 * Uw10**2), -9999
                                 )))))))))))))
    
    return da


@numba.njit
def kaw_T(
    water_temp_c: xr.DataArray,
    kaw_20: xr.DataArray,
    theta: xr.DataArray
) -> xr.DataArray:
    """Calculate the temperature adjusted wind oxygen reaeration velocity (m/d)

    Args:
        water_temp_c: Water temperature in Celsius
        kaw_20: Wind oxygen reaeration velocity at 20 degrees Celsius
        theta: Arrhenius coefficient
    """
    return arrhenius_correction(water_temp_c, kaw_20, theta)


@numba.njit
def ka_T(
    kah_T: xr.DataArray,
    kaw_T: xr.DataArray,
    depth: xr.DataArray
) -> xr.DataArray:
    """Compute the oxygen reaeration rate, adjusted for temperature (1/d)

    Args:
        kah_T: Oxygen reaeration rate adjusted for temperature (1/d)
        kaw_T: Wind oxygen reaeration velocity adjusted for temperature (m/d)
        depth: Average water depth in cell (m)
    """
    return kaw_T / depth + kah_T

def SOD_tc(
    SOD_20: xr.DataArray,
    t_water_C: xr.DataArray,
    theta: xr.DataArray,
    DOX: xr.DataArray,
    KsSOD: xr.DataArray,
    use_DOX: xr.DataArray
) -> xr.DataArray:
    """Compute the sediment oxygen demand corrected by temperature and dissolved oxygen concentration

    Args:
        SOD_20: Sediment oxygen demand at 20 degrees celsius (mg-O2/m2)
        t_water_C: Water temperature in degrees C
        theta: Arrhenius coefficient
        use_DOX: Option to consider DOX concentration in water in calculation of sediment oxygen demand
    """
    SOD_tc = arrhenius_correction(t_water_C, SOD_20, theta)

    da: xr.DataArray = xr.where(use_DOX == True, SOD_tc * DOX / (DOX + KsSOD), SOD_tc)

    return da