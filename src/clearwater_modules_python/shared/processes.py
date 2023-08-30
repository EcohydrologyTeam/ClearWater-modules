"""Process functions used by one or more modules"""
import numba


@numba.njit
def celsius_to_kelvin(tempc: float) -> float:
    return tempc + 273.16


@numba.njit
def kelvin_to_celsius(tempk: float) -> float:
    return tempk - 273.16


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
def mf_d_esat_dT(
    TwaterK: float,
    a1: float,
    a2: float,
    a3: float,
    a4: float,
    a5: float,
    a6: float,
) -> float:
    """
    Compute the derivative of function computing saturation vapor pressure 
    as a function of water temperature (Kelvin)

    Fitting parameters for vapor pressure:
    Brutsaert (1982) Evaporation into the Atmosphere, p42
    """
    return (
        a1 +
        (2.0 * a2 * TwaterK) +
        (3.0 * a3 * TwaterK**2.0) +
        (4.0 * a4 * TwaterK**3.0) +
        (5.0 * a5 * TwaterK**4.0) +
        (6.0 * a6 * TwaterK**5.0)
    )


# -----------------------------------------------------------------------------------
# Define functions to be used in the latent heat formulation
# -----------------------------------------------------------------------------------

@numba.njit
def mf_q_longwave_down(
    TairK: float,
    emissivity_air: float,
    cloudiness: float,
    stefan_boltzmann: float,
) -> float:
    """
    Compute downwelling longwave radiation (W/m2)

    Parameters:
        TairK (float):              Air temperature (Kelvin)
        emissivity_air (float):     Emissivity of air (unitless)
        cloudiness (float):         Cloudiness (fraction)

    Returns:
        Downwelling longwave radiation (W/m2, float)
    """

    return (1.0 + 0.17 * cloudiness**2) * emissivity_air * stefan_boltzmann * TairK**4.0


@numba.njit
def mf_q_longwave_up(
    TwaterK: float,
    emissivity_water: float,
    stefan_boltzmann: float,
) -> float:
    """
    Compute upwelling longwave radiation (W/m2) as a function of water temperature (Kelvin)
    """

    # logger.debug(f'mf_q_longwave_up({TwaterK:.2f})')

    return emissivity_water * stefan_boltzmann * TwaterK**4.0


@numba.njit
def mf_esat_mb(
    TwaterK: float,
    a0: float,
    a1: float,
    a2: float,
    a3: float,
    a4: float,
    a5: float,
    a6: float,
) -> float:
    """
    Compute the saturation vapor pressure as a function of water temperature (Kelvin)

    Fitting parameters for vapor pressure are defined in:
    Brutsaert (1982) Evaporation into the Atmosphere, p42.
    """

    return (
        a0 +
        TwaterK *
        (
            a1 +
            TwaterK *
            (
                a2 +
                TwaterK * (
                    a3 +
                    TwaterK *
                    (a4 + TwaterK * (a5 + TwaterK * a6))
                )
            )
        )
    )


# Temperature conversion functions

@numba.njit
def ri_number(
    gravity: float,
    density_air: float,
    density_air_sat: float,
    wind_speed: float,
) -> float:
    """Calculates the Richardson Number.

    Args:
        gravity: Gravity (m/s2)
        density_air: Density of air (kg/m3)
        density_air_sat: Saturation density of air (kg/m3)
        wind_speed: Wind speed (m/s)
    """
    return (
        gravity *
        (density_air - density_air_sat) *
        2.0 / (density_air * (wind_speed**2.0))
    )


@numba.njit
def ri_function(ri_number: float) -> float:
    """Calculates the Richardson Function from the Richardson Number.
    Richardson Number:
        Unstable: 0.01 >= ri_function
        Stable: 0.01 <= ri_function < 2
        Neutral: -0.01 <  ri_function < 0.01
    """

    # Set bounds
    if ri_number > 2.0:
        ri_number = 2.0
    elif ri_number < -1.0:
        ri_number = -1.0

    if ri_number < 0.0:
        if ri_number >= - 0.01:
            # neutral
            return 1.0
        else:
            # unstable
            return (1.0 - 22.0 * ri_number)**0.80
    else:
        if ri_number <= 0.01:
            # neutral
            return 1.0
        else:
            # stable
            return (1.0 + 34.0 * ri_number)**(-0.80)


@numba.njit
def mf_latent_heat_vaporization(TwaterK: float) -> float:
    """
    Compute the latent heat of vaporization (W/m2) as a function of water temperature (Kelvin)
    """

    return 2499999 - 2385.74 * TwaterK


@numba.njit
def mf_density_water(water_temp_c: float) -> float:
    """
    Compute density of water (kg/m3) as a function of water temperature (Celsius)
    """

    # logger.debug(f'mf_density_water({water_temp_c:.2f})')

    return (
        999.973 *
        (1.0 - (
            (
                (water_temp_c - 3.9863) *
                (water_temp_c - 3.9863) *
                (water_temp_c + 288.9414)
            ) /
            (
                508929.2 *
                (water_temp_c + 68.12963)
            )
        ))
    )


@numba.njit
def mf_density_air_sat(TwaterK: float, esat_mb: float, pressure_mb: float) -> float:
    """
    Compute the density of saturated air at water surface temperature.

    Parameters:
        TwaterK (float):        Water temperature (Kelvin)
        esat_mb (float):        Saturation vapor pressure in millibars
        pressure_mb (float):    Air pressure in millibars

    Returns:
        Density of saturated air at water surface temperature (kg/m3, float)
    """

    # logger.debug(f'mf_density_air_sat({TwaterK:.2f}, {esat_mb:.2f}, {pressure_mb:.2f})')

    mixing_ratio_sat = 0.622 * esat_mb / (pressure_mb - esat_mb)
    return 0.348 * (pressure_mb / TwaterK) * (1.0 + mixing_ratio_sat) / (1.0 + 1.61 * mixing_ratio_sat)


@numba.njit
def mf_cp_water(water_temp_c: float) -> float:
    """
    Compute the specific heat of water (J/kg/K) as a function of water temperature (Celsius).
    This is used in computing the source/sink term.
    """

    if water_temp_c <= 0.0:
        return 4218.0
    elif water_temp_c <= 5.0:
        return 4202.0
    elif water_temp_c <= 10.0:
        return 4192.0
    elif water_temp_c <= 15.0:
        return 4186.0
    elif water_temp_c <= 20.0:
        return 4182.0
    elif water_temp_c <= 25.0:
        return 4180.0
    else:
        return 4178.0
