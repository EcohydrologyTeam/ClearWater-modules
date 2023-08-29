import numba
"""
=======================================================================================
ClearWater Modules: Water quality equations
=======================================================================================

Developed by:
* Dr. Todd E. Steissberg (ERDC-EL)
* Dr. Billy E. Johnson (ERDC-EL, LimnoTech)
* Dr. Zhonglong Zhang (Portland State University)
* Mr. Mark Jensen (HEC)

This module computes the water quality of a single computational cell. The algorithms 
and structure of this program were adapted from the Fortran 95 version of this module, 
developed by:
* Dr. Billy E. Johnson (ERDC-EL)
* Dr. Zhonglong Zhang (Portland State University)
* Mr. Mark Jensen (HEC)

Version 1.0

Initial Version: April 10, 2021
Last Revision Date: April 11, 2021
"""


def ArrheniusCorrection(TwaterC: float, rc20: float, theta: float):
    """
    Computes an adjusted kinetics reaction rate coefficient for the specified water 
    temperature using the van't Hoff form of the Arrhenius equation

    Parameters
    ----------
    TwaterC : float
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
    return rc20 * theta**(TwaterC - 20.0)


class TempCorrection:
    """
    Temperature correction class
    """

    def __init__(self, rc20: float, theta: float):
        self.rc20 = rc20
        self.theta = theta

    def arrhenius_correction(self, TwaterC: float):
        return ArrheniusCorrection(TwaterC, self.rc20, self.theta)

    # Functions to set the pathways dictionary


def set_pathways_float(pathways, value: float, name: str, full_name: str, units: str = "", absolute_min: float = -1.0e3, absolute_max: float = 1.0e3, expected_min: float = -1.0e3, expected_max: float = 1.0e3, description: str = ""):
    """
    Set the pathway for a floating point value
    """
    pathways[name] = {"value": value, "name": name, "full_name": full_name, "units": units,
                      "absolute_min": absolute_min, "absolute_max": absolute_max, "expected_min": expected_min,
                      "expected_max": expected_max, "description": description}


def set_pathways_int(pathways, value: int, name: str, full_name: str, units: str = "", absolute_min: int = -1e3, absolute_max: int = 1e3, expected_min: int = -1e3, expected_max: int = 1e3, description: str = ""):
    """
    Set the pathway for an integer value
    """
    pathways[name] = {"value": value, "name": name, "full_name": full_name, "units": units,
                      "absolute_min": absolute_min, "absolute_max": absolute_max, "expected_min": expected_min,
                      "expected_max": expected_max, "description": description}


def set_pathways_bool(pathways, value: bool, name: str, full_name: str, description: str = ""):
    """
    Set the pathway for a boolean value
    """
    pathways[name] = {"value": value, "name": name,
                      "full_name": full_name, "description": description}


@numba.njit
def mf_d_esat_dT(TwaterK: float) -> float:
    """
    Compute the derivative of function computing saturation vapor pressure 
    as a function of water temperature (Kelvin)

    Fitting parameters for vapor pressure:
    Brutsaert (1982) Evaporation into the Atmosphere, p42
    """

    # logger.debug(f'mf_d_esat_dT({TwaterK})')

    return a1 + 2.0*a2*TwaterK + 3.0*a3*TwaterK**2.0 + 4.0*a4*TwaterK**3.0 + 5.0*a5*TwaterK**4.0 + 6.0*a6*TwaterK**5.0


# -----------------------------------------------------------------------------------
# Define functions to be used in the latent heat formulation
# -----------------------------------------------------------------------------------

@numba.njit
def mf_q_longwave_down(TairK: float, emissivity_air: float, cloudiness: float, stefan_boltzmann: float) -> float:
    """
    Compute downwelling longwave radiation (W/m2)

    Parameters:
        TairK (float):              Air temperature (Kelvin)
        emissivity_air (float):     Emissivity of air (unitless)
        cloudiness (float):         Cloudiness (fraction)

    Returns:
        Downwelling longwave radiation (W/m2, float)
    """

    # logger.debug(f'mf_q_longwave_down({TairK:.2f}, {emissivity_air:.2f}, {cloudiness:.2f})')

    return (1.0 + 0.17 * cloudiness**2) * emissivity_air * stefan_boltzmann * TairK**4.0


@numba.njit
def mf_q_longwave_up(TwaterK: float, emissivity_water: float, stefan_boltzmann: float) -> float:
    """
    Compute upwelling longwave radiation (W/m2) as a function of water temperature (Kelvin)
    """

    # logger.debug(f'mf_q_longwave_up({TwaterK:.2f})')

    return emissivity_water * stefan_boltzmann * TwaterK**4.0


@numba.njit
def mf_esat_mb(TwaterK: float, a0: float, a1: float, a2: float, a3: float, a4: float, a5: float, a6: float) -> float:
    """
    Compute the saturation vapor pressure as a function of water temperature (Kelvin)

    Fitting parameters for vapor pressure are defined in:
    Brutsaert (1982) Evaporation into the Atmosphere, p42.
    """

    # logger.debug(f'mf_esat_mb({TwaterK:.2f})')

    return a0 + TwaterK*(a1 + TwaterK*(a2 + TwaterK * (a3 + TwaterK * (a4 + TwaterK*(a5 + TwaterK*a6)))))


# Temperature conversion functions

@numba.njit
def celsius_to_kelvin(tempc: float) -> float:
    return tempc + 273.16


@numba.njit
def kelvin_to_celsius(tempk: float) -> float:
    return tempk - 273.16


@numba.njit
def RichardsonNumber(wind_speed: float, density_air_sat: float, density_air: float, gravity: float) -> list:
    """
    Compute the Richardson Number. This is used in latent and sensible heat flux 
    computations to correct for atmospheric stability.

    Richardson Number:
        0.01 >= Ri_fxn        -> unstable
        0.01 <= Ri_fxn <  2   -> stable
        -0.01 <  Ri_fxn < 0.01 -> neutral

    Parameters:
        wind_speed (float):         Wind speed (m/s)
        density_air_sat (float):    Saturation density of air (kg/m3)
        density_air (float):        Density of air (kg/m3)

    Returns:
        Richardson Number and Richardson Function (list)
    """

    # logger.debug(f'RichardsonNumber({wind_speed:.2f}, {density_air_sat:.2f}, {density_air:.2f})')

    Ri_fxn: float = 0.0
    Ri_No: float = gravity * (density_air - density_air_sat) * \
        2.0 / (density_air * (wind_speed**2.0))

    # Set bounds
    if (Ri_No > 2.0):
        Ri_No = 2.0
    if (Ri_No < -1.0):
        Ri_No = -1.0

    if (Ri_No < 0.0):
        if (Ri_No >= - 0.01):
            # neutral
            Ri_fxn = 1.0
        else:
            # unstable
            Ri_fxn = (1.0 - 22.0 * Ri_No)**0.80
    else:
        if (Ri_No <= 0.01):
            # neutral
            Ri_fxn = 1.0
        else:
            # stable
            Ri_fxn = (1.0 + 34.0 * Ri_No)**(-0.80)
    return (Ri_No, Ri_fxn)


@numba.njit
def mf_latent_heat_vaporization(TwaterK: float) -> float:
    """
    Compute the latent heat of vaporization (W/m2) as a function of water temperature (Kelvin)
    """

    # logger.debug(f'mf_latent_heat_vaporization({TwaterK:.2f})')

    return 2499999 - 2385.74 * TwaterK


@numba.njit
def mf_density_water(TwaterC: float) -> float:
    """
    Compute density of water (kg/m3) as a function of water temperature (Celsius)
    """

    # logger.debug(f'mf_density_water({TwaterC:.2f})')

    return 999.973 * (1.0 -
                      (((TwaterC - 3.9863) * (TwaterC - 3.9863) * (TwaterC + 288.9414)) /
                       (508929.2 * (TwaterC + 68.12963))))


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
def mf_Cp_water(TwaterC: float) -> float:
    """
    Compute the specific heat of water (J/kg/K) as a function of water temperature (Celsius).
    This is used in computing the source/sink term.
    """

    # logger.debug(f'mf_esat_mb({TwaterC:.2f})')

    if TwaterC <= 0.0:
        Cp_water = 4218.0
    elif (TwaterC <= 5.0):
        Cp_water = 4202.0
    elif (TwaterC <= 10.0):
        Cp_water = 4192.0
    elif (TwaterC <= 15.0):
        Cp_water = 4186.0
    elif (TwaterC <= 20.0):
        Cp_water = 4182.0
    elif (TwaterC <= 25.0):
        Cp_water = 4180.0
    else:
        Cp_water = 4178.0
    return Cp_water


def print_pathways(pathways):
    """
    Print the values of the pathway variables computed by TSM
    """
    print('Pathways:')
    for key in pathways.keys():
        p = pathways[key]
        print("%20s%12.3f %s" % (p['name'], p['value'], p['units']))
