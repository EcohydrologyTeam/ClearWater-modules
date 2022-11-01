'''
=======================================================================================
Water Temperature Simulation Module (TSM)
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

Initial Version: April 9, 2021
'''

import logging
import os
import numba
from collections import namedtuple
import types
import typing

'''
Initialize class

Parameters:
    TwaterC (float):        Water temperature (degrees Celsius)
    surface_area (float):   Cell surface area (m2)
    volume (float):         Cell volume (m3)
'''


# Initialize debugging log file
LOG_LEVEL = logging.DEBUG
LOG_FORMAT = '%(levelname)-8s %(asctime)20s: %(message)s'
src_path = os.path.dirname(os.path.realpath(__file__))
log_filename = os.path.join(src_path, 'TSM.log')
logging.basicConfig(filename=log_filename, level=LOG_LEVEL,
                    format=LOG_FORMAT, filemode='w')
logger = logging.getLogger()

# Constants
# Stefan Boltzman constant (W/m2/K)
stefan_boltzmann = 5.67E-8

# Specific heat of air (J/kg/K)
Cp_air = 1005.0

# Emissivity of water (unitless)
emissivity_water = 0.97

# Gravitational acceleration (m/s2)
gravity = -9.806

# Vapor pressure fitting parameters
# Defined in Brutsaert (1982) Evaporation into the Atmosphere, p42
a0 = 6984.505294
a1 = -188.903931
a2 = 2.133357675
a3 = -1.288580973E-2
a4 = 4.393587233E-5
a5 = -8.023923082E-8
a6 = 6.136820929E-11

# Initial values
pb = 1600.0
Cps = 1673.0
h2 = 0.1
alphas = 0.0432
Richardson_option = True

# Initial dictionary that will hold the pathways to return from this function
pathways = {}


class Parameter:
    def __init__(self, value: typing.Any, name: str, full_name: str, units: str = "", 
            absolute_min: float = -1.0e3, absolute_max: float = 1.0e3, 
            expected_min: float = -1.0e3, expected_max: float = 1.0e3, 
            description: str = ""):
        self.value = value
        self.name = name
        self.full_name = full_name
        self.units = units
        self.absolute_min = absolute_min
        self.absolute_max = absolute_max
        self.expected_min = expected_min
        self.expected_max = expected_max
        self.description = description


def energy_budget_method(TwaterC: float, surface_area: float, volume: float, TairC: float, q_solar: float, pressure_mb: float, eair_mb: float, cloudiness: float, wind_speed: float, wind_a: float, wind_b: float, wind_c: float, wind_kh_kw: float, use_SedTemp: bool = False, TsedC: float = 0.0, num_iterations: int = 10, tolerance: float = 0.01, time_step: float = 60):
    '''
    Compute water temperature kinetics using the energy budget method

    Parameters:

        TwaterC (float):        Water temperature entering cell (degrees Celsius)
        surface_area (float):   Surface area of cell face (m2?)
        volume (float):         Volume of cell (m^3???)
          
        TairC (float):          Air temperature (degrees Celsius)
        q_solar (float):        Solar radiation (W/m2)
        pressure_mb (float):    Air pressure (mb)
        eair_mb (float):        Vapor pressure (mb)
        cloudiness (float):     Cloudiness (fraction)
        wind_speed (float):     Wind speed (m/s)
        wind_a (float):         "a" coefficient of the wind function
        wind_b (float):         "b" coefficient of the wind function
        wind_c (float):         "c" coefficient of the wind function
        wind_kh_kw (float):     Diffusivity ratio (unitless)
        use_SedTemp (bool):     Compute surface temperature (on/off)
        TsedC (float):          Sediment temperature (degrees Celsius)
        num_iterations (int)    Number of iterations to use in Newton-Raphson algorithm
        tolerance (float)       Tolerance for stopping iteration
        time_step (float):      Time step for one calculation of temperature change due to heat flux (s) 

        Returns:
            Dictionary of pathway variables
    '''

    # logger.debug(
    #     f'energy_budget_method({TairC:.2f}, {q_solar:.2f}, {pressure_mb:.2f}, {eair_mb:.2f}, {cloudiness:.2f}, {wind_speed:.2f}, {wind_a:.2f}, {wind_b:.2f}, {wind_c:.2f}, {wind_kh_kw:.2f}, use_SedTemp={use_SedTemp}, TsedC={TsedC:.2f}, num_iterations={num_iterations}, tolerance={tolerance})')

    # Temperature
    TwaterK: float = celsius_to_kelvin(TwaterC)
    TairK: float = celsius_to_kelvin(TairC)

    # Pressure
    # Saturated vapor pressure computed from water temperature at previous time step (mb)

    # ----------------------------------------------------------------------------------------------
    # Wind function, stability and flux partitioning
    # ----------------------------------------------------------------------------------------------

    # Mixing ratio (unitless)
    mixing_ratio_air = 0.622 * eair_mb / (pressure_mb - eair_mb)

    # Compute density of air (kg/m3)
    density_air = 0.348 * (pressure_mb / TairK) * (1.0 +
                                                   mixing_ratio_air) / (1.0 + 1.61 * mixing_ratio_air)

    # ----------------------------------------------------------------------------------------------
    #  Computations that are not a function of water temperature
    #  Note: solar radiation comes in directly from the interface
    # ----------------------------------------------------------------------------------------------

    # Downwelling (atmospheric) longwave radiation

    # Emissivity of air (unitless)
    emissivity_air = 0.00000937 * TairK**2.0

    # Atmospheric (downwelling) longwave radiation (W/m2)
    q_longwave_down = mf_q_longwave_down(
        TairK, emissivity_air, cloudiness, stefan_boltzmann)

    pathways["q_longwave_down"] = q_longwave_down

    # Wind function for latent and sensible heat (unitless)
    wind_function = wind_a / 1000000.0 + wind_b / 1000000.0 * wind_speed**wind_c

    # ____________________________________________________________________________________
    #  Energy balance computations that are functions of water temperature
    # ____________________________________________________________________________________
    # Get water temperature from previous time step plus 1/2 AD change??? DLL can not do it.

    # Physical values that are functions of water temperature

    # Latent heat of vaporization (J/kg)
    Lv = mf_latent_heat_vaporization(TwaterK)

    # Density (kg/m3)
    density_water = mf_density_water(TwaterC)

    # Specific heat of water (J/kg/K)
    Cp_water = mf_Cp_water(TwaterC)

    # Saturated vapor pressure computed from water temperature at previous time step (mb)
    esat_mb = mf_esat_mb(TwaterK)


    # ------------------------------------------------------------------------
    # Upwelling (back or water surface) longwave radiation (W/m2)
    q_longwave_up = mf_q_longwave_up(TwaterK)

    # ------------------------------------------------------------------------
    # Surface fluxes
    # ------------------------------------------------------------------------

    # ------------------------------------------------------------------------
    # Compute Richardson number and check stability

    # Richardson's number (unitless)
    Ri_No: float = 0.0

    # Richardson's stablility function (unitless)
    Ri_fxn: float = 1.0

    if (wind_speed > 0.0 and Richardson_option):
        # Density of air computed at water surface temperature (kg/m3)
        density_air_sat = mf_density_air_sat(TwaterK, esat_mb, pressure_mb)
        (Ri_No, Ri_fxn) = RichardsonNumber(
            wind_speed, density_air_sat, density_air)

    # ------------------------------------------------------------------------
    # Latent heat flux (W/m2)
    q_latent = Ri_fxn * (0.622 / pressure_mb) * Lv * \
        density_water * wind_function * (esat_mb - eair_mb)

    # ------------------------------------------------------------------------
    # Sensible heat flux
    q_sensible = wind_kh_kw * Ri_fxn * Cp_air * \
        density_water * wind_function * (TairK - TwaterK)

    # ------------------------------------------------------------------------
    # Compute sediment heat flux and temperature change
    q_sediment: float = 0.0
    dTsedCdt: float = 0.0
    if (use_SedTemp):
        q_sediment = pb * Cps * alphas / 0.5 / h2 * (TsedC - TwaterC) / 86400.0 # 86400 converts the sediment thermal diffusivity from units of m^2/d to m^2/s 
        dTsedCdt = alphas / (0.5 * h2 * h2) * (TwaterC - TsedC)

    # ------------------------------------------------------------------------
    # Net heat flux
    q_net = q_sensible - q_latent - q_longwave_up + \
        q_longwave_down + q_solar + q_sediment

    # ------------------------------------------------------------------------
    # Compute water temperature change
    dTwaterCdt = q_net * surface_area / \
        (volume * density_water * Cp_water) * time_step
    TwaterC += dTwaterCdt

    # ------------------------------------------------------------------------------------
    # Difference between air and water temperature (Celsius or Kelvins)
    Ta_Tw = TairC - TwaterC

    # Difference between saturated and current vapor pressure (mb)
    Esat_Eair = esat_mb - eair_mb

    # logger.debug('Computed energy balance values:')
    # logger.debug('q_net = %.2f' % q_net)
    # logger.debug('q_sensible = %.2f' % q_sensible)
    # logger.debug('q_latent = %.2f' % q_latent)
    # logger.debug('q_longwave_up = %.2f' % q_longwave_up)
    # logger.debug('q_longwave_down = %.2f' % q_longwave_down)
    # logger.debug('q_solar = %.2f' % q_solar)
    # logger.debug('q_sediment = %.2f' % q_sediment)
    # logger.debug('wind_function = %.2f' % wind_function)

    # Equilibrium temperature for current met conditions (C)
    set_pathways_float(q_net, 'q_net', 'Net Solar Radiation', "W/m2")
    # TODO: check why q_sensible was computed twice
    set_pathways_float(q_sensible, 'q_sensible', 'Sensible Radiation', 'W/m2')
    # TODO: check why q_sediment was computed twice
    set_pathways_float(q_sediment, 'q_sediment', 'Sediment Heat Flux', 'W/m2')
    set_pathways_float(q_latent, 'q_latent', 'Latent Heat', 'W/m2')
    set_pathways_float(q_longwave_up, 'q_longwave_up',
                       'Upwelling Longwave Radiation', 'W/m2')
    set_pathways_float(q_longwave_down, 'q_longwave_down',
                       'Downwelling Longwave Radiation', 'W/m2')
    set_pathways_float(
        Ta_Tw, 'Ta_Tw', 'Difference between Air and Water Temperature', 'degC')
    set_pathways_float(Esat_Eair, 'Esat_Eair',
                       'Difference between Saturation and Air Vapor Pressure', 'mb')
    set_pathways_float(Ri_No, 'Ri_No', 'Richardson Number', '')
    set_pathways_float(Ri_fxn, 'Ri_fxn', 'Richardson Function', '')
    set_pathways_float(dTsedCdt, 'dTsedCdt',
                       'Sediment Temperature Rate of Change', '')
    set_pathways_float(dTwaterCdt, 'dTwaterCdt',
                       'Water Temperature Rate of Change', 'degC')
    set_pathways_float(TwaterC, 'TwaterC', 'Water Temperature', 'degC')
    
    ##Check ri_no
    set_pathways_float(density_air, 'Density Air', 'Water Temperature', 'degC')
    set_pathways_float(density_air_sat, 'Density Sat', 'Water Temperature', 'degC')
    set_pathways_float(wind_speed, 'Wind Speed', 'Water Temperature', 'degC')
    
    return TwaterC


# Functions to set the pathways dictionary

def set_pathways_float(value: float, name: str, full_name: str, units: str = "", absolute_min: float = -1.0e3, absolute_max: float = 1.0e3, expected_min: float = -1.0e3, expected_max: float = 1.0e3, description: str = ""):
    '''
    Set the pathway for a floating point value
    '''
    pathways[name] = {"value": value, "name": name, "full_name": full_name, "units": units,
                      "absolute_min": absolute_min, "absolute_max": absolute_max, "expected_min": expected_min,
                      "expected_max": expected_max, "description": description}


def set_pathways_int(value: int, name: str, full_name: str, units: str = "", absolute_min: int = -1e3, absolute_max: int = 1e3, expected_min: int = -1e3, expected_max: int = 1e3, description: str = ""):
    '''
    Set the pathway for an integer value
    '''
    pathways[name] = {"value": value, "name": name, "full_name": full_name, "units": units,
                      "absolute_min": absolute_min, "absolute_max": absolute_max, "expected_min": expected_min,
                      "expected_max": expected_max, "description": description}


def set_pathways_bool(value: bool, name: str, full_name: str, description: str = ""):
    '''
    Set the pathway for a boolean value
    '''
    pathways[name] = {"value": value, "name": name,
                      "full_name": full_name, "description": description}


@numba.njit
def mf_d_esat_dT(TwaterK: float) -> float:
    '''
    Compute the derivative of function computing saturation vapor pressure 
    as a function of water temperature (Kelvin)

    Fitting parameters for vapor pressure:
    Brutsaert (1982) Evaporation into the Atmosphere, p42
    '''

    # logger.debug(f'mf_d_esat_dT({TwaterK})')

    return a1 + 2.0*a2*TwaterK + 3.0*a3*TwaterK**2.0 + 4.0*a4*TwaterK**3.0 + 5.0*a5*TwaterK**4.0 + 6.0*a6*TwaterK**5.0


# -----------------------------------------------------------------------------------
# Define functions to be used in the latent heat formulation
# -----------------------------------------------------------------------------------

@numba.njit
def mf_q_longwave_down(TairK: float, emissivity_air: float, cloudiness: float, stefan_boltzmann: float) -> float:
    '''
    Compute downwelling longwave radiation (W/m2)

    Parameters:
        TairK (float):              Air temperature (Kelvin)
        emissivity_air (float):     Emissivity of air (unitless)
        cloudiness (float):         Cloudiness (fraction)

    Returns:
        Downwelling longwave radiation (W/m2, float)
    '''

    # logger.debug(f'mf_q_longwave_down({TairK:.2f}, {emissivity_air:.2f}, {cloudiness:.2f})')

    return (1.0 + 0.17 * cloudiness**2) * emissivity_air * stefan_boltzmann * TairK**4.0


@numba.njit
def mf_q_longwave_up(TwaterK: float) -> float:
    '''
    Compute upwelling longwave radiation (W/m2) as a function of water temperature (Kelvin)
    '''

    # logger.debug(f'mf_q_longwave_up({TwaterK:.2f})')

    return emissivity_water * stefan_boltzmann * TwaterK**4.0


@numba.njit
def mf_esat_mb(TwaterK: float) -> float:
    '''
    Compute the saturation vapor pressure as a function of water temperature (Kelvin)

    Fitting parameters for vapor pressure are defined in:
    Brutsaert (1982) Evaporation into the Atmosphere, p42.
    '''

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
def RichardsonNumber(wind_speed: float, density_air_sat: float, density_air: float) -> list:
    '''
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
    '''

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
    '''
    Compute the latent heat of vaporization (W/m2) as a function of water temperature (Kelvin)
    '''

    # logger.debug(f'mf_latent_heat_vaporization({TwaterK:.2f})')

    return 2499999 - 2385.74 * TwaterK


@numba.njit
def mf_density_water(TwaterC: float) -> float:
    '''
    Compute density of water (kg/m3) as a function of water temperature (Celsius)
    '''

    # logger.debug(f'mf_density_water({TwaterC:.2f})')

    return 999.973 * (1.0 -
                      (((TwaterC - 3.9863) * (TwaterC - 3.9863) * (TwaterC + 288.9414)) /
                       (508929.2 * (TwaterC + 68.12963))))


@numba.njit
def mf_density_air_sat(TwaterK: float, esat_mb: float, pressure_mb: float) -> float:
    '''
    Compute the density of saturated air at water surface temperature.

    Parameters:
        TwaterK (float):        Water temperature (Kelvin)
        esat_mb (float):        Saturation vapor pressure in millibars
        pressure_mb (float):    Air pressure in millibars

    Returns:
        Density of saturated air at water surface temperature (kg/m3, float)
    '''

    # logger.debug(f'mf_density_air_sat({TwaterK:.2f}, {esat_mb:.2f}, {pressure_mb:.2f})')

    mixing_ratio_sat = 0.622 * esat_mb / (pressure_mb - esat_mb)
    return 0.348 * (pressure_mb / TwaterK) * (1.0 + mixing_ratio_sat) / (1.0 + 1.61 * mixing_ratio_sat)


@numba.njit
def mf_Cp_water(TwaterC: float) -> float:
    '''
    Compute the specific heat of water (J/kg/K) as a function of water temperature (Celsius).
    This is used in computing the source/sink term.
    '''

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


def print_pathways():
    '''
    Print the values of the pathway variables computed by TSM
    '''
    print('Pathways:')
    for key in pathways.keys():
        p = pathways[key]
        print("%20s%12.3f %s" % (p['name'], p['value'], p['units']))
