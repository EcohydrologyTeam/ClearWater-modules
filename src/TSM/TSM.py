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
# %%

import logging
import os
import sys
import numba
from collections import namedtuple, OrderedDict
import types
import typing

src_path = os.path.dirname(
            os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(src_path))
from src import water_quality_functions as wqf

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

class Temperature :
    
    def __init__(self, global_module_choices, global_vars, met_constant_changes, temperature_constant_changes):
        
        self.global_module_choices = global_module_choices
        self.global_vars = global_vars
        self.met_constant_changes = met_constant_changes
        self.temperature_constant_changes = temperature_constant_changes
        self.pathways = OrderedDict()

        self.met_constant = OrderedDict()
        self.met_constant = {
            'TairC' : 20,
            'q_solar' : 400,
            'TsedC' : 5.0,
            'eair_mb' : 1.0,
            'pressure_mb' : 1013.0,
            'cloudiness' : 0.1,
            'wind_speed' : 3.0,
            'wind_a' : 0.3,
            'wind_b' : 1.5,
            'wind_c' : 1.0,
            'wind_kh_kw' : 1.0

        }
        
        for key in self.met_constant_changes.keys() :
            if key in self.met_constant:
                self.met_constant[key] = self.met_constant_changes[key]

        self.temperature_constant = OrderedDict()
        self.temperature_constant = {
            'stefan_boltzmann' : 5.67E-8,
            'Cp_air' : 1005.0,
            'emissivity_water' : 0.97,
            'gravity' : -9.806,
            'a0' : 6984.505294,
            'a1' : -188.903931,
            'a2' : 2.133357675,
            'a3' : -1.288580973E-2,
            'a4' : 4.393587233E-5,
            'a5' : -8.023923082E-8,
            'a6' : 6.136820929E-11,
            'pb' : 1600.0,
            'Cps' : 1673.0,
            'h2' : 0.1,
            'alphas' : 0.0432,
            'Richardson_option' : True

        }

        for key in self.temperature_constant_changes.keys() :
            if key in self.temperature_constant:
                self.temperature_constant[key] = self.temperature_constant_changes[key]
    
    def energy_budget_method(self):
        '''
        Compute water temperature kinetics using the energy budget method

        Parameters:
            
            Global module choices (T/F)
            use_SedTemp (bool):         Compute surface temperature (on/off)

            Global Variables
            TwaterC (float):            Water temperature entering cell (degrees C)
            surface_area (float):       Surface area of cell face (m^2?)
            volume (float):             Volume of cell (m^3???)
            
            Meteorological Constants
            TairC (float):              Air temperature (degrees Celsius)
            q_solar (float):            Solar radiation (W/m^2)
            TsedC (float):              Sediment temperature (degrees Celsius)
            pressure_mb (float):        Air pressure (mb)
            eair_mb (float):            Vapor pressure (mb)
            cloudiness (float):         Cloudiness (%)
            wind_speed (float):         Wind speed (m/s)
            wind_a (float):             "a" coefficient of the wind function 
            wind_b (float):             "b" coefficient of the wind function
            wind_c (float):             "c" coefficient of the wind function
            wind_kh_kw (float):         Diffusivity ratio (unitless)
            
            Temperature Constants
            stefan_boltzmann (float):   Constant relating emitted radiation to temperature of matter (W/(m^2*K^4))             
            Cp_air (float):             Specific heat capacity of air (J/kg*C)
            emissivity_water (float):   Emissivity of water (unitless)
            gravity (float):            Acceleration of gravity (m/s^2)
            a0 (float):                 Saturation vapor pressure constant 0 (mb)
            a1 (float):                 Saturation vapor pressure constant 1 
            a2 (float):                 Saturation vapor pressure constant 2
            a3 (float):                 Saturation vapor pressure constant 3
            a4 (float):                 Saturation vapor pressure constant 4
            a5 (float):                 Saturation vapor pressure constant 5
            a6 (float):                 Saturation vapor pressure constant 6
            pb (float):                 Sediment bulk density (kg/m^3)
            Cps (float):                Specific heat capacity of sediments (J/Kg*C)
            h2 (float):                 Sediment active layer thickness (m)
            alphas (float):             Sediment thermal diffusivity (m^2/s)
            Richardson_option (bool):   Richardson option

            Not Currently Used
            num_iterations (int)    Number of iterations to use in Newton-Raphson algorithm
            tolerance (float)       Tolerance for stopping iteration
            time_step (float):      Time step for one calculation of temperature change due to heat flux (s) 

            Returns:
                Dictionary of pathway variables
        '''

        # logger.debug(
        #     f'energy_budget_method({TairC:.2f}, {q_solar:.2f}, {pressure_mb:.2f}, {eair_mb:.2f}, {cloudiness:.2f}, {wind_speed:.2f}, {wind_a:.2f}, {wind_b:.2f}, {wind_c:.2f}, {wind_kh_kw:.2f}, use_SedTemp={use_SedTemp}, TsedC={TsedC:.2f}, num_iterations={num_iterations}, tolerance={tolerance})')

        TwaterC = self.global_vars['TwaterC']
        surface_area = self.global_vars['surface_area']
        volume = self.global_vars['volume']

        pathways=OrderedDict()

        TairC = self.met_constant['TairC']
        q_solar = self.met_constant['q_solar']
        TsedC = self.met_constant['TsedC']
        pressure_mb = self.met_constant['pressure_mb']
        eair_mb = self.met_constant['eair_mb']
        cloudiness = self.met_constant['cloudiness']
        wind_speed = self.met_constant['wind_speed']
        wind_a = self.met_constant['wind_a']
        wind_b = self.met_constant['wind_b']
        wind_c = self.met_constant['wind_c']
        wind_kh_kw = self.met_constant['wind_kh_kw']

        
        stefan_boltzmann = self.temperature_constant['stefan_boltzmann']
        Cp_air = self.temperature_constant['Cp_air']
        emissivity_water = self.temperature_constant['emissivity_water']
        gravity = self.temperature_constant['gravity']
        a0 = self.temperature_constant['a0']
        a1 = self.temperature_constant['a1']
        a2 = self.temperature_constant['a2']
        a3 = self.temperature_constant['a3']
        a4 = self.temperature_constant['a4']
        a5 = self.temperature_constant['a5']
        a6 = self.temperature_constant['a6']
        pb = self.temperature_constant['pb']
        Cps = self.temperature_constant['Cps']
        h2 = self.temperature_constant['h2']
        alphas = self.temperature_constant['alphas']
        Richardson_option = self.temperature_constant['Richardson_option']
        

        # Temperature
        TwaterK: float = wqf.celsius_to_kelvin(TwaterC)
        TairK: float = wqf.celsius_to_kelvin(TairC)

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
        q_longwave_down = wqf.mf_q_longwave_down(
            TairK, emissivity_air, cloudiness, stefan_boltzmann)

        #pathways["q_longwave_down"] = q_longwave_down

        # Wind function for latent and sensible heat (unitless)
        wind_function = wind_a / 1000000.0 + wind_b / 1000000.0 * wind_speed**wind_c

        # ____________________________________________________________________________________
        #  Energy balance computations that are functions of water temperature
        # ____________________________________________________________________________________
        # Get water temperature from previous time step plus 1/2 AD change??? DLL can not do it.

        # Physical values that are functions of water temperature

        # Latent heat of vaporization (J/kg)
        Lv = wqf.mf_latent_heat_vaporization(TwaterK)

        # Density (kg/m3)
        density_water = wqf.mf_density_water(TwaterC)

        # Specific heat of water (J/kg/K)
        Cp_water = wqf.mf_Cp_water(TwaterC)

        # Saturated vapor pressure computed from water temperature at previous time step (mb)
        esat_mb = wqf.mf_esat_mb(TwaterK, a0, a1, a2, a3, a4, a5, a6)


        # ------------------------------------------------------------------------
        # Upwelling (back or water surface) longwave radiation (W/m2)
        q_longwave_up = wqf.mf_q_longwave_up(TwaterK, emissivity_water, stefan_boltzmann)

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
            density_air_sat = wqf.mf_density_air_sat(TwaterK, esat_mb, pressure_mb)
            (Ri_No, Ri_fxn) = wqf.RichardsonNumber(
                wind_speed, density_air_sat, density_air, gravity)

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
        if (self.global_module_choices['use_SedTemp']):
            q_sediment = pb * Cps * alphas / 0.5 / h2 * (TsedC - TwaterC) / 86400.0 # 86400 converts the sediment thermal diffusivity from units of m^2/d to m^2/s 
            dTsedCdt = alphas / (0.5 * h2 * h2) * (TwaterC - TsedC)

        # ------------------------------------------------------------------------
        # Net heat flux
        q_net = q_sensible - q_latent - q_longwave_up + \
            q_longwave_down + q_solar + q_sediment

        # ------------------------------------------------------------------------
        # Compute water temperature change
        dTwaterCdt = q_net * surface_area / \
            (volume * density_water * Cp_water) 
        #TwaterC += dTwaterCdt

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
        wqf.set_pathways_float(pathways, q_net, 'q_net', 'Net Solar Radiation', "W/m2")
        # TODO: check why q_sensible was computed twice
        wqf.set_pathways_float(pathways, q_sensible, 'q_sensible', 'Sensible Radiation', 'W/m2')
        # TODO: check why q_sediment was computed twice
        wqf.set_pathways_float(pathways, q_sediment, 'q_sediment', 'Sediment Heat Flux', 'W/m2')
        wqf.set_pathways_float(pathways, q_latent, 'q_latent', 'Latent Heat', 'W/m2')
        wqf.set_pathways_float(pathways, q_longwave_up, 'q_longwave_up',
                        'Upwelling Longwave Radiation', 'W/m2')
        wqf.set_pathways_float(pathways, q_longwave_down, 'q_longwave_down',
                        'Downwelling Longwave Radiation', 'W/m2')
        wqf.set_pathways_float(
            pathways, Ta_Tw, 'Ta_Tw', 'Difference between Air and Water Temperature', 'degC')
        wqf.set_pathways_float(pathways, Esat_Eair, 'Esat_Eair',
                        'Difference between Saturation and Air Vapor Pressure', 'mb')
        wqf.set_pathways_float(pathways, Ri_No, 'Ri_No', 'Richardson Number', '')
        wqf.set_pathways_float(pathways, Ri_fxn, 'Ri_fxn', 'Richardson Function', '')
        wqf.set_pathways_float(pathways, dTsedCdt, 'dTsedCdt',
                        'Sediment Temperature Rate of Change', '')
        wqf.set_pathways_float(pathways, dTwaterCdt, 'dTwaterCdt',
                        'Water Temperature Rate of Change', 'degC')
        wqf.set_pathways_float(pathways, TwaterC, 'TwaterC', 'Water Temperature', 'degC')
        
        ##Check ri_no
        wqf.set_pathways_float(pathways, density_air, 'Density Air', 'Water Temperature', 'degC')
        if Richardson_option == True:
            wqf.set_pathways_float(pathways, density_air_sat, 'Density Sat', 'Water Temperature', 'degC')

        wqf.set_pathways_float(pathways, wind_speed, 'Wind Speed', 'Water Temperature', 'degC')
        
        return dTwaterCdt


# %%
