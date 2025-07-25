import functools
from processes.base import Process
from datetime import datetime
from variables import VariableRegistry
import xarray as xr
import numpy as np
from utils import constants, conversions

from custom_types import ArrayLike

# Vapor pressure fitting parameters
# Defined in Brutsaert (1982) Evaporation into the Atmosphere, p42
__A0 = 6984.505294
__A1 = -188.903931
__A2 = 2.133357675
__A3 = -1.288580973E-2
__A4 = 4.393587233E-5
__A5 = -8.023923082E-8
__A6 = 6.136820929E-11


class Temperature(Process):
    """
    Temperature process.
    """
    variables = []

    def __init__(self, wind_a:float, wind_b:float, wind_c:float, diffusivity_ratio:float=1.0) -> None:
        self.wind_a = wind_a
        self.wind_b = wind_b
        self.wind_c = wind_c
        self.diffusivity_ratio = diffusivity_ratio

    def run(self, time_step:datetime, variables:VariableRegistry) -> None:
        """
        Run the temperature process.
        """

        #pull out variables from the registry
        water_temperature = variables.get("water_temperature").sel(time=time_step)
        surface_area = variables.get("surface_area")
        volume = variables.get("volume").sel(time=time_step)

        cloudiness = variables.get("cloudiness").sel(time=time_step)
        air_temperature = variables.get("air_temperature").sel(time=time_step)
        solar_flux = variables.get("solar_radiation").sel(time=time_step)
        wind_speed = variables.get("wind_speed").sel(time=time_step)
        
        #we do not want to perform temperature calculations on dry cells
        #water_temperature = water_temperature.where(volume > 0)
        
        #compute the new water temperature
        updated_water_temperature = self.temperature_change(water_temperature, surface_area, volume, cloudiness, air_temperature, solar_flux, wind_speed)

        #change the temperature in the registry
        #TODO: remove temporary increase to testing
        water_temperature *= 1.01
        #TODO: Ask Anthony/Sarah if there is a better way to update an array in place
        #water_temperature += updated_water_temperature.fillna(0)

    #TODO: implement
    def flux_solar(self): 
        return 0
        raise NotImplementedError()   

    def flux_atmospheric_longwave(self,water_temperature: ArrayLike) -> xr.DataArray:
        """
        Compute the atmospheric longwave flux in of the grid in (W/m^2)

        Parameters:
            water_temperature [xr.DataArray]: Water temperature in units of celsius
        Returns:
            flux_atmospheric_longwave [xr.DataArray]: Atmospheric longwave flux in units of W/m^2
        """
        return (
            -0.97  # upwelling is a negative flow from the water to the atmosphere
            * constants.STEFAN_BOLTZMANN_CONSTANT
            * conversions.celsius_to_kelvin(water_temperature) ** 4
        )


    def flux_upwelling_longwave(self,
        air_temperature: ArrayLike,
        cloudiness: ArrayLike,
    ) -> ArrayLike:
        """
        Compute the upwelling longwave flux in of the grid in (W/m^2)

        Parameters:
            air_temperature [xr.DataArray]: Air temperature in units of celsius
        Returns:
            flux_upwelling_longwave [xr.DataArray]: Upwelling longwave flux in units of W/m^2
        """

        return (
            9.37e-6  # note this was 0.937E-5 in original equation
            * (1 + 0.17 * cloudiness**2)
            * constants.STEFAN_BOLTZMANN_CONSTANT
            # This equation is for air temperature in Kelvin
            * conversions.celsius_to_kelvin(air_temperature) ** 6
        )


    def flux_latent_heat(self,
        atmospheric_pressure: ArrayLike,
        water_temperature: ArrayLike,
        wind_speed: ArrayLike,
        atmospheric_vapor_pressure: ArrayLike,
    ) -> xr.DataArray:
        """

        Parameters:
            atmospheric_pressure (ArrayLike): atmospheric pressure scaled or grid rectified in units of millibars
            water_temperature (ArrayLike): _description_
            wind_speed (ArrayLike): _description_
            wind_function (ArrayLike): _description_

        Returns:
            xr.DataArray: _description_
        """
        return (
            0.622
            / atmospheric_pressure
            * self.latent_heat_vaporization(water_temperature)  
            * self.density_water(water_temperature)  
            * self.wind_function(wind_speed)
            * (self.saturation_vapor_pressure(water_temperature) - atmospheric_vapor_pressure)  # need to implement saturation vapor pressure
        )


    def flux_sensible(self, water_temperature: ArrayLike, air_temperature: ArrayLike, wind_speed: ArrayLike) -> ArrayLike:
        """Compute the sensible heat flux in of the grid in (W/m^2)
        
        Sensible heat describes the flux of heat through molecular or turbulent
        transfer between the air and water surface
        """
        water_temperature_kelvin = conversions.celsius_to_kelvin(water_temperature)
        air_temperature_kelvin = conversions.celsius_to_kelvin(air_temperature)

        return (
            self.diffusivity_ratio * #diffusivity is a class parameter
            constants.SPECIFIC_HEAT_AIR *
            self.density_water(water_temperature) *
            self.wind_function(wind_speed)[1] *
            (air_temperature_kelvin - water_temperature_kelvin)
        )

    #TODO: implement
    def flux_sediment(self): 
        return 0
        raise NotImplementedError()

    def flux_net(self, water_temperature:ArrayLike, cloudiness:ArrayLike, air_temperature:ArrayLike, solar_flux:ArrayLike, wind_speed:ArrayLike,) -> ArrayLike:
        """
        Compute the net heatflux in of the grid in (W/m^2)
        """
        return (
            self.flux_sensible(water_temperature, air_temperature, wind_speed)
            + solar_flux #provided as direct input
            + self.flux_sediment()
            + self.flux_atmospheric_longwave(water_temperature)  # shouldn't the value handle the need for the negative?
            + self.flux_upwelling_longwave(cloudiness, air_temperature)  # upwelling flux is typically a loss of energy back to atmosphere
            - self.flux_latent_heat(water_temperature, atmospheric_pressure, water_temperature, wind_speed, atmospheric_vapor_pressure)
        )


    def water_specific_heat(self, temperature: ArrayLike) -> ArrayLike:
        """Approximate the heat capacity of water as a function of water temperature (Celsius)
        Parameters: 
            temperature - Water temperature in units of Celsius
        Returns: 
            DataArray/Float with value for the heat capacity of water in units of J/kg/K
        """
        return np.select(
            condlist=[
                temperature <= 0.0,
                temperature <= 5.0,
                temperature <= 10.0,
                temperature <= 15.0,
                temperature <= 20.0,
                temperature <= 25.0,
            ],
            choicelist=[
                4218.0,
                4202.0,
                4192.0,
                4186.0,
                4182.0,
                4180.0,
            ],
            default=4178.0,
        )

    def temperature_change(
        self,
        water_temperature: ArrayLike,
        surface_area: ArrayLike,
        volume: ArrayLike,
        cloudiness: ArrayLike,
        air_temperature: ArrayLike,
        solar_flux: ArrayLike,
        wind_speed: ArrayLike,
    ) -> ArrayLike:
        """
        Compute the change in temperature of the grid in (C)
        Parameters:
            water_temperature - Water temperature in units of Celsius
            surface_area - Surface area of the grid in units of m^2
            volume - Volume of the grid in units of m^3
        Returns:
            DataArray/Float with value for the change in temperature in units of Celsius
        """
        return (
            self.flux_net(water_temperature, cloudiness, air_temperature, solar_flux)  # computed though net flux
            * surface_area  # passed in from Riverine
            / (
                volume
                * self.density_water(water_temperature)
                * self.water_specific_heat(water_temperature)
            )
        )

    @functools.lru_cache(maxsize=1)
    def density_water(self, temperature: ArrayLike) -> ArrayLike:
        """Compute the density of water (kg/m3) as a function of water temperature (Celsius)
        Parameters:
            temperature - Water temperature in units of Celsius
        Returns:
            DataArray/Float with value for the density of water in units of kg/m3
        """
        return 999.973 * (
            1.0
            - (
                ((temperature - 3.9863) ** 2 * (temperature + 288.9414))
                / (508929.2 * (temperature + 68.12963))
            )
        )

    @functools.lru_cache(maxsize=1)
    def latent_heat_vaporization(
        self,
        water_temperature: ArrayLike
    ) -> ArrayLike:
        """
        Compute the latent heat of vaporization (J/kg) as a function of water temperature (Kelvin)
        Parameters:
            water_temperature - Water temperature in units of Celsius
        Returns:
            DataArray/Float with value for the latent heat of vaporization in units of J/kg
        """
        return 2499999 - 2385.74 * conversions.celsius_to_kelvin(water_temperature)
    
    @functools.lru_cache(maxsize=1)
    def saturation_vapor_pressure(self, water_temperature: ArrayLike) -> ArrayLike:
        """
        Compute the saturation vapor pressure (mb) as a function of water temperature (Kelvin)
        Parameters:
            temperature - Water temperature in units of Celsius
        Returns:
            DataArray/Float with value for the saturation vapor pressure in units of mb
        """
        water_temperature_kelvin = conversions.celsius_to_kelvin(water_temperature)
        return (
            __A0 +
            water_temperature_kelvin *
            (
                __A1 +
                water_temperature_kelvin *
                (
                    __A2 +
                    water_temperature_kelvin * (
                        __A3 +
                        water_temperature_kelvin *
                        (__A4 + water_temperature_kelvin * (__A5 + water_temperature_kelvin * __A6))
                    )
                )
            )
        )

    #TODO: this needs the richardson function
    @functools.lru_cache(maxsize=1)
    def wind_function(self, wind_speed: ArrayLike) -> ArrayLike:
        return self.wind_a/1_000_000 + self.wind_b/1_000_000 * wind_speed ** self.wind_c


    def richardson_number(self, wind_speed: ArrayLike, density_air_sat: ArrayLike, density_air: ArrayLike) -> tuple[float,float]:
        """
        Compute the Richardson Number. This is used in latent and sensible heat flux 
        computations to correct for atmospheric stability.
        
        Richardson Number:
            0.01 >= richardson_function        -> unstable
            0.01 <= richardson_function <  2   -> stable
            -0.01 <  richardson_function < 0.01 -> neutral
        
        Parameters
        ----------
        wind_speed : double
            Wind speed (m/s)
        density_air_sat : double
            Saturation density of air (kg/m3)
        density_air : double
            Density of air (kg/m3)
        
        Returns
        ----------
        list
            Richardson Number and Richardson Function
        """

        richardson_function: float = 0.0
        richardson_number: float = -1 * constants.GRAVITY * (density_air - density_air_sat) * 2.0 / (density_air * (wind_speed**2.0))

        # Set bounds
        if (richardson_number > 2.0):
            richardson_number = 2.0
        if (richardson_number < -1.0):
            richardson_number = -1.0

        if (richardson_number < 0.0):
            if (richardson_number >= - 0.01):
                # neutral
                richardson_function = 1.0
            else:
                # unstable
                richardson_function = (1.0 - 22.0 * richardson_number)**0.80
        else:
            if (richardson_number <= 0.01):
                # neutral
                richardson_function = 1.0
            else:
                # stable
                richardson_function = (1.0 + 34.0 * richardson_number)**(-0.80)
        return (richardson_number, richardson_function)