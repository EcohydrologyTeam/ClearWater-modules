from .base import Process, ProcessFactory
from datetime import datetime, timedelta
from clearwater_data.variables import VariableRegistry
import xarray as xr
import numpy as np
from clearwater_modules_v2.utils import constants, conversions

from clearwater_data.custom_types import ArrayLike

# References:
# https://erdc-library.erdc.dren.mil/server/api/core/bitstreams/81b728f8-87a7-4ef8-e053-411ac80adeb3/content


class Temperature(Process):
    """
    The temperature process simulates an energy balance to update water temperature.

    This is an implementation of the methodology outlined in Water Temperature Simulation Module
    developed by the U.S. Army Engineer Research and Development Center (ERDC) and presented in the
    Aquatic Nutrient Simulation Modules (NSMs) Developed for Hydrologic and Hydraulic Models report.
    """

    # Vapor pressure fitting parameters
    # Defined in Brutsaert (1982) Evaporation into the Atmosphere, p42
    __A0 = 6984.505294
    __A1 = -188.903931
    __A2 = 2.133357675
    __A3 = -1.288580973e-2
    __A4 = 4.393587233e-5
    __A5 = -8.023923082e-8
    __A6 = 6.136820929e-11

    variables = [
        "water_temperature",
        "wetted_surface_area",
        "volume",
        "cloudiness",
        "air_temperature",
        "solar_radiation",
        "wind_speed",
        "atmospheric_pressure",
        "atmospheric_vapor_pressure",
        "sediment_temperature",
        "sediment_thickness",
    ]

    def __init__(
        self,
        wind_a: float,
        wind_b: float,
        wind_c: float,
        sediment_density: ArrayLike = 1600.0,
        sediment_specific_heat: float = 1000.0,
        air_diffusivity_ratio: float = 1.0,
        sediment_diffusivity: float = 0.0061,
        time_step: timedelta = timedelta(minutes=5),
        use_sediment_temperature: bool = True,
    ) -> None:
        """
        Initialize the temperature process.

        Parameters:
            wind_a (float): Wind function parameter
            wind_b (float): Wind function parameter
            wind_c (float): Wind function parameter
            sediment_density (ArrayLike): Sediment density in units of kg/m^3
            sediment_specific_heat (float): Sediment specific heat in units of J/kg/C
            air_diffusivity_ratio (float): Air diffusivity ratio
            sediment_diffusivity (float): Sediment diffusivity in units of m^2/s
            time_step_frequency (timedelta): Time step frequency
        """
        # TODO: We should get Billy to both explain and provide guidance on wind parameters
        self.wind_a = wind_a
        self.wind_b = wind_b
        self.wind_c = wind_c
        self.sediment_density = sediment_density
        self.sediment_specific_heat = sediment_specific_heat
        self.air_diffusivity_ratio = air_diffusivity_ratio
        self.sediment_diffusivity = sediment_diffusivity
        self.use_sediment_temperature = use_sediment_temperature

        #V1 of the coupling had timestep 1 be skipped and started processing on timestep 2.
        self.__skip_first_time_step = True


        Process.__init__(self, time_step)

    @ProcessFactory.register("temperature")
    @staticmethod
    def from_config(config: dict, variable_registry: VariableRegistry) -> "Temperature":
        return Temperature(**config)

    def run(self, time: datetime, registry: VariableRegistry) -> None:
        """
        Run the temperature process.
        """

        if self.__skip_first_time_step:
            self.__skip_first_time_step = False
            return

        # pull out variables from the registry
        water_temperature = registry.get_at_time("water_temperature", time)
        surface_area = registry.get_at_time("wetted_surface_area", time)
        volume = registry.get_at_time("volume", time)

        cloudiness = registry.get_at_time("cloudiness", time)
        air_temperature = registry.get_at_time("air_temperature", time)
        solar_flux = registry.get_at_time("solar_radiation", time)
        wind_speed = registry.get_at_time("wind_speed", time)
        atmospheric_pressure = registry.get_at_time("atmospheric_pressure", time)
        atmospheric_vapor_pressure = registry.get_at_time(
            "atmospheric_vapor_pressure", time
        )

        # TODO: We should make the get method handle time selected
        # time independent .... for now
        sediment_temperature = registry.get_at_time("sediment_temperature", time)
        sediment_thickness = registry.get_at_time("sediment_thickness", time)

        # compute the new water temperature
        delta_water_temperature = self.temperature_change(
            water_temperature=water_temperature,
            surface_area=surface_area,
            volume=volume,
            cloudiness=cloudiness,
            air_temperature=air_temperature,
            solar_flux=solar_flux,
            wind_speed=wind_speed,
            sediment_temperature=sediment_temperature,
            sediment_thickness=sediment_thickness,
            atmospheric_pressure=atmospheric_pressure,
            atmospheric_vapor_pressure=atmospheric_vapor_pressure,
        )

        # we only want to update the temperature in cells that have water
        delta_water_temperature = xr.where(volume > 0, delta_water_temperature, 0)
        updated_water_temperature = water_temperature + delta_water_temperature

        # write changes back to registry
        registry.set_at_time("water_temperature", time, updated_water_temperature)

    #### Energy balance calculations ####

    def flux_upwelling_longwave(self, water_temperature: ArrayLike) -> xr.DataArray:
        """
        Compute the atmospheric longwave flux in of the grid in (W/m^2)

        Parameters:
            water_temperature [xr.DataArray]: Water temperature in units of celsius
        Returns:
            flux_upwelling_longwave [xr.DataArray]: Upwelling longwave flux in units of W/m^2
        """
        flux = (
            -0.97  # upwelling is a negative flow from the water to the atmosphere
            * constants.STEFAN_BOLTZMANN
            * conversions.celsius_to_kelvin(water_temperature) ** 4
        )
        return flux

    def flux_atmospheric_longwave(
        self,
        air_temperature: ArrayLike,
        cloudiness: ArrayLike,
    ) -> ArrayLike:
        """
        Compute the atmospheric longwave flux in of the grid in (W/m^2)

        Parameters:
            air_temperature [xr.DataArray]: Air temperature in units of celsius (C)
        Returns:
            flux_atmospheric_longwave [xr.DataArray]: Atmospheric longwave flux in units of W/m^2
        """

        # TODO: convert this to log statements, but we cannot assume these will be float convertible
        # print(f"    Longwave down terms:")
        # print(f"      cloudiness_term: {float(1.0 + 0.17 * cloudiness**2)}")
        ##print(f"        cloudiness_frac: {float(cloudiness)}")
        # print(
        #    f"      emissivity_air: {float(9.37e-6 * conversions.celsius_to_kelvin(air_temperature) ** 2)}"
        # )
        # print(f"      stefan_boltzmann: {float(constants.STEFAN_BOLTZMANN)}")
        # print(
        #    f"      air_temp_term: {float(conversions.celsius_to_kelvin(air_temperature) ** 4)}"
        # )
        # print(
        #    f"        air_temp_k: {float(conversions.celsius_to_kelvin(air_temperature))}"
        # )

        flux = (
            # TODO: Should this change as a function of temperature?
            # This is emissivity of air, which in our simply model is
            # a function of air temperature
            (
                9.37e-6 * conversions.celsius_to_kelvin(air_temperature) ** 2
            )  # note this was 0.937E-5 in original equation
            * (1 + 0.17 * cloudiness**2)
            * constants.STEFAN_BOLTZMANN
            # This equation is for air temperature in Kelvin
            # note the original equation air_temperature raised to the 6th power
            # but the STEFAN_BOLTZMANN constant is only to the 4th power.
            # So I believe this is a typo in the original equation as the deminsions don't work out otherwise.
            * conversions.celsius_to_kelvin(air_temperature) ** 4
        )
        return flux

    def flux_latent_heat(
        self,
        atmospheric_pressure: ArrayLike,
        water_temperature: ArrayLike,
        wind_speed: ArrayLike,
        atmospheric_vapor_pressure: ArrayLike,
        richardson_function: ArrayLike,
    ) -> xr.DataArray:
        """
        Returns latent heat flux in of the grid in units of (W/m^2)

        Parameters:
            atmospheric_pressure (ArrayLike): atmospheric pressure scaled or grid rectified in units of millibars
            water_temperature (ArrayLike): water temperature in units of celsius
            wind_speed (ArrayLike): wind speed in units of m/s
            atmospheric_vapor_pressure (ArrayLike): atmospheric vapor pressure in units of millibars

        Returns:
            xr.DataArray: latent heat flux in units of W/m^2
        """
        # TODO: we could consider keeping this for log statements, but we cannot assume they will be float convertible
        # print(f"    Latent heat terms:")
        # print(f"      atmospheric pressure: {float(atmospheric_pressure)}")
        # print(
        #    f"      latent_heat_vaporization: {float(self.latent_heat_vaporization(water_temperature))}"
        # )
        # print(f"        water_temperature: {float(water_temperature)}")
        # print(f"      water_density: {float(self.water_density(water_temperature))}")
        # print(
        #    f"      wind_function: {float(self.wind_function(wind_speed, richardson_function))}"
        # )
        # print(f"        wind_speed: {float(wind_speed)}")
        # print(
        #    f"      saturation_vapor_pressure: {float(self.saturation_vapor_pressure(water_temperature))}"
        # )
        # print(f"      atmospheric_vapor_pressure: {float(atmospheric_vapor_pressure)}")

        flux = (
            -0.622
            / atmospheric_pressure
            * self.latent_heat_vaporization(water_temperature)
            * self.water_density(water_temperature)
            * self.wind_function(wind_speed, richardson_function)
            * (
                self.saturation_vapor_pressure(water_temperature)
                - atmospheric_vapor_pressure
            )
        )
        return flux

    def flux_sensible(
        self,
        water_temperature: ArrayLike,
        air_temperature: ArrayLike,
        wind_speed: ArrayLike,
        richardson_function: ArrayLike,
    ) -> ArrayLike:
        """Compute the sensible heat flux in of the grid in (W/m^2)

        Sensible heat describes the flux of heat through molecular or turbulent
        transfer between the air and water surface
        """
        water_temperature_kelvin = conversions.celsius_to_kelvin(water_temperature)
        air_temperature_kelvin = conversions.celsius_to_kelvin(air_temperature)

        flux = (
            self.air_diffusivity_ratio
            * constants.AIR_SPECIFIC_HEAT
            * self.water_density(water_temperature)
            * self.wind_function(wind_speed, richardson_function)
            * (air_temperature_kelvin - water_temperature_kelvin)
        )
        return flux

    def flux_sediment(
        self,
        water_temperature: ArrayLike,
        sediment_temperature: ArrayLike,
        sediment_thickness: ArrayLike,
    ) -> ArrayLike:
        # optional flag to disable sediment temperature
        if not self.use_sediment_temperature:
            return 0.0

        flux = (
            self.sediment_density
            * self.sediment_specific_heat
            * self.sediment_diffusivity
            / 0.5  # TODO: determine why we need this 0.5
            / sediment_thickness
            * (sediment_temperature - water_temperature)
            / 86400.0  # convert days to seconds
        )
        return flux

    def flux_net(
        self,
        water_temperature: ArrayLike,
        cloudiness: ArrayLike,
        air_temperature: ArrayLike,
        solar_flux: ArrayLike,
        wind_speed: ArrayLike,
        atmospheric_pressure: ArrayLike,
        atmospheric_vapor_pressure: ArrayLike,
        sediment_temperature: ArrayLike,
        sediment_thickness: ArrayLike,
    ) -> ArrayLike:
        """
        Compute the net heatflux in of the grid in (W/m^2)
        """
        mixing_ratio_air = self.mixing_ratio_air(
            atmospheric_vapor_pressure, atmospheric_pressure
        )
        density_air = self.density_air(
            atmospheric_pressure, air_temperature, mixing_ratio_air
        )
        density_air_sat = self.density_air_sat(water_temperature, atmospheric_pressure)

        _, richardson_function = self.richardson_number(
            wind_speed,
            density_air_sat=density_air_sat,
            density_air=density_air,
        )
        sensible = self.flux_sensible(
            water_temperature, air_temperature, wind_speed, richardson_function
        )
        latent = self.flux_latent_heat(
            water_temperature=water_temperature,
            atmospheric_pressure=atmospheric_pressure,
            wind_speed=wind_speed,
            atmospheric_vapor_pressure=atmospheric_vapor_pressure,
            richardson_function=richardson_function,
        )
        sediment = self.flux_sediment(
            water_temperature, sediment_temperature, sediment_thickness
        )
        atmospheric = self.flux_atmospheric_longwave(air_temperature, cloudiness)
        upwelling = self.flux_upwelling_longwave(water_temperature)

        flux = (
            sensible
            + solar_flux  # provided as direct input
            + sediment
            + atmospheric
            + upwelling
            + latent
        )
        # print(f"    sensible: {float(sensible)}")
        # print(f"    solar: {float(solar_flux)}")
        # print(f"    sediment: {float(sediment)}")
        # print(f"    longwave: {float(longwave)}")
        # print(f"    upwelling: {float(upwelling)}")
        # print(f"    latent: {float(latent)}")
        # print(f"    net flux: {float(flux)}")
        return flux

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
        sediment_temperature: ArrayLike,
        sediment_thickness: ArrayLike,
        atmospheric_pressure: ArrayLike,
        atmospheric_vapor_pressure: ArrayLike,
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
        # flux * surface area = energy applied to the system
        # volume * density * specific heat = heat capacity of the system
        # energy / heat capacity = change in temperature
        return (
            self.flux_net(
                water_temperature=water_temperature,
                cloudiness=cloudiness,
                air_temperature=air_temperature,
                solar_flux=solar_flux,
                wind_speed=wind_speed,
                sediment_temperature=sediment_temperature,
                sediment_thickness=sediment_thickness,
                atmospheric_pressure=atmospheric_pressure,
                atmospheric_vapor_pressure=atmospheric_vapor_pressure,
            )
            * surface_area
            * self.time_step_seconds
            / (
                volume
                * self.water_density(water_temperature)
                * self.water_specific_heat(water_temperature)
            )
        )

    # @functools.lru_cache(maxsize=1)
    def water_density(self, temperature: ArrayLike) -> ArrayLike:
        """Compute the density of water (kg/m3) as a function of water temperature (Celsius)
        Parameters:
            temperature - Water temperature in units of Celsius
        Returns:
            DataArray/Float with value for the density of water in units of kg/m3
        """
        # TODO: verify if this equation is correct for both fresh and salt water
        return 999.973 * (
            1.0
            - (
                ((temperature - 3.9863) ** 2 * (temperature + 288.9414))
                / (508929.2 * (temperature + 68.12963))
            )
        )

    # @functools.lru_cache(maxsize=1)
    def latent_heat_vaporization(self, water_temperature: ArrayLike) -> ArrayLike:
        """
        Compute the latent heat of vaporization (J/kg) as a function of water temperature (Kelvin)
        Parameters:
            water_temperature - Water temperature in units of Celsius
        Returns:
            DataArray/Float with value for the latent heat of vaporization in units of J/kg
        """
        return 2499999 - 2385.74 * conversions.celsius_to_kelvin(water_temperature)

    # @functools.lru_cache(maxsize=1)
    def saturation_vapor_pressure(self, water_temperature: ArrayLike) -> ArrayLike:
        """
        Compute the saturation vapor pressure (mb) as a function of water temperature (celsius)
        Parameters:
            water_temperature - Water temperature in units of Celsius
        Returns:
            DataArray/Float with value for the saturation vapor pressure in units of mb
        """
        water_temperature_kelvin = conversions.celsius_to_kelvin(water_temperature)
        saturation_vapor_pressure = self.__A0 + water_temperature_kelvin * (
            self.__A1
            + water_temperature_kelvin
            * (
                self.__A2
                + water_temperature_kelvin
                * (
                    self.__A3
                    + water_temperature_kelvin
                    * (
                        self.__A4
                        + water_temperature_kelvin
                        * (self.__A5 + water_temperature_kelvin * self.__A6)
                    )
                )
            )
        )
        return saturation_vapor_pressure

    # TODO: this needs the richardson function
    # @functools.lru_cache(maxsize=1)
    def wind_function(
        self, wind_speed: ArrayLike, richardson_function: ArrayLike
    ) -> ArrayLike:
        # print(f"    Wind Function terms:")
        # print(f"      richardson_function: {float(richardson_function)}")
        # print(f"      wind_a: {float(self.wind_a)}")
        # print(f"      wind_b: {float(self.wind_b)}")
        # print(f"      wind_c: {float(self.wind_c)}")
        # print(f"      wind_speed: {float(wind_speed)}")

        return richardson_function * (
            (self.wind_a / 1_000_000)
            + (self.wind_b / 1_000_000) * (wind_speed**self.wind_c)
        )

    def mixing_ratio_air(
        self, atmospheric_vapor_pressure: ArrayLike, atmospheric_pressure: ArrayLike
    ) -> ArrayLike:
        """Calculate air mixing ratio (unitless).

        Args:
            atmospheric_vapor_pressure: Atmospheric vapour pressure of air (mb)
            atmospheric_pressure: Atmospheric pressure (mb)
        """
        # TODO: what if atmospheric_pressure == atmospheric_vapor_pressure?
        if atmospheric_pressure == atmospheric_vapor_pressure:
            return 0.0
        mixing_ratio = (
            0.622
            * atmospheric_vapor_pressure
            / (atmospheric_pressure - atmospheric_vapor_pressure)
        )
        return mixing_ratio

    def density_air(
        self,
        atmospheric_pressure: ArrayLike,
        air_temperature: ArrayLike,
        mixing_ratio_air: ArrayLike,
    ) -> ArrayLike:
        """Calculate air density (kg/m^3).

        Args:
            atmospheric_pressure: Atmospheric pressure (mb)
            air_temperature: Air temperature (Celsius)
            mixing_ratio_air: Air mixing ratio (unitless)
        """
        air_temperature_kelvin = conversions.celsius_to_kelvin(air_temperature)
        return (
            0.348
            * (atmospheric_pressure / air_temperature_kelvin)
            * (1.0 + mixing_ratio_air)
            / (1.0 + 1.61 * mixing_ratio_air)
        )

    def density_air_sat(
        self, water_temperature: ArrayLike, atmospheric_pressure: ArrayLike
    ) -> ArrayLike:
        """
        Compute the density of saturated air at water surface temperature.

        Parameters:
            water_temperature (float): Water temperature (Celsius)
            atmospheric_pressure (float): Atmospheric pressure (millibars)

        Returns:
            Density of saturated air at water surface temperature (kg/m3, float)
        """
        water_temperature_kelvin = conversions.celsius_to_kelvin(water_temperature)
        saturation_vapor_pressure = self.saturation_vapor_pressure(water_temperature)
        mixing_ratio_sat = (
            0.622
            * saturation_vapor_pressure
            / (atmospheric_pressure - saturation_vapor_pressure)
        )

        return (
            0.348
            * (atmospheric_pressure / water_temperature_kelvin)
            * (1.0 + mixing_ratio_sat)
            / (1.0 + 1.61 * mixing_ratio_sat)
        )

    def richardson_number(
        self, wind_speed: ArrayLike, density_air_sat: ArrayLike, density_air: ArrayLike
    ) -> tuple[ArrayLike, ArrayLike]:
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

        richardson_number: ArrayLike = (
            # -1 #TODO: check original equation to see if this multiplication by negative one is needed (not in v1 of code)
            constants.GRAVITY
            * (density_air - density_air_sat)
            * 2.0
            / (density_air * (wind_speed**2.0))
        )
        # print(f"    Richardson Number: {float(richardson_number)}")
        # print(f"      gravity: {float(constants.GRAVITY)}")
        # print(f"      density_air: {float(density_air)}")
        # print(f"      density_air_sat: {float(density_air_sat)}")
        # print(f"      wind_speed: {float(wind_speed)}")

        # Set bounds for richardson number
        # print(f"    Richardson Number before bounds: {richardson_number}")

        richardson_number = xr.where(richardson_number > 2.0, 2.0, richardson_number)
        richardson_number = xr.where(richardson_number < -1.0, -1.0, richardson_number)

        # Calculate richardson function
        # TODO: can we find a more efficient way to calculate this?
        # four where clauses is a little rough
        richardson_function: ArrayLike = 1.0

        ### CASE 1 ###
        # neutral rn < 0
        richardson_function = xr.where(
            (richardson_number < 0.0) & (richardson_number >= -0.01),
            1.0,
            richardson_function,
        )
        # if (richardson_number < 0.0) & (richardson_number >= -0.01):
        #    print('Case 1 == True')
        # else:
        #    print('Case 1 == False')

        ### CASE 2 ###
        # unstable
        richardson_function = xr.where(
            (richardson_number < 0.0) & (richardson_number < -0.01),
            (1.0 - 22.0 * richardson_number) ** 0.80,
            richardson_function,
        )
        # if (richardson_number < 0.0) & (richardson_number < -0.01):
        #    print('Case 2 == True')
        # else:
        #    print('Case 2 == False')

        ### CASE 3 ###
        # neutral rn > 0
        richardson_function = xr.where(
            (richardson_number >= 0.0) & (richardson_number <= 0.01),
            1.0,
            richardson_function,
        )
        # if (richardson_number >= 0.0) & (richardson_number <= 0.01):
        #    print('Case 3 == True')
        # else:
        #    print('Case 3 == False')

        ### CASE 4 ###
        # stable
        richardson_function = xr.where(
            (richardson_number >= 0.0) & (richardson_number > 0.01),
            (1.0 + 34.0 * richardson_number) ** (-0.80),
            richardson_function,
        )

        return (richardson_number, richardson_function)
