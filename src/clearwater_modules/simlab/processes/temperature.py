import xsimlab as xs
from .riverine import Riverine
import xarray as xr
import numpy as np
import constants
import conversions

ArrayLike = xr.DataArray | np.ndarray | float


@xs.process
class Temperature:
    """
    This process is a placeholder for temperature-related calculations.
    It is intended to replace the TSM module of the old ClearwaterModules.
    """

    volume = xs.foreign(
        other_process_cls=Riverine,
        var_name="volume",
        intent="in",
    )
    wetted_surface_area = xs.foreign(
        other_process_cls=Riverine,
        var_name="wetted_surface_area",
        intent="in",
    )
    water_temperature = xs.global_ref(
        name="water_temperature",
        intent="in",
    )
    #water_temperature = xs.global_ref(
    #    other_process_cls=Riverine,
    #    var_name="temperature",
    #    intent="out", #also used as an input
    #)
    #this is what should be calculated by this process via energy balancre
    #temperature = xs.variable(
    #    dims=("time", "nface"),
    #    intent="out",
    #    description="Temperature of the riverine mesh at each time step and face",
    #)
    air_temperature = xs.variable(
        dims=("time"),
        intent="in",
    )

    #air_pressure = xs.variable(
    #    dims=("time","air_pressure"),
    #    intent="in",
    #)
    #cloudiness = xs.variable(
    #    dims=("time","percent_cloud_cover"),
    #    intent="in",
    #)
    solar_radiation = xs.variable(
        dims=("time"),
        intent="in",
    )

    #   wind_speed = xs.variable(
    #       dims=("time","wind_speed"),
    #       intent="in",
    #   )
    #   #This one might be something that is calculated
    #   wind_function = xs.variable(
    #       dims=("time","wind_function"),
    #       intent="in",
    #   )
    
    """Inputs we need 

    cell temperature - Riverine
    cell volume - Riverine
    cell surface area (wetted) - Riverine
    solar radiation ????
    air temperature ????
    air pressure ???? - i'm not sure about this one
    """


    def run_step(self, dt):
        # Placeholder for temperature calculations
        pass

        # Goals
        # We want replicate the temperature processes (TSM) to be a vectorized calculation
        # It should take a subset (slice) of xarray.DataArray and return an
        # xarray.DataArray of the same size and dimensions

        # looking at the TSM module
        # the new water temperature is the previous water temperature plus the change in water temperature
        # the change in water temperature is a function of the net flux times the surface area / volumetric heat capacity heat of the water

        """
        

        Technical References:
            Much of the temperature module's methodology is outlined in this report
                https://erdc-library.erdc.dren.mil/server/api/core/bitstreams/81b728f8-87a7-4ef8-e053-411ac80adeb3/content
        """
        pass


    def flux_solar(self): 
        raise NotImplementedError()   


    def flux_atmospheric_longwave(self,water_temperature: ArrayLike) -> xr.DataArray:
        """
        Compute the atmospheric longwave flux in of the grid in (W/m^2)

        Inputs: water_temperature [xr.DataArray] Water temperature in units of celsius
        Outputs: flux_atmospheric_longwave [xr.DataArray] Atmospheric longwave flux in units of W/m^2
        """
        return (
            -0.97  # upwelling is a negative flow from the water to the atmosphere
            * constants.STEFAN_BOLTZMANN_CONSTANT
            * conversions.celcius_to_kelvin(water_temperature) ** 4
        )


    def flux_upwelling_longwave(self,
        air_temperature: xr.DataArray,
        cloudiness: ArrayLike,
    ) -> xr.DataArray:
        """
        Compute the upwelling longwave flux in of the grid in (W/m^2)

        Inputs: air_temperature [xr.DataArray] Air temperature in units of celsius
        Outputs: flux_upwelling_longwave [xr.DataArray] Upwelling longwave flux in units of W/m^2
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
        wind_function: ArrayLike,  ##not sure we need this as an argument
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
            * latent_heat_of_varporization(water_temperature)  # need to implement
            * water_density(water_temperature)  # also need to implement
            * wind_function  # need to implement
            * (esat(water_temperature) - eair(water_temperature))  # need to implement both
        )


    def flux_sensible(self):
        raise NotImplementedError()

    def flux_sediment(self): 
        raise NotImplementedError()

    def flux_net(self):
        """
        Compute the net heatflux in of the grid in (W/m^2)
        """
        return (
            self.flux_sensible()
            + self.flux_solar()
            + self.flux_sediment()
            + self.flux_atmospheric_longwave()  # shouldn't the value handle the need for the negative?
            + self.flux_upwelling_longwave()  # upwelling flux is typically a loss of energy back to atmosphere
            - self.flux_latent_heat()
        )


    def water_specific_heat_capacity(self, temperature: xr.Dataset) -> xr.Dataset:
        """Approximate the heat capacity of water as a function of water temperature (Celsius)
        Inputs: water_temperature [xr.Dataset] Water temperature in units of Celsius
        Outputs: water_specific_heat_capacity [xr.Dataset] Heat capacity of water in units of J/kg/K
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
        water_temperature: xr.Dataset,
        surface_area: xr.Dataset,
        volume: xr.Dataset,
    ) -> xr.Dataset:
        """
        Compute the change in temperature of the grid in (C)
        """
        return (
            self.flux_net()  # computed though net flux
            * surface_area  # passed in from Riverine
            / (
                volume
                * self.density_water(water_temperature)
                * self.water_specific_heat_capacity(water_temperature)
            )
        )
