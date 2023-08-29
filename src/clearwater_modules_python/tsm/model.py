import numba
import numpy as np
from enum import Enum
from clearwater_modules_python.tsm import (
    constants,
    equations,
)
import clearwater_modules_python.shared_equations as shared_equations
from typing import (
    TypedDict,
    Protocol,
    runtime_checkable,
    Optional,
)


class EnergyBalanceInputs(TypedDict):
    """A dictionary of energy balance inputs.

    Args:
        water_temp_c: Water temperature entering cell (degrees C)
        surface_area: Surface area of cell face (m^2?)
        volume: Volume of cell (m^3???)
    """
    water_temp_c: float
    surface_area: float
    volume: float


class EnergyBalanceOutputs(TypedDict):
    """A dictionary of energy balance outputs.

    Args:
    """
    # TODO: copy these from the original code.
    ...


@runtime_checkable
class Model(Protocol):
    """A protocol class for our model"""

    def __init__(self, **kwargs) -> None:
        ...

    def increment(self) -> None:
        ...


@runtime_checkable
class Process(Protocol):
    """A protocol class for our model"""

    def __init__(
        self,
        **kwargs,
    ) -> None:
        ...

    def run(
        self,
        variables: TypedDict,
    ) -> TypedDict:
        """Runs the process using configured constants and input variables.

        Args:
            variables: A typed dictionary instance of input variables.
        Returns:
            A typed dictionary instance of output variables.

        """
        ...
        # TODO: figure out how best to do this vectorized on xarray.
        # Lets use apply_ufunc(parallelized=True)


class EnergyBudget(Process):
    """"""

    def __init__(
        self,
        meteo_constants: Optional[dict[str, float]] = None,
        temp_constants: Optional[dict[str, float]] = None,
        use_sed_temp: bool = False,
    ) -> None:
        self.__meteo_constants: constants.Meteorological = constants.DEFAULT_METEOROLOGICAL
        self.__temp_constants: constants.Temperature = constants.DEFAULT_TEMPERATURE

        if meteo_constants is None:
            meteo_constants = {}
        if temp_constants is None:
            temp_constants = {}

        # set default values
        for key, value in self.__meteo_constants.items():
            self.__meteo_constants[key] = meteo_constants.get(
                key,
                value,
            )

        for key, value in self.__temp_constants.items():
            self.__temp_constants[key] = temp_constants.get(
                key,
                value,
            )

        self.use_sed_temp = use_sed_temp

    @property
    def met_constants(self) -> constants.Meteorological:
        return self.__meteo_constants

    @property
    def temp_constants(self) -> constants.Temperature:
        return self.__temp_constants

    def run(
        self,
        variables: EnergyBalanceInputs,
    ) -> float:

        # Temperature
        water_temp_k: float = shared_equations.celsius_to_kelvin(
            variables['water_temp_c'],
        )
        air_temp_k: float = shared_equations.celsius_to_kelvin(
            self.met_constants['air_temp_c'],
        )

        # Pressure
        # Saturated vapor pressure computed from water temperature at previous time step (mb)

        # ----------------------------------------------------------------------------------------------
        # Wind function, stability and flux partitioning
        # ----------------------------------------------------------------------------------------------

        mixing_ratio_air: float = equations.mixing_ratio_air(
            self.met_constants['eair_mb'],
            self.met_constants['pressure_mb'],
        )

        density_air: float = equations.density_air(
            self.met_constants['pressure_mb'],
            air_temp_k,
            mixing_ratio_air,
        )

        # ----------------------------------------------------------------------------------------------
        #  Computations that are not a function of water temperature
        #  Note: solar radiation comes in directly from the interface
        # ----------------------------------------------------------------------------------------------

        emissivity_air: float = equations.emissivity_air(
            air_temp_k,
        )

        # Atmospheric (downwelling) longwave radiation (W/m2)
        q_longwave_down: float = shared_equations.mf_q_longwave_down(
            air_temp_k,
            emissivity_air,
            self.met_constants['cloudiness'],
            self.temp_constants['stephan_boltzmann'],
        )

        # Wind function for latent and sensible heat (unitless)
        wind_function: float = equations.wind_function(
            self.met_constants['wind_a'],
            self.met_constants['wind_b'],
            self.met_constants['wind_c'],
            self.met_constants['wind_speed'],
        )

        # ____________________________________________________________________________________
        #  Energy balance computations that are functions of water temperature
        # ____________________________________________________________________________________
        # Get water temperature from previous time step plus 1/2 AD change??? DLL can not do it.

        # Physical values that are functions of water temperature

        # Latent heat of vaporization (J/kg)
        lv: float = shared_equations.mf_latent_heat_vaporization(water_temp_k)

        # Density (kg/m3)
        density_water: float = shared_equations.mf_density_water(
            variables['water_temp_c'],
        )

        # Specific heat of water (J/kg/K)
        cp_water: float = shared_equations.mf_cp_water(
            variables['water_temp_c'],
        )

        # Saturated vapor pressure computed from water temperature at previous time step (mb)
        esat_mb: float = shared_equations.mf_esat_mb(
            water_temp_k,
            self.temp_constants['a0'],
            self.temp_constants['a1'],
            self.temp_constants['a2'],
            self.temp_constants['a3'],
            self.temp_constants['a4'],
            self.temp_constants['a5'],
            self.temp_constants['a6'],
        )

        # ------------------------------------------------------------------------
        # Upwelling (back or water surface) longwave radiation (W/m2)
        q_longwave_up: float = shared_equations.mf_q_longwave_up(
            water_temp_k,
            self.temp_constants['emissivity_water'],
            self.temp_constants['stephan_boltzmann'],
        )

        # ------------------------------------------------------------------------
        # Surface fluxes
        # ------------------------------------------------------------------------

        # ------------------------------------------------------------------------
        # Compute Richardson number and check stability

        # Richardson's number (unitless)
        ri_number: float = 0.0

        # Richardson's stablility function (unitless)
        ri_function: float = 1.0

        if (self.met_constants['wind_speed'] > 0.0 and self.temp_constants['richardson_option']):
            # Density of air computed at water surface temperature (kg/m3)
            density_air_sat: float = shared_equations.mf_density_air_sat(
                water_temp_k,
                esat_mb,
                self.met_constants['pressure_mb'],
            )
            ri_number, ri_function = shared_equations.RichardsonNumber(
                self.met_constants['wind_speed'],
                density_air_sat,
                density_air,
                self.temp_constants['gravity'],
            )

        # ------------------------------------------------------------------------
        # Latent heat flux (W/m2)
        q_latent: float = equations.q_latent(
            ri_function,
            self.met_constants['pressure_mb'],
            density_water,
            lv,
            wind_function,
            esat_mb,
            self.met_constants['eair_mb'],
        )

        # ------------------------------------------------------------------------
        # Sensible heat flux
        q_sensible: float = equations.q_sensible(
            self.met_constants['wind_kh_kw'],
            ri_function,
            self.temp_constants['Cp_air'],
            density_water,
            wind_function,
            air_temp_k,
            water_temp_k,
        )
        # ------------------------------------------------------------------------
        # Compute sediment heat flux and temperature change
        q_sediment: float = 0.0
        dTdt_sediment_c: float = 0.0
        if self.use_sed_temp:
            q_sediment: float = equations.q_sediment(
                self.temp_constants['pb'],
                self.temp_constants['Cps'],
                self.temp_constants['alphas'],
                self.temp_constants['h2'],
                self.met_constants['sed_temp_c'],
                variables['water_temp_c'],
            )
            dTdt_sediment_c: float = equations.dTdt_sediment_c(
                self.temp_constants['alphas'],
                self.temp_constants['h2'],
                variables['water_temp_c'],
                self.met_constants['sed_temp_c'],
            )

        # ------------------------------------------------------------------------
        # Net heat flux
        q_net: float = equations.q_net(
            q_sensible,
            q_latent,
            q_longwave_down,
            q_longwave_up,
            self.met_constants['q_solar'],
            q_sediment,
        )

        # ------------------------------------------------------------------------
        # Compute water temperature change
        dTdt_water_c: float = equations.dTdt_water_c(
            q_net,
            variables['surface_area'],
            variables['volume'],
            density_water,
            cp_water,
        )
        # TODO: why was the following commented out here -> water_temp_c += dTdt_water_c

        # ------------------------------------------------------------------------------------
        # Difference between air and water temperature (Celsius or Kelvins)
        ta_tw: float = self.met_constants['air_temp_c'] - \
            variables['water_temp_c']

        # Difference between saturated and current vapor pressure (mb)
        esat_eair: float = esat_mb - self.met_constants['eair_mb']

        return dTdt_water_c


if __name__ == '__main__':
    print(constants.Temperature)
    print(constants.Meteorological)
    print(constants.DEFAULT_TEMPERATURE)
    print(constants.DEFAULT_METEOROLOGICAL)
    model_const = EnergyBudget(
        meteo_constants={'air_temp_c': 20},
        temp_constants={'richardson_option': False},
    )
    print(model_const.met_constants)
    assert model_const.met_constants['air_temp_c'] == 20
    assert model_const.temp_constants['richardson_option'] is False
