"""Water Temperature Simulation Model (TSM) module."""
import numba
import numpy as np
from enum import Enum
from clearwater_modules_python.tsm import (
    parameters,
    processes,
)
from clearwater_modules_python.base import (
    ParametersDict,
    Variable,
    Model,
)
import clearwater_modules_python.shared.processes as shared_processes
from typing import (
    TypedDict,
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


class EnergyBudget(Model):
    """"""

    def __init__(
        self,
        meteo_parameters: Optional[dict[str, float]] = None,
        temp_parameters: Optional[dict[str, float]] = None,
        use_sed_temp: bool = False,
    ) -> None:
        self.__meteo_parameters: parameters.Meteorological = parameters.DEFAULT_METEOROLOGICAL
        self.__temp_parameters: parameters.Temperature = parameters.DEFAULT_TEMPERATURE

        if meteo_parameters is None:
            meteo_parameters = {}
        if temp_parameters is None:
            temp_parameters = {}

        # set default values
        for key, value in self.__meteo_parameters.items():
            self.__meteo_parameters[key] = meteo_parameters.get(
                key,
                value,
            )

        for key, value in self.__temp_parameters.items():
            self.__temp_parameters[key] = temp_parameters.get(
                key,
                value,
            )

        self.use_sed_temp = use_sed_temp

    @property
    def met_parameters(self) -> parameters.Meteorological:
        return self.__meteo_parameters

    @property
    def temp_parameters(self) -> parameters.Temperature:
        return self.__temp_parameters


if __name__ == '__main__':
    print(parameters.Temperature)
    print(parameters.Meteorological)
    print(parameters.DEFAULT_TEMPERATURE)
    print(parameters.DEFAULT_METEOROLOGICAL)
    model_const = EnergyBudget(
        meteo_parameters={'air_temp_c': 20},
        temp_parameters={'richardson_option': False},
    )
    print(model_const.met_parameters)
    assert model_const.met_parameters['air_temp_c'] == 20
    assert model_const.temp_parameters['richardson_option'] is False
