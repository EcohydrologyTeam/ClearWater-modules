"""Water Temperature Simulation Model (TSM) module."""
import numba
import numpy as np
from enum import Enum
from clearwater_modules_python.tsm import (
    constants,
    equations,
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


class EnergyBalanceOutputs(TypedDict):
    """A dictionary of energy balance outputs.

    Args:
    """
    # TODO: copy these from the original code.
    ...


class EnergyBudget:
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
