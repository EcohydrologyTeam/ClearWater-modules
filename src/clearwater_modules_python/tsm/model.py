import numba
import numpy as np
import xarray as xr
from enum import Enum
import constants
import processes
from typing import (
    TypedDict,
    Protocol,
    runtime_checkable,
    Optional,
)


class EnergyBalanceInputs(TypedDict):
    """A dictionary of energy balance inputs.

    Args:
        Twater_C: Water temperature entering cell (degrees C)
        surface_area: Surface area of cell face (m^2?)
        volume: Volume of cell (m^3???)
    """
    Twater_C: float
    surface_area: float
    surface_volume: float


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


class EnergyBudget(Process):
    """"""

    def __init__(
        self,
        meteo_constants: Optional[dict[str, float]] = None,
        temp_constants: Optional[dict[str, float]] = None,
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

    @property
    def met_constants(self) -> constants.Meteorological:
        return self.__meteo_constants

    @property
    def temp_constants(self) -> constants.Temperature:
        return self.__temp_constants

    @property
    def Twater_K(self) -> float:
        return self.__meteo_constants['Twater_C'] + 273.15

    @property
    def Tair_K(self) -> float:
        return self.__meteo_constants['Tair_C'] + 273.15

    def run(
        self,
        variables: EnergyBalanceInputs,
    ) -> None:
        pass


if __name__ == '__main__':
    print(constants.Temperature)
    print(constants.Meteorological)
    print(constants.DEFAULT_TEMPERATURE)
    print(constants.DEFAULT_METEOROLOGICAL)
    model_const = EnergyBudget(
        meteo_constants={'Tair_C': 20},
        temp_constants={'richardson_option': False},
    )
    print(model_const.met_constants)
    assert model_const.met_constants['Tair_C'] == 20
    assert model_const.temp_constants['richardson_option'] is False
