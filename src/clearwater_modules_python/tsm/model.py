"""Water Temperature Simulation Model (TSM) module."""
import numba
import numpy as np
from enum import Enum
from clearwater_modules_python.tsm import (
    parameters,
    processes,
)
from clearwater_modules_python import base
import clearwater_modules_python.shared.processes as shared_processes
from typing import (
    Optional,
)


class EnergyBudget(base.Model):
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
        super().__init__(
            initial_state_values={},
            static_variable_values={},
            track_dynamic_variables=True,
            hotstart_dataset=None,
            time_dim=None,
        )

    @property
    def met_parameters(self) -> parameters.Meteorological:
        return self.__meteo_parameters

    @property
    def temp_parameters(self) -> parameters.Temperature:
        return self.__temp_parameters


