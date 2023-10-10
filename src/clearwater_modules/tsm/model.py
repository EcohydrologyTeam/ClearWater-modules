"""Water Nutrient Simulation Model (NSM1) module."""
import xarray as xr
import numpy as np
from enum import Enum
from clearwater_modules.tsm import (
    constants,
)
from clearwater_modules import base
import clearwater_modules.shared.processes as shared_processes
from typing import (
    Optional,
)


class EnergyBudget(base.Model):
    """"""
    _variables: list[base.Variable] = []

    def __init__(
        self,
        initial_state_values: Optional[base.InitialVariablesDict] = None,
        meteo_parameters: Optional[dict[str, float]] = None,
        temp_parameters: Optional[dict[str, float]] = None,
        track_dynamic_variables: bool = True,
        hotstart_dataset: Optional[xr.Dataset] = None,
        time_dim: Optional[str] = None,
    ) -> None:
        self.__meteo_parameters: constants.Meteorological = constants.DEFAULT_METEOROLOGICAL
        self.__temp_parameters: constants.Temperature = constants.DEFAULT_TEMPERATURE

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

        static_variable_values = {
            **self.__meteo_parameters, **self.__temp_parameters}

        # TODO: make sure this feature works -> test it, but post demo
        # static_variable_values['use_sed_temp'] = use_sed_temp

        super().__init__(
            initial_state_values=initial_state_values,
            static_variable_values=static_variable_values,
            track_dynamic_variables=track_dynamic_variables,
            hotstart_dataset=hotstart_dataset,
            time_dim=time_dim,
        )

    @property
    def met_parameters(self) -> constants.Meteorological:
        return self.__meteo_parameters

    @property
    def temp_parameters(self) -> constants.Temperature:
        return self.__temp_parameters
