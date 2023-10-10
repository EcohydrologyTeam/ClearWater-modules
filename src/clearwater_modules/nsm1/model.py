"""Water Nutrient Simulation Model 1 (NSM1) module."""
import xarray as xr
import numpy as np
from enum import Enum
from clearwater_modules.nsm1 import (
    constants,
)
from clearwater_modules import base
import clearwater_modules.shared.processes as shared_processes
from typing import (
    Optional,
)


class NutrientBudget(base.Model):
    """"""
    _variables: list[base.Variable] = []

    def __init__(
        self,
        initial_state_values: Optional[base.InitialVariablesDict] = None,
        algae_parameters: Optional[dict[str, float]] = None,
        track_dynamic_variables: bool = True,
        hotstart_dataset: Optional[xr.Dataset] = None,
        time_dim: Optional[str] = None,
    ) -> None:
        self.__algae_parameters: constants.Algae = constants.DEFAULT_Algae

        if algae_parameters is None:
            algae_parameters_parameters = {}

        # set default values
        for key, value in self.__algae_parameters.items():
            self.__algae_parameters[key] = algae_parameters.get(
                key,
                value,
            )

        static_variable_values = {
            **self.__algae_parameters}

        super().__init__(
            initial_state_values=initial_state_values,
            static_variable_values=static_variable_values,
            track_dynamic_variables=track_dynamic_variables,
            hotstart_dataset=hotstart_dataset,
            time_dim=time_dim,
        )

    @property
    def algae_parameters(self) -> constants.algae:
        return self.__algae_parameters
