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
        alkalinity_parameters: Optional[dict[str, float]] = None,
        balgae_parameters: Optional[dict[str, float]] = None,
        carbon_parameters: Optional[dict[str, float]] = None,
        CBOD_parameters: Optional[dict[str, float]] = None,
        DOX_parameters: Optional[dict[str, float]] = None,
        nitrogen_parameters: Optional[dict[str, float]] = None,
        POM_parameters: Optional[dict[str, float]] = None,           
        track_dynamic_variables: bool = True,
        hotstart_dataset: Optional[xr.Dataset] = None,
        time_dim: Optional[str] = None,
    ) -> None:
        self.__algae_parameters: constants.Algae = constants.DEFAULT_ALGAE
        self.__alkalinity_parameters: constants.Alkalinity = constants.DEFAULT_ALKALINITY
        self.__balgae_parameters: constants.Balgae = constants.DEFAULT_BALGAE
        self.__carbon_parameters: constants.Carbon = constants.DEFAULT_CARBON
        self.__CBOD_parameters: constants.CBOD = constants.DEFAULT_CBOD
        self.__DOX_parameters: constants.DOX = constants.DEFAULT_DOX
        self.__nitrogen_parameters: constants.DOX = constants.DEFAULT_DOX
        self.__POM_parameters: constants.DOX = constants.DEFAULT_DOX
        

        if algae_parameters is None:
            algae_parameters = {}
        if alkalinity_parameters is None:
            alkalinity_parameters = {}
        if balgae_parameters is None:
            balgae_parameters = {}
        if carbon_parameters is None:
            carbon_parameters = {}
        if CBOD_parameters is None:
            CBOD_parameters = {}
        if DOX_parameters is None:
            DOX_parameters = {}




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
