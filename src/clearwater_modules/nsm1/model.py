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
        updateable_static_variables: Optional[list[str]] = None,
        algae_parameters: Optional[dict[str, float]] = None,
        alkalinity_parameters: Optional[dict[str, float]] = None,
        balgae_parameters: Optional[dict[str, float]] = None,
        carbon_parameters: Optional[dict[str, float]] = None,
        CBOD_parameters: Optional[dict[str, float]] = None,
        DOX_parameters: Optional[dict[str, float]] = None,
        nitrogen_parameters: Optional[dict[str, float]] = None,
        POM_parameters: Optional[dict[str, float]] = None,
        N2_parameters: Optional[dict[str, float]] = None,
        phosphorus_parameters: Optional[dict[str, float]] = None,
        pathogen_parameters: Optional[dict[str, float]] = None, 
        global_parameters: Optional[dict[str, float]] = None, 
        global_vars: Optional[dict[str, float]] = None,             
        track_dynamic_variables: bool = True,
        hotstart_dataset: Optional[xr.Dataset] = None,
        time_dim: Optional[str] = None,
    ) -> None:
        self.__algae_parameters: constants.AlgaeStaticVariables = constants.DEFAULT_ALGAE
        self.__alkalinity_parameters: constants.AlkalinityStaticVariables = constants.DEFAULT_ALKALINITY
        self.__balgae_parameters: constants.BalgaeStaticVariables = constants.DEFAULT_BALGAE
        self.__carbon_parameters: constants.CarbonStaticVariables = constants.DEFAULT_CARBON
        self.__CBOD_parameters: constants.CBODStaticVariables = constants.DEFAULT_CBOD
        self.__DOX_parameters: constants.DOXStaticVariables = constants.DEFAULT_DOX
        self.__nitrogen_parameters: constants.NitrogenStaticVariables = constants.DEFAULT_NITROGEN
        self.__POM_parameters: constants.POMStaticVariables = constants.DEFAULT_POM
        self.__N2_parameters: constants.N2StaticVariables = constants.DEFAULT_N2
        self.__phosphorus_parameters: constants.PhosphorusStaticVariables = constants.DEFAULT_PHOSPHORUS
        self.__pathogen_parameters: constants.PathogenStaticVariables = constants.DEFAULT_PATHOGEN
        self.__global_parameters: constants.PathogenStaticVariables = constants.DEFAULT_GLOBALPARAMETERS
        self.__global_vars: constants.PathogenStaticVariables = constants.DEFAULT_GLOBALVARS
        

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
        if nitrogen_parameters is None:
            nitrogen_parameters = {}
        if POM_parameters is None:
            POM_parameters = {}
        if N2_parameters is None:
            N2_parameters = {}
        if phosphorus_parameters is None:
            phosphorus_parameters = {}
        if pathogen_parameters is None:
            pathogen_parameters = {}
        if global_parameters is None:
            global_parameters = {}
        if global_vars is None:
            global_vars = {}
            
        # set default values
        for key, value in self.__algae_parameters.items():
            self.__algae_parameters[key] = algae_parameters.get(
                key,
                value,
            )
        for key, value in self.__alkalinity_parameters.items():
            self.__alkalinity_parameters[key] = alkalinity_parameters.get(
                key,
                value,
            )
        for key, value in self.__balgae_parameters.items():
            self.__balgae_parameters[key] = balgae_parameters.get(
                key,
                value,
            )
        for key, value in self.__carbon_parameters.items():
            self.__carbon_parameters[key] = carbon_parameters.get(
                key,
                value,
            )
        for key, value in self.__CBOD_parameters.items():
            self.__CBOD_parameters[key] = CBOD_parameters.get(
                key,
                value,
            )
        for key, value in self.__DOX_parameters.items():
            self.__DOX_parameters[key] = DOX_parameters.get(
                key,
                value,
            )
        for key, value in self.__nitrogen_parameters.items():
            self.__nitrogen_parameters[key] = nitrogen_parameters.get(
                key,
                value,
            )
        for key, value in self.__POM_parameters.items():
            self.__POM_parameters[key] = POM_parameters.get(
                key,
                value,
            )
        for key, value in self.__N2_parameters.items():
            self.__N2_parameters[key] = N2_parameters.get(
                key,
                value,
            )
        for key, value in self.__phosphorus_parameters.items():
            self.__phosphorus_parameters[key] = phosphorus_parameters.get(
                key,
                value,
            )
        for key, value in self.__pathogen_parameters.items():
            self.__pathogen_parameters[key] = pathogen_parameters.get(
                key,
                value,
            )
            
        for key, value in self.__global_parameters.items():
            self.__global_parameters[key] = global_parameters.get(
                key,
                value,
            )
            
        for key, value in self.__global_vars.items():
            self.__global_vars[key] = global_vars.get(
                key,
                value,
            )

        static_variable_values = {
            **self.__algae_parameters,
            **self.__alkalinity_parameters,
            **self.__balgae_parameters,
            **self.__carbon_parameters,
            **self.__CBOD_parameters,
            **self.__DOX_parameters,
            **self.__nitrogen_parameters,
            **self.__POM_parameters,
            **self.__N2_parameters,
            **self.__phosphorus_parameters,
            **self.__pathogen_parameters,
            **self.__global_parameters,
            **self.__global_vars}
        
        # TODO: make sure this feature works -> test it, but post demo
        #static_variable_values['use_sed_temp'] = use_sed_temp

        super().__init__(
            initial_state_values=initial_state_values,
            static_variable_values=static_variable_values,
            updateable_static_variables=updateable_static_variables,            
            track_dynamic_variables=track_dynamic_variables,
            hotstart_dataset=hotstart_dataset,
            time_dim=time_dim,
        )

    @property
    def algae_parameters(self) -> constants.AlgaeStaticVariables:
        return self.__algae_parameters

    @property
    def alkalinity_parameters(self) -> constants.AlkalinityStaticVariables:
        return self.__alkalinity_parameters
    
    @property
    def balgae_parameters(self) -> constants.BalgaeStaticVariables:
        return self.__balgae_parameters
    
    @property
    def carbon_parameters(self) -> constants.CarbonStaticVariables:
        return self.__carbon_parameters
    
    @property
    def CBOD_parameters(self) -> constants.CBODStaticVariables:
        return self.__CBOD_parameters
    
    @property
    def DOX_parameters(self) -> constants.DOXStaticVariables:
        return self.__DOX_parameters
    
    @property
    def nitrogen_parameters(self) -> constants.NitrogenStaticVariables:
        return self.__nitrogen_parameters
    
    @property
    def POM_parameters(self) -> constants.POMStaticVariables:
        return self.__POM_parameters
    
    @property
    def N2_parameters(self) -> constants.N2StaticVariables:
        return self.__N2_parameters
    
    @property
    def phosphorus_parameters(self) -> constants.PhosphorusStaticVariables:
        return self.__phosphorus_parameters
    
    @property
    def pathogen_parameters(self) -> constants.PathogenStaticVariables:
        return self.__pathogen_parameters
    
    @property
    def global_parameters(self) -> constants.GlobalParameters:
        return self.__global_parameters
    
    @property
    def global_vars(self) -> constants.GlobalVars:
        return self.__global_vars
