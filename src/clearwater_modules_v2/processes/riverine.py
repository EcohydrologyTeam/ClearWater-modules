from processes.base import Process
from datetime import datetime
from variables import VariableRegistry
import xarray as xr
import numpy as np
import clearwater_riverine as cwr


class Riverine(Process):
    """
    Riverine process.
    """
    variables = []

    def __init__(self, riverine_instance:cwr.ClearwaterRiverine) -> None:
        self.riverine_instance = riverine_instance

    def init_process(self, variables:VariableRegistry) -> None:
        """
        Initialize the riverine process.
        """
        #register the water temperature, volume, and surface area to the registry
        variables.register("water_temperature", self.riverine_instance.mesh.temperature.copy(deep=False))
        variables.register("volume", self.riverine_instance.mesh.volume.copy(deep=False))
        variables.register("surface_area", self.riverine_instance.mesh.faces_surface_area.copy(deep=False))

    def run(self, time_step:datetime, variables:VariableRegistry) -> None:
        """
        Run the riverine process.
        """
        #run next riverine time step
        self.riverine_instance.update()
        #temperature updates are handled by other modules
        #which access the Riverine data directly.