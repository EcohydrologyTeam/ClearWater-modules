from processes.base import Process
from datetime import datetime, timedelta
from variables import VariableRegistry
import clearwater_riverine as cwr


class Riverine(Process):
    """
    Riverine process.
    """

    variables = []

    def __init__(
        self,
        riverine_instance: cwr.ClearwaterRiverine,
        time_step_frequency: timedelta = timedelta(seconds=30),
    ) -> None:
        self.riverine_instance = riverine_instance
        Process.__init__(self, time_step_frequency)

    def init_process(self, variables: VariableRegistry) -> None:
        """
        Initialize the riverine process.
        """
        # register the water temperature, volume, and surface area to the registry
        variables.register(
            "water_temperature",
            self.riverine_instance.mesh.temperature.copy(deep=False),
        )
        variables.register(
            "volume", self.riverine_instance.mesh.volume.copy(deep=False)
        )
        variables.register(
            "surface_area",
            self.riverine_instance.mesh.faces_surface_area.copy(deep=False),
        )

    def run(self, time_step: datetime, variables: VariableRegistry) -> None:
        """
        Run the riverine process.
        """
        # run next riverine time step
        self.riverine_instance.update()
        # temperature updates are handled by other modules
        # which access the Riverine data directly.
        water_temperature = variables.get("water_temperature")
        prt = 1
