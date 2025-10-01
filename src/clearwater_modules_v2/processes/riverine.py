from processes.base import Process
from datetime import datetime, timedelta
from clearwater_data.variables import VariableRegistry, DataArrayVariable, FloatVariable
import clearwater_riverine as cwr
import clearwater_riverine.utilities as cwr_utils
from pathlib import Path
from processes.nutrients.floating_algae import FloatingAlgae
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from model import Model


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

    def from_file_path(
        configuration_path: str | Path,
        start_datetime: str,
        end_datetime: str,
        time_step_frequency: timedelta = timedelta(seconds=30),
    ) -> "Riverine":
        return Riverine(
            cwr.ClearwaterRiverine(
                config_filepath=configuration_path,
                datetime_range=(start_datetime, end_datetime),
            ),
            time_step_frequency,
        )

    def init_process(self, model: "Model", registry: VariableRegistry) -> None:
        """
        Initialize the riverine process.
        """
        # TODO: Ideally Riverine will register these to the registry as part of it's initialization

        # register the water temperature, volume, and surface area to the registry
        registry.register(
            "water_temperature",
            DataArrayVariable(self.riverine_instance.mesh.temperature.copy(deep=False)),
        )
        registry.register(
            "volume",
            DataArrayVariable(self.riverine_instance.mesh.volume.copy(deep=False)),
        )

        # 'wetted_surface_area' is not calculated by default
        # We may need to specifically call the calculate_wetted_surface_area function
        if "wetted_surface_area" not in self.riverine_instance.mesh:
            cwr_utils.calculate_wetted_surface_area(self.riverine_instance.mesh)

        registry.register(
            "surface_area",
            DataArrayVariable(
                self.riverine_instance.mesh.wetted_surface_area.copy(deep=False)
            ),
        )

        # TODO: update once Riverine can register variables to the registry
        if model.has_process(FloatingAlgae):
            registry.register(
                "algae_floating",
                DataArrayVariable(self.riverine_instance.mesh.Ap.copy(deep=False)),
            )
            registry.register(
                "ammonium",
                DataArrayVariable(self.riverine_instance.mesh.NH4.copy(deep=False)),
            )
            registry.register(
                "nitrate",
                DataArrayVariable(self.riverine_instance.mesh.NO3.copy(deep=False)),
            )
            registry.register(
                "phosphorus_total_inorganic",
                DataArrayVariable(self.riverine_instance.mesh.TIP.copy(deep=False)),
            )
            # TODO: Follow-up with Sarah about having riverine compute depth
            registry.register(
                "depth",
                FloatVariable(1.5),
            )

        # The riverine model use current time_step as the start point and
        # updates the model at the time_step + delta time.
        # The rest of modules uses a definition of time_step as the time
        # to be updated. This boolean allows us to skip the first time_step
        self.__skip_first_time_step = True

    def run(self, time_step: datetime, registry: VariableRegistry) -> None:
        """
        Run the riverine process.
        """
        if self.__skip_first_time_step:
            self.__skip_first_time_step = False
            return

        # run the next time step
        self.riverine_instance.update()

        # previous couplings with riverine model required passing data arrays
        # to the model. Now we are using the registry to access the data.
        # ClearWater-Modules update the memory directly through the registry.
