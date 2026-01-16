from .base import Process, ProcessFactory
from datetime import datetime, timedelta
from clearwater_data.variables import VariableRegistry, DataArrayVariable, FloatVariable
import clearwater_riverine as cwr
import clearwater_riverine.utilities as cwr_utils
from pathlib import Path
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
        time_step: timedelta = timedelta(seconds=30),
    ) -> None:
        self.riverine_instance = riverine_instance
        Process.__init__(self, time_step)

    @staticmethod
    def from_file_path(
        configuration_path: str | Path,
        start_datetime: str | datetime,
        end_datetime: str | datetime,
        time_step: timedelta = timedelta(seconds=30),
    ) -> "Riverine":
        # TODO: This will be removed once Riverine is updated to use datetime objects
        # or the pandas.to_datetime function is used to convert the start and end datetimes
        if isinstance(start_datetime, datetime):
            start_datetime = start_datetime.strftime("%m-%d-%y %H:%M:%S")
        if isinstance(end_datetime, datetime):
            end_datetime = end_datetime.strftime("%m-%d-%y %H:%M:%S")

        return Riverine(
            cwr.ClearwaterRiverine(
                config_filepath=configuration_path,
                datetime_range=(start_datetime, end_datetime),
            ),
            time_step,
        )

    @ProcessFactory.register("riverine")
    @staticmethod
    def from_config(config: dict) -> "Riverine":
        return Riverine.from_file_path(**config)

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
        if model.has_process("FloatingAlgae"):
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
            registry.register(
                "oxygen_dissolved",
                DataArrayVariable(self.riverine_instance.mesh.DOX.copy(deep=False)),
            )
            registry.register(
                "depth",
                DataArrayVariable(self.riverine_instance.mesh.depth.copy(deep=False)),
            )

        # The riverine model use current time_step as the start point and
        # updates the model at the time_step + delta time.
        # The rest of modules uses a definition of time_step as the time
        # to be updated. This boolean allows us to skip the first time_step
        self.__skip_first_time_step = True

    def run(self, time: datetime, registry: VariableRegistry) -> None:
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
