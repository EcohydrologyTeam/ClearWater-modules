from processes.base import Process, ProcessFactory
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
        variable_registry: VariableRegistry,
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
                start_datetime=start_datetime,
                end_datetime=end_datetime,
                variable_registry=variable_registry,
            ),
            time_step,
        )

    @ProcessFactory.register("riverine")
    @staticmethod
    def from_config(config: dict, variable_registry: VariableRegistry) -> "Riverine":
        return Riverine.from_file_path(**config, variable_registry=variable_registry)

    def init_process(self, model: "Model", registry: VariableRegistry) -> None:
        """
        Initialize the riverine process.
        """

        def check_variable_in_registry(variable_name: str) -> None:
            if variable_name not in registry:
                raise KeyError(
                    f"Variable '{variable_name}' not found in registry. Did models initialize correctly."
                )

        # register the water temperature, volume, and surface area to the registry
        # verify required varables have initialized correctly
        check_variable_in_registry("water_temperature")
        check_variable_in_registry("volume")
        try:
            check_variable_in_registry("wetted_surface_area")
        except KeyError:
            cwr_utils.calculate_wetted_surface_area(self.riverine_instance)

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

            # TODO: replace this with depth calculation
            registry.register(
                "depth",
                DataArrayVariable(
                    self.riverine_instance.mesh.wetted_surface_area.copy(deep=False)
                ),
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
