import xsimlab as xs

import clearwater_riverine as cwr

import pickle
import os


@xs.process
class Riverine:
    riverine_path = xs.variable(
        dims=("path"),
        intent="in",
        description="Path to the configuration file for the riverine model",
    )
    # This should eventually be a model parameter and not specific to the process
    datetime_range = xs.variable(
        dims=("indicies"),
        intent="in",
        description="Datetime range for the simulation",
    )
    water_temperature = xs.variable(
        dims=("nface"),
        intent="out",
        global_name="water_temperature",
        description="Water temperature of the grid",
    )
    volume = xs.variable(
        dims=("nface"),
        intent="out",
        description="Volume of the cell",
    )
    wetted_surface_area = xs.variable(
        dims=("nface"),
        intent="out",
        description="Wetted surface area of the cell",
    )
    # volume (volume)
    # wetted_surface_area (faces_surface_area)
    # Sarah noted we could use surface area for a stub

    def initialize(self):
        config_filepath = self.riverine_path[0]
        print(config_filepath)
        start_index, end_index = self.datetime_range
        print(start_index, end_index)

        # This is old for testing purposes
        # The use of the pickle speeds up  interactions by avoiding model initialization
        # each time it runs
        fp = "riverine.plk"
        if os.path.exists(fp):
            with open(fp, "rb") as f:
                self._transport = cwr.ClearwaterRiverine(
                    config_filepath=config_filepath, datetime_range=(960, 1920)
                )
        else:
            self._transport = cwr.ClearwaterRiverine(
                config_filepath=config_filepath, datetime_range=(960, 1920)
            )
            with open(fp, "wb") as f:
                pickle.dump(self._transport, f)

    @xs.runtime(args="step_delta")
    def run_step(self, dt):
        self._transport.update()
        self.water_temperature = self._transport.mesh.temperature.values
        self.volume = self._transport.mesh.volume.values
        self.wetted_surface_area = self._transport.mesh.faces_surface_area.values
