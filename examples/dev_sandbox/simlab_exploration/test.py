import numpy as np
import matplotlib.pyplot as plt
import xsimlab as xs
import xarray as xr

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
    temperature = xs.variable(
        dims=("time", "nface"),
        intent="out",
        description="Temperature of the riverine mesh at each time step and face",
    )
    volume = xs.variable(
        dims=("time", "nface"),
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

    def finalize_step(self):
        self.temperature = self._transport.mesh.temperature.values
        self.volume = self._transport.mesh.volume.values
        self.wetted_surface_area = self._transport.mesh.faces_surface_area.values
        # self.temperature = self._transport.mesh[crosswalk["temperature"]]


@xs.process
class ArbitraryProcess:
    """
    This is an arbitrary process that does nothing but is used to test the model.
    """

    calls = xs.variable(
        dims=(),
        intent="out",
        description="Just a number indicating home many times the process has been called",
    )

    def initialize(self) -> None:
        self._number_of_calls = 0

    @xs.runtime(args="step_delta")
    def run_step(self, dt):
        self._number_of_calls += 1

    def finalize_step(self):
        self.calls = self._number_of_calls


model = xs.Model({"transport": Riverine, "arbitrary": ArbitraryProcess})

config_path = r"C:\Users\ptomasula\Repositories\ClearWater-modules\examples\dev_sandbox\data_temp\sumwere_creek_coarse_p48\demo_config.yml"
path = xr.DataArray(
    [config_path],
    dims=("path"),
)

start_index = int(
    8 * 60 * (60 / 30)
)  # start at 8:00 am on the first day of the simulation (30 second model)
end_index = start_index + int(8 * 60 * (60 / 30))

input_dataset = xs.create_setup(
    model=model,
    clocks={"step": np.arange(start_index, end_index, 1)},
    input_vars={
        "transport__riverine_path": path,  # path,
        "transport__datetime_range": ("indicies", [start_index, end_index]),
        #'transport__mesh': (mesh.values),
    },
    output_vars={
        "arbitrary__calls": None,
        "transport__temperature": None,
        "transport__volume": None,
        "transport__wetted_surface_area": None,
    },
)


output_dataset = input_dataset.xsimlab.run(model=model)

print(output_dataset)
