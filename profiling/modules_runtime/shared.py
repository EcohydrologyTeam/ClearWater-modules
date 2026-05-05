import datetime

from clearwater_modules_v2.model import Model
from clearwater_modules_v2.processes import Temperature
from clearwater_data.variables import VariableRegistry, DataArrayVariable, FloatVariable

import xarray as xr
import pandas as pd
import numpy as np

# This will set up the model for execution.
# We don't want to profile this, so there is apples to apples comparison.

START_TIME = "2026-01-1 00:00:00"
END_TIME = "2026-01-2 00:00:00"
TIME_STEP = "30s"


def make_data_array_variable(
    start_time: datetime.datetime, end_time: datetime.datetime, x: int, y: int
) -> DataArrayVariable:
    return DataArrayVariable(
        data_array=xr.DataArray(
            data=20.0,
            dims=["x", "y", "time"],
            coords={
                "x": np.arange(x),
                "y": np.arange(y),
                "time": pd.date_range(start=start_time, end=end_time, freq="30s"),
            },
        ),
        time_dimension="time",
        space_dimension=["x", "y"],
    )


def make_timeseries(
    start_time: datetime.datetime, end_time: datetime.datetime
) -> DataArrayVariable:
    return DataArrayVariable(
        data_array=xr.DataArray(
            data=20,
            dims=["time"],
            coords={"time": pd.date_range(start_time, end_time, freq="30s")},
        )
    )


def init_variable_registry(
    grid_size: tuple(int, int),
    start_time: datetime.datetime,
    end_time: datetime.datetime,
    time_step,
) -> VariableRegistry:
    registry = VariableRegistry()

    # water_temperature
    # gridded timeseries
    registry.register(
        "water_temperature",
        make_data_array_variable(start_time, end_time, grid_size[0], grid_size[1]),
    )

    # wetted_surface_area
    # gridded timeseries
    registry.register(
        "wetted_surface_area",
        make_data_array_variable(start_time, end_time, grid_size[0], grid_size[1]),
    )

    # volume
    # gridded timeseries
    registry.register(
        "volume",
        make_data_array_variable(start_time, end_time, grid_size[0], grid_size[1]),
    )

    # cloudiness
    registry.register("cloudiness", FloatVariable(0.1))

    # air_temperature
    # timeseries, but not gridded
    registry.register("air_temperature", make_timeseries(start_time, end_time))

    # solar_radiation
    # timeseries, but not gridded
    registry.register("solar_radiation", make_timeseries(start_time, end_time))

    # wind_speed
    registry.register("wind_speed", FloatVariable(3.0))

    # atmospheric_pressure
    registry.register("atmospheric_pressure", FloatVariable(1013.0))

    # atmospheric_vapor_pressure
    registry.register("atmospheric_vapor_pressure", FloatVariable(1013.0))

    # sediment_temperature
    registry.register("sediment_temperature", FloatVariable(20.0))

    # sediment_thickness
    registry.register("sediment_thickness", FloatVariable(0.1))

    return registry


def init_model(grid_x: int, grid_y: int) -> Model:
    temperature = Temperature(1.3, 1.5, 3.0)

    model = Model(
        start_time=pd.to_datetime(START_TIME),
        end_time=pd.to_datetime(END_TIME),
        time_step=pd.to_timedelta(TIME_STEP),
        processes=(temperature,),
        output_variables=[],
        variable_data_sources={},
        simulation_directory=None,
        variable_registry=init_variable_registry(
            (grid_x, grid_y), START_TIME, END_TIME, TIME_STEP
        ),
    )
    model.init_model()
    return model


def main():
    import sys

    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <x> <y>")
        sys.exit(1)

    x, y = int(sys.argv[1]), int(sys.argv[2])
    print(f"Model size {x} x {y}")
    model = init_model(x, y)
    model.run()


if __name__ == "__main__":
    main()
