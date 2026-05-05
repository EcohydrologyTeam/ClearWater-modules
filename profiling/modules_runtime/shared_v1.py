import numpy as np
import pandas as pd
import xarray as xr

from clearwater_modules.tsm import EnergyBudget

START_TIME = "2026-01-1 00:00:00"
END_TIME = "2026-01-2 00:00:00"
TIME_STEP = "30s"
ITERS = len(pd.date_range(start=START_TIME, end=END_TIME, freq=TIME_STEP))


def make_data_array_variable(grid_x: int, grid_y: int) -> xr.DataArray:
    """Create NxM xarray."""
    return xr.DataArray(
            data=20.0,
            dims=["x", "y", "time"],
            coords={
                "x": np.arange(grid_x),
                "y": np.arange(grid_y),
                "time": pd.date_range(start=START_TIME, end=END_TIME, freq=TIME_STEP),
            },
        )

def init_state_values(grid_x:int, grid_y: int):
    """Initialize state variables."""
    return {
        'water_temp_c': make_data_array_variable(grid_x, grid_y),
        'surface_area': make_data_array_variable(grid_x, grid_y),
        'volume': make_data_array_variable(grid_x, grid_y),
    }


def init_model(grid_x: int, grid_y: int):
    """Initialize Energy Budget model."""
    model = EnergyBudget(
        time_steps=ITERS,
        initial_state_values=init_state_values(grid_x, grid_y),
        meteo_parameters={'wind_c': 1.0},
        track_dynamic_variables=False
    )
    return model


def run_model(model: EnergyBudget):
    """Run for all timesteps."""
    for _ in range(ITERS):
        model.increment_timestep()


def main():
    import sys
    
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <x> <y>")
        sys.exit(1)
    
    x, y = int(sys.argv[1]), int(sys.argv[2])
    print(f"Model size {x} x {y}")
    model = init_model(x, y)
    run_model(model)


if __name__ == "__main__":
    main()