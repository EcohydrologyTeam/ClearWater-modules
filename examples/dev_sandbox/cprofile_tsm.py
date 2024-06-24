import clearwater_modules
import numpy as np
import xarray as xr
import cProfile
from pathlib import Path


def create_data_array(
    gridsize: int,
    value: int | float  
):
    data = np.full(gridsize, value)
    data_array = xr.DataArray(data)
    return data_array

def run_performance_test(
        iters: int,
        gridsize: int
):
    """Log the performance for different numbers of iterations."""

    if gridsize == 1:
        state_i = {
            'water_temp_c': 40.0,
            'surface_area': 1.0,
            'volume': 1.0,
        }
        meteo_parameters = {'wind_c': 1.0}
    
    else:
        state_i = {
            'water_temp_c': create_data_array(gridsize, 40.0),
            'surface_area': create_data_array(gridsize, 1.0),
            'volume': create_data_array(gridsize, 1.0),
        }
        meteo_parameters = {
            'wind_c': 1.0
        }
    
        # instantiate the TSM module
        tsm = clearwater_modules.tsm.EnergyBudget(
            time_steps=iters,
            initial_state_values=state_i,
            meteo_parameters=meteo_parameters,
            track_dynamic_variables=False,
        )

        for _ in range(iters):
            tsm.increment_timestep()


if __name__ == '__main__':
    iterations_list = [300000, 400000, 500000] # [1000, 10000, 100000, 150000, 200000]
    gridsize_list = [10000, 100000]
    cwd = Path.cwd()
    print(cwd)
    root = Path('./examples/dev_sandbox/profiling/')
    # iterations_list = [1000]
    # gridsize_list = [1]
    # detailed_profile = True
    for iteration in iterations_list:
        for gridsize in gridsize_list:
            cprofile_file = root / f'{iteration}_iters_{gridsize}_gridsize.prof'
            profiler = cProfile.Profile()
            profiler.enable()
            run_performance_test(iteration, gridsize)
            profiler.disable()
            profiler.dump_stats(cprofile_file)
