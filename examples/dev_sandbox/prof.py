"""A script to allow for debugging of the TSM module."""
import clearwater_modules
import time 
import sys
import xarray as xr
import numpy as np


def main(iters: int, baseline: bool):
    ti = time.time()
    # define starting state values
    if baseline:
        state_i = {
            'water_temp_c': 40.0,
            'surface_area': 1.0,
            'volume': 1.0,
        }
    else:
        state_i = {
            'water_temp_c': xr.DataArray(
                np.full(10, 40),
                dims='cell',
                coords={'cell': np.arange(10)}),
            'surface_area': xr.DataArray(
                np.full(10, 1.0),
                dims='cell',
                coords={'cell': np.arange(10)}),
            'volume': xr.DataArray(
                np.full(10, 1.0),
                dims='cell',
                coords={'cell': np.arange(10)}),
        }

    # instantiate the TSM module
    tsm = clearwater_modules.tsm.EnergyBudget(
        time_steps=iters,
        initial_state_values=state_i,
        meteo_parameters={'wind_c': 1.0},
        updateable_static_variables=['wind_c']
    )
    print(tsm.static_variable_values)
    t2 = time.time()
    for _ in range(iters):
        tsm.increment_timestep()
    print(f'Increment timestep speed (average of {iters}): {(time.time() - t2) / 100}')
    print(f'Run time: {time.time() - ti}')

if __name__ == '__main__':
    if len(sys.argv) > 1:
        try:
            iters = int(sys.argv[1])
            print(f'Running {iters} iterations.')
        except ValueError:
            raise ValueError('Argument must be an integer # of iterations.')
    else:
        print('No argument given, defaulting to 100 iteration.')
        iters = 100
            
    main(iters=iters, baseline=True)
