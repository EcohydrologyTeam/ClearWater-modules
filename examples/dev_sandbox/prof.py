"""A script to allow for debugging of the TSM module."""
import clearwater_modules
import time 
import sys

def main(iters: int):
    ti = time.time()
    # define starting state values
    state_i = {
        'water_temp_c': 40.0,
        'surface_area': 1.0,
        'volume': 1.0,
    }

    # instantiate the TSM module
    tsm = clearwater_modules.tsm.EnergyBudget(
        initial_state_values=state_i,
        meteo_parameters={'wind_c': 1.0},
    )
    print(tsm.static_variable_values)
    t2 = time.time()
    for _ in range(iters):
        tsm.increment_timestep()
    print(f'Increment timestep speed (average of {iters}): {(time.time() - t2) / iters}')
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
            
    main(iters=iters)
