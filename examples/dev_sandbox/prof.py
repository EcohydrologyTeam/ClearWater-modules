"""A script to allow for debugging of the TSM module."""
import clearwater_modules
import time 

def main():
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
    for _ in range(100):
        tsm.increment_timestep()
    print(f'Increment timestep speed (average of 100): {(time.time() - t2) / 100}')
    print(f'Run time: {time.time() - ti}')
if __name__ == '__main__':
    main()
