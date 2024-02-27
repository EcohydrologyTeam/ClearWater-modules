import clearwater_modules
import time
import sys
import logging
import numpy as np
import xarray as xr

from typing import (
    List
)

def create_data_array(
    gridsize: int,
    value: int | float  
):
    data = np.full(gridsize, value)
    data_array = xr.DataArray(data)
    return data_array

def run_performance_test(
        iterations_list: List[int],
        gridsize_list: List[int],
        log_file: str,
        detailed_profile: bool
):
    """Log the performance for different numbers of iterations."""
    # set up logger
    log_format = '%(asctime)s,%(levelname)s,%(message)s'
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format=log_format,
    )
    if detailed_profile:
        logging.info(f"timestep,MB,time_per_iter")
    else:
        logging.info(f"iters,gridsize,avg_increment_time,total_increment_time")

    for iters in iterations_list:
        for gridsize in gridsize_list:
            # define starting state values
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
            
            ti = time.time()

            # instantiate the TSM module
            tsm = clearwater_modules.tsm.EnergyBudget(
                time_steps=iters,
                initial_state_values=state_i,
                meteo_parameters=meteo_parameters,
            )

            t2 = time.time()
            curr_time = t2
            for _ in range(iters):
                tsm.increment_timestep()

                # detailed profiling
                if detailed_profile:
                    if (_ % 100 == 0) and (_ != 0):
                        current_time = (time.time() - curr_time) / 100
                        logging.info(f"{_},{tsm.dataset.nbytes * 0.000001}, {current_time}")
                        curr_time = time.time()
            
            # wrap up
            avg_increment_time = (time.time() - t2) / iters
            total_run_time = time.time() - ti
            if not detailed_profile:
                logging.info(f"{iters},{gridsize},{avg_increment_time},{total_run_time}")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python performance_profiling_tsm.py <log_file>")
        sys.exit(1)

    log_file = sys.argv[1]
    iterations_list = [1, 10, 100, 1000, 10000, 100000]
    gridsize_list = [1, 1000, 10000]
    # iterations_list = [10000]
    # gridsize_list = [10000]
    detailed_profile = False
    run_performance_test(iterations_list, gridsize_list, log_file, detailed_profile)
