"""A script to allow for debugging of the TSM module."""
import clearwater_modules
import numpy as np


def main():
    # define starting state values
    state_i = {
        'water_temp_c': 20.0,
        'surface_area': 1.0,
        'volume': 1.0,
    }

    # instantiate the TSM module
    tsm = clearwater_modules.tsm.EnergyBudget(
        initial_state_values=state_i,
    )

    input_dataset = tsm.dataset.isel(time_step=0).copy()

    inputs = map(
        lambda x: clearwater_modules.utils._prep_inputs(input_dataset, x),
        tsm.computation_order,
    )
    dims = input_dataset.dims

    for name, func, arrays in inputs:
        print(name)
        array: np.ndarray = func(*arrays)
        input_dataset[name] = (dims, array)
        print(array)
        print()
    return input_dataset


if __name__ == '__main__':
    main()
