"""Xarray utility functions"""
import xarray as xr
import numpy as np
from clearwater_modules.shared.types import (
    Variable,
)
from typing import Callable
import clearwater_modules.sorter as sorter
import numba


def validate_arrays(array: xr.DataArray, *args: xr.DataArray) -> None:
    """Validate that all DataArrays have the same dimensions."""
    for arg in args:
        if not isinstance(arg, xr.DataArray):
            raise TypeError(
                'All arguments must be of type xarray.DataArray.'
            )
        if tuple(arg.dims) != tuple(array.dims):
            print(tuple(arg.dims), tuple(array.dims))
            raise ValueError(
                'All DataArrays must have the same dimensions.'
            )
        # TODO: do we want to verify coords?
        # if arg.coords != array.coords:
        #    raise ValueError(
        #        'All DataArrays must have the same coordinates.'
        #    )


def _prep_inputs(
    input_dataset: xr.Dataset,
    var: Variable,
) -> tuple[str, Callable, list[np.ndarray]]:
    """Prepare inputs for computation. This is used to speed up computation.

    Returns:
        A tuple with (
            name:str, 
            function:callable, 
            args:tuple[str], 
            arrays:list[np.ndarray]
        )
    """
    args: list[str] = sorter.get_process_args(var.process)
    return (
        var.name,
        var.process,
        [input_dataset[name].values for name in args],
    )


def iter_computations(
    input_dataset: xr.Dataset,
    compute_order: list[Variable],
) -> xr.Dataset:
    """Iterate over the computation order."""
    inputs = map(lambda x: _prep_inputs(input_dataset, x), compute_order)
    dims = input_dataset.dims

    for name, func, arrays in inputs:
        array: np.ndarray = func(*arrays)
        input_dataset[name] = (dims, array)

    return input_dataset
