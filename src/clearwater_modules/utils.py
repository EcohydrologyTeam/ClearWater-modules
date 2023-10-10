"""Xarray utility functions"""
import xarray as xr
from clearwater_modules.shared.types import (
    Variable,
)
import clearwater_modules.sorter as sorter
import numba


def validate_arrays(array: xr.DataArray, *args: xr.DataArray) -> None:
    """Validate that all DataArrays have the same dimensions."""
    for arg in args:
        if not isinstance(arg, xr.DataArray):
            raise TypeError(
                'All arguments must be of type xarray.DataArray.'
            )
        if arg.dims != array.dims:
            raise ValueError(
                'All DataArrays must have the same dimensions.'
            )
        # TODO: do we want to verify coords?
        # if arg.coords != array.coords:
        #    raise ValueError(
        #        'All DataArrays must have the same coordinates.'
        #    )


@numba.jit(forceobj=True)
def iter_computations(
    input_array: xr.DataArray,
    compute_order: list[Variable],
) -> xr.DataArray:
    """Iterate over the computation order."""
    for var in compute_order:
        input_vars: list[str] = sorter.get_process_args(var.process)
        input_array[var.name] = xr.apply_ufunc(
            var.process,
            *[input_array[name] for name in input_vars],
        )
    return input_array
