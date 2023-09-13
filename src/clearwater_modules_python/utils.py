"""Xarray utility functions"""
import xarray as xr

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
        if arg.coords != array.coords:
            raise ValueError(
                'All DataArrays must have the same coordinates.'
            )
