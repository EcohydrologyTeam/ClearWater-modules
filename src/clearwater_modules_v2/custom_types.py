"""Defines custom types for ClearWater modules"""

from typing import Union
import xarray as xr
import numpy as np

ArrayLike = Union[xr.DataArray, np.ndarray, float]