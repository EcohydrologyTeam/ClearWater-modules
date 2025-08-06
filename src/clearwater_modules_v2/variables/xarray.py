from variables.base import Variable
import xarray as xr
from datetime import datetime


class DataArrayVariable(Variable):
    def __init__(self, data_array: xr.DataArray) -> None:
        self.data_array = data_array

    # TODO : lazy loading implementation

    def get(self) -> xr.DataArray:
        return self.data_array

    def get_at_time(self, time: datetime) -> xr.DataArray:
        # if data a has time dimension, return the value at that time
        if "time" in self.data_array.dims:
            return self.data_array.sel(time=time)
        # otherwise return the value
        return self.data_array
