from variables.base import Variable
import xarray as xr
from datetime import datetime


class DataArrayVariable(Variable):
    def __init__(self, data_array: xr.DataArray):
        self.data_array = data_array

    def get(self) -> xr.DataArray:
        return self.data_array

    def get_at_time(self, time: datetime) -> xr.DataArray:
        return self.data_array.sel(time=time)
