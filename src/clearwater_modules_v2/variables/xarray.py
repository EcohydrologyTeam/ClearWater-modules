from variables.base import Variable
import xarray as xr
from datetime import datetime


class DataArrayVariable(Variable):
    def __init__(
        self, data_array: xr.DataArray, time_dimension: str | None = "time"
    ) -> None:
        self.data_array = data_array
        self.__time_dimension = time_dimension

    @property
    def time_dimension(self) -> str | None:
        return self.__time_dimension

    def get(self) -> xr.DataArray:
        return self.data_array

    def get_at_time(self, time: datetime) -> xr.DataArray:
        # if data a has time dimension, return the value at that time
        if "time" in self.data_array.dims:
            return self.data_array.sel({self.time_dimension: time})
        # otherwise return the value
        return self.data_array
