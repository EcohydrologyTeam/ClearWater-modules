from data_io.base import DataSource
from custom_types import ArrayLike
from pathlib import Path
import xarray as xr
from variables import FloatVariable


class FloatDataSource(DataSource):
    def __init__(self, **kwargs) -> None:
        self.value = kwargs.pop("value")

    def read(self, parameter_name: str) -> ArrayLike:
        return FloatVariable(self.value)
