from data.base import DataProvider
from custom_types import ArrayLike
from pathlib import Path
import xarray as xr


class FloatDataProvider(DataProvider):
    def __init__(self, **kwargs) -> None:
        self.value = kwargs.pop("value")

    def read(self) -> ArrayLike:
        return xr.DataArray(self.value, dims=["time"])

    def write_to_store(
        self, store_path: Path, variable_name: str, field_name: str | None = None
    ) -> None:
        store_path = store_path / variable_name
        self.__data.to_zarr(store_path)
