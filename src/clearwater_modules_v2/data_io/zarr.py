from custom_types import ArrayLike
from datetime import datetime, timedelta
from pathlib import Path
import xarray as xr
import dask.array as da
import pandas as pd
from variables.xarray import DataArrayVariable


class ZarrDataSource:
    def __init__(self, **kwargs) -> None:
        self.store_path: Path = kwargs.pop("store_path")
        self.__dataset = xr.open_zarr(self.store_path)

    def read(self, parameter_name: str) -> DataArrayVariable:
        return DataArrayVariable(self.__dataset[parameter_name])


class ZarrDataStore:
    def __init__(self, **kwargs) -> None:
        self.store_path: Path = kwargs.pop("store_path")
        self.start_date: datetime = kwargs.pop("start_date")
        self.end_date: datetime = kwargs.pop("end_date")
        self.time_step: timedelta = kwargs.pop("time_step")
        self.variables: list[str] = kwargs.pop("variables")

        self.__init_zarr_store()

    def __init_zarr_store(self) -> None:
        time = pd.date_range(self.start_date, self.end_date, freq=self.time_step)
        template_dataset = xr.Dataset(
            {
                v: (("time", "scalar"), da.empty((time.shape[0], 1), dtype="float"))
                for v in self.variables
            },
            coords={"time": time, "scalar": [1]},
        )

        # write the template out to generate zarr
        template_dataset.to_zarr(self.store_path, mode="w", compute=False)

    def write(self, data: ArrayLike, parameter_name: str) -> None:
        data.to_zarr(self.store_path, mode="a")


class ChunkedZarrDataStore:
    def __init__(self, **kwargs) -> None:
        self.store_path = kwargs.pop("store_path")

    def write_chunk(
        self,
        data: ArrayLike,
        parameter_name: str,
        start_time: datetime,
        end_time: datetime,
    ) -> None:
        data.to_zarr(self.store_path, mode="a")


# This should be the init method for the ChunkedZarrDataStore
"""
def init_data_store(
    root_directory: Path,
    start_time: datetime,
    end_time: datetime,
    time_step: timedelta,
    variables: list[str],
) -> Path:
    # define template
    time = pd.date_range(start_time, end_time, freq=time_step)
    template_dataset = xr.Dataset(
        {
            v: (("time", "scalar"), da.empty((time.shape[0], 1), dtype="float"))
            for v in variables
        },
        coords={"time": time, "scalar": [1]},
    )

    # write the template out to generate zarr
    template_dataset.to_zarr(root_directory / ZARR_NAME, mode="w", compute=False)

    return Path(root_directory / ZARR_NAME)
"""
