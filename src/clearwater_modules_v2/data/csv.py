from data.base import DataProvider
from custom_types import ArrayLike
from pathlib import Path
from datetime import datetime, timedelta
import pandas as pd


class CSVDataProvider(DataProvider):
    def __init__(self, **kwargs) -> None:
        self.file_path = kwargs.pop("file_path")
        self.time_field = kwargs.pop("time_field", None)
        self.interpolation_method = kwargs.pop("interpolation_method", "linear")
        self.__data: ArrayLike | None = None

    def read(self) -> ArrayLike:
        df = pd.read_csv(self.file_path)

        # rename time field
        if self.time_field is not None:
            df = df.rename(columns={self.time_field: "time"})
        else:
            df = df.rename(index={0: "time"})

        # convert time to datetime
        df["time"] = pd.to_datetime(df["time"])
        df = df.set_index("time")

        return df.to_xarray()

    def write_to_store(
        self,
        store_path: Path,
        start_time: datetime,
        end_time: datetime,
        time_step: timedelta,
        variable_name: str,
        field_name: str | None = None,
    ) -> None:
        # check if data has been loaded
        if self.__data is None:
            self.__data = self.read()

        # rename the field to the variable
        self.__data = self.__data.rename({field_name: variable_name})

        # extract the array
        array_out = self.__data[variable_name]

        # bound the data to the model time range
        array_out = array_out.sel(time=slice(start_time, end_time))

        # align time dimension to model specification
        # TODO: Do we need to handle downsampling?
        array_out = array_out.resample(time=time_step).interpolate(
            self.interpolation_method
        )

        # add scalar dimension to align with zarr input template
        array_out = array_out.expand_dims({"scalar": 1})

        # write this to the variable specific data array in our zarr store
        array_out.to_zarr(store_path, mode="a")
