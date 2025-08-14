"""Intermediate zarr data for model inputs"""

import xarray as xr
from pathlib import Path
from datetime import datetime, timedelta
import pandas as pd
from dask import array as da


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
    template_dataset.to_zarr(
        root_directory / "model_data.zarr", mode="w", compute=False
    )

    return Path(root_directory / "model_data.zarr")
