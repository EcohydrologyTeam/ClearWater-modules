import numpy as np
import xarray as xr
import pandas as pd
from pathlib import Path
import os
import pickle

# import clearwater_riverine as cwr
from datetime import datetime, timedelta

from utils import conversions

from processes import Riverine
from processes import Temperature
from variables import VariableRegistry


import clearwater_riverine as cwr

# init a clearwater instance
data_root = Path(
    r"C:\Users\ptomasula\Repositories\ClearWater-modules\data_temp\sumwere_creek_coarse_p48"
)
config_path = data_root / "demo_config.yml"

# maybe p38 for a coarse model?

START_DATETIME = "05-14-2022 00:00:00"
END_DATETIME = "05-15-2022 00:00:00"

fp = "riverine.plk"
if os.path.exists(fp):
    with open(fp, "rb") as f:
        transport = cwr.ClearwaterRiverine(
            config_filepath=config_path, datetime_range=(START_DATETIME, END_DATETIME)
        )
else:
    transport = cwr.ClearwaterRiverine(
        config_filepath=config_path, datetime_range=(START_DATETIME, END_DATETIME)
    )
    with open(fp, "wb") as f:
        pickle.dump(transport, f)

riverine_process = Riverine(transport)

### Prepare inputs
# define input variables
variables = VariableRegistry()

hourly = xr.date_range(START_DATETIME, END_DATETIME, freq="h")
minute_5 = xr.date_range(START_DATETIME, END_DATETIME, freq="5min")

### solar radiation
solar = pd.read_csv(data_root / "cwr_boundary_conditions_q_Solar_p28.csv")
# rename columns
solar = solar.rename(columns={"Datetime": "time", "q_Solar": "solar"})
# convert to datetime
solar["time"] = pd.to_datetime(solar["time"])
solar = solar.set_index("time")
# resample to 5 minute data
solar = solar.resample("5min").interpolate()

# air temperature
air_temperature = pd.read_csv(data_root / "cwr_boundary_conditions_TairC_p28.csv")
# rename columns
air_temperature = air_temperature.rename(
    columns={"Datetime": "time", "TairC": "air_temperature"}
)
# convert to datetime
air_temperature["time"] = pd.to_datetime(air_temperature["time"])
air_temperature = air_temperature.set_index("time")
# air is already 5 minute data

# other meteorological variables
cloudiness = pd.DataFrame(
    {"time": minute_5, "cloudiness": np.ones(minute_5.size) * 0.1}
)
cloudiness = cloudiness.set_index("time")

wind_speed = pd.DataFrame({"time": minute_5, "wind_speed": np.ones(minute_5.size) * 3})
wind_speed = wind_speed.set_index("time")

atmospheric_pressure = pd.DataFrame(
    {"time": minute_5, "atmospheric_pressure": np.ones(minute_5.size) * 1013.0}
)
atmospheric_pressure = atmospheric_pressure.set_index("time")

atmospheric_vapor_pressure = pd.DataFrame(
    {
        "time": minute_5,
        "atmospheric_vapor_pressure": np.ones(minute_5.size) * 1013.0
        - np.random.normal(0, 1, minute_5.size),
    }
)
atmospheric_vapor_pressure = atmospheric_vapor_pressure.set_index("time")

inputs = xr.Dataset(
    {
        "solar_radiation": xr.DataArray(
            solar["solar"].values,
            coords={"time": solar.index},
            dims=("time"),
        ),
        "air_temperature": xr.DataArray(
            air_temperature["air_temperature"].values,
            coords={"time": air_temperature.index},
            dims=("time"),
        ),
        "cloudiness": xr.DataArray(
            cloudiness["cloudiness"].values,
            coords={"time": cloudiness.index},
            dims=("time"),
        ),
        "wind_speed": xr.DataArray(
            wind_speed["wind_speed"].values,
            coords={"time": wind_speed.index},
            dims=("time"),
        ),
        "atmospheric_pressure": xr.DataArray(
            atmospheric_pressure["atmospheric_pressure"].values,
            coords={"time": atmospheric_pressure.index},
            dims=("time"),
        ),
        "atmospheric_vapor_pressure": xr.DataArray(
            atmospheric_vapor_pressure["atmospheric_vapor_pressure"].values,
            coords={"time": atmospheric_vapor_pressure.index},
            dims=("time"),
        ),
    },
    coords={"time": minute_5},
)

variables.register("solar_radiation", inputs["solar_radiation"])
variables.register("air_temperature", inputs["air_temperature"])
variables.register("cloudiness", inputs["cloudiness"])
variables.register("wind_speed", inputs["wind_speed"])
variables.register("atmospheric_pressure", inputs["atmospheric_pressure"])
variables.register("atmospheric_vapor_pressure", inputs["atmospheric_vapor_pressure"])
variables.register("sediment_temperature", 20.0)  # degrees celsius
variables.register("sediment_thickness", 0.1)  # meters

## define our model
from model import Model

model = Model(
    # processes=[riverine_process],
    processes=(riverine_process, Temperature(wind_a=0.3, wind_b=1.5, wind_c=3.0)),
    variables=variables,
    start_time=datetime(2022, 5, 14),
    end_time=datetime(2022, 5, 15),
    time_step=timedelta(seconds=30),
)

model.run()

prt = 1
