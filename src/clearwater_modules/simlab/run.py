import numpy as np
import xsimlab as xs
import xarray as xr
import pandas as pd

from pathlib import Path


### This is just a development script, so far from official documentation
#But I'll start getting this into something that is more representative of
#a notebook based workflow.

###########################################################
# Step 1 - define Riverine model inputs
###########################################################

#model pathway
data_root = Path(r"C:\Users\ptomasula\Repositories\ClearWater-modules\data_temp\sumwere_creek_coarse_p48")
config_path = data_root / "demo_config.yml"
riverine_path = xr.DataArray(
    [config_path],
    dims=("path"),
)

#model date range
start_index = int(
    8 * 60 * (60 / 30)
)  # start at 8:00 am on the first day of the simulation (30 second model)
end_index = start_index + int(8 * 60 * (60 / 30))

###########################################################
# Step 2 - Define ClearWater Modules Inputs
###########################################################

### Import note
# By default, xarray should auto align coordinates (e.g. time) by filling them in with NaNs.
# But xarray-simlab does not appear to support that behavior, and instead rejects the inputs and fails.
# To get around this, we need to generate an input xarray.Dataset object and then pass or inputs from that array to the xsimlab model.

hourly = xr.date_range("2025-01-01", "2025-01-02", freq="h")
minute_5 = xr.date_range("2025-01-01", "2025-01-02", freq="5min")

### temperature inputs
#TODO: these should be loaded in from time series but I have question
#why are these datetime indexed, but the riverrine inputs are integer indexed?
#how did you resolve the two?

#solar_radiation_path = data_root / "cwr_boundary_conditions_q_Solar_p28.csv"
#solar_radiation = pd.read_csv(solar_radiation_path)
#solar_radiation = solar_radiation.iloc[start_index:end_index]  # filter rows to match simulation period
#solar_radiation = xr.DataArray(
#    solar_radiation.values,
#    dims=("time", "solar_radiation"),
#)
# placeholder, but long term we need to load this from a file
solar = pd.DataFrame({'time': hourly, 'solar': np.random.rand(hourly.size)})

#air_temperature_path = data_root / "cwr_boundary_conditions_TairC_p28.csv"
#air_temperature = pd.read_csv(air_temperature_path)
#air_temperature = air_temperature.iloc[start_index:end_index]  # filter rows to match simulation period
# TODO: air temperature is sampled at a different frequency than the other datasets
# ideally we find a way where the frequency doesn't matter, but I'm not sure how best to handle that
#air_temperature = air_temperature.resample("1h").mean()
air_temperature = pd.DataFrame({'time': minute_5, 'air_temperature': np.random.rand(minute_5.size)})

water_temperature = pd.DataFrame({'time': minute_5, 'water_temperature': np.random.rand(minute_5.size)})

###########################################################
## Step 3 - Consolidate inputs into single aligned dataset
###########################################################
inputs = xr.Dataset(
    {
        "solar_radiation": xr.DataArray(
            solar["solar"].values,
            coords={"time": solar["time"]},
            dims=("time"),
        ),
        "air_temperature": xr.DataArray(
            air_temperature["air_temperature"].values,
            coords={"time": air_temperature["time"]},
            dims=("time"),
        ),
    },
    coords={"time": minute_5}
)

###########################################################
## Step 4 - Define Model Steps (Processes)
###########################################################
#TODO: move this back up top, it's just that the import of Riverine is really slow and this is killing iterations
import processes

model = xs.Model(
    {
        "transport": processes.Riverine,
        "temperature": processes.Temperature,
    }
)

### map inputs to model
input_dataset = xs.create_setup(
    model=model,
    clocks={"time": minute_5},
    main_clock="time",
    input_vars={
        #transport model (riverine)
        "transport__riverine_path": riverine_path,  # path,
        "transport__datetime_range": ("indicies", [start_index, end_index]),
        #'transport__mesh': (mesh.values),
        #temperature model
        "temperature__solar_radiation": inputs["solar_radiation"],
        "temperature__air_temperature": inputs["air_temperature"],
    },
    output_vars={
        "transport__riverine_path": None,
        #"water_temperature": None,
        #"transport__volume": None,
        #"transport__wetted_surface_area": None,
    },
)


###########################################################
## Step 5 - Execute Model
###########################################################

output_dataset = input_dataset.xsimlab.run(model=model)

print(output_dataset)
