import numpy as np
import xsimlab as xs
import xarray as xr
import pandas as pd

import processes

from pathlib import Path


data_root = Path(r"C:\Users\ptomasula\Repositories\ClearWater-modules\data_temp\sumwere_creek_coarse_p48")

### riverine inputs
config_path = data_root / "demo_config.yml"
riverine_path = xr.DataArray(
    [config_path],
    dims=("path"),
)

start_index = int(
    8 * 60 * (60 / 30)
)  # start at 8:00 am on the first day of the simulation (30 second model)
end_index = start_index + int(8 * 60 * (60 / 30))


### temperature inputs



#solar_radiation = xr.DataArray(
#    [np.arange(start_index,end_index,1),np.zeros(end_index - start_index)],
#    dims=("time", "solar_radiation"),
#)

#air_temperature = xr.DataArray(
#    [np.arange(start_index,end_index,1),np.zeros(end_index - start_index)],
#    dims=("time", "air_temperature"),
#)

#TODO: these should be loaded in from time series but I have question
#why are these datetime indexed, but the riverrine inputs are integer indexed?
#how did you resolve the two?


solar_radiation_path = data_root / "cwr_boundary_conditions_q_Solar_p28.csv"
solar_radiation = pd.read_csv(solar_radiation_path)
solar_radiation = solar_radiation.iloc[start_index:end_index]  # filter rows to match simulation period
solar_radiation = xr.DataArray(
    solar_radiation.values,
    dims=("time", "solar_radiation"),
) 


air_temperature_path = data_root / "cwr_boundary_conditions_TairC_p28.csv"
air_temperature = pd.read_csv(air_temperature_path)
air_temperature = air_temperature.iloc[start_index:end_index]  # filter rows to match simulation period
# TODO: air temperature is sampled at a different frequency than the other datasets
# ideally we find a way where the frequency doesn't matter, but I'm not sure how best to handle that
#air_temperature = air_temperature.resample("1h").mean()

air_temperature = xr.DataArray(
    air_temperature.values,
    dims=("time", "air_temperature"),
) 


### define model structure

model = xs.Model(
    {
        "transport": processes.Riverine,
        "temperature": processes.Temperature,
    }
)

### define model inputs
input_dataset = xs.create_setup(
    model=model,
    clocks={"step": np.arange(start_index, end_index, 1)},
    main_clock="step",
    input_vars={
        #transport model (riverine)
        "transport__riverine_path": riverine_path,  # path,
        "transport__datetime_range": ("indicies", [start_index, end_index]),
        #'transport__mesh': (mesh.values),
        #temperature model
        "temperature__solar_radiation": solar_radiation,
        "temperature__air_temperature": air_temperature,
    },
    output_vars={
        #"water_temperature": None,
        "transport__volume": None,
        "transport__wetted_surface_area": None,
    },
)

output_dataset = input_dataset.xsimlab.run(model=model)

print(output_dataset)
