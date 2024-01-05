# memory profiling

# import packages
import sys
import numpy as np
import pandas as pd
import holoviews as hv
import xarray as xr

from clearwater_modules.tsm.model import EnergyBudget



# SET PATHS
root = r'C:\Users\sjordan\OneDrive - LimnoTech\Documents\GitHub\ClearWater-riverine\examples\dev_sandbox\data\sumwere_test_cases\plan28_testTSM'


def interpolate_to_model_timestep(meteo_df, model_dates, col_name):
    merged_df = pd.merge_asof(
        pd.DataFrame({'time': model_dates}),
        meteo_df,
        left_on='time',
        right_on='Datetime')
    merged_df[col_name.lower()] = merged_df[col_name].interpolate(method='linear')
    merged_df.drop(
        columns=['Datetime', col_name],
        inplace=True)
    merged_df.rename(
        columns={'time': 'Datetime'},
        inplace=True,
    )
    return merged_df

def run_n_timesteps_profiled(
    time_steps: int,
    reaction: EnergyBudget,
    transport: xr.Dataset,
    meteo_params,
    concentration_update = None
):
    for i in range(1, time_steps):
        if i % 1000 == 0:
            print(i)
        # increment the timestep of the reaction model, using meteo parameters + output from transport model
        reaction.increment_timestep()


def initialize_clearwater_riverine():
    transport_model = xr.open_zarr(
        r'W:\2ERDC12 - Clearwater\Clearwater_testing_TSM\plan28_testTSM_pwrPlnt_May2022\full_test_output\mesh_output_full_2023_12_20.zarr',
    )
    return transport_model

def define_meteo_params(transport_model):
    
    xarray_time_index = pd.DatetimeIndex(
        transport_model.time.values
    )

    q_solar = pd.read_csv(
        f'{root}/cwr_boundary_conditions_q_solar_p28.csv', 
        parse_dates=['Datetime'])
    q_solar.dropna(axis=0, inplace=True)

    q_solar_interp = interpolate_to_model_timestep(
        q_solar,
        xarray_time_index,
        'q_Solar'
    )

    air_temp_c = pd.read_csv(
        f'{root}/cwr_boundary_conditions_TairC_p28.csv', 
        parse_dates=['Datetime'])
    air_temp_c.dropna(axis=0, inplace=True)

    air_temp_c_interp = interpolate_to_model_timestep(
        air_temp_c,
        xarray_time_index,
        'TairC'
    )

    air_temp_c_interp['air_temp_c'] = (air_temp_c_interp.tairc - 32)* (5/9)


    q_solar_array = q_solar_interp.q_solar.to_numpy()
    air_temp_array = air_temp_c_interp.air_temp_c.to_numpy()


    # for each individual timestep
    all_meteo_params = {
        'q_solar': q_solar_array,
        'air_temp_c': air_temp_array,
    }


    # for initial conditions
    initial_meteo_params = {
        'air_temp_c': air_temp_array[0],
        'q_solar': q_solar_array[0],
    }
    return all_meteo_params, initial_meteo_params

def initialize_clearwater_modules(track_dynamic, initial_state_values, initial_meteo_params):
    # updateable static variable = unexpected keyword argument
    reaction_model = EnergyBudget(
        initial_state_values,
        time_dim='seconds',
        meteo_parameters= initial_meteo_params,
        track_dynamic_variables = track_dynamic,
        use_sed_temp = True,
        updateable_static_variables=['air_temp_c', 'q_solar']
    )
    return reaction_model


# Clearwater-Riverine
# Paths
def main(iters, track_dynamic):

    transport_model = initialize_clearwater_riverine()

    initial_state_values = {
        'water_temp_c': transport_model['concentration'].isel(time=0, nface=slice(0, transport_model.nreal + 1))* 0 + 100,
        'volume': transport_model['volume'].isel(time=0, nface=slice(0, transport_model.nreal + 1))* 0 + 8,
        'surface_area': transport_model['faces_surface_area'].isel(nface=slice(0, transport_model.nreal + 1))* 0 + 4,
    }


    all_meteo_params, initial_meteo_params = define_meteo_params(transport_model)

    reaction_model = initialize_clearwater_modules(track_dynamic, initial_state_values, initial_meteo_params)

    run_n_timesteps_profiled(
        iters,
        reaction_model,
        transport_model,
        all_meteo_params,
    )


if __name__ == '__main__':
    if len(sys.argv) > 1:
        try:
            iters = int(sys.argv[1])
            print(f'Running {iters} iterations.')
        except ValueError:
            raise ValueError('Argument must be an integer # of iterations.')
        try:
            track_dynamic = sys.argv[2]
            print(track_dynamic)
            if track_dynamic == True:
                print('Tracking dynamic variables.')
            else:
                print('Not tracking dynamic variables.')
        except ValueError:
            track_dynamic = False
            print('Defaulting to not tracking dynamic variables.')

    else:
        print('No argument given, defaulting to 100 iteration and not tracking dynamic variables.')
        iters = 100
        track_dynamic = False

            
    main(iters=iters, track_dynamic=track_dynamic)
