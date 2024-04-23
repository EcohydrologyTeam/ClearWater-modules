"""A script to allow for debugging of the NSM module."""
import time
import sys
import clearwater_modules as cwm
from clearwater_modules.nsm1.model import NutrientBudget




initial_state_values = {
    'Ap': 1,
    'Ab': 1,
    'NH4': 1,
    'NO3': 1,
    'OrgN': 1,
    'N2': 1,
    'TIP': 1,
    'OrgP': 1,
    'POC': 1,
    'DOC': 1,
    'DIC': 1,
    'POM': 1,
    'CBOD': 1,
    'DOX': 1,
    'PX': 1,
    'Alk': 1
}

algae_parameters = {
    'AWd': 100,
    'AWc': 40,
    'AWn': 7.2,
    'AWp': 1,
    'AWa': 1000,
    'KL': 10,
    'KsN': 0.04,
    'KsP': 0.0012,
    'mu_max_20': 1,
    'kdp_20': 0.15,
    'krp_20': 0.2,
    'vsap': 0.15,
    'growth_rate_option': 1,
    'light_limitation_option': 1,
    'lambda0': .5,
    'lambda1': .5,
    'lambda2': .5,
    'lambdas': .5,
    'lambdam': .5, 
    'Fr_PAR': .5  
}

balgae_parameters = {
    'BWd': 100,
    'BWc': 40,
    'BWn': 7.2,
    'BWp': 1,
    'BWa': 3500,

    'KLb': 10,
    'KsNb': 0.25,
    'KsPb': 0.125,
    'Ksb': 10,
    'mub_max_20': 0.4,
    'krb_20': 0.2,
    'kdb_20': 0.3,
    'b_growth_rate_option': 1,
    'b_light_limitation_option': 1,
    'Fw': 0.9,
    'Fb': 0.9
}

nitrogen_parameters = {
    'KNR': 0.6,
    'knit_20': 0.1,
    'kon_20': 0.1,
    'kdnit_20': 0.002,
    'rnh4_20': 0,
    'vno3_20': 0,
    'KsOxdn': 0.1,
    'PN': 0.5,
    'PNb': 0.5
}

phosphorus_parameters = {
    'kop_20': 0.1,
    'rpo4_20': 0
}

POM_parameters = {
    'kpom_20': 0.1
}

CBOD_parameters = {
    'KsOxbod': 0.5,
    'kbod_20': 0.12,
    'ksbod_20': 0
}

carbon_parameters = {
    'f_pocp': 0.9,
    'kdoc_20': 0.01,
    'f_pocb': 0.9,
    'kpoc_20': 0.005,
    'K_sOxmc': 1,
    'pCO2': 383,
    'FCO2': 0.2
}

pathogen_parameters = {
    'kdx': 0.8,
    'apx': 1,
    'vx': 1
}

alkalinity_parameters = {
    'r_alkaa': 1,
    'r_alkan': 1,
    'r_alkn': 1,
    'r_alkden': 1,
    'r_alkba': 1,
    'r_alkbn': 1 
}

global_parameters = {
    'use_NH4': True,
    'use_NO3': True, 
    'use_OrgN': True,
    'use_TIP': True,  
    'use_SedFlux': False,
    'use_DOX': True,
    'use_Algae': True,
    'use_Balgae': True,
    'use_OrgP': True,
    'use_POC': True,
    'use_DOC': True,
    'use_DIC': True,
    'use_N2': True,
    'use_Pathogen': True,
    'use_Alk': True,
    'use_POM': True 
}


global_vars = {
    'vson': 0.01,
    'vsoc': 0.01,
    'vsop': 999,
    'vs': 999,
    'SOD_20': 999,
    'SOD_theta': 999,
    'vb': 0.01,
    'fcom': 0.4,
    'kaw_20_user': 999,
    'kah_20_user': 999,
    'hydraulic_reaeration_option': 2,
    'wind_reaeration_option': 2,    
    'timestep': 86400,
    'depth': 1,
    'TwaterC': 20,
    'theta': 1.047,
    'velocity': 1,
    'flow': 2,
    'topwidth': 1,
    'slope': 2,
    'shear_velocity': 4,
    'pressure_atm': 2,
    'wind_speed': 4,
    'q_solar': 4,
    'Solid': 1, 
}

DOX_parameters = {
    'DOX': 6.5,
}
N2_parameters = {}

def main(iters: int):
    ti = time.time()
    # define starting state values
    nsm_model = NutrientBudget(
        initial_state_values=initial_state_values,  # mandatory
        algae_parameters=algae_parameters,
        alkalinity_parameters=alkalinity_parameters,
        balgae_parameters=balgae_parameters,
        carbon_parameters=carbon_parameters,
        CBOD_parameters=CBOD_parameters,
        DOX_parameters=DOX_parameters,
        nitrogen_parameters=nitrogen_parameters,
        POM_parameters=POM_parameters,
        N2_parameters=N2_parameters,
        phosphorus_parameters=phosphorus_parameters,
        pathogen_parameters=pathogen_parameters,
        global_parameters=global_parameters,
        global_vars=global_vars,  
        track_dynamic_variables=True,  # default is true
        hotstart_dataset=None,  # default is None
        time_dim='year',  # default is "timestep"
    )
    # print(nsm_model.get_variable_names())
    print(nsm_model.dynamic_variables_names)

    for _ in range(iters):
        nsm_model.increment_timestep()


if __name__ == '__main__':
    if len(sys.argv) > 1:
        try:
            iters = int(sys.argv[1])
            print(f'Running {iters} iterations.')
        except ValueError:
            raise ValueError('Argument must be an integer # of iterations.')
    else:
        print('No argument given, defaulting to 100 iteration.')
        iters = 100
            
    main(iters=iters)
