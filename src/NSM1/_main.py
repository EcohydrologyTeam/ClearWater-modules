import time
from collections import OrderedDict

st=time.time()

from _algae import Algae

class create_global_constant_changes :
    def __init__(self) :
        pass
    def assign (self):
        self.global_constant_changes = OrderedDict()
        self.global_constant_changes = {
         #    'use_Algae': True,
         #   'use_NH4': True,
         #   'use_NO3': True,
         #   'use_TIP': True,
         #   'use_POC': False,
         #   'use_DOC': False,

        }
        

        return self.global_constant_changes

class create_globals_vars :
    def __init__(self) :
        pass
    def assign (self):
        self.globals_vars = OrderedDict()
        self.globals_vars = {
            'Ap': 100,
            'NH4':100,
            'NO3':100,
            'TIP':100,
            'TwaterC':20,
            'depth':1,

            'Ab':100,

            'lambda': 1,
            'fdp': 0.5,
            'PAR': 100
        }

        return self.globals_vars

class create_algae_constant_changes:
    def __init__(self) :
        pass
    def assign (self):
        self.algae_constant_changes = OrderedDict()
        self.algae_constant_changes = {
        
        #   'AWd':100,              
        #   'AWc':40,               
        #   'AWn':7.2,              
        #   'AWp': 1,               
        #   'AWa':1000,        

        #   'KL':10,
        #   'KsN':0.04,
        #   'KsP':0.0012,
        #   'PN':0.5,
        #   'Fpocp':0.9,
        #   'mu_max': 1,
        #   'kdp':0.15,
        #   'krp': 0.2,
        #   'vsap':0.15,

        #   'growth_rate_option':1,
        #   'light_limitation_option': 1

    }

        return self.algae_constant_changes


class run_script:

    def __init__(self):
        pass

global_vars = create_globals_vars().assign()
algae_constant_changes = create_algae_constant_changes().assign()
global_constant_changes = create_global_constant_changes().assign()

Algae(global_vars, algae_constant_changes, global_constant_changes).Calculations()


et = time.time()
elapsed_time = et - st 
print('Execution time:', elapsed_time, 'seconds')

