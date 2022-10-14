import time
from collections import OrderedDict

st=time.time()

from _algae import Algae
# from _nitrogen import Nitrogen
# from _benthic_algae import BenthicAlgae
# from _carbon import Carbon
# from _cbod import CBOD
# from _dox import DOX
# from _n2 import N2
# from _pathogen import Pathogen
# from _phosphorus import Phosphorus
# from _pom import POM
# from _sed_flux import SedFlux
# from _alkalinity import Alkalinity

#Variables to return
output_variables = OrderedDict()
output_variables = {

} 

#Two classes for global variables (true/false) modules to use and user defined variables
class variables :
    def __init__(self) :
        pass
    
    #True/False module use, user defined
    def assign_global_module_choices (self):
        self.global_module_choices = OrderedDict()
        self.global_module_choices = {
            'use_Algae': True,
            'use_NH4': True,
            'use_NO3': True,
            'use_TIP': True,
            'use_POC': False,
            'use_DOC': False,

            'use_BAlgae': False,
            'use_OrgN' : True,
            'use_OrgP' : True,

            'use_SedFlux' : False,
            'use_DOX': True,

            'use_DIC': False,
            'use_N2' : False,
            'use_Pathogen' : False,
            'use_Alk' : False,
            'use_POM2' : False

        }
        
        return self.global_module_choices

    #User-defined global variables
    def assign_global_vars (self):
        self.globals_vars = OrderedDict()
        self.globals_vars = {
            'Ap': 100,
            'NH4':100,
            'NO3':100,
            'TIP':100,
            'TwaterC':20,
            'depth':1,

            'Ab':100,

            'DOX' : 100,

            'lambda': 1,
            'fdp': 0.5,
            'PAR': 100
        }

        return self.globals_vars

    #Algae Module Optional Changes
    def assign_algae_constant_changes (self):
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

    #Nitrogen Module Optional Changes
    def assign_nitorgen_constant_changes (self):
        self.nitrogen_constant_changes = OrderedDict()
        self.nitrogen_constant_changes = {
        
        #   'vson' : 0.01,
        #   'KNR' : 0,
        #   'knit' : 0.1,
        #   'kon'  : 0.1,
        #   'kdnit' : 0.002,
        #   'rnh4'  : 0,
        #   'KsOxdn' : 0.1
        
    #TODO should these go in this function or come out of Algae/BAlgae?

        #   'PN' : 0.5 (ALSO IN ALGAE CONSTANT_CHANGE)
        #   'PNB' : 0.5 (ALSO IN BENTHIC ALGAE CONSTANT_CHANGE)
        #   'Fw' : 0.9 (ALSO IN BENTHIC ALGAE CONSTANT_CHANGE)
        #   'Fb' : 0.9 (ALSO IN BNETHIC ALGAE CONSTANT_CHANGE)

        }
        return self.nitrogen_constant_changes

global_vars = variables().assign_global_vars()
algae_constant_changes = variables().assign_algae_constant_changes()
global_module_choices = variables().assign_global_module_choices()
nitrogen_constant_changes = variables().assign_nitorgen_constant_changes()

if global_module_choices['use_Algae'] :
    dApdt, algae_rates = Algae(global_vars, algae_constant_changes, global_module_choices).Calculations()

if global_module_choices['use_BAlgae'] :
# Call Benthic Algea 
# TEMP here until integrate Balgae module to get the variables for nitrogen
    Balgae_rates = {
        'rnb' : 100,
        'AbGrowth' : 100,
        'AbDeath' : 100,
        'AbRespiration' : 100,
    }

if global_module_choices['use_OrgP'] or global_module_choices['use_TIP']:
    pass

if global_module_choices['use_POC'] or global_module_choices['use_DOC'] or global_module_choices['use_DIC'] :
    pass

if global_module_choices['use_SedFlux'] :
    #TEMP here unitl integrate Sediment Flux module to get the variables for nitrogen
    sed_Flux = {
        'JNH4': 100,
        'JNO3': 100
    }

if global_module_choices['use_NH4'] or global_module_choices['use_NO3'] or global_module_choices['use_OrgN'] :
    #Nitrogen(algae_rates, Balgae_rates, sed_Flux, global_module_choices, global_vars, nitrogen_constant_changes).Calculations
    pass 

et = time.time()
elapsed_time = et - st 
print('Execution time:', elapsed_time, 'seconds')

