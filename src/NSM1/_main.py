import time
from collections import OrderedDict
from numba.typed import Dict
from numba import types

st=time.time()

from _ApDeath import ApDeath
from _ApRespiration import ApRespiration
from _ApGrowth import ApGrowth

from _algae import Algae
from _nitrogen import Nitrogen
from _benthic_algae import BenthicAlgae
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
output_variables = Dict.empty(key_type=types.unicode_type, value_type=types.float64)
output_variables = {

} 

#True/False module use, user defined

global_module_choices = Dict.empty(key_type=types.unicode_type, value_type=types.float64)
global_module_choices = {
    'use_Algae': True,
    'use_NH4': True,
    'use_NO3': True,
    'use_TIP': True,
    'use_POC': False,
    'use_DOC': False,

    'use_BAlgae': True,
    'use_OrgN' : True,
    'use_OrgP' : True,

    'use_SedFlux' : False,
    'use_DOX': True,

    'use_DIC': False,
    'use_N2' : False,
    'use_Pathogen' : False,
    'use_Alk' : False,
    'use_POM2' : False,
}
        
#User-defined global variables
global_vars = Dict.empty(key_type=types.unicode_type, value_type=types.float64)
global_vars = {
    'Ap': 100,
    'NH4':100,
    'NO3': 100,
    'TIP':100,
    'TwaterC':20,
    'depth':1,

    'Ab':100,

    'DOX' : 100,
    'OrgN' : 100,
    'vson':1,

    'lambda': 1,
    'fdp': 0.5,
    'PAR': 100
}

#Algae Module Optional Changes 
#TODO potentially make three columes which distinguishes the module it is for
algae_constant_changes = Dict.empty(key_type=types.unicode_type, value_type=types.float64)
algae_constant_changes = {
        
#   'AWd':100,              
#   'AWc':40,               
#   'AWn':7.2,              
#   'AWp': 1,               
#   'AWa':1000,        

    'KL':10,
    'KsN':0.04,
    'KsP':0.0012,
    'mu_max': 1,
    'kdp':0.15,
    'krp': 0.2,
    'vsap':0.15,

    'growth_rate_option':1,
    'light_limitation_option': 1
}

#Benthic algae module optional changes
Balgae_constant_changes = Dict.empty(key_type=types.unicode_type, value_type=types.float64)
Balgae_constant_changes = {
#    'BWd': 100,       
#    'BWc': 40,      
#    'Bwn' : 7.2,      
#    'BWp' : 1,      
#    'BWa' : 3500,       

#    'KLb': 10,       
#    'KsNb' : 0.25,       
#    'KsPb' : 0.125,    
#    'Ksb' : 10,       

#    'mub_max' : 0.4,   
#    'krb' : 0.2,     
#    'kdb': 0.3,   
            
#    'b_growth_rate_option' : 1,     
#    'b_light_limitation_option' : 1

}

#Nitrogen Module Optional Changes
nitrogen_constant_changes = Dict.empty(key_type=types.unicode_type, value_type=types.float64)
nitrogen_constant_changes = {
        
#   'vson' : 0.01,
#   'KNR' : 0,
#   'knit' : 0.1,
#   'kon'  : 0.1,
#   'kdnit' : 0.002,
#   'rnh4'  : 0,
#   'KsOxdn' : 0.1
        
#   'PN' : 0.5 
#   'PNB' : 0.5 
#   'Fw' : 0.9 
#   'Fb' : 0.9 

}



#Create some univresal variables
if global_module_choices['use_Algae'] or global_module_choices['use_NH4'] or global_module_choices['use_NO3'] or global_module_choices['use_OrgN'] :
    #Call growth rate calculation
    ApGrowth_rate = ApGrowth(global_vars['Ap'], global_vars['TwaterC'], global_vars['PAR'], global_vars['lambda'], global_vars['depth'], global_vars['NH4'], 
        global_vars['NO3'], global_vars['TIP'], global_vars['fdp'], algae_constant_changes['mu_max'], algae_constant_changes['KsN'], algae_constant_changes['KsP'], 
        algae_constant_changes['KL'], global_module_choices['use_NH4'], global_module_choices['use_NO3'], global_module_choices['use_TIP'], algae_constant_changes['light_limitation_option'], 
        algae_constant_changes['growth_rate_option']).Calculation()
   
    #Call death rate calculation
    ApDeath_rate = ApDeath(algae_constant_changes['kdp'], global_vars['TwaterC'], global_vars['Ap']).Calculation()
    
    #Call respiration calculation
    ApRespiration_rate=ApRespiration(algae_constant_changes['krp'], global_vars['TwaterC'], global_vars['Ap']).Calculation()

#TODO should I send extra variables in a dictionary or each individual variable seperatly, or each individual variable in a dicitonary
#Call Algae module
if global_module_choices['use_Algae'] :
    output_variables['dApdt'] = Algae(global_vars['Ap'], global_vars['depth'], algae_constant_changes['vsap'], ApGrowth_rate, ApDeath_rate, ApRespiration_rate).Calculations()

#Call Benthic Algae module
if global_module_choices['use_BAlgae'] :
# Call Benthic Algea 
# TEMP here until integrate Balgae module to get the variables for nitrogen
    dAbdt, Balgae_pathways = BenthicAlgae(global_module_choices, global_vars, Balgae_constant_changes).Calculations()
else :
    Balgae_pathways = {}

#Call Phosphorus module
if global_module_choices['use_OrgP'] or global_module_choices['use_TIP']:
    pass

#Call carbon module
if global_module_choices['use_POC'] or global_module_choices['use_DOC'] or global_module_choices['use_DIC'] :
    pass

#Call Sediment Flux module
if global_module_choices['use_SedFlux'] :
    #TEMP here unitl integrate Sediment Flux module to get the variables for nitrogen
    sedFlux_pathways = {
        'JNH4': 100,
        'JNO3': 100
    }
else : 
    sedFlux_pathways ={}

#Call Nitrogen module
if global_module_choices['use_NH4'] or global_module_choices['use_NO3'] or global_module_choices['use_OrgN'] :

    output_variables['DIN'], output_variables['TON'], output_variables['TKN'], output_variables['TN']= Nitrogen(global_module_choices, global_vars, Balgae_pathways, sedFlux_pathways, nitrogen_constant_changes).Calculations()

#TODO create for alkalinity, DOX, CBOD, N2, Pathogen, and POM
 
et = time.time()
elapsed_time = et - st 

print('Execution time:', elapsed_time, 'seconds')

