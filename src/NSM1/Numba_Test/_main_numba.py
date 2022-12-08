
import time
from collections import OrderedDict
from numba import types, typed, njit

from _algae_numba import Calculations
#from _nitrogen import Nitrogen
#from _benthic_algae import BenthicAlgae
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

'''
@njit
def create_dict_float(items):
    dictionary=typed.Dict.empty(key_type=types.unicode_type, value_type=types.float64)
    print (items)
    for k,v in items:
        dictionary[k] = v
    #    print(k)

    return dictionary

nb_gv = create_dict_float(tuple(global_vars.items()))
'''
st=time.time()

@njit(cache = True)
def creat_dicts() : 
    #True/False module use, user defined
    global_module_choices =typed.Dict.empty(key_type=types.unicode_type, value_type=types.boolean)
    global_module_choices['use_Algae'] = True
    global_module_choices['use_NH4'] = True
    global_module_choices['use_NO3'] = True
    global_module_choices['use_TIP'] = True
    global_module_choices['use_POC'] = True
    global_module_choices['use_DOC'] = True
    global_module_choices['use_OrgN'] = True
    global_module_choices['use_OrgP'] = True
    global_module_choices['use_SedFlux'] = False
    global_module_choices['use_DOX'] = True
    global_module_choices['use_DIC'] = False
    global_module_choices['use_N2'] = False
    global_module_choices['use_Pathogen'] = False
    global_module_choices['use_Alk'] = False
    global_module_choices['use_POM2'] = False
            
    #User-defined global variables
    global_vars =typed.Dict.empty(key_type=types.unicode_type, value_type=types.float64)
    global_vars['Ap'] = 100.0
    global_vars['NH4'] = 100.0
    global_vars['NO3'] = 100.0
    global_vars['TIP'] = 100.0
    global_vars['TwaterC'] = 20.0
    global_vars['depth'] = 1.0

    global_vars['Ab'] = 100.0
    global_vars['DOX'] = 100.0
    global_vars['OrgN'] = 100.0
    global_vars['vson'] = 0.01
    global_vars['lambda'] = 1.0
    global_vars['fdp'] = 0.5
    global_vars['PAR'] = 100.0

    #Algae Module Optional Changes 
    algae_constant_changes = typed.Dict.empty(key_type=types.unicode_type, value_type=types.float64)
    algae_constant_changes['AWd'] = 100
    algae_constant_changes['AWc'] = 40
    algae_constant_changes['AWn'] = 7.2
    algae_constant_changes['AWp'] = 1
    algae_constant_changes['AWa'] = 1000

    algae_constant_changes['KL'] = 10
    algae_constant_changes['KsN'] = 0.04
    algae_constant_changes['KsP'] = 0.0012
    algae_constant_changes['mu_max'] = 1
    algae_constant_changes['kdp'] = 0.15
    algae_constant_changes['krp'] = 0.2
    algae_constant_changes['vsap'] = 0.15
    algae_constant_changes['growth_rate_option'] = 1
    algae_constant_changes['light_limitation_option'] = 1
    return global_module_choices, global_vars, algae_constant_changes

#Benthic algae module optional changes
Balgae_constant_changes = OrderedDict()
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
nitrogen_constant_changes = OrderedDict()
nitrogen_constant_changes = {
        
#   'vson' : 0.01,
#   'KNR' : 0,
#   'knit' : 0.1,
#   'kon'  : 0.1,
#   'kdnit' : 0.002,
#   'rnh4'  : 0,
#   'KsOxdn' : 0.1
#   'vno3' : 0
        
#   'PN' : 0.5 
#   'PNB' : 0.5 
#   'Fw' : 0.9 
#   'Fb' : 0.9 

}

#TODO should I send extra variables in a dictionary or each individual variable seperatly, or each individual variable in a dicitonary
#Call Algae module


global_module_choices, global_vars, algae_constant_changes = creat_dicts()

if global_module_choices['use_Algae'] :
    #output_variables['dApdt'], algae_pathways = 
    Calculations(global_module_choices, global_vars, algae_constant_changes)   
else:
    algae_pathways={}

et = time.time()
elapsed_time = et - st 

print('Execution time:', elapsed_time, 'seconds')

'''
#Call Benthic Algae module
if global_module_choices['use_BAlgae'] :
# Call Benthic Algea 
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
        'JNH4': 2,
        'JNO3': 2
    }
else : 
    sedFlux_pathways ={}

#Call Nitrogen module
if global_module_choices['use_NH4'] or global_module_choices['use_NO3'] or global_module_choices['use_OrgN'] :
    output_variables['DIN'], output_variables['TON'], output_variables['TKN'], output_variables['TN'], output_variables['dOrgNdt'], output_variables['dNH4dt'], output_variables['dNO3dt']= Nitrogen(global_module_choices, global_vars, algae_pathways, Balgae_pathways, sedFlux_pathways, nitrogen_constant_changes).Calculations()
  
#TODO create for alkalinity, DOX, CBOD, N2, Pathogen, and POM
''' 


