
import time
from collections import OrderedDict

st=time.time()

from _algae import Algae
from _nitrogen import Nitrogen
from _benthic_algae import BenthicAlgae
from _phosphorus import Phosphorus
# from _carbon import Carbon
from _cbod import CBOD
# from _dox import DOX
# from _n2 import N2
# from _pathogen import Pathogen
# from _pom import POM
# from _sed_flux import SedFlux
from _alkalinity import Alkalinity

#Variables to return
output_variables = OrderedDict()
output_variables = {

} 

#True/False module use, user defined
global_module_choices =OrderedDict
global_module_choices = {
    'use_Algae' : True,
    'use_BAlgae': True,
    'use_NH4' : True,
    'use_NO3' : True,
    'use_TIP' : True,
    'use_POC' : True,
    'use_DOC' : True,
    'use_OrgN' : True,
    'use_OrgP' : True,
    'use_SedFlux' : False,
    'use_DOX' : True,
    'use_DIC' : False,
    'use_N2' : False,
    'use_Pathogen' : False,
    'use_Alk' : True,
    'use_POM2' : False,
    'use_CBOD': True,
}

#User-defined global variables
global_vars =OrderedDict()
global_vars = {
    #Algae
    'Ap' : 100.0,
    'NH4' : 100.0,
    'NO3' : 100.0,
    'TIP' : 100.0,
    'TwaterC' : 25.0,
    'depth' : 1.0,

    #Benthic algae
    'Ab' : 100.0,

    #Nitrogen
    'DOX' : 100.0,
    'OrgN' : 100.0,

    #Phosphrous
    'OrgP' : 100,

    #CBOD
    'CBOD' : 100,

    #Parameters
    'lambda' : 1.0,
    'fdp' : 0.5,
    'PAR' : 100.0,
    'vs' : 1,
    'vson' : 0.01,
    'vsop' : 0.01,
    'vb' : 0.0025, 
    'h2' : 0.1,

    #SedFlux
    'Salinity': 100,
    'dt' : 0.1,
    'POM2' : 100,
    'TsedC' : 100,

    #Alkalinity
    'Alk' : 100,
    'pH' : 7,
    'DIC' : 100,

}

#User-defined global variables only used in SedFlux
global_var_sedflux = OrderedDict()
global_var_sedflux = {
    'NH41' : 100,
    'NO31' : 100,
    'TIP1' : 100,
    'CH41' : 100,
    'SO41' : 100,
    'TH2S1' : 100,
    'DIC1' : 100,
    
    'NH42' : 100,
    'NO32' : 100,
    'TIP2' : 100,
    'CH42' : 100,
    'SO42' : 100,
    'TH2S2' : 100,
    'DIC2' : 100,

    'POC2_1' : 100,
    'PON2_1' : 100,
    'POP2_1' : 100,

    'POC2_2' : 100,
    'PON2_2' : 100,
    'POP2_2' : 100,

    'POC2_3' : 100,
    'PON2_3' : 100,
    'POP2_3' : 100,

    # 't' : 1,

}
#Algae Module Optional Changes 
algae_constant_changes=OrderedDict()
algae_constant_changes = {
#    'AWd': 100,
#   'AWc' : 40,
#    'AWn' : 7.2,
#    'AWp' : 1,
#    'AWa' : 1000,

#    'KL' : 10,
#    'KsN' : 0.04,
#    'KsP' : 0.0012,
#    'mu_max' : 1,
#    'kdp' : 0.15,
#    'krp' : 0.2,
#    'vsap' : 0.15,
#    'growth_rate_option' : 1,
#    'light_limitation_option' : 1
}

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

#   'Fw' : 0.9 
#   'Fb' : 0.9 

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


}

#Phosphorus Module Optional Changes
phosphorous_constant_changes = OrderedDict()
phosphorous_constant_changes = {
        
#   'kop' : 0.1,
#   'rpo4' : 0

}

CBOD_constant_changes = OrderedDict()
CBOD_constant_changes= {
    # 'kbod' : 0.12,
    # 'ksbod' : 0,
    # 'KsOxbod' : 0.5,
}

alkalinity_constant_changes=OrderedDict()
alkalinity_constant_changes = {
#    'ralkca' : 14/106/12/1000, #translating algal and balgal growht into Alk if NH2 is the N source (eq/ug-chla)
#    'ralkcn' : 18/106/12/1000, #ratio translating algal and balgal growth into Alk if NO3 is the N source (eq/ug-Chla)
#    'ralkn' : 2/14/1000, #nitrification
#    'ralkden' : 4/14/1000, # denitrification
#    'pH_solution' : 2,
#    'imax' : 13, # maximum iteration number for computing pH
#    'es' : 0.003, # maximum relative error for computing pH
}

#Call Algae module
if global_module_choices['use_Algae'] :
    output_variables['dApdt'], algae_pathways = Algae(global_module_choices, global_vars, algae_constant_changes).Calculations()
else:
    algae_pathways={}

#Call Benthic Algae module
if global_module_choices['use_BAlgae'] :
# Call Benthic Algea 
    output_variables ['dAbdt'], balgae_pathways = BenthicAlgae(global_module_choices, global_vars, Balgae_constant_changes).Calculations()
else :
    balgae_pathways = {}

#Call CBOD module
if global_module_choices['use_CBOD'] :
    output_variables['dCBODdt'] = CBOD(global_vars['CBOD'],global_vars['TwaterC'], global_vars['DOX'], global_module_choices['use_DOX'], CBOD_constant_changes).Calculation()

#Call Sediment Flux module
if global_module_choices['use_SedFlux'] :
    #TEMP here unitl integrate Sediment Flux module to get the variables for nitrogen
    sedFlux_pathways = {
        'JNH4': 2,
        'JNO3': 2,
        'JDIP' : 2,
    }
else : 
    sedFlux_pathways ={}

#Call Phosphorus module
if global_module_choices['use_OrgP'] or global_module_choices['use_TIP']:
    output_variables['dOrgPdt'], output_variables['dTIPdt'], output_variables['TOP'], output_variables['TP'] \
        = Phosphorus(global_module_choices, global_vars, algae_pathways, balgae_pathways, sedFlux_pathways, phosphorous_constant_changes).Calculation()

#Call Nitrogen module
if global_module_choices['use_NH4'] or global_module_choices['use_NO3'] or global_module_choices['use_OrgN'] :
    output_variables['DIN'], output_variables['TON'], output_variables['TKN'], output_variables['TN'], output_variables['dOrgNdt'], output_variables['dNH4dt'], output_variables['dNO3dt'], nitrogen_pathways = Nitrogen(global_module_choices, global_vars, algae_pathways, balgae_pathways, sedFlux_pathways, nitrogen_constant_changes).Calculations()
else:
    nitrogen_pathways = {}

#Call Alkalinity module and pH
if global_module_choices['use_Alk'] :
    output_variables['dAlkdt'], output_variables['pH'] = Alkalinity(global_module_choices, global_vars, algae_pathways, balgae_pathways, nitrogen_pathways, alkalinity_constant_changes).Calculations_Alk()
else:


#TODO create for carbon, DOX, N2, Pathogen, and POM

et = time.time()
elapsed_time = et - st 

print('Execution time:', elapsed_time, 'seconds')
