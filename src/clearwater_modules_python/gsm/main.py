
import time
from collections import OrderedDict
from general_constituents import GeneralConstituent

output_variables = OrderedDict()
output_variables = {
}

global_vars = OrderedDict ()
global_vars = {
    "GC": 100,
    "TwaterC": 20,
    "depth": 1.0,
}

gsm_constant_changes = OrderedDict()
gsm_constant_changes = {
    "k_rc20": 1,
    "k_theta" : 1.047,

    "rgc_rc20" : 1,
    "rgc_theta" : 1.047,
    "order" : 1,
    "release" : True,
    "settling" : True,
    "settling_rate" : 0.002,
}

output_variables['dGCdt'] = GeneralConstituent(global_vars, gsm_constant_changes).Calculations()