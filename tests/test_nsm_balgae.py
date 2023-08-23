from collections import OrderedDict

import unittest
from numba.typed import Dict
from numba import types
import os
import sys

from _benthic_algae import BenthicAlgae

class Test_Balgae(unittest.TestCase):

    def setUp(self) :
        self.global_module_choices = Dict.empty(key_type=types.unicode_type, value_type=types.float64)
        self.global_module_choices = {
            'use_Algae': True,
            'use_NH4': True,
            'use_NO3': True,
            'use_TIP': True,
            'use_POC': True,
            'use_DOC': True,

            'use_BAlgae': True,
            'use_OrgN' : True,
            'use_OrgP' : True,

            'use_SedFlux' : True,
            'use_DOX': True,

            'use_DIC': False,
            'use_N2' : False,
            'use_Pathogen' : False,
            'use_Alk' : False,
            'use_POM2' : False,
        }

        self.global_vars = Dict.empty(key_type=types.unicode_type, value_type=types.float64)
        self.global_vars = {
            'Ap': 100,
            'NH4':100,
            'NO3': 100,
            'TIP':100,
            'TwaterC':20,
            'depth':1,

            'Ab':100,

            'DOX' : 100,
            'OrgN' : 200,
            'vson': 0.01,

            'lambda': 1,
            'fdp': 0.5,
            'PAR': 100
        }
                
        self.Balgae_constant_changes = Dict.empty(key_type=types.unicode_type, value_type=types.float64)
        self.Balgae_constant_changes = {
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

#Original test
    def test_original (self):

        Change,path=BenthicAlgae(self.global_module_choices,self.global_vars, self.Balgae_constant_changes).Calculations()
        self.assertAlmostEqual(Change, -47.1515, 3)
#Change Ab
    def test_Ab (self):
        self.global_vars['Ab'] = 200
        Change,path=BenthicAlgae(self.global_module_choices,self.global_vars, self.Balgae_constant_changes).Calculations()
        self.assertAlmostEqual(Change, -97.0159, 3)

#Change mub_max
    def test_mub_max (self):
        self.Balgae_constant_changes['mub_max'] = 0.2
        Change,path=BenthicAlgae(self.global_module_choices,self.global_vars, self.Balgae_constant_changes).Calculations()
        self.assertAlmostEqual(Change, -48.57576, 3)

#Change growth_rate_option
    def test_GRO (self):
        self.Balgae_constant_changes['b_growth_rate_option'] = 2
        Change,path=BenthicAlgae(self.global_module_choices,self.global_vars, self.Balgae_constant_changes).Calculations()
        self.assertAlmostEqual(Change, -47.14797, 3)

#Change light_limitation_option
    def test_LLO (self):
        self.Balgae_constant_changes['b_light_limitation_option'] = 2
        Change,path=BenthicAlgae(self.global_module_choices,self.global_vars, self.Balgae_constant_changes).Calculations()
        self.assertAlmostEqual(Change, -46.50408888, 3)

if __name__ == '__main__':
    unittest.main()
