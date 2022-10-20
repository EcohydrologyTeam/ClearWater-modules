from collections import OrderedDict

import unittest
from numba.typed import Dict
from numba import types
import os
import sys

from _algae import Algae

class Test_algae(unittest.TestCase):

    def setUp(self) :
        self.global_module_choices = Dict.empty(key_type=types.unicode_type, value_type=types.float64)
        self.global_module_choices = {
            'use_Algae': True,
            'use_NH4': True,
            'use_NO3': True,
            'use_TIP': True,
            'use_POC': False,
            'use_DOC': False,

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
        
        self.algae_constant_changes = Dict.empty(key_type=types.unicode_type, value_type=types.float64)
        self.algae_constant_changes = {
                
        #   'AWd':100,              
        #   'AWc':40,               
        #   'AWn':7.2,              
        #   'AWp': 1,               
        #   'AWa':1000,        

        #    'KL':10,
        #    'KsN':0.04,
        #    'KsP':0.0012,
        #    'mu_max': 1,
        #    'kdp':0.30,
        #    'krp': 0.2,
        #    'vsap':0.15,

        #    'growth_rate_option':1,
        #    'light_limitation_option': 1
        }
  
#Change NH4 and KL
    def test_NH4_KL (self):

        self.global_vars['NH4']=200
        self.algae_constant_changes['KL']=20

        Change,path=Algae(self.global_module_choices,self.global_vars, self.algae_constant_changes).Calculations()

        self.assertAlmostEqual(Change, 24.8049, 3)

#Change NO3_lambda
    def test_NO3_lambda (self) :

        self.global_vars['NO3']= 200
        self.global_vars['lambda'] = 2

        Change,path = Algae(self.global_module_choices,self.global_vars, self.algae_constant_changes).Calculations()

        self.assertAlmostEqual(Change, 27.0906, 3)
 
#Change Ksn and Depth
    def test_Ksn_Depth (self) :

        self.global_vars['depth']= 2
        self.algae_constant_changes['KsN'] = 0.08

        Change, path = Algae(self.global_module_choices,self.global_vars, self.algae_constant_changes).Calculations()

        self.assertAlmostEqual(Change, 34.570, 3)

#Change mu_max and PAR
    def test_mu_max_PAR (self) :

        self.global_vars['PAR']= 200
        self.algae_constant_changes['mu_max'] = 2

        Change,path = Algae(self.global_module_choices,self.global_vars, self.algae_constant_changes).Calculations()

        self.assertAlmostEqual(Change, 134.229, 3)

#Change Ap
    def test_Ap(self) :

        self.global_vars['Ap']= 200

        Change,path = Algae(self.global_module_choices,self.global_vars, self.algae_constant_changes).Calculations()

        self.assertAlmostEqual(Change, 70.93267, 3)

#Change light_limitation and growth_rate options
    def test_light_growth_option(self) :

        self.algae_constant_changes['light_limitation_option']= 2
        self.algae_constant_changes['growth_rate_option'] =2

        Change,path = Algae(self.global_module_choices,self.global_vars, self.algae_constant_changes).Calculations()

        self.assertAlmostEqual(Change, 48.4313, 3)

#Change light_limitation and growth_rate options and KL + fdp
    def test_LLO_GRO_KL_fdp(self) :

        self.algae_constant_changes['light_limitation_option']= 2
        self.algae_constant_changes['growth_rate_option'] =2
        self.algae_constant_changes['KL'] = 20
        self.global_vars['fdp'] = 0.75

        Change,path = Algae(self.global_module_choices,self.global_vars, self.algae_constant_changes).Calculations()

        self.assertAlmostEqual(Change, 44.28346, 3)

#Change light_limitation and growth_rate options and KsP + lambda
    def test_LLO_GRO_KsP_lambda(self) :

        self.algae_constant_changes['light_limitation_option']= 2
        self.algae_constant_changes['growth_rate_option'] =2
        self.algae_constant_changes['KsP'] = 0.0024
        self.global_vars['lambda'] = 4

        Change,path = Algae(self.global_module_choices,self.global_vars, self.algae_constant_changes).Calculations()

        self.assertAlmostEqual(Change, 20.38781, 3)

#Change light_limitation and growth_rate options and TIP + depth
    def test_LLO_GRO_TIP_depth(self) :

        self.algae_constant_changes['light_limitation_option']= 2
        self.algae_constant_changes['growth_rate_option'] =2
        self.global_vars['TIP'] = 200
        self.global_vars['depth'] = 5

        Change,path = Algae(self.global_module_choices,self.global_vars, self.algae_constant_changes).Calculations()

        self.assertAlmostEqual(Change, 20.60617, 3)   

#Change light_limitation and growth_rate options and PAR
    def test_LLO_GRO_PAR(self) :

        self.algae_constant_changes['light_limitation_option']= 2
        self.algae_constant_changes['growth_rate_option'] =2
        self.global_vars['PAR'] = 50

        Change,path = Algae(self.global_module_choices,self.global_vars, self.algae_constant_changes).Calculations()

        self.assertAlmostEqual(Change, 44.2834564, 3)   

#Change light_limitation and growth_rate options and Ap + TwaterC
    def test_LLO_GRO_Ap_TwaterC(self) :

        self.algae_constant_changes['light_limitation_option']= 2
        self.algae_constant_changes['growth_rate_option'] =2
        self.global_vars['Ap'] = 200
        self.global_vars['TwaterC'] = 10

        Change,path = Algae(self.global_module_choices,self.global_vars, self.algae_constant_changes).Calculations()

        self.assertAlmostEqual(Change, 50.14322, 3)   

#Change light_limitation and growth_rate options and KL
    def test_LLO_GRO_KL(self) :

        self.algae_constant_changes['light_limitation_option']= 3
        self.algae_constant_changes['growth_rate_option'] =3
        self.algae_constant_changes['KL'] = 5

        Change,path = Algae(self.global_module_choices,self.global_vars, self.algae_constant_changes).Calculations()

        self.assertAlmostEqual(Change, -49.8267, 3)   

#Change light_limitation and growth_rate options and krp + lambda
    def test_LLO_GRO_krp_lambda(self) :

        self.algae_constant_changes['light_limitation_option']= 3
        self.algae_constant_changes['growth_rate_option'] =3
        self.algae_constant_changes['krp'] = 0.4
        self.global_vars['lambda'] = 0.5

        Change,path = Algae(self.global_module_choices,self.global_vars, self.algae_constant_changes).Calculations()

        self.assertAlmostEqual(Change, -68.7626, 3)   

#Change light_limitation and growth_rate options and depth
    def test_LLO_GRO_depth(self) :

        self.algae_constant_changes['light_limitation_option']= 3
        self.algae_constant_changes['growth_rate_option'] =3
        self.global_vars['depth'] = 0.5

        Change,path = Algae(self.global_module_choices,self.global_vars, self.algae_constant_changes).Calculations()

        self.assertAlmostEqual(Change, -63.7626, 3)   

#Change light_limitation and growth_rate options and vsap + PAR
    def test_LLO_GRO_vsap_PAR(self) :

        self.algae_constant_changes['light_limitation_option']= 3
        self.algae_constant_changes['growth_rate_option'] =3
        self.global_vars['PAR'] = 200
        self.algae_constant_changes['vsap'] = 0.3

        Change,path = Algae(self.global_module_choices,self.global_vars, self.algae_constant_changes).Calculations()

        self.assertAlmostEqual(Change, -64.826684, 3)   

#Change light_limitation and growth_rate options and vsap + PAR
    def test_LLO_GRO_kdp_Ap(self) :

        self.algae_constant_changes['light_limitation_option']= 3
        self.algae_constant_changes['growth_rate_option'] =3
        self.global_vars['Ap'] = 200
        self.algae_constant_changes['kdp'] = 0.3

        Change,path = Algae(self.global_module_choices,self.global_vars, self.algae_constant_changes).Calculations()

        self.assertAlmostEqual(Change, -116.298, 3)  

if __name__ == '__main__':
    unittest.main()
