from collections import OrderedDict

import unittest
from numba.typed import Dict
from numba import types

from _phosphorus import Phosphorus

class Test_phosphorus(unittest.TestCase):

#Define inital array values
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

            'use_SedFlux' : False,
            'use_DOX': True,

            'use_DIC': False,
            'use_N2' : False,
            'use_Pathogen' : False,
            'use_Alk' : False,
            'use_POM2' : False,
        }

        self.global_vars = Dict.empty(key_type=types.unicode_type, value_type=types.float64)
        self.global_vars = {
            #Algae
            'Ap' : 100.0,
            'NH4' : 100.0,
            'NO3' : 100.0,
            'TIP' : 100.0,
            'TwaterC' : 20.0,
            'depth' : 1.0,

            #Benthic algae
            'Ab' : 100.0,

            #Nitrogen
            'DOX' : 100.0,
            'OrgN' : 100.0,
            'vson' : 0.01,

            #Phosphrous
            'OrgP' : 100,
            'vs' : 1,
            'vsop' : 0.01,

            #Parameters
            'lambda' : 1.0,
            'fdp' : 0.5,
            'PAR' : 100.0
        }
        
        self.algae_pathways = Dict.empty(key_type=types.unicode_type, value_type=types.float64)
        self.algae_pathways={
            'ApGrowth': 10,
            'ApDeath' : 20,
            'ApRespiration': 30,
            'rpa' : 0.5
        } 

        self.Balgae_pathways= Dict.empty(key_type=types.unicode_type, value_type=types.float64)
        self.Balgae_pathways = {
            'AbGrowth' : 30,
            'AbDeath' : 20,
            'AbRespiration' : 10,
            'rpb' : 0.25,
            'Fw' : 0.9,
            'Fb' : 0.9,
        }

        self.sedFlux_pathways = Dict.empty(key_type=types.unicode_type, value_type=types.float64)
        self.sedFlux_pathways = {
            'JDIP': 2,
        }

        self.P_constant_changes = Dict.empty(key_type=types.unicode_type, value_type=types.float64)
        self.P_constant_changes = {     
        #   'kop' : 0.1,
        #   'rpo4' : 0
        }
  
#Original test
    def test_phosphorus (self):

        dOrgPdt, dTIPdt, TOP, TP = Phosphorus(self.global_module_choices, self.global_vars, self.algae_pathways, self.Balgae_pathways, self.sedFlux_pathways, self.P_constant_changes).Calculation()
        self.assertAlmostEqual(dOrgPdt, 3.05, 3)
        self.assertAlmostEqual(dTIPdt, -34.25, 3)
        self.assertAlmostEqual(TOP, 150, 3)
        self.assertAlmostEqual(TP, 250, 3)

#Change OrgP
    def test_OrgP (self):
        
        self.global_vars['OrgP'] = 200
        dOrgPdt, dTIPdt, TOP, TP = Phosphorus(self.global_module_choices, self.global_vars, self.algae_pathways, self.Balgae_pathways, self.sedFlux_pathways, self.P_constant_changes).Calculation()
        self.assertAlmostEqual(dOrgPdt, -7.95, 3)
        self.assertAlmostEqual(dTIPdt, -24.25, 3)
        self.assertAlmostEqual(TOP, 250, 3)
        self.assertAlmostEqual(TP, 350, 3)

#Change rpo4
    def test_rpo4 (self):
        self.P_constant_changes['rpo4'] = 0.4
        dOrgPdt, dTIPdt, TOP, TP = Phosphorus(self.global_module_choices, self.global_vars, self.algae_pathways, self.Balgae_pathways, self.sedFlux_pathways, self.P_constant_changes).Calculation()
        self.assertAlmostEqual(dOrgPdt, 3.05, 3)
        self.assertAlmostEqual(dTIPdt, -33.85, 3)
        self.assertAlmostEqual(TOP, 150, 3)
        self.assertAlmostEqual(TP, 250, 3)

if __name__ == '__main__':
    unittest.main()
