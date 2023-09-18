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
            'use_TIP': True,
            'use_BAlgae': True,
            'use_OrgP' : True,
            'use_SedFlux' : False,

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

#Change rpo4
    def test_rpo4 (self):
        self.P_constant_changes['rpo4'] = 0.5
        dOrgPdt, dTIPdt, TOP, TP = Phosphorus(self.global_module_choices, self.global_vars, self.algae_pathways, self.Balgae_pathways, self.sedFlux_pathways, self.P_constant_changes).Calculation()
        self.assertAlmostEqual(dOrgPdt, 3.05, 3)
        self.assertAlmostEqual(dTIPdt, -33.75, 3)
        self.assertAlmostEqual(TOP, 150, 3)
        self.assertAlmostEqual(TP, 250, 3)

#Change kop
    def test_kop (self):
        self.P_constant_changes['kop'] = 0.3
        dOrgPdt, dTIPdt, TOP, TP = Phosphorus(self.global_module_choices, self.global_vars, self.algae_pathways, self.Balgae_pathways, self.sedFlux_pathways, self.P_constant_changes).Calculation()
        self.assertAlmostEqual(dOrgPdt, -16.95, 3)
        self.assertAlmostEqual(dTIPdt, -14.25, 3)
        self.assertAlmostEqual(TOP, 150, 3)
        self.assertAlmostEqual(TP, 250, 3)

#Change vsop
    def test_vsop (self):
        self.global_vars['vsop'] = 0.1
        dOrgPdt, dTIPdt, TOP, TP = Phosphorus(self.global_module_choices, self.global_vars, self.algae_pathways, self.Balgae_pathways, self.sedFlux_pathways, self.P_constant_changes).Calculation()
        self.assertAlmostEqual(dOrgPdt, -5.95, 3)
        self.assertAlmostEqual(dTIPdt, -34.25, 3)
        self.assertAlmostEqual(TOP, 150, 3)
        self.assertAlmostEqual(TP, 250, 3)

#Change TIP
    def test_TIP (self):
        self.global_vars['TIP'] = 50
        dOrgPdt, dTIPdt, TOP, TP = Phosphorus(self.global_module_choices, self.global_vars, self.algae_pathways, self.Balgae_pathways, self.sedFlux_pathways, self.P_constant_changes).Calculation()
        self.assertAlmostEqual(dOrgPdt, 3.05, 3)
        self.assertAlmostEqual(dTIPdt, -9.25, 3)
        self.assertAlmostEqual(TOP, 150, 3)
        self.assertAlmostEqual(TP, 200, 3)

#Change OrgP
    def test_OrgP (self):
        self.global_vars['OrgP'] = 200
        dOrgPdt, dTIPdt, TOP, TP = Phosphorus(self.global_module_choices, self.global_vars, self.algae_pathways, self.Balgae_pathways, self.sedFlux_pathways, self.P_constant_changes).Calculation()
        self.assertAlmostEqual(dOrgPdt, -7.95, 3)
        self.assertAlmostEqual(dTIPdt, -24.25, 3)
        self.assertAlmostEqual(TOP, 250, 3)
        self.assertAlmostEqual(TP, 350, 3)

#Change depth
    def test_depth (self):
        self.global_vars['depth'] = 5
        dOrgPdt, dTIPdt, TOP, TP = Phosphorus(self.global_module_choices, self.global_vars, self.algae_pathways, self.Balgae_pathways, self.sedFlux_pathways, self.P_constant_changes).Calculation()
        self.assertAlmostEqual(dOrgPdt, 0.61, 3)
        self.assertAlmostEqual(dTIPdt, 11.15, 3)
        self.assertAlmostEqual(TOP, 150, 3)
        self.assertAlmostEqual(TP, 250, 3)

#Change TwaterC
    def test_TwaterC (self):
        self.global_vars['TwaterC'] = 35
        dOrgPdt, dTIPdt, TOP, TP = Phosphorus(self.global_module_choices, self.global_vars, self.algae_pathways, self.Balgae_pathways, self.sedFlux_pathways, self.P_constant_changes).Calculation()
        self.assertAlmostEqual(dOrgPdt, -6.8659, 3)
        self.assertAlmostEqual(dTIPdt, -24.334, 3)
        self.assertAlmostEqual(TOP, 150, 3)
        self.assertAlmostEqual(TP, 250, 3)

#Change vs
    def test_vs (self):
        self.global_vars['vs'] = 4
        dOrgPdt, dTIPdt, TOP, TP = Phosphorus(self.global_module_choices, self.global_vars, self.algae_pathways, self.Balgae_pathways, self.sedFlux_pathways, self.P_constant_changes).Calculation()
        self.assertAlmostEqual(dOrgPdt, 3.05, 3)
        self.assertAlmostEqual(dTIPdt, -184.25, 3)
        self.assertAlmostEqual(TOP, 150, 3)
        self.assertAlmostEqual(TP, 250, 3)

#Change fdp
    def test_fdp (self):
        self.global_vars['fdp'] = 0.25
        dOrgPdt, dTIPdt, TOP, TP = Phosphorus(self.global_module_choices, self.global_vars, self.algae_pathways, self.Balgae_pathways, self.sedFlux_pathways, self.P_constant_changes).Calculation()
        self.assertAlmostEqual(dOrgPdt, 3.05, 3)
        self.assertAlmostEqual(dTIPdt, -59.25, 3)
        self.assertAlmostEqual(TOP, 150, 3)
        self.assertAlmostEqual(TP, 250, 3)

#Change Ap
    def test_Ap(self):
        self.global_vars['Ap'] = 50
        dOrgPdt, dTIPdt, TOP, TP = Phosphorus(self.global_module_choices, self.global_vars, self.algae_pathways, self.Balgae_pathways, self.sedFlux_pathways, self.P_constant_changes).Calculation()
        self.assertAlmostEqual(dOrgPdt, 3.05, 3)
        self.assertAlmostEqual(dTIPdt, -34.25, 3)
        self.assertAlmostEqual(TOP, 125, 3)
        self.assertAlmostEqual(TP, 225, 3)

#Change rpa
    def test_rpa(self):
        self.algae_pathways['rpa'] = 0.25
        dOrgPdt, dTIPdt, TOP, TP = Phosphorus(self.global_module_choices, self.global_vars, self.algae_pathways, self.Balgae_pathways, self.sedFlux_pathways, self.P_constant_changes).Calculation()
        self.assertAlmostEqual(dOrgPdt, -1.95, 3)
        self.assertAlmostEqual(dTIPdt, -39.25, 3)
        self.assertAlmostEqual(TOP, 125, 3)
        self.assertAlmostEqual(TP, 225, 3)

#Change ApDeath
    def test_ApDeath(self):
        self.algae_pathways['ApDeath'] = 70
        dOrgPdt, dTIPdt, TOP, TP = Phosphorus(self.global_module_choices, self.global_vars, self.algae_pathways, self.Balgae_pathways, self.sedFlux_pathways, self.P_constant_changes).Calculation()
        self.assertAlmostEqual(dOrgPdt, 28.05, 3)
        self.assertAlmostEqual(dTIPdt, -34.25, 3)
        self.assertAlmostEqual(TOP, 150, 3)
        self.assertAlmostEqual(TP, 250, 3)

#Change ApRespiration
    def test_ApRespiration(self):
        self.algae_pathways['ApRespiration'] = 15
        dOrgPdt, dTIPdt, TOP, TP = Phosphorus(self.global_module_choices, self.global_vars, self.algae_pathways, self.Balgae_pathways, self.sedFlux_pathways, self.P_constant_changes).Calculation()
        self.assertAlmostEqual(dOrgPdt, 3.05, 3)
        self.assertAlmostEqual(dTIPdt, -41.75, 3)
        self.assertAlmostEqual(TOP, 150, 3)
        self.assertAlmostEqual(TP, 250, 3)

#Change ApGrowth
    def test_ApGrowth(self):
        self.algae_pathways['ApGrowth'] = 30
        dOrgPdt, dTIPdt, TOP, TP = Phosphorus(self.global_module_choices, self.global_vars, self.algae_pathways, self.Balgae_pathways, self.sedFlux_pathways, self.P_constant_changes).Calculation()
        self.assertAlmostEqual(dOrgPdt, 3.05, 3)
        self.assertAlmostEqual(dTIPdt, -44.25, 3)
        self.assertAlmostEqual(TOP, 150, 3)
        self.assertAlmostEqual(TP, 250, 3)

#Change rpb
    def test_rpb(self):
        self.Balgae_pathways['rpb'] = 0.75
        dOrgPdt, dTIPdt, TOP, TP = Phosphorus(self.global_module_choices, self.global_vars, self.algae_pathways, self.Balgae_pathways, self.sedFlux_pathways, self.P_constant_changes).Calculation()
        self.assertAlmostEqual(dOrgPdt, 11.15, 3)
        self.assertAlmostEqual(dTIPdt, -42.75, 3)
        self.assertAlmostEqual(TOP, 150, 3)
        self.assertAlmostEqual(TP, 250, 3)

#Change Fw
    def test_Fw(self):
        self.Balgae_pathways['Fw'] = 0.5
        dOrgPdt, dTIPdt, TOP, TP = Phosphorus(self.global_module_choices, self.global_vars, self.algae_pathways, self.Balgae_pathways, self.sedFlux_pathways, self.P_constant_changes).Calculation()
        self.assertAlmostEqual(dOrgPdt, 1.25, 3)
        self.assertAlmostEqual(dTIPdt, -34.25, 3)
        self.assertAlmostEqual(TOP, 150, 3)
        self.assertAlmostEqual(TP, 250, 3)

#Change Fb
    def test_Fb(self):
        self.Balgae_pathways['Fb'] = 0.3
        dOrgPdt, dTIPdt, TOP, TP = Phosphorus(self.global_module_choices, self.global_vars, self.algae_pathways, self.Balgae_pathways, self.sedFlux_pathways, self.P_constant_changes).Calculation()
        self.assertAlmostEqual(dOrgPdt, 0.35, 3)
        self.assertAlmostEqual(dTIPdt, -29.75, 3)
        self.assertAlmostEqual(TOP, 150, 3)
        self.assertAlmostEqual(TP, 250, 3)

#Change AbGrowth
    def test_AbGrowth(self):
        self.Balgae_pathways['AbGrowth'] = 20
        dOrgPdt, dTIPdt, TOP, TP = Phosphorus(self.global_module_choices, self.global_vars, self.algae_pathways, self.Balgae_pathways, self.sedFlux_pathways, self.P_constant_changes).Calculation()
        self.assertAlmostEqual(dOrgPdt, 3.05, 3)
        self.assertAlmostEqual(dTIPdt, -32.0, 3)
        self.assertAlmostEqual(TOP, 150, 3)
        self.assertAlmostEqual(TP, 250, 3)

#Change AbRespiration
    def test_AbRespiration(self):
        self.Balgae_pathways['AbRespiration'] = 50
        dOrgPdt, dTIPdt, TOP, TP = Phosphorus(self.global_module_choices, self.global_vars, self.algae_pathways, self.Balgae_pathways, self.sedFlux_pathways, self.P_constant_changes).Calculation()
        self.assertAlmostEqual(dOrgPdt, 3.05, 3)
        self.assertAlmostEqual(dTIPdt, -24.25, 3)
        self.assertAlmostEqual(TOP, 150, 3)
        self.assertAlmostEqual(TP, 250, 3)

#Change AbDeath
    def test_AbDeath(self):
        self.Balgae_pathways['AbDeath'] = 10
        dOrgPdt, dTIPdt, TOP, TP = Phosphorus(self.global_module_choices, self.global_vars, self.algae_pathways, self.Balgae_pathways, self.sedFlux_pathways, self.P_constant_changes).Calculation()
        self.assertAlmostEqual(dOrgPdt, 1.025, 3)
        self.assertAlmostEqual(dTIPdt, -34.25, 3)
        self.assertAlmostEqual(TOP, 150, 3)
        self.assertAlmostEqual(TP, 250, 3)

#Change JDIP and use_SedFlux
    def test_JDIP_use_SedFlux(self):
        self.global_module_choices['use_SedFlux'] = True
        self.sedFlux_pathways['JDIP'] = 5
        dOrgPdt, dTIPdt, TOP, TP = Phosphorus(self.global_module_choices, self.global_vars, self.algae_pathways, self.Balgae_pathways, self.sedFlux_pathways, self.P_constant_changes).Calculation()
        self.assertAlmostEqual(dOrgPdt, 3.05, 3)
        self.assertAlmostEqual(dTIPdt, -29.25, 3)
        self.assertAlmostEqual(TOP, 150, 3)
        self.assertAlmostEqual(TP, 250, 3)

#Change use_OrgP
    def test_use_OrgP(self):
        self.global_module_choices['use_OrgP'] = False

        dOrgPdt, dTIPdt, TOP, TP = Phosphorus(self.global_module_choices, self.global_vars, self.algae_pathways, self.Balgae_pathways, self.sedFlux_pathways, self.P_constant_changes).Calculation()
        self.assertAlmostEqual(dOrgPdt, 0, 3)
        self.assertAlmostEqual(dTIPdt, -44.25, 3)
        self.assertAlmostEqual(TOP, 50, 3)
        self.assertAlmostEqual(TP, 150, 3)

#Change use_TIP
    def test_use_OrgP(self):
        self.global_module_choices['use_TIP'] = False

        dOrgPdt, dTIPdt, TOP, TP = Phosphorus(self.global_module_choices, self.global_vars, self.algae_pathways, self.Balgae_pathways, self.sedFlux_pathways, self.P_constant_changes).Calculation()
        self.assertAlmostEqual(dOrgPdt, 3.05, 3)
        self.assertAlmostEqual(dTIPdt, 0, 3)
        self.assertAlmostEqual(TOP, 150, 3)
        self.assertAlmostEqual(TP, 150, 3)

#Change use_Algae
    def test_use_Algae(self):
        self.global_module_choices['use_Algae'] = False

        dOrgPdt, dTIPdt, TOP, TP = Phosphorus(self.global_module_choices, self.global_vars, self.algae_pathways, self.Balgae_pathways, self.sedFlux_pathways, self.P_constant_changes).Calculation()
        self.assertAlmostEqual(dOrgPdt, -6.95, 3)
        self.assertAlmostEqual(dTIPdt, -44.25, 3)
        self.assertAlmostEqual(TOP, 100, 3)
        self.assertAlmostEqual(TP, 200, 3)        

#Change use_BAlgae
    def test_use_BAlgae(self):
        self.global_module_choices['use_BAlgae'] = False

        dOrgPdt, dTIPdt, TOP, TP = Phosphorus(self.global_module_choices, self.global_vars, self.algae_pathways, self.Balgae_pathways, self.sedFlux_pathways, self.P_constant_changes).Calculation()
        self.assertAlmostEqual(dOrgPdt, -1, 3)
        self.assertAlmostEqual(dTIPdt, -30, 3)
        self.assertAlmostEqual(TOP, 150, 3)
        self.assertAlmostEqual(TP, 250, 3)    

if __name__ == '__main__':
    unittest.main()
