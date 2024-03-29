from collections import OrderedDict

import unittest
from numba.typed import Dict
from numba import types

from _nitrogen import Nitrogen

class Test_nitrogen(unittest.TestCase):

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
            'OrgN' : 100,
            'vson': 0.01,

            'lambda': 1,
            'fdp': 0.5,
            'PAR': 100
        }
        
        self.algae_pathways = Dict.empty(key_type=types.unicode_type, value_type=types.float64)
        self.algae_pathways={
            'ApGrowth': 10,
            'ApDeath' : 20,
            'ApRespiration': 30,
            'rna' : 0.5
        } 

        self.Balgae_pathways= Dict.empty(key_type=types.unicode_type, value_type=types.float64)
        self.Balgae_pathways = {
            'AbGrowth' : 30,
            'AbDeath' : 20,
            'AbRespiration' : 10,
            'rnb' : 0.25,
            'Fw' : 0.9,
            'Fb' : 0.9,
        }

        self.sedFlux_pathways = Dict.empty(key_type=types.unicode_type, value_type=types.float64)
        self.sedFlux_pathways = {
            'JNH4': 2,
            'JNO3': 2
        }

        self.nitrogen_constant_changes = Dict.empty(key_type=types.unicode_type, value_type=types.float64)
        self.nitrogen_constant_changes = {
                
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

#Original test
    def test_nitrogen1 (self):

        DIN, TON, TKN, TN, dOrgNdt, dNH4dt, dNO3dt= Nitrogen(self.global_module_choices, self.global_vars, self.algae_pathways, 
            self.Balgae_pathways, self.sedFlux_pathways, self.nitrogen_constant_changes).Calculations()
        
        self.assertAlmostEqual(DIN, 200, 3)
        self.assertAlmostEqual(TON, 150, 3)
        self.assertAlmostEqual(TKN, 250, 3)
        self.assertAlmostEqual(TN, 350, 3)

        self.assertAlmostEqual(dNH4dt, 11.625, 3)
        self.assertAlmostEqual(dNO3dt, 4.1248, 3)
        self.assertAlmostEqual(dOrgNdt, 3.05, 3)

#Change DOX and OrgN
    def test_change_DOX_OrgN (self):

        self.global_vars['DOX'] = 200
        self.global_vars['OrgN'] = 200

        DIN, TON, TKN, TN, dOrgNdt, dNH4dt, dNO3dt= Nitrogen(self.global_module_choices, self.global_vars, self.algae_pathways, 
            self.Balgae_pathways, self.sedFlux_pathways, self.nitrogen_constant_changes).Calculations()
        
        self.assertAlmostEqual(DIN, 200, 3)
        self.assertAlmostEqual(TON, 250, 3)
        self.assertAlmostEqual(TKN, 350, 3)
        self.assertAlmostEqual(TN, 450, 3)

        self.assertAlmostEqual(dNH4dt, 21.625, 3)
        self.assertAlmostEqual(dNO3dt, 4.1249, 3)
        self.assertAlmostEqual(dOrgNdt, -7.95, 3)

#Change knit, kon, kdnit
    def test_change_knit_kon_kdnit (self):

        self.nitrogen_constant_changes['knit'] = 0.2
        self.nitrogen_constant_changes['kon'] = 0.2
        self.nitrogen_constant_changes['kdnit'] = 0.004

        DIN, TON, TKN, TN, dOrgNdt, dNH4dt, dNO3dt= Nitrogen(self.global_module_choices, self.global_vars, self.algae_pathways, 
            self.Balgae_pathways, self.sedFlux_pathways, self.nitrogen_constant_changes).Calculations()
        
        self.assertAlmostEqual(DIN, 200, 3)
        self.assertAlmostEqual(TON, 150, 3)
        self.assertAlmostEqual(TKN, 250, 3)
        self.assertAlmostEqual(TN, 350, 3)

        self.assertAlmostEqual(dNH4dt, 11.625, 3)
        self.assertAlmostEqual(dNO3dt, 14.1246, 3)
        self.assertAlmostEqual(dOrgNdt, -6.95, 3)    

#Change ApGrowth and ApDeath
    def test_change_ApGrowth_ApDeath (self):

        self.algae_pathways['ApGrowth'] = 30
        self.algae_pathways['ApDeath'] = 30

        DIN, TON, TKN, TN, dOrgNdt, dNH4dt, dNO3dt= Nitrogen(self.global_module_choices, self.global_vars, self.algae_pathways, 
            self.Balgae_pathways, self.sedFlux_pathways, self.nitrogen_constant_changes).Calculations()
        
        self.assertAlmostEqual(DIN, 200, 3)
        self.assertAlmostEqual(TON, 150, 3)
        self.assertAlmostEqual(TKN, 250, 3)
        self.assertAlmostEqual(TN, 350, 3)

        self.assertAlmostEqual(dNH4dt, 6.625, 3)
        self.assertAlmostEqual(dNO3dt, -0.8752, 3)
        self.assertAlmostEqual(dOrgNdt, 8.05, 3) 

#Change depth
    def test_change_depth (self):

        self.global_vars['depth'] = 5

        DIN, TON, TKN, TN, dOrgNdt, dNH4dt, dNO3dt= Nitrogen(self.global_module_choices, self.global_vars, self.algae_pathways, 
            self.Balgae_pathways, self.sedFlux_pathways, self.nitrogen_constant_changes).Calculations()
        
        self.assertAlmostEqual(DIN, 200, 3)
        self.assertAlmostEqual(TON, 150, 3)
        self.assertAlmostEqual(TKN, 250, 3)
        self.assertAlmostEqual(TN, 350, 3)

        self.assertAlmostEqual(dNH4dt, 14.325, 3)
        self.assertAlmostEqual(dNO3dt, 6.8248, 3)
        self.assertAlmostEqual(dOrgNdt, 0.61, 3) 

#Change rna
    def test_change_rna (self):

        self.algae_pathways['rna'] = 0.25

        DIN, TON, TKN, TN, dOrgNdt, dNH4dt, dNO3dt= Nitrogen(self.global_module_choices, self.global_vars, self.algae_pathways, 
            self.Balgae_pathways, self.sedFlux_pathways, self.nitrogen_constant_changes).Calculations()
        
        self.assertAlmostEqual(DIN, 200, 3)
        self.assertAlmostEqual(TON, 125, 3)
        self.assertAlmostEqual(TKN, 225, 3)
        self.assertAlmostEqual(TN, 325, 3)

        self.assertAlmostEqual(dNH4dt, 5.375, 3)
        self.assertAlmostEqual(dNO3dt, 5.3748, 3)
        self.assertAlmostEqual(dOrgNdt, -1.95, 3)

#Change rnb
    def test_change_rnb (self):

        self.Balgae_pathways['rnb'] = 0.5

        DIN, TON, TKN, TN, dOrgNdt, dNH4dt, dNO3dt= Nitrogen(self.global_module_choices, self.global_vars, self.algae_pathways, 
            self.Balgae_pathways, self.sedFlux_pathways, self.nitrogen_constant_changes).Calculations()
        
        self.assertAlmostEqual(DIN, 200, 3)
        self.assertAlmostEqual(TON, 150, 3)
        self.assertAlmostEqual(TKN, 250, 3)
        self.assertAlmostEqual(TN, 350, 3)

        self.assertAlmostEqual(dNH4dt, 10.75, 3)
        self.assertAlmostEqual(dNO3dt, 0.7498, 3)
        self.assertAlmostEqual(dOrgNdt, 7.1, 3)  

#Change KNR, vno3, vson
    def test_change_KNR_vno3_vson (self):

        self.nitrogen_constant_changes['KRN'] = 0.3
        self.nitrogen_constant_changes['vno3'] = 1
        self.global_vars['vson'] = 0.02

        DIN, TON, TKN, TN, dOrgNdt, dNH4dt, dNO3dt= Nitrogen(self.global_module_choices, self.global_vars, self.algae_pathways, 
            self.Balgae_pathways, self.sedFlux_pathways, self.nitrogen_constant_changes).Calculations()
        
        self.assertAlmostEqual(DIN, 200, 3)
        self.assertAlmostEqual(TON, 150, 3)
        self.assertAlmostEqual(TKN, 250, 3)
        self.assertAlmostEqual(TN, 350, 3)

        self.assertAlmostEqual(dNH4dt, 11.625, 3)
        self.assertAlmostEqual(dNO3dt, -95.8752, 3)
        self.assertAlmostEqual(dOrgNdt, 2.05, 3)  

#Change Fb
    def test_change_Fb (self):

        self.Balgae_pathways['Fb'] = 0.75

        DIN, TON, TKN, TN, dOrgNdt, dNH4dt, dNO3dt= Nitrogen(self.global_module_choices, self.global_vars, self.algae_pathways, 
            self.Balgae_pathways, self.sedFlux_pathways, self.nitrogen_constant_changes).Calculations()
        
        self.assertAlmostEqual(DIN, 200, 3)
        self.assertAlmostEqual(TON, 150, 3)
        self.assertAlmostEqual(TKN, 250, 3)
        self.assertAlmostEqual(TN, 350, 3)

        self.assertAlmostEqual(dNH4dt, 12.1875, 3)
        self.assertAlmostEqual(dNO3dt, 4.6873, 3)
        self.assertAlmostEqual(dOrgNdt, 2.375, 3) 

#Change Fw, rnh4, KsOxdn
    def test_change_Fw_rnh4_KsOxdn (self):

        self.Balgae_pathways['Fw'] = 0.5
        self.nitrogen_constant_changes['rnh4'] = 1
        self.nitrogen_constant_changes['KsOxdn'] = 0.2

        DIN, TON, TKN, TN, dOrgNdt, dNH4dt, dNO3dt= Nitrogen(self.global_module_choices, self.global_vars, self.algae_pathways, 
            self.Balgae_pathways, self.sedFlux_pathways, self.nitrogen_constant_changes).Calculations()
        
        self.assertAlmostEqual(DIN, 200, 3)
        self.assertAlmostEqual(TON, 150, 3)
        self.assertAlmostEqual(TKN, 250, 3)
        self.assertAlmostEqual(TN, 350, 3)

        self.assertAlmostEqual(dNH4dt, 12.625, 3)
        self.assertAlmostEqual(dNO3dt, 4.1246, 3)
        self.assertAlmostEqual(dOrgNdt, 1.25, 3) 

#Change AbGrowth, AbDeath
    def test_change_AbGrowth_AbDeath (self):

        self.Balgae_pathways['AbGrowth'] = 10
        self.Balgae_pathways['AbDeath'] = 30

        DIN, TON, TKN, TN, dOrgNdt, dNH4dt, dNO3dt= Nitrogen(self.global_module_choices, self.global_vars, self.algae_pathways, 
            self.Balgae_pathways, self.sedFlux_pathways, self.nitrogen_constant_changes).Calculations()
        
        self.assertAlmostEqual(DIN, 200, 3)
        self.assertAlmostEqual(TON, 150, 3)
        self.assertAlmostEqual(TKN, 250, 3)
        self.assertAlmostEqual(TN, 350, 3)

        self.assertAlmostEqual(dNH4dt, 13.875, 3)
        self.assertAlmostEqual(dNO3dt, 6.3748, 3)
        self.assertAlmostEqual(dOrgNdt, 5.075, 3) 

#Change ApRespiration
    def test_change_ApRespiration (self):

        self.algae_pathways['ApRespiration'] = 15

        DIN, TON, TKN, TN, dOrgNdt, dNH4dt, dNO3dt= Nitrogen(self.global_module_choices, self.global_vars, self.algae_pathways, 
            self.Balgae_pathways, self.sedFlux_pathways, self.nitrogen_constant_changes).Calculations()
        
        self.assertAlmostEqual(DIN, 200, 3)
        self.assertAlmostEqual(TON, 150, 3)
        self.assertAlmostEqual(TKN, 250, 3)
        self.assertAlmostEqual(TN, 350, 3)

        self.assertAlmostEqual(dNH4dt, 4.125, 3)
        self.assertAlmostEqual(dNO3dt, 4.1248, 3)
        self.assertAlmostEqual(dOrgNdt, 3.05, 3) 

#Change Ap and AbRespiration
    def test_change_Ap_AbRespiration (self):

        self.Balgae_pathways['AbRespiration'] = 20
        self.global_vars['Ap'] = 200

        DIN, TON, TKN, TN, dOrgNdt, dNH4dt, dNO3dt= Nitrogen(self.global_module_choices, self.global_vars, self.algae_pathways, 
            self.Balgae_pathways, self.sedFlux_pathways, self.nitrogen_constant_changes).Calculations()
        
        self.assertAlmostEqual(DIN, 200, 3)
        self.assertAlmostEqual(TON, 200, 3)
        self.assertAlmostEqual(TKN, 300, 3)
        self.assertAlmostEqual(TN, 400, 3)

        self.assertAlmostEqual(dNH4dt, 14.125, 3)
        self.assertAlmostEqual(dNO3dt, 4.1248, 3)
        self.assertAlmostEqual(dOrgNdt, 3.05, 3) 

#Change NH4
    def test_change_NH4 (self):

        self.global_vars['NH4'] = 200

        DIN, TON, TKN, TN, dOrgNdt, dNH4dt, dNO3dt= Nitrogen(self.global_module_choices, self.global_vars, self.algae_pathways, 
            self.Balgae_pathways, self.sedFlux_pathways, self.nitrogen_constant_changes).Calculations()
        
        self.assertAlmostEqual(DIN, 300, 3)
        self.assertAlmostEqual(TON, 150, 3)
        self.assertAlmostEqual(TKN, 350, 3)
        self.assertAlmostEqual(TN, 450, 3)

        self.assertAlmostEqual(dNH4dt, -0.3333, 3)
        self.assertAlmostEqual(dNO3dt, 16.08313, 3)
        self.assertAlmostEqual(dOrgNdt, 3.05, 3) 

#Change NO3
    def test_change_NO3 (self):

        self.global_vars['NO3'] = 200

        DIN, TON, TKN, TN, dOrgNdt, dNH4dt, dNO3dt= Nitrogen(self.global_module_choices, self.global_vars, self.algae_pathways, 
            self.Balgae_pathways, self.sedFlux_pathways, self.nitrogen_constant_changes).Calculations()
        
        self.assertAlmostEqual(DIN, 300, 3)
        self.assertAlmostEqual(TON, 150, 3)
        self.assertAlmostEqual(TKN, 250, 3)
        self.assertAlmostEqual(TN, 450, 3)

        self.assertAlmostEqual(dNH4dt, 13.5833, 3)
        self.assertAlmostEqual(dNO3dt, 2.16627, 3)
        self.assertAlmostEqual(dOrgNdt, 3.05, 3) 

#Change TwaterC
    def test_change_TwaterC (self):

        self.global_vars['TwaterC'] = 10

        DIN, TON, TKN, TN, dOrgNdt, dNH4dt, dNO3dt= Nitrogen(self.global_module_choices, self.global_vars, self.algae_pathways, 
            self.Balgae_pathways, self.sedFlux_pathways, self.nitrogen_constant_changes).Calculations()
        
        self.assertAlmostEqual(DIN, 200, 3)
        self.assertAlmostEqual(TON, 150, 3)
        self.assertAlmostEqual(TKN, 250, 3)
        self.assertAlmostEqual(TN, 350, 3)

        self.assertAlmostEqual(dNH4dt, 12.01719, 3)
        self.assertAlmostEqual(dNO3dt, -1.369915, 3)
        self.assertAlmostEqual(dOrgNdt, 8.15269, 3) 

#Change PN
    def test_change_PN (self):

        self.nitrogen_constant_changes['PN'] = 0.75

        DIN, TON, TKN, TN, dOrgNdt, dNH4dt, dNO3dt= Nitrogen(self.global_module_choices, self.global_vars, self.algae_pathways, 
            self.Balgae_pathways, self.sedFlux_pathways, self.nitrogen_constant_changes).Calculations()
        
        self.assertAlmostEqual(DIN, 200, 3)
        self.assertAlmostEqual(TON, 150, 3)
        self.assertAlmostEqual(TKN, 250, 3)
        self.assertAlmostEqual(TN, 350, 3)

        self.assertAlmostEqual(dNH4dt, 10.375, 3)
        self.assertAlmostEqual(dNO3dt, 5.3748, 3)
        self.assertAlmostEqual(dOrgNdt, 3.05, 3) 

#Change PNb
    def test_change_PNb (self):

        self.nitrogen_constant_changes['PNb'] = 0.25

        DIN, TON, TKN, TN, dOrgNdt, dNH4dt, dNO3dt= Nitrogen(self.global_module_choices, self.global_vars, self.algae_pathways, 
            self.Balgae_pathways, self.sedFlux_pathways, self.nitrogen_constant_changes).Calculations()
        
        self.assertAlmostEqual(DIN, 200, 3)
        self.assertAlmostEqual(TON, 150, 3)
        self.assertAlmostEqual(TKN, 250, 3)
        self.assertAlmostEqual(TN, 350, 3)

        self.assertAlmostEqual(dNH4dt, 13.3125, 3)
        self.assertAlmostEqual(dNO3dt, 2.43730, 3)
        self.assertAlmostEqual(dOrgNdt, 3.05, 3) 

#Change use_SedFlux
    def test_change_use_SedFlux(self):

        self.global_module_choices['use_SedFlux'] = True

        DIN, TON, TKN, TN, dOrgNdt, dNH4dt, dNO3dt= Nitrogen(self.global_module_choices, self.global_vars, self.algae_pathways, 
            self.Balgae_pathways, self.sedFlux_pathways, self.nitrogen_constant_changes).Calculations()
        
        self.assertAlmostEqual(DIN, 200, 3)
        self.assertAlmostEqual(TON, 150, 3)
        self.assertAlmostEqual(TKN, 250, 3)
        self.assertAlmostEqual(TN, 350, 3)

        self.assertAlmostEqual(dNH4dt, 13.625, 3)
        self.assertAlmostEqual(dNO3dt, 2.12480, 3)
        self.assertAlmostEqual(dOrgNdt, 3.05, 3) 

#Change JNH4 and JNO3
    def test_change_JNH4_JNO3(self):

        self.global_module_choices['use_SedFlux'] = True
        self.sedFlux_pathways['JNH4'] = 4
        self.sedFlux_pathways['JNO3'] = 4

        DIN, TON, TKN, TN, dOrgNdt, dNH4dt, dNO3dt= Nitrogen(self.global_module_choices, self.global_vars, self.algae_pathways, 
            self.Balgae_pathways, self.sedFlux_pathways, self.nitrogen_constant_changes).Calculations()
        
        self.assertAlmostEqual(DIN, 200, 3)
        self.assertAlmostEqual(TON, 150, 3)
        self.assertAlmostEqual(TKN, 250, 3)
        self.assertAlmostEqual(TN, 350, 3)

        self.assertAlmostEqual(dNH4dt, 15.625, 3)
        self.assertAlmostEqual(dNO3dt, 0.12480, 3)
        self.assertAlmostEqual(dOrgNdt, 3.05, 3) 

#Change use_NO4
    def test_use_NH4(self):

        self.global_module_choices['use_NH4'] = False

        DIN, TON, TKN, TN, dOrgNdt, dNH4dt, dNO3dt= Nitrogen(self.global_module_choices, self.global_vars, self.algae_pathways, 
            self.Balgae_pathways, self.sedFlux_pathways, self.nitrogen_constant_changes).Calculations()

        self.assertAlmostEqual(dNH4dt, 0, 3)
        self.assertAlmostEqual(dNO3dt, -11.75, 3)
        self.assertAlmostEqual(dOrgNdt, 3.05, 3)

#Change use_NO3
    def test_use_NO3(self):

        self.global_module_choices['use_NO3'] = False

        DIN, TON, TKN, TN, dOrgNdt, dNH4dt, dNO3dt= Nitrogen(self.global_module_choices, self.global_vars, self.algae_pathways, 
            self.Balgae_pathways, self.sedFlux_pathways, self.nitrogen_constant_changes).Calculations()

        self.assertAlmostEqual(dNH4dt, 5.75, 3)
        self.assertAlmostEqual(dNO3dt, 0, 3)
        self.assertAlmostEqual(dOrgNdt, 3.05, 3) 

#Change use_OrgN
    def test_use_OrgN(self):

        self.global_module_choices['use_OrgN'] = False

        DIN, TON, TKN, TN, dOrgNdt, dNH4dt, dNO3dt= Nitrogen(self.global_module_choices, self.global_vars, self.algae_pathways, 
            self.Balgae_pathways, self.sedFlux_pathways, self.nitrogen_constant_changes).Calculations()

        self.assertAlmostEqual(dNH4dt, 1.625, 3)
        self.assertAlmostEqual(dNO3dt, 4.1248, 3)
        self.assertAlmostEqual(dOrgNdt, 0, 3) 

#Change use_DOX
    def test_use_DOX(self):

        self.global_module_choices['use_DOX'] = False

        DIN, TON, TKN, TN, dOrgNdt, dNH4dt, dNO3dt= Nitrogen(self.global_module_choices, self.global_vars, self.algae_pathways, 
            self.Balgae_pathways, self.sedFlux_pathways, self.nitrogen_constant_changes).Calculations()

        self.assertAlmostEqual(dNH4dt, 11.625, 3)
        self.assertAlmostEqual(dNO3dt, 4.125, 3)
        self.assertAlmostEqual(dOrgNdt, 3.05, 3) 

#Change use_Algae
    def test_use_Algae(self):

        self.global_module_choices['use_Algae'] = False

        DIN, TON, TKN, TN, dOrgNdt, dNH4dt, dNO3dt= Nitrogen(self.global_module_choices, self.global_vars, self.algae_pathways, 
            self.Balgae_pathways, self.sedFlux_pathways, self.nitrogen_constant_changes).Calculations()

        self.assertAlmostEqual(dNH4dt, -0.875, 3)
        self.assertAlmostEqual(dNO3dt, 6.6248, 3)
        self.assertAlmostEqual(dOrgNdt, -6.95, 3) 

#Change use_BAlgae
    def test_use_BAlgae(self):

        self.global_module_choices['use_BAlgae'] = False

        DIN, TON, TKN, TN, dOrgNdt, dNH4dt, dNO3dt= Nitrogen(self.global_module_choices, self.global_vars, self.algae_pathways, 
            self.Balgae_pathways, self.sedFlux_pathways, self.nitrogen_constant_changes).Calculations()

        self.assertAlmostEqual(dNH4dt, 12.5, 3)
        self.assertAlmostEqual(dNO3dt, 7.4998, 3)
        self.assertAlmostEqual(dOrgNdt, -1, 3) 
                
if __name__ == '__main__':
    unittest.main()
