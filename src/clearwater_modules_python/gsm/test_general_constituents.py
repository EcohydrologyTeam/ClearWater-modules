'''
=========================================================================================
Unit Tests - Water Temperature Simulation Module (TSM): Temperature Energy Budget Module
=========================================================================================

Developed by:
Dr. Todd E. Steissberg (ERDC-EL)

Date: April 11, 2021
'''

import unittest
from numba.typed import Dict
from numba import types
from general_constituents import GeneralConstituent

class Test_GSM (unittest.TestCase):
    def setUp(self):
        self.global_vars = Dict.empty(key_type=types.unicode_type, value_type=types.float64)
        self.global_vars = {
            "GC": 10,
            "TwaterC": 25,
            "depth": 1.0,
        }

        self.gsm_constant_changes = Dict.empty(key_type=types.unicode_type, value_type=types.float64)
        self.gsm_constant_changes = {
            "k_rc20": 0.5,
            "k_theta" : 1.047,

            "rgc_rc20" : 0.5,
            "rgc_theta" : 1.047,
            "order" : 1,
            "release" : True,
            "settling" : True,
            "settling_rate" : 0.0002,
        }

    def test_TwaterC (self):
        self.global_vars['TwaterC'] = 35
        Change = GeneralConstituent(self.global_vars, self.gsm_constant_changes).Calculations()
        self.assertAlmostEqual(Change, -8.964160948, places=4)
    
    def test_depth (self):
        self.global_vars['depth'] = 5
        Change = GeneralConstituent(self.global_vars, self.gsm_constant_changes).Calculations()
        self.assertAlmostEqual(Change, -6.165349003, places=4)   

    def test_k_rc20_1 (self):
        self.gsm_constant_changes['k_rc20'] = 0.7
        Change = GeneralConstituent(self.global_vars, self.gsm_constant_changes).Calculations()
        self.assertAlmostEqual(Change, -8.179993575, places=4) 

    def test_k_rc20_2 (self):
        self.gsm_constant_changes['k_rc20'] = 0.7
        self.gsm_constant_changes['order'] =2
        Change = GeneralConstituent(self.global_vars, self.gsm_constant_changes).Calculations()
        self.assertAlmostEqual(Change, -87.44362361, places=4)     

    def test_k_rc20_0 (self):
        self.gsm_constant_changes['k_rc20'] = 0.2
        self.gsm_constant_changes['order'] =0
        Change = GeneralConstituent(self.global_vars, self.gsm_constant_changes).Calculations()
        self.assertAlmostEqual(Change, 0.375445857, places=4)     

    def test_k_theta (self):
        self.gsm_constant_changes['k_theta'] = 1.1
        Change = GeneralConstituent(self.global_vars, self.gsm_constant_changes).Calculations()
        self.assertAlmostEqual(Change, -7.425473571, places=4)  

    def test_rgc_rc20 (self):
        self.gsm_constant_changes['rgc_rc20'] = 0.2
        Change = GeneralConstituent(self.global_vars, self.gsm_constant_changes).Calculations()
        self.assertAlmostEqual(Change, -6.041133717, places=4)        

    def test_order0 (self):
        self.gsm_constant_changes['order'] = 0
        Change = GeneralConstituent(self.global_vars, self.gsm_constant_changes).Calculations()
        self.assertAlmostEqual(Change, -0.002, places=4)  

    def test_order2 (self):
        self.gsm_constant_changes['order'] = 2
        Change = GeneralConstituent(self.global_vars, self.gsm_constant_changes).Calculations()
        self.assertAlmostEqual(Change, -62.28056646, places=4)   

    def test_release (self):
        self.gsm_constant_changes['release'] = False
        Change = GeneralConstituent(self.global_vars, self.gsm_constant_changes).Calculations()
        self.assertAlmostEqual(Change, -6.292764289, places=4)   

    def test_settling (self):
        self.gsm_constant_changes['settling'] = False
        Change = GeneralConstituent(self.global_vars, self.gsm_constant_changes).Calculations()
        self.assertAlmostEqual(Change, -5.66168786, places=4)   

    def test_settling_rate (self):
        self.gsm_constant_changes['settling_rate'] = 0.005
        Change = GeneralConstituent(self.global_vars, self.gsm_constant_changes).Calculations()
        self.assertAlmostEqual(Change, -5.71168786, places=4)   
'''            

def call_GeneralConstituentKinetics(params):
    return GeneralConstituentKinetics(
        params["GC"],
        params["TwaterC"],
        params["order"],
        params["k_rc20"],
        params["k_theta"],
        params["rgc_rc20"],
        params["rgc_theta"],
        release = params["release"],
        settling = params["settling"],
        depth = params["depth"],
        settling_rate = params["settling_rate"]
    )


class Test_GeneralConstituentKinetics(unittest.TestCase):
    def setUp(self):
        self.params = {
            "GC": 10.0,              # General constituent concentration
            "TwaterC": 25.0,         # Water temperature in degrees Celsius
            "order": 1,              # Compute 1st order kinetics
            "k_rc20": 0.5,           # Reaction rate at 20 degrees Celsius, decay
            "k_theta": 1.047,        # Arrhenius temperature correction factor, decay
            "rgc_rc20": 0.5,         # Reaction rate at 20 degrees Celsius, settling
            "rgc_theta": 1.047,      # Arrhenius temperature correction factor, settling
            "release": True,         # Turn suspension on
            "settling": True,        # Turn settling on
            "depth": 1.0,            # Bed depth
            "settling_rate": 0.0002  # Settling rate
        }


    def test_GC(self):
        test_params = self.params.copy()

        test_params["GC"] = 0.0
        self.assertAlmostEqual(call_GeneralConstituentKinetics(test_params), 0.6291, places=4)

        test_params["GC"] = 10.0
        self.assertAlmostEqual(call_GeneralConstituentKinetics(test_params), -5.6637, places=4)

        test_params["GC"] = 20.0
        self.assertAlmostEqual(call_GeneralConstituentKinetics(test_params), -11.9565, places=4)

        test_params["GC"] = 30.0
        self.assertAlmostEqual(call_GeneralConstituentKinetics(test_params), -18.2492, places=4)


    def test_TwaterC(self):
        test_params = self.params.copy()
        test_params["TwaterC"] = 0.0
        self.assertAlmostEqual(call_GeneralConstituentKinetics(test_params), -1.7979, places=4)
        test_params["TwaterC"] = 10.0
        self.assertAlmostEqual(call_GeneralConstituentKinetics(test_params), -2.8448, places=4)
        test_params["TwaterC"] = 20.0
        self.assertAlmostEqual(call_GeneralConstituentKinetics(test_params), -4.5020, places=4)
        test_params["TwaterC"] = 30.0
        self.assertAlmostEqual(call_GeneralConstituentKinetics(test_params), -7.1253, places=4)
'''
if __name__ == '__main__':
    unittest.main()

