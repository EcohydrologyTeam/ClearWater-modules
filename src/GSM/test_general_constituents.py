'''
=========================================================================================
Unit Tests - Water Temperature Simulation Module (TSM): Temperature Energy Budget Module
=========================================================================================

Developed by:
Dr. Todd E. Steissberg (ERDC-EL)

Date: April 11, 2021
'''

# %%
import unittest
from general_constituents import GeneralConstituentKinetics

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