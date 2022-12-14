# %%
import unittest
import os
import sys
from numba.typed import Dict
from numba import types

# %%
## Code for getting to correct file path for local module import ##
base_path = os.path.dirname(
    os.path.dirname(
        os.path.dirname(
            os.path.realpath(__file__))))
tsm_path = os.path.join(base_path, 'src', 'TSM')
sys.path.append(tsm_path)
####################################################################

from TSM import Temperature
# %%
class Test_TSM_energy_budget_method(unittest.TestCase):
    # Locate module path and add to system path

    def setUp(self):
        
        self.global_module_choices = Dict.empty(key_type=types.unicode_type, value_type=types.float64)
        self.global_module_choices = {
            'use_SedTemp' : True
        }
        
        self.global_vars = Dict.empty(key_type=types.unicode_type, value_type=types.float64)
        self.global_vars = {
            'TwaterC' : 20,
            'surface_area' : 1,
            'volume' : 1

        }
        
        self.temperature_constant_changes = Dict.empty(key_type=types.unicode_type, value_type=types.float64)
        self.met_constant_changes = Dict.empty(key_type=types.unicode_type, value_type=types.float64)

    
    def test_TwaterC (self):

        self.global_vars['TwaterC'] = 20

        Change = Temperature(self.global_module_choices,self.global_vars,self.temperature_constant_changes,self.met_constant_changes).energy_budget_method()
        #Change = Change * 60

        self.assertAlmostEqual(Change, -5.39*10**-5, 3)


if __name__ == '__main__':
    unittest.main()
# %%
