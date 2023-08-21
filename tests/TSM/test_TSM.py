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

#Define initial array values
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

    
    def test_TwaterC_20 (self):

        self.global_vars['TwaterC'] = 20
        
        dTwaterCdt = Temperature(self.global_module_choices,self.global_vars,self.met_constant_changes,self.temperature_constant_changes).energy_budget_method()
        TwaterC = self.global_vars['TwaterC'] + dTwaterCdt * 60
        self.assertAlmostEqual(TwaterC, 19.997, 3)

    def test_TwaterC_40 (self):

        self.global_vars['TwaterC'] = 40
        
        dTwaterCdt = Temperature(self.global_module_choices,self.global_vars,self.met_constant_changes,self.temperature_constant_changes).energy_budget_method()
        TwaterC = self.global_vars['TwaterC'] + dTwaterCdt * 60
        self.assertAlmostEqual(TwaterC, 39.964, 3)

    def test_surface_area_2 (self):

        self.global_vars['surface_area'] = 2
        
        dTwaterCdt = Temperature(self.global_module_choices,self.global_vars,self.met_constant_changes,self.temperature_constant_changes).energy_budget_method()
        TwaterC = self.global_vars['TwaterC'] + dTwaterCdt * 60
        self.assertAlmostEqual(TwaterC, 19.994, 3)

    def test_surface_area_4 (self):

        self.global_vars['surface_area'] = 4
        
        dTwaterCdt = Temperature(self.global_module_choices,self.global_vars,self.met_constant_changes,self.temperature_constant_changes).energy_budget_method()
        TwaterC = self.global_vars['TwaterC'] + dTwaterCdt * 60
        self.assertAlmostEqual(TwaterC, 19.987, 3)

    def test_volume_2 (self):

        self.global_vars['volume'] = 2
        
        dTwaterCdt = Temperature(self.global_module_choices,self.global_vars,self.met_constant_changes,self.temperature_constant_changes).energy_budget_method()
        TwaterC = self.global_vars['TwaterC'] + dTwaterCdt * 60
        self.assertAlmostEqual(TwaterC, 19.998, 3)

    def test_volume_4 (self):

        self.global_vars['volume'] = 4
        
        dTwaterCdt = Temperature(self.global_module_choices,self.global_vars,self.met_constant_changes,self.temperature_constant_changes).energy_budget_method()
        TwaterC = self.global_vars['TwaterC'] + dTwaterCdt * 60
        self.assertAlmostEqual(TwaterC, 19.999, 3)

    def test_TairC_30 (self):

        self.met_constant_changes['TairC'] = 30
        
        dTwaterCdt = Temperature(self.global_module_choices,self.global_vars,self.met_constant_changes,self.temperature_constant_changes).energy_budget_method()
        TwaterC = self.global_vars['TwaterC'] + dTwaterCdt * 60
        self.assertAlmostEqual(TwaterC, 20, 3)

    def test_TairC_40 (self):

        self.met_constant_changes['TairC'] = 40
        
        dTwaterCdt = Temperature(self.global_module_choices,self.global_vars,self.met_constant_changes,self.temperature_constant_changes).energy_budget_method()
        TwaterC = self.global_vars['TwaterC'] + dTwaterCdt * 60
        self.assertAlmostEqual(TwaterC, 20.001, 3)







    def test_TsedC_10 (self):

        self.met_constant_changes['TsedC'] = 10
        
        dTwaterCdt = Temperature(self.global_module_choices,self.global_vars,self.met_constant_changes,self.temperature_constant_changes).energy_budget_method()
        TwaterC = self.global_vars['TwaterC'] + dTwaterCdt * 60
        self.assertAlmostEqual(TwaterC, 19.999, 3)

    def test_TsedC_15 (self):

        self.met_constant_changes['TsedC'] = 15
        
        dTwaterCdt = Temperature(self.global_module_choices,self.global_vars,self.met_constant_changes,self.temperature_constant_changes).energy_budget_method()
        TwaterC = self.global_vars['TwaterC'] + dTwaterCdt * 60
        self.assertAlmostEqual(TwaterC, 20.001, 3)

    def test_q_solar_450 (self):

        self.met_constant_changes['q_solar'] = 450
        
        dTwaterCdt = Temperature(self.global_module_choices,self.global_vars,self.met_constant_changes,self.temperature_constant_changes).energy_budget_method()
        TwaterC = self.global_vars['TwaterC'] + dTwaterCdt * 60
        self.assertAlmostEqual(TwaterC, 19.997, 3)

    def test_q_solar_350 (self):

        self.met_constant_changes['q_solar'] = 350
        
        dTwaterCdt = Temperature(self.global_module_choices,self.global_vars,self.met_constant_changes,self.temperature_constant_changes).energy_budget_method()
        TwaterC = self.global_vars['TwaterC'] + dTwaterCdt * 60
        self.assertAlmostEqual(TwaterC, 19.996, 3)

    def test_wind_kh_kw_0_5 (self):

        self.met_constant_changes['wind_kh_kw'] = .5
        
        dTwaterCdt = Temperature(self.global_module_choices,self.global_vars,self.met_constant_changes,self.temperature_constant_changes).energy_budget_method()
        TwaterC = self.global_vars['TwaterC'] + dTwaterCdt * 60
        self.assertAlmostEqual(TwaterC, 19.997, 3)

    def test_wind_kh_kw_1_5 (self):

        self.met_constant_changes['wind_kh_kw'] = 1.5
        
        dTwaterCdt = Temperature(self.global_module_choices,self.global_vars,self.met_constant_changes,self.temperature_constant_changes).energy_budget_method()
        TwaterC = self.global_vars['TwaterC'] + dTwaterCdt * 60
        self.assertAlmostEqual(TwaterC, 19.997, 3)

    def test_eair_mb_2 (self):

        self.met_constant_changes['eair_mb'] = 2
        
        dTwaterCdt = Temperature(self.global_module_choices,self.global_vars,self.met_constant_changes,self.temperature_constant_changes).energy_budget_method()
        TwaterC = self.global_vars['TwaterC'] + dTwaterCdt * 60
        self.assertAlmostEqual(TwaterC, 19.997, 3)

    def test_eair_mb_5 (self):

        self.met_constant_changes['eair_mb'] = 5
        
        dTwaterCdt = Temperature(self.global_module_choices,self.global_vars,self.met_constant_changes,self.temperature_constant_changes).energy_budget_method()
        TwaterC = self.global_vars['TwaterC'] + dTwaterCdt * 60
        self.assertAlmostEqual(TwaterC, 19.997, 3)




    def test_pressure_mb_970 (self):

        self.met_constant_changes['pressure_mb'] = 970
        
        dTwaterCdt = Temperature(self.global_module_choices,self.global_vars,self.met_constant_changes,self.temperature_constant_changes).energy_budget_method()
        TwaterC = self.global_vars['TwaterC'] + dTwaterCdt * 60
        self.assertAlmostEqual(TwaterC, 19.997, 3)

    def test_pressure_mb_1050 (self):

        self.met_constant_changes['pressure_mb'] = 1050
        
        dTwaterCdt = Temperature(self.global_module_choices,self.global_vars,self.met_constant_changes,self.temperature_constant_changes).energy_budget_method()
        TwaterC = self.global_vars['TwaterC'] + dTwaterCdt * 60
        self.assertAlmostEqual(TwaterC, 19.997, 3)

    def test_cloudiness_0 (self):

        self.met_constant_changes['cloudiness'] = 0
        
        dTwaterCdt = Temperature(self.global_module_choices,self.global_vars,self.met_constant_changes,self.temperature_constant_changes).energy_budget_method()
        TwaterC = self.global_vars['TwaterC'] + dTwaterCdt * 60
        self.assertAlmostEqual(TwaterC, 19.997, 3)

    def test_cloudiness_0_5 (self):

        self.met_constant_changes['cloudiness'] = .5
        
        dTwaterCdt = Temperature(self.global_module_choices,self.global_vars,self.met_constant_changes,self.temperature_constant_changes).energy_budget_method()
        TwaterC = self.global_vars['TwaterC'] + dTwaterCdt * 60
        self.assertAlmostEqual(TwaterC, 19.997, 3)

    def test_wind_speed_5 (self):

        self.met_constant_changes['wind_speed'] = 5
        
        dTwaterCdt = Temperature(self.global_module_choices,self.global_vars,self.met_constant_changes,self.temperature_constant_changes).energy_budget_method()
        TwaterC = self.global_vars['TwaterC'] + dTwaterCdt * 60
        self.assertAlmostEqual(TwaterC, 19.997, 3)

    def test_wind_speed_30 (self):

        self.met_constant_changes['wind_speed'] = 30
        
        dTwaterCdt = Temperature(self.global_module_choices,self.global_vars,self.met_constant_changes,self.temperature_constant_changes).energy_budget_method()
        TwaterC = self.global_vars['TwaterC'] + dTwaterCdt * 60
        self.assertAlmostEqual(TwaterC, 19.983, 3)

    def test_wind_a_1En7 (self):

        self.met_constant_changes['wind_a'] = 1*10**-7
        
        dTwaterCdt = Temperature(self.global_module_choices,self.global_vars,self.met_constant_changes,self.temperature_constant_changes).energy_budget_method()
        TwaterC = self.global_vars['TwaterC'] + dTwaterCdt * 60
        self.assertAlmostEqual(TwaterC, 19.997, 3)

    def test_wind_a_7En7 (self):

        self.met_constant_changes['wind_a'] = 7*10**-7
        
        dTwaterCdt = Temperature(self.global_module_choices,self.global_vars,self.met_constant_changes,self.temperature_constant_changes).energy_budget_method()
        TwaterC = self.global_vars['TwaterC'] + dTwaterCdt * 60
        self.assertAlmostEqual(TwaterC, 19.997, 3)

    def test_wind_b_1En6 (self):

        self.met_constant_changes['wind_b'] = 1*10**-6
        
        dTwaterCdt = Temperature(self.global_module_choices,self.global_vars,self.met_constant_changes,self.temperature_constant_changes).energy_budget_method()
        TwaterC = self.global_vars['TwaterC'] + dTwaterCdt * 60
        self.assertAlmostEqual(TwaterC, 19.997, 3)

    def test_wind_b_2En6 (self):

        self.met_constant_changes['wind_b'] = 2*10**-6
        
        dTwaterCdt = Temperature(self.global_module_choices,self.global_vars,self.met_constant_changes,self.temperature_constant_changes).energy_budget_method()
        TwaterC = self.global_vars['TwaterC'] + dTwaterCdt * 60
        self.assertAlmostEqual(TwaterC, 19.996, 3)

    def test_wind_c_0_5 (self):

        self.met_constant_changes['wind_c'] = .5
        
        dTwaterCdt = Temperature(self.global_module_choices,self.global_vars,self.met_constant_changes,self.temperature_constant_changes).energy_budget_method()
        TwaterC = self.global_vars['TwaterC'] + dTwaterCdt * 60
        self.assertAlmostEqual(TwaterC, 19.998, 3)

    def test_wind_c_3 (self):

        self.met_constant_changes['wind_c'] = 3
        
        dTwaterCdt = Temperature(self.global_module_choices,self.global_vars,self.met_constant_changes,self.temperature_constant_changes).energy_budget_method()
        TwaterC = self.global_vars['TwaterC'] + dTwaterCdt * 60
        self.assertAlmostEqual(TwaterC, 19.980, 3)

if __name__ == '__main__':
    unittest.main()
