import unittest
import os
import sys


## Code for getting to correct file path for local module import ##
base_path = os.path.dirname(
    os.path.dirname(
        os.path.dirname(
            os.path.realpath(r'C:\Users\b2edhijm\Documents\Projects\ERDC Steiss\DLL\Github Practice\GitTest1'))))
tsm_path = os.path.join(base_path, 'src', 'TSM')
tsm_path=r'C:\Users\b2edhijm\Documents\Projects\ERDC Steiss\DLL\Github Practice\GitTest1'
sys.path.append(tsm_path)
####################################################################

import TSM

class Test_TSM_energy_budget_method(unittest.TestCase):
    # Locate module path and add to system path

    def setUp(self):
        self.TwaterC = 20.0
        self.surface_area = 1.0
        self.volume = 1.0
        self.TairC = 20.0
        self.TsedC = 5.0
        self.q_solar = 400.0
        self.wind_kh_kw = 1.0
        self.eair_mb = 1.0
        self.pressure_mb = 1013.0
        self.cloudiness = 0.1
        self.wind_speed = 3.0
        self.wind_a = 0.3
        self.wind_b = 1.5
        self.wind_c = 1.0
        self.use_SedTemp = True
        self.num_iterations = 10
        self.tolerance = 0.01
        self.time_step = 60

    def test_TwaterC_eq_20degC(self):
        print('Test if TwaterC = 20 degC yields TwaterC == 19.997 degC')
        self.TwaterC = 20
        TSM.energy_budget_method(self.TwaterC, self.surface_area, self.volume,
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance, self.time_step)
        TSM.print_pathways()
        TwaterC = TSM.pathways['TwaterC']['value']
        self.assertAlmostEqual(TwaterC, 19.997, 3)

    def test_TwaterC_eq_40degC(self):
        print('Test if TwaterC = 40 degC yields TwaterC == 39.964 degC')
        self.TwaterC = 40
        TSM.energy_budget_method(self.TwaterC, self.surface_area, self.volume,
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance, self.time_step)
        TSM.print_pathways()
        TwaterC = TSM.pathways['TwaterC']['value']
        self.assertAlmostEqual(TwaterC, 39.964, 3)
        
    def test_surface_area_eq_2m2(self):
        print('Test if surface_area = 2 m^2 yields TwaterC == 19.994 degC')
        self.surface_area = 2
        TSM.energy_budget_method(self.TwaterC, self.surface_area, self.volume,
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance, self.time_step)
        TSM.print_pathways()
        TwaterC = TSM.pathways['TwaterC']['value']
        self.assertAlmostEqual(TwaterC, 19.994, 3)

    def test_surface_area_eq_4m2(self):
        print('Test if surface_area = 4 m^2 yields TwaterC == 19.987 degC')
        self.surface_area = 4
        TSM.energy_budget_method(self.TwaterC, self.surface_area, self.volume,
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance, self.time_step)
        TSM.print_pathways()
        TwaterC = TSM.pathways['TwaterC']['value']
        self.assertAlmostEqual(TwaterC, 19.987, 2)
        
    def test_volume_eq_2m3(self):
        print('Test if volume = 2 m^3 yields TwaterC == 19.998 degC')
        self.volume = 2
        TSM.energy_budget_method(self.TwaterC, self.surface_area, self.volume,
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance, self.time_step)
        TSM.print_pathways()
        TwaterC = TSM.pathways['TwaterC']['value']
        self.assertAlmostEqual(TwaterC, 19.998, 2)

    def test_volume_eq_4m3(self):
        print('Test if volume = 4 m^3 yields TwaterC == 19.999 degC')
        self.volume = 4
        TSM.energy_budget_method(self.TwaterC, self.surface_area, self.volume,
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance, self.time_step)
        TSM.print_pathways()
        TwaterC = TSM.pathways['TwaterC']['value']
        self.assertAlmostEqual(TwaterC, 19.999, 2)
        
    def test_TairC_eq_30degC(self):
        print('Test if TairC = 30 degC yields TwaterC == 20.000 degC')
        self.TairC = 30
        TSM.energy_budget_method(self.TwaterC, self.surface_area, self.volume,
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance, self.time_step)
        TSM.print_pathways()
        TwaterC = TSM.pathways['TwaterC']['value']
        self.assertAlmostEqual(TwaterC, 20.000, 3)

    def test_TairC_eq_40degC(self):
        print('Test if TairC = 40 degC yields TwaterC == 20.001 degC')
        self.TairC = 40
        TSM.energy_budget_method(self.TwaterC, self.surface_area, self.volume,
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance, self.time_step)
        TSM.print_pathways()
        TwaterC = TSM.pathways['TwaterC']['value']
        self.assertAlmostEqual(TwaterC, 20.001, 3)
        
    def test_TsedC_eq_10degC(self):
        print('Test if TsedC = 10 degC yields TwaterC == 19.999 degC')
        self.TsedC = 10
        TSM.energy_budget_method(self.TwaterC, self.surface_area, self.volume,
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance, self.time_step)
        TSM.print_pathways()
        TwaterC = TSM.pathways['TwaterC']['value']
        self.assertAlmostEqual(TwaterC, 19.999, 3)

    def test_TsedC_eq_15degC(self):
        print('Test if TsedC = 15 degC yields TwaterC == 20.001 degC')
        self.TsedC = 15
        TSM.energy_budget_method(self.TwaterC, self.surface_area, self.volume,
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance, self.time_step)
        TSM.print_pathways()
        TwaterC = TSM.pathways['TwaterC']['value']
        self.assertAlmostEqual(TwaterC, 20.001, 3)
        
    def test_q_solar_eq_450w_m2(self):
        print('Test if q_solar = 450 w/m^2 yields TwaterC == 19.997 degC')
        self.q_solar = 450
        TSM.energy_budget_method(self.TwaterC, self.surface_area, self.volume,
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance, self.time_step)
        TSM.print_pathways()
        TwaterC = TSM.pathways['TwaterC']['value']
        self.assertAlmostEqual(TwaterC, 19.997, 3)

    def test_q_solar_eq_350w_m2(self):
        print('Test if q_solar = 350 w/m^2 yields TwaterC == 19.996 degC')
        self.q_solar = 350
        TSM.energy_budget_method(self.TwaterC, self.surface_area, self.volume,
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance, self.time_step)
        TSM.print_pathways()
        TwaterC = TSM.pathways['TwaterC']['value']
        self.assertAlmostEqual(TwaterC, 19.996, 3)
        
    def test_wind_kh_kw_eq_p5(self):
        print('Test if wind_kh_kw = .5 yields TwaterC == 19.997 degC')
        self.wind_kh_kw = .5
        TSM.energy_budget_method(self.TwaterC, self.surface_area, self.volume,
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance, self.time_step)
        TSM.print_pathways()
        TwaterC = TSM.pathways['TwaterC']['value']
        self.assertAlmostEqual(TwaterC, 19.997, 3)
        
    def test_wind_kh_kw_eq_1(self):
        print('Test if wind_kh_kw = 1.5 yields TwaterC == 19.997 degC')
        self.wind_kh_kw = 1
        TSM.energy_budget_method(self.TwaterC, self.surface_area, self.volume,
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance, self.time_step)
        TSM.print_pathways()
        TwaterC = TSM.pathways['TwaterC']['value']
        self.assertAlmostEqual(TwaterC, 19.997, 3)
        
    def test_eair_mb_eq_2(self):
        print('Test if eair_mb = 2 mb yields TwaterC == 19.997 degC')
        self.eair_mb = 2
        TSM.energy_budget_method(self.TwaterC, self.surface_area, self.volume,
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance, self.time_step)
        TSM.print_pathways()
        TwaterC = TSM.pathways['TwaterC']['value']
        self.assertAlmostEqual(TwaterC, 19.997, 3)

    def test_eair_mb_eq_5(self):
        print('Test if eair_mb = 5 mb yields TwaterC == 19.997 degC')
        self.eair_mb = 5
        TSM.energy_budget_method(self.TwaterC, self.surface_area, self.volume,
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance, self.time_step)
        TSM.print_pathways()
        TwaterC = TSM.pathways['TwaterC']['value']
        self.assertAlmostEqual(TwaterC, 19.997, 3)
        
    def test_pressure_mb_eq_970(self):
        print('Test if pressure_mb = 970 mb yields TwaterC == 19.997 degC')
        self.pressure_mb = 970
        TSM.energy_budget_method(self.TwaterC, self.surface_area, self.volume,
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance, self.time_step)
        TSM.print_pathways()
        TwaterC = TSM.pathways['TwaterC']['value']
        self.assertAlmostEqual(TwaterC, 19.997, 3)
        
    def test_pressure_mb_eq_1050(self):
        print('Test if pressure_mb = 1050 mb yields TwaterC == 19.997 degC')
        self.pressure_mb = 1050
        TSM.energy_budget_method(self.TwaterC, self.surface_area, self.volume,
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance, self.time_step)
        TSM.print_pathways()
        TwaterC = TSM.pathways['TwaterC']['value']
        self.assertAlmostEqual(TwaterC, 19.997, 3)

    def test_cloudiness_eq_0(self):
        print('Test if cloudiness = 0 yields TwaterC == 19.997 degC')
        self.cloudiness = 0
        TSM.energy_budget_method(self.TwaterC, self.surface_area, self.volume,
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance, self.time_step)
        TSM.print_pathways()
        TwaterC = TSM.pathways['TwaterC']['value']
        self.assertAlmostEqual(TwaterC, 19.997, 3)

    def test_cloudiness_eq_p5(self):
        print('Test if cloudiness = 0 yields TwaterC == 19.997 degC')
        self.cloudiness = .5
        TSM.energy_budget_method(self.TwaterC, self.surface_area, self.volume,
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance, self.time_step)
        TSM.print_pathways()
        TwaterC = TSM.pathways['TwaterC']['value']
        self.assertAlmostEqual(TwaterC, 19.997, 3)
        
    def test_wind_speed_eq_5m_s(self):
        print('Test if wind_speed = 5 m/s yields TwaterC == 19.996 degC')
        self.wind_speed = 5
        TSM.energy_budget_method(self.TwaterC, self.surface_area, self.volume,
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance, self.time_step)
        TSM.print_pathways()
        TwaterC = TSM.pathways['TwaterC']['value']
        self.assertAlmostEqual(TwaterC, 19.996, 3)

    def test_wind_speed_eq_30m_s(self):
        print('Test if wind_speed = 30 m/s yields TwaterC == 19.983 degC')
        self.wind_speed = 30
        TSM.energy_budget_method(self.TwaterC, self.surface_area, self.volume,
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance, self.time_step)
        TSM.print_pathways()
        TwaterC = TSM.pathways['TwaterC']['value']
        self.assertAlmostEqual(TwaterC, 19.983, 3)
        
    def test_wind_a_eq_1En7(self):
        print('Test if wind_a = .1 yields TwaterC == 19.997 degC')
        self.wind_a = .1
        TSM.energy_budget_method(self.TwaterC, self.surface_area, self.volume,
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance, self.time_step)
        TSM.print_pathways()
        TwaterC = TSM.pathways['TwaterC']['value']
        self.assertAlmostEqual(TwaterC, 19.997, 3)

    def test_wind_a_eq_7En7(self):
        print('Test if wind_a = .7 yields TwaterC == 19.997 degC')
        self.wind_a = .7
        TSM.energy_budget_method(self.TwaterC, self.surface_area, self.volume,
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance, self.time_step)
        TSM.print_pathways()
        TwaterC = TSM.pathways['TwaterC']['value']
        self.assertAlmostEqual(TwaterC, 19.997, 3)
        
    def test_wind_b_eq_1En6(self):
        print('Test if wind_b = 1 yields TwaterC == 19.997 degC')
        self.wind_b = 1
        TSM.energy_budget_method(self.TwaterC, self.surface_area, self.volume,
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance, self.time_step)
        TSM.print_pathways()
        TwaterC = TSM.pathways['TwaterC']['value']
        self.assertAlmostEqual(TwaterC, 19.997, 3)
        
    def test_wind_b_eq_2En6(self):
        print('Test if wind_b = 2 yields TwaterC == 19.996 degC')
        self.wind_b = 2
        TSM.energy_budget_method(self.TwaterC, self.surface_area, self.volume,
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance, self.time_step)
        TSM.print_pathways()
        TwaterC = TSM.pathways['TwaterC']['value']
        self.assertAlmostEqual(TwaterC, 19.996, 3)

    def test_wind_c_eq_p5(self):
        print('Test if wind_c = .5 yields TwaterC == 19.998 degC')
        self.wind_c = .5
        TSM.energy_budget_method(self.TwaterC, self.surface_area, self.volume,
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance, self.time_step)
        TSM.print_pathways()
        TwaterC = TSM.pathways['TwaterC']['value']
        self.assertAlmostEqual(TwaterC, 19.998, 3)
        
    def test_wind_c_eq_3(self):
        print('Test if wind_c = 3 yields TwaterC == 19.980 degC')
        self.wind_c = 3
        TSM.energy_budget_method(self.TwaterC, self.surface_area, self.volume,
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance, self.time_step)
        TSM.print_pathways()
        TwaterC = TSM.pathways['TwaterC']['value']
        self.assertAlmostEqual(TwaterC, 19.980, 3)
        
    def test_time_step_eq_1000(self):
        print('Test if time_step = 1000 yields TwaterC == 19.946 degC')
        self.time_step = 1000
        TSM.energy_budget_method(self.TwaterC, self.surface_area, self.volume,
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance, self.time_step)
        TSM.print_pathways()
        TwaterC = TSM.pathways['TwaterC']['value']
        self.assertAlmostEqual(TwaterC, 19.946, 3)
        
    def test_time_step_eq_10000(self):
        print('Test if time_step = 10000 yields TwaterC == 19.461 degC')
        self.time_step = 10000
        TSM.energy_budget_method(self.TwaterC, self.surface_area, self.volume,
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance, self.time_step)
        TSM.print_pathways()
        TwaterC = TSM.pathways['TwaterC']['value']
        self.assertAlmostEqual(TwaterC, 19.461, 3)

if __name__ == '__main__':
    unittest.main()