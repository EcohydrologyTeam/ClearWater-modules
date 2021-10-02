import unittest
import os
import sys

base_path = os.path.dirname(
    os.path.dirname(
        os.path.dirname(
            os.path.realpath(__file__))))
print(f'base_path={base_path}')
tsm_path = os.path.join(base_path, 'src', 'TSM')
sys.path.append(tsm_path)

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

    def test_wind_speed_eq_3mps(self):
        print('Test if wind_speed = 3.0 m/s yields TeqC == 24.65 degC')
        self.wind_speed = 3.0
        tsm = TSM.TSM(self.TwaterC, self.surface_area, self.volume)
        tsm.energy_budget_method(
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance)
        tsm.print_pathways()
        TeqC = tsm.pathways['TeqC']['value']
        # self.assertAlmostEqual(TeqC, 24.65, 2)

    def test_wind_speed_eq_30mps(self):
        print('Test if wind_speed = 30.0 m/s yields TeqC == 24.65 degC')
        self.wind_speed = 30.0
        tsm = TSM.TSM(self.TwaterC, self.surface_area, self.volume)
        tsm.energy_budget_method(
            self.TairC, self.q_solar, self.pressure_mb, self.eair_mb, self.cloudiness, self.wind_speed, self.wind_a,
            self.wind_b, self.wind_c, self.wind_kh_kw, self.use_SedTemp, self.TsedC, self.num_iterations, self.
            tolerance)
        tsm.print_pathways()
        TeqC = tsm.pathways['TeqC']['value']
        # self.assertAlmostEqual(TeqC, 5.74, 2)


class Test_TSM_equilibrium_temperature_method(unittest.TestCase):
    def setUp(self):
        self.TwaterC = 20.0
        self.TeqC = 15.0
        self.KT = 0.1
        self.TsedC = 5.0
        self.surface_area = 1.0
        self.volume = 1.0
        self.use_SedTemp = True

    def test_KT_eq_1(self):
        print('Test if KT = 0.1 yields q_net == -401.52')
        self.KT = 1.0
        tsm = TSM.TSM(self.TwaterC, self.surface_area, self.volume)
        tsm.equilibrium_temperature_method(self.TeqC, self.KT, self.use_SedTemp, self.TsedC)
        tsm.print_pathways()
        q_net = tsm.pathways['q_net']['value']
        # self.assertAlmostEqual(q_net, -401.52, 2)

    def test_KT_eq_100(self):
        self.KT = 100.0
        print('Test if KT = 0.9 yields q_net == -401.52')
        tsm = TSM.TSM(self.TwaterC, self.surface_area, self.volume)
        tsm.equilibrium_temperature_method(self.TeqC, self.KT, self.use_SedTemp, self.TsedC)
        tsm.print_pathways()
        q_net = tsm.pathways['q_net']['value']
        # self.assertAlmostEqual(q_net, -401.52, 2)


if __name__ == '__main__':
    unittest.main()
