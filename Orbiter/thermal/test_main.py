from unittest import TestCase

from Orbiter.thermal.main import *

# Constants
boltzman = 5.67 * 10 ** (-8)

# Planet list [Sun distance, radius, albedo factor, radiating temperature, *closest approach during gravity assist]
planets_list_test = {'Planet_Test': [2.3762 * 10 ** 8, 69720 / 2, 0.3, 60, 150],
                     'Uranus': [2872500000, 51118 / 2, 0.51, 58.2],
                     'Venus': [108200000, 12104 / 2, 0.65, 227, 200]}


# Unit Tests Functions
class Test(TestCase):
    def test_solar_intensity(self, r=156, planet='Planet_Test', planet_list=planets_list_test):
        J_solar_true = 543.4512698
        J_solar_numerical = solar_intensity(r, planet, planet_list)
        self.assertAlmostEqual(J_solar_numerical, J_solar_true)

    def test_albedo(self, J=10 ** 2, r=156, planet="Planet_Test", planet_list=planets_list_test):
        J_albedo_true = 29.73328906
        J_albedo_numerical = albedo(J, r, planet, planet_list)
        self.assertAlmostEqual(J_albedo_numerical, J_albedo_true)

    def test_power_absorbed_solar_true(self, r=156, A_rec=2.56, alpha=0.02, epsilon=0.76, planet="Planet_Test",
                                       planet_list=planets_list_test, solar=True):
        P_abs_true = 37.51488367
        P_abs_numerical = power_absorbed(r, A_rec, alpha, epsilon, planet, planet_list, solar)
        self.assertAlmostEqual(P_abs_numerical, P_abs_true)

    def test_power_absorbed_solar_false(self, r=156, A_rec=2.56, alpha=0.02, epsilon=0.76, planet="Planet_Test",
                                        planet_list=planets_list_test, solar=False):
        P_abs_true = 9.690178656
        P_abs_numerical = power_absorbed(r, A_rec, alpha, epsilon, planet, planet_list, solar)
        self.assertAlmostEqual(P_abs_numerical, P_abs_true)

    def test_power_emitted(self, A_emit=11.56, epsilon=0.76, T_opt=200):
        P_emit_true = 797.029632
        P_emit_numerical = power_emitted(A_emit, epsilon, T_opt)
        self.assertAlmostEqual(P_emit_numerical, P_emit_true)

    def test_power_dissipated(self, power_em=200, power_abs=150):
        P_dis_true = -50
        P_dis_numerical = power_dissipated(power_em, power_abs)
        self.assertAlmostEqual(P_dis_numerical, P_dis_true)

    def test_distance_rtg(self, n_rtg=2, p_rtg=500, p_diss=-70, A_rec=5.9, alpha=0.02):
        d_rtg_true = 0.3662579427
        d_rtg_numerical = distance_rtg(n_rtg, p_rtg, p_diss, A_rec, alpha)
        self.assertAlmostEqual(d_rtg_numerical, d_rtg_true)

    def test_louvres_area(self, p_diss=-70, A_rec=5.73, alpha=0.02, d_rtg_uranus=2.5, A_rtg=7.63, p_rtg_tot=11370,
                          n_rtg=5, A_single_l=0.2):
        n_louvre_true = 30
        n_louvre_numerical = louvres_area(p_diss, A_rec, alpha, d_rtg_uranus, A_rtg, p_rtg_tot, n_rtg, A_single_l)
        self.assertEqual(n_louvre_numerical, n_louvre_true)

    def test_power_phases(self, A_rec=5.73, A_emit=11.56, n_rtg=5, planet_list=planets_list_test, r_orbit=156, alpha=0.02,
                          epsilon=0.76, p_rtg_tot=11370, A_single_l=0.2, T_operational = 200, A_rtg=7.63):
        areas_true = [('Planet_Test', 191), ('Venus', 191)]
        d_rtg_true = 0.8082663753
        d_rtg_numerical, areas_numerical = power_phases(A_rec, A_emit, n_rtg, planet_list, r_orbit, alpha, epsilon, p_rtg_tot, A_single_l, T_operational, A_rtg)
        self.assertEqual(areas_numerical[0], areas_true[0])
        self.assertEqual(areas_numerical[1], areas_true[1])
        self.assertAlmostEqual(d_rtg_numerical, d_rtg_true, 5)

