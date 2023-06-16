from unittest import TestCase
from atmospheric_vehicle.glider.structure.main import *

tau_yield = 508.2 * 10 ** 6 * 1.1
sigma_yield = 924 * 10 ** 6 * 1.1
E = 120 * 10 ** 9
alpha = 9.1 * 10 ** -6
c_w_root = 0.785  # [m], root chord, given from Luigi
c_w_tip = 0.713  # [m], tip chord, given from Luigi
b_w = 6  # [m], wing span, given from Luigi
c_w = calculate_c_w(c_w_root, c_w_tip, b_w)[0]  # array of chord values that will be used for later estimations
b_range = calculate_c_w(c_w_root, c_w_tip, b_w)[1]
db = calculate_c_w(c_w_root, c_w_tip, b_w)[2]
A = 0.28
B = 0.08
C = 0.68
a_2 = c_w * A
t = c_w * B
a_3 = c_w * C
L_distr = 224.482
T_space = 2.73
T_01 = 53
T_1 = 76
T_20 = 193

Xcentr_tot = calculate_tau_xz_torque(tau_yield, c_w, b_w, t, a_2, a_3, L_distr)[1]
Ty = calculate_tau_xz_torque(tau_yield, c_w, b_w, t, a_2, a_3, L_distr)[2]
Am = calculate_tau_xz_torque(tau_yield, c_w, b_w, t, a_2, a_3, L_distr)[3]

t_t_bend = calculate_sigma_bend_xz(sigma_yield, c_w_root, c_w_tip, b_w, L_distr, b_range, t, a_2, a_3)[0]
Mx = calculate_sigma_bend_xz(sigma_yield, c_w_root, c_w_tip, b_w, L_distr, b_range, t, a_2, a_3)[1]
Ixx = calculate_sigma_bend_xz(sigma_yield, c_w_root, c_w_tip, b_w, L_distr, b_range, t, a_2, a_3)[2]

max_delta = calculate_sigma_thermal(T_space, T_01, T_1, T_20, sigma_yield, alpha, E)[1]

t_t_I, S_shear = calculate_sigma_xz_buckling(E, c_w, db, L_distr, t)

a = calculate_hoop_stress_wing(sigma_yield, c_w, b_w)[1]

class GliderStructureTests(TestCase):

    def test_calculate_c_w(self):
        self.assertEqual(c_w[0], c_w_root)
        self.assertEqual(c_w[-1], c_w_tip)
        self.assertEqual(c_w[15], 0.749)

    def test_calculate_tau_xz_torque(self):
    #Xcentr_tot correctness checked on analytical calculations, found a missing /2 for L_res
        L_res = L_distr * (c_w_root + c_w_tip) * b_w / 4
        self.assertAlmostEqual(abs(Ty[0]), abs(L_res * (Xcentr_tot[0] - t[0] / 2 - a_2[0] / 2)))
        self.assertAlmostEqual(Am[0], t[0] * (np.pi / 8 * t[0] + a_2[0] + a_3[0] / 2))
    #Am and Ty already correct, hence t_t correct as well

    def test_calculate_sigma_bend_xz(self):
        self.assertAlmostEqual(Mx[0], -2 / 3 * (c_w_root - c_w_tip) / b_w * L_distr * b_range[0]**3 +
                               c_w_root / 2 * L_distr * b_range[0]**2)
        self.assertAlmostEqual(Ixx[0], t_t_bend[0] * (1 / 2 * t[0]**2 * a_2[0] + 1 / 6 * t[0]**3 + (np.pi / 16 + 1 / 12)
                                                      * t[0]**3 + 1 / 8 * t[0]**2 * a_3[0] + 1 / 24 * t[0]**3))
    #shear in the xz correctness checked on analytical calculations, found a sign error in equation for moment equivalence
    def test_calculate_sigma_thermal(self):
        self.assertAlmostEqual(max_delta, 117)

    def test_calculate_sigma_xz_buckling(self): # found error in using thin-walled approx: didn't neglect one 2nd order term
        self.assertAlmostEqual(S_shear[0], (c_w[0] + c_w[1]) * db / 2)
        self.assertAlmostEqual(t_t_I[0], (c_w[0] + c_w[1]) * db / 2 * 12 / (np.pi ** 2 * E * t[1]))

    #shear in yz shall be redone, hence unit testing will be done later (if there is time)

    def test_calculate_hoop_stress_wing(self):
        self.assertAlmostEqual(a, 0.6875 * np.pi / 180, 6)

    # mass function checked analytically, everything correct
