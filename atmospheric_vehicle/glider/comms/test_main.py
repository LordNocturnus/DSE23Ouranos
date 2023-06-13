from unittest import TestCase

import numpy as np

from atmospheric_vehicle.glider.comms.main import *


class Test(TestCase):
    def test_unit_to_db(self, unit=73.5):
        x_dB_true = 18.66287339
        x_dB_numerical = unit_to_db(unit)
        self.assertAlmostEqual(x_dB_numerical, x_dB_true)

    def test_db_to_unit(self, db=-273):
        x_unit_true = 5.01187234 * 10 ** -28
        x_unit_numerical = db_to_unit(db)
        self.assertAlmostEqual(x_unit_numerical, x_unit_true)

    def test_l_fs(self, freq=2.0 * 10 ** 6, r=35.2 * 10 ** 9):
        FSL_true = 189.3932254
        FSL_numerical = l_fs(freq, r)
        self.assertAlmostEqual(FSL_numerical, FSL_true)

    def test_s(self, r=2):
        s_true = 4 * np.pi
        s_numerical = s(r)
        self.assertAlmostEqual(s_numerical, s_true)

    def test_p_rx(self, p_tx=5, g_tx_ao=10.6, l_tx=2.3, l_fs=3.9, l_m=12.9, g_rx=23.7, l_rx=9.6):
        p_rx_true = 10.6
        p_rx_numerical = p_rx(p_tx, g_tx_ao, l_tx, l_fs, l_m, g_rx, l_rx)
        self.assertAlmostEqual(p_rx_numerical, p_rx_true)

    def test_e_b(self, p_rx=10.6, bit_rate=123):
        E_b_true = 0.086178862
        E_b_numerical = e_b(p_rx, bit_rate)
        self.assertAlmostEqual(E_b_numerical, E_b_true)

    def test_n_0(self, t_s=2.36 * 10 ** 16):
        n_0_true = 3.26 * 10 ** -7
        n_0_numerical = n_0(t_s)
        self.assertAlmostEqual(n_0_numerical, n_0_true)

    def test_energy_per_bit_criteria_true(self, e_b=50, n_0=2):
        ratio_true = 25
        test_true = True
        ratio_numerical, test_numerical = energy_per_bit_criteria(e_b, n_0)
        self.assertEqual(test_numerical, test_true)
        self.assertAlmostEqual(ratio_numerical, ratio_true)

    def test_energy_per_bit_criteria_true(self, e_b=2, n_0=5):
        ratio_true = 0.4
        test_true = False
        ratio_numerical, test_numerical = energy_per_bit_criteria(e_b, n_0)
        self.assertEqual(test_numerical, test_true)
        self.assertAlmostEqual(ratio_numerical, ratio_true)
