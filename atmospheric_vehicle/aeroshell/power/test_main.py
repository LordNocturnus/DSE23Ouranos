from unittest import TestCase
from atmospheric_vehicle.aeroshell.power.main import *

class Test(TestCase):
    def test_f_tot_p_req(self, adcs_p_req=100, pr_p_req=25):
        tot_p_req_true = 125
        tot_p_req_numerical =f_tot_p_req(adcs_p_req, pr_p_req)
        self.assertAlmostEqual(tot_p_req_numerical, tot_p_req_true)

    def test_f_t(self, days=5):
        t_true = 120
        t_numerical = f_t(days)
        self.assertAlmostEqual(t_numerical, t_true)

    def test_f_m_bat(self, tot_p_req=125, t=120, dod=0.9, n_bat=0.7, n_cab=0.3, spec_energy=128):
        m_bat_true = 620.03968253968
        m_bat_numerical = f_m_bat(tot_p_req, t, dod, n_bat, n_cab, spec_energy)
        self.assertAlmostEqual(m_bat_numerical, m_bat_true)

    def test_f_v_bat(self, tot_p_req=125, t=120, dod=0.9, n_bat=0.7, n_cab=0.3, vol_energy=293):
        v_bat_true = 270.87057803781
        v_bat_numerical = f_v_bat(tot_p_req, t, dod, n_bat, n_cab, vol_energy)
        self.assertAlmostEqual(v_bat_numerical, v_bat_true)
