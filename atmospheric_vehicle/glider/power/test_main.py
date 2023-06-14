from unittest import TestCase
from atmospheric_vehicle.glider.power.main import *

class Test(TestCase):
    def test_f_tot_p_req(self,dh_p_req=670, ttc_p_req=52, adcs_p_req=63, pl_p_req=723, tm_p_req=1.23):
        tot_p_req_true = 1509.23
        tot_p_req_numerical = f_tot_p_req(dh_p_req, ttc_p_req, adcs_p_req, pl_p_req, tm_p_req)
        self.assertAlmostEquals(tot_p_req_numerical, tot_p_req_true)

    def test_f_t(self, days=3):
        t_true = 72
        t_numerical = f_t(days)
        self.assertAlmostEquals(t_numerical, t_true)

    def test_f_m_bat(self, tot_p_req=1509.23, t=72, dod=0.53, n_bat=0.95, n_cab=0.67, spec_energy=123.2):
        m_bat_true = 2614.586054
        m_bat_numerical = f_m_bat(tot_p_req, t, dod, n_bat, n_cab, spec_energy)
        self.assertAlmostEquals(m_bat_numerical, m_bat_true,6)


    def test_f_v_bat(self, tot_p_req=1509.23, t=72, dod=0.53, n_bat=0.95, n_cab=0.67,vol_energy=23.5):
        v_bat_true=13707.10646
        v_bat_numerical=f_v_bat(tot_p_req, t, dod, n_bat, n_cab,vol_energy)
        self.assertAlmostEquals(v_bat_numerical, v_bat_true,5)

