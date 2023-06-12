from unittest import TestCase
from Orbiter.propulsion.main import *


class Test(TestCase):
    def test_mass_prop(self, m_orbiter=200, m_combined=350, dV_orbiter=223, dV_combined=576, g=9.81, I_sp=325):
        m_prop_orbiter_true = 14.4896946
        m_prop_comb_true = 69.30414003
        m_prop_true_true = 83.79383463
        m_prop_orbiter_numerical, m_prop_comb_numerical, m_prop_true = mass_prop(m_orbiter, dV_orbiter, m_combined,
                                                                                 dV_combined, g, I_sp)
        self.assertAlmostEqual(m_prop_orbiter_numerical, m_prop_orbiter_true)
        self.assertAlmostEqual(m_prop_comb_numerical, m_prop_comb_true)
        self.assertAlmostEqual(m_prop_true, m_prop_true_true)

    def test_burntimecombined(self, T=100, m_f=150, dV_combined=576, m_orbiter=200, m_combined=350, dV_orbiter=223,
                              g=9.81, I_sp=325):
        tb_true = 0.3740701354
        tb_numerical = burntimecombined(T, m_f, dV_combined, m_orbiter, m_combined, dV_orbiter, g, I_sp)
        self.assertAlmostEqual(tb_numerical, tb_true)

    def test_burntimeorbiter(self,T=100, m_f=150, dV_orbiter=223, m_orbiter=200, m_combined=350, dV_combined=576, g=9.81, I_sp=325):
        tb_true = 0.1018922275
        tb_numerical = burntimeorbiter(T, m_f, dV_orbiter,m_orbiter, m_combined, dV_combined, g, I_sp)
        self.assertAlmostEqual(tb_numerical, tb_true)
