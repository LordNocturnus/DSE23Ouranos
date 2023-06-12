from unittest import TestCase
from Orbiter.power.main import *


class Test(TestCase):
    def test_powerdecay(self, P_0=2561, tau=2.5, t=30.1):
        P_decay_true = 0.6081468067
        P_decay_numerical = powerdecay(P_0, tau, t)
        self.assertAlmostEquals(P_decay_numerical, P_decay_true)

    def test_number_rtg(self, P_req=156, P_start=2561, tau=2.5, t=30.1):
        M_true = 301.7847167
        N_true = 302
        M_numerical, N_numerical = numberRTG(P_req, tau, t, P_start)
        self.assertAlmostEquals(M_numerical, M_true)
        self.assertAlmostEquals(N_numerical, N_true)

    def test_mass_rtg(self, m_RTG=2.5, tau=2.5, t=30.1, P_start=2561, P_req=156):
        mass_true = 755
        mass_numerical = massRTG(m_RTG, tau, t, P_start, P_req)
        self.assertAlmostEquals(mass_numerical, mass_true)

    def test_cost_rtg(self, costRTG1=2.0* 10**3, P_req = 156, P_start=2561, tau=2.5, t=30.1):
        cost_true = 604000
        cost_numerical = costRTG(costRTG1, P_req, P_start, tau, t)