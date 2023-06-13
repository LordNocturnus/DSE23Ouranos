from unittest import TestCase
from atmospheric_model.flight_path.DeploymentTrajectory import *


class Test(TestCase):
    def test_d_v_dt(self,rho=1.56, V=56, gamma=5/8*np.pi, C_D=0.032, m=100, g=6.2, S=12.37):
        dvdt_true = -15.41061617
        dvdt_numerical = dV_dt(rho, V, gamma, C_D, m, g, S)
        self.assertAlmostEqual(dvdt_numerical, dvdt_true)


    def test_dgamma_dt(self, rho=1.56, V=56, gamma=5/8*np.pi, C_L=0.72, m=100, g=6.2, S=12.37):
        dgammadt_true = 3.932684043
        dgammadt_numerical = dgamma_dt(rho, V, gamma, C_L, m, g, S)
        self.assertAlmostEqual(dgammadt_numerical, dgammadt_true)

    def test_v_deploy(self, V_0=25, dVdt=2.3, dt=0.1):
        v_deploy_true = 25.23
        v_deploy_numerical = v_deploy(V_0, dVdt, dt)
        self.assertAlmostEqual(v_deploy_numerical, v_deploy_true)


    def test_gamma_deploy(self, gamma_0=np.pi/3, dgammadt=0.234, dt=0.1):
        gamma_deploy_test = 1.070597551
        gamma_deploy_numerical = gamma_deploy(gamma_0, dgammadt, dt)
        self.assertAlmostEqual(gamma_deploy_numerical, gamma_deploy_test)
