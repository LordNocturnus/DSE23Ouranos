from unittest import TestCase
from atmospheric_model.flight_path.GlideTrajectory import *

class Test(TestCase):
    def test_velocity(self, rho=1.56, m=100, g=6.2, S=12.37, C_L=0.72):
        velocity_true = 9.447076787
        velocity_numerical = velocity(rho, m, g, S, C_L)
        self.assertAlmostEqual(velocity_numerical, velocity_true)

    def test_range(self, h=1236, C_L_C_D=12.3):
        R_true = 15202.8
        R_numerical = range(h, C_L_C_D)
        self.assertAlmostEqual(R_numerical, R_true)

    def test_time_of_flight(self, delta_h=100, rho=1.56, S=12.37, m=100, g=6.2, C_L=0.72, C_D=0.032):
        t_true = 238.168912
        t_numerical = time_of_flight(delta_h, rho, S, m, g, C_L, C_D)
        self.assertAlmostEqual(t_numerical, t_true)