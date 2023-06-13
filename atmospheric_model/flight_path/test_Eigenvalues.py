from unittest import TestCase
from atmospheric_model.flight_path.Eigenvalues import *

class Test(TestCase):
    # def test_eigenvalue_short_period_v_constant_stable(self):
    #
    # def test_eigenvalue_short_period_v_constant_unstable(self):
    #
    # def test_eigenvalue_short_period_v_constant_gamma_constant_stable(self):
    #
    # def test_eigenvalue_short_period_v_constant_gamma_constant_unstable(self):
    #
    #
    # def test_eigenvalue_phugoid_q_dot_zero_alpha_zero_stable(self):
    #
    #
    # def test_eigenvalue_phugoid_q_dot_zero_alpha_zero_unstable(self):
    #
    #
    # def test_eigenvalue_phugoid_q_dot_zero_alpha_dot_zero_stable(self):
    #
    # def test_eigenvalue_phugoid_q_dot_zero_alpha_dot_zero_unstable(self):
    #
    # def test_eigenvalue_heavily_damped_aperiodic_roll_stable(self):
    # def test_eigenvalue_heavily_damped_aperiodic_roll_unstable(self):
    #
    # def test_eigenvalue_dutch_roll_phi_zero_stable(self):
    # def test_eigenvalue_dutch_roll_phi_zero_unstable(self):
    # def test_eigenvalue_dutch_roll_phi_zero_yaw_only_stable(self):
    #
    # def test_eigenvalue_dutch_roll_phi_zero_yaw_only_unstable(self):
    #
    # def test_eigenvalue_aperiodic_spiral_stable(self):
    #
    # def test_eigenvalue_aperiodic_spiral_unstable(self):
    #
    # def test_eigenvalue_dutch_roll_plus_aperiodic_spiral_motion_stable(self):
    #
    # def test_eigenvalue_dutch_roll_plus_aperiodic_spiral_motion_unstable(self):

    def test_symmetric_real_eigenvalues_analysis_real(self, lambda1=-0.237, V=50, c=2.3):
        T_half_true = 0.1345348958
        tau_true = 0.194092827
        T_half_numerical, tau_numerical = symmetric_real_eigenvalues_analysis(lambda1,V, c)
        self.assertAlmostEqual(T_half_numerical, T_half_true)
        self.assertAlmostEqual(tau_numerical, tau_true)

    def test_symmetric_real_eigenvalues_analysis_complex(self, lambda1=2+3j, V=50, c=2.3):
        output_true = False
        output_numerical = symmetric_real_eigenvalues_analysis(lambda1,V, c)
        self.assertAlmostEqual(output_numerical, output_true)

    # def test_symmetric_complex_eigenvalues_analysis_real(self):
    #
    # def test_symmetric_complex_eigenvalues_analysis_complex(self):
    # def test_asymmetric_real_eigenvalues_analysis_real(self):
    #
    # def test_asymmetric_real_eigenvalues_analysis_complex(self):
    # def test_asymmetric_complex_eigenvalues_analysis_real(self):
    # def test_asymmetric_complex_eigenvalues_analysis_complex(self):
