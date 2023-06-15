from unittest import TestCase
from atmospheric_model.flight_path.Eigenvalues import *


class Test(TestCase):
    def test_eigenvalue_short_period_v_constant(self, mu_c=1, K_Y_2=1, C_Z_alpha_dot=1, C_Z_alpha=1, C_Z_q=1
                                                ,C_m_alpha_dot=1, C_m_q=1, C_m_alpha=1):
        lambda1_true = 3.302775638
        lambda1_numerical = eigenvalue_short_period_V_constant(mu_c, K_Y_2, C_Z_alpha_dot, C_Z_alpha, C_Z_q, C_m_alpha_dot, C_m_q, C_m_alpha)
        self.assertAlmostEqual(lambda1_numerical, lambda1_true)


    def test_eigenvalue_short_period_v_constant_gamma_constant(self,mu_c=1, K_Y_2=1, C_m_alpha_dot=1,
                                                      C_m_q=1, C_m_alpha=1):
        lambda1_true = -0.5 + 0.5j
        lambda1_numerical = eigenvalue_short_period_V_constant_gamma_constant(mu_c, K_Y_2, C_m_alpha_dot,
                                                      C_m_q, C_m_alpha)
        self.assertAlmostEqual(lambda1_numerical, lambda1_true)


    def test_eigenvalue_phugoid_q_dot_zero_alpha_zero(self, mu_c=1, C_X_u=1, C_Z_u=1, C_Z_0=1):
        lambda1_true = 0.25 - 0.4330127019j
        lambda1_numerical = eigenvalue_phugoid_q_dot_zero_alpha_zero(mu_c, C_X_u, C_Z_u, C_Z_0)
        self.assertAlmostEqual(lambda1_numerical, lambda1_true)



    def test_eigenvalue_phugoid_q_dot_zero_alpha_dot_zero(self, mu_c=2, C_Z_alpha=1, C_m_q=1, C_m_alpha=-1,
                                                 C_X_u=1, C_m_u=2, C_X_alpha=1, C_Z_u=1,
                                                 C_Z_0=-1):
        lambda1_true = 0.7898979486 + 0j
        lambda1_numerical = eigenvalue_phugoid_q_dot_zero_alpha_dot_zero(mu_c, C_Z_alpha, C_m_q, C_m_alpha,
                                                 C_X_u, C_m_u, C_X_alpha, C_Z_u, C_Z_0)
        self.assertAlmostEqual(lambda1_numerical, lambda1_true)


    def test_eigenvalue_heavily_damped_aperiodic_roll(self, mu_b=1, K_X_2=1, C_l_p=1):
        lambda1_true = 0.25
        lambda1_numerical = eigenvalue_heavily_damped_aperiodic_roll(mu_b, K_X_2, C_l_p)
        self.assertAlmostEqual(lambda1_numerical, lambda1_true)

    def test_eigenvalue_dutch_roll_phi_zero(self, mu_b=1, K_Z_2=1, C_n_r=1, C_Y_beta=1,C_n_beta=1):
        lambda1_true = 0.375 + 0.6959705454j
        lambda1_numerical = eigenvalue_dutch_roll_phi_zero(mu_b, K_Z_2, C_n_r, C_Y_beta,C_n_beta)
        self.assertAlmostEqual(lambda1_numerical, lambda1_true)

    def test_eigenvalue_dutch_roll_phi_zero_yaw_only(self, mu_b=1, K_Z_2=1, C_n_r=1, C_n_beta=1):
        lamda1_true = 0.125 - 0.6959705454j
        lamda1_numerical = eigenvalue_dutch_roll_phi_zero_yaw_only(mu_b, K_Z_2, C_n_r, C_n_beta)
        self.assertAlmostEqual(lamda1_numerical, lamda1_true)


    def test_eigenvalue_aperiodic_spiral(self, C_L=1, C_l_beta=2, C_n_r=1, C_n_beta=1, C_l_r=1,
                                C_l_p=1, C_Y_beta=1, mu_b=1, C_n_p=1):
        lambda1_true = -0.5
        lambda1_numerical = eigenvalue_aperiodic_spiral(C_L, C_l_beta, C_n_r, C_n_beta, C_l_r,
                                C_l_p, C_Y_beta, mu_b, C_n_p)
        self.assertAlmostEqual(lambda1_numerical, lambda1_true)


    def test_eigenvalue_dutch_roll_plus_aperiodic_spiral_motion(self, mu_b=1, K_X_2=2, K_Z_2=1, K_XZ=1, C_l_r=1,
                                                       C_n_p=1, C_n_r=1, C_l_p=1, C_l_beta=1,C_n_beta=1):
        lambdas_true = [0.625 + 1.053268722j, 0.625 - 1.053268722j, 0]
        lambdas_numerical = eigenvalue_dutch_roll_plus_aperiodic_spiral_motion(mu_b, K_X_2, K_Z_2, K_XZ, C_l_r,
                                                       C_n_p, C_n_r, C_l_p, C_l_beta, C_n_beta)

        self.assertAlmostEqual(lambdas_numerical[0], lambdas_true[0])
        self.assertAlmostEqual(lambdas_numerical[1], lambdas_true[1])
        self.assertAlmostEqual(lambdas_numerical[2], lambdas_true[2])


    def test_symmetric_real_eigenvalues_analysis_real(self, lambda1=-0.237, V=50, c=2.3):
        T_half_true = 0.1345348958
        tau_true = 0.194092827
        T_half_numerical, tau_numerical = symmetric_real_eigenvalues_analysis(lambda1, V, c)
        self.assertAlmostEqual(T_half_numerical, T_half_true)
        self.assertAlmostEqual(tau_numerical, tau_true)

    def test_symmetric_real_eigenvalues_analysis_complex(self, lambda1=2 + 3j, V=50, c=2.3):
        output_true = False
        output_numerical = symmetric_real_eigenvalues_analysis(lambda1, V, c)
        self.assertAlmostEqual(output_numerical, output_true)

    def test_symmetric_complex_eigenvalues_analysis_real(self, lambda1=2 + 0j, V=50, c=2.3):
        output_true = False
        output_numerical = symmetric_complex_eigenvalues_analysis(lambda1, V, c)
        self.assertAlmostEqual(output_numerical, output_true)

    def test_symmetric_complex_eigenvalues_analysis_complex(self, lambda1=2 + 3j, V=50, c=2.3):
        T_half_true = -0.01594238515
        P_true = 0.09634217471
        P_numerical, T_half_numerical = symmetric_complex_eigenvalues_analysis(lambda1, V, c)
        self.assertAlmostEqual(T_half_numerical, T_half_true)
        self.assertAlmostEqual(P_numerical, P_true)

    def test_asymmetric_real_eigenvalues_analysis_real(self, lambda1=2 + 3j, V=50, b=20):
        output_true = False
        output_numerical = asymmetric_real_eigenvalues_analysis(lambda1, V, b)
        self.assertAlmostEqual(output_numerical, output_true)

    def test_asymmetric_real_eigenvalues_analysis_complex(self, lambda1=2, V=50, b=20):
        T_half_true = -0.1386294361
        T_half_numerical = asymmetric_real_eigenvalues_analysis(lambda1, V, b)
        self.assertAlmostEqual(T_half_numerical, T_half_true)

    def test_asymmetric_complex_eigenvalues_analysis_real(self, lambda1=1.0, V=50, b=20):
        output_true = False
        output_numerical = asymmetric_complex_eigenvalues_analysis(lambda1, V, b)
        self.assertAlmostEqual(output_true, output_numerical)

    def test_asymmetric_complex_eigenvalues_analysis_complex(self, lambda1=1.0 + 5j, V=50, b=20):
        P_true = 0.5026548246
        T_half_true = -0.2772588722
        P_numerical, T_half_numerical = asymmetric_complex_eigenvalues_analysis(lambda1, V, b)
        self.assertAlmostEqual(P_numerical, P_true)
        self.assertAlmostEqual(T_half_numerical, T_half_true)
