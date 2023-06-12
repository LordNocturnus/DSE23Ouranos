import cmath
from Derivatives import *
from GliderVariables import *


# Equations for Eigenvalue Calculataions
def eigenvalue_short_period_V_constant():
    A = 2 * mu_c * K_Y_2 * (2 * mu_c - C_Z_alpha_dot)
    B = -2 * mu_c * K_Y_2 * C_Z_alpha - (2 * mu_c + C_Z_q) * C_m_alpha_dot - (2 * mu_c - C_Z_alpha_dot) * C_m_q
    C = C_Z_alpha * C_m_q - (2 * mu_c + C_Z_q) * C_m_alpha
    lambda1 = -B / (2 * A) + cmath.sqrt(B * B - 4 * A * C) / (2 * A)
    if lambda1.real < 0:
        print("The short period with constant velocity is stable (lambda = " +str(lambda1) + " )")
    else:
        print("The short period with constant velocity is not stable (lambda = " +str(lambda1) + " )")
    return lambda1


def eigenvalue_short_period_V_constant_gamma_constant():
    A = 2 * mu_c * K_Y_2
    B = C_m_alpha_dot + C_m_q
    C = C_m_alpha
    lambda1 = -B / (2 * A) + cmath.sqrt(B * B - 4 * A * C) / (2 * A)
    if lambda1.real < 0:
        print("The short period with constant velocity and constant flight path angle is stable (lambda = " +str(lambda1) + " )")
    else:
        print("The short period with constant and constant flight path angle velocity is not stable (lambda = " +str(lambda1) + " )")
    return lambda1


def eigenvalue_phugoid_q_dot_zero_alpha_zero():
    A = -4 * mu_c * mu_c
    B = 2 * mu_c * C_X_u
    C = - C_Z_u * C_Z_0
    lambda1 = -B / (2 * A) + cmath.sqrt(B * B - 4 * A * C) / (2 * A)
    if lambda1.real < 0:
        print("The phugoid with q_dot and angle of attack zero is stable (lambda = " +str(lambda1) + " )")
    else:
        print("The phugoid with q_dot and angle of attack zero is not stable (lambda = " +str(lambda1) + " )")
    return lambda1


def eigenvalue_phugoid_q_dot_zero_alpha_dot_zero():
    A = 2 * mu_c * (C_Z_alpha * C_m_q - 2 * mu_c * C_m_alpha)
    B = 2 * mu_c * (C_X_u * C_m_alpha - C_m_u * C_X_alpha) + C_m_q * (C_Z_u * C_X_alpha - C_X_u * C_Z_alpha)
    C = C_Z_0 * (C_m_u * C_Z_alpha - C_Z_u * C_m_alpha)
    lambda1 = -B / (2 * A) + cmath.sqrt(B * B - 4 * A * C) / (2 * A)
    if lambda1.real < 0:
        print("The phugoid with q_dot and alpha_dot zero is stable (lambda = " +str(lambda1) + " )")
    else:
        print("The phugoid with q_dot and alpha_dot zero is not stable (lambda = " +str(lambda1) + " )")
    return lambda1


def eigenvalue_heavily_damped_aperiodic_roll():
    A = 4 * mu_b * K_X_2
    lambda1 = C_l_p / A
    if lambda1 < 0:
        print("The aperiodic roll is stable (lambda = " +str(lambda1) + " )")
    else:
        print("The aperiodic roll is not stable (lambda = " +str(lambda1) + " )")
    return lambda1


def eigenvalue_dutch_roll_phi_zero():
    A = 8 * mu_b * mu_b * K_Z_2
    B = -2 * mu_b * (C_n_r + 2 * K_Z_2 * C_Y_beta)
    C = 4 * mu_b * C_n_beta + C_Y_beta * C_n_r
    lambda1 = -B / (2 * A) + cmath.sqrt(B * B - 4 * A * C) / (2 * A)
    if lambda1.real < 0:
        print("The dutch roll with phi zero is stable (lambda = " +str(lambda1) + " )")
    else:
        print("The dutch roll with phi zero is not stable (lambda = " +str(lambda1) + " )")
    return lambda1


def eigenvalue_dutch_roll_phi_zero_yaw_only():
    A = -2 * mu_b * K_Z_2
    B = 0.5 * C_n_r
    C = -C_n_beta
    lambda1 = -B / (2 * A) + cmath.sqrt(B * B - 4 * A * C) / (2 * A)
    if lambda1.real < 0:
        print("The dutch roll with phi zero and yaw only is stable (lambda = " +str(lambda1) + " )")
    else:
        print("The dutch roll with phi zero and yaw only is not stable (lambda = " +str(lambda1) + " )")
    return lambda1


def eigenvalue_aperiodic_spiral():
    A = 2 * C_L_opt * (C_l_beta * C_n_r - C_n_beta * C_l_r)
    B = C_l_p * (C_Y_beta * C_n_r + 4 * mu_b * C_n_beta)
    C = C_n_p * (C_Y_beta * C_l_r + 4 * mu_b * C_l_beta)
    lambda1 = A / (B - C)
    if lambda1 < 0:
        print("The aperiodic spiral is stable (lambda = " +str(lambda1) + " )")
    else:
        print("The aperiodic spiral is not stable (lambda = " +str(lambda1) + " )")
    return lambda1


def eigenvalue_dutch_roll_plus_aperiodic_spiral_motion():
    A = 4 * mu_b * mu_b * (K_X_2 * K_Z_2 - K_XZ * K_XZ)
    B = -mu_b * ((C_l_r + C_n_p) * K_XZ + C_n_r * K_X_2 + C_l_p * K_Z_2)
    C = 2 * mu_b * (C_l_beta * K_XZ + C_n_beta * K_X_2) + 0.25 * (C_l_p * C_n_r - C_n_p * C_l_r)
    D = 0.5 * (C_l_beta * C_n_p - C_n_beta * C_l_p)
    lambdas = np.roots([A, B, C, D])
    lambda1 = lambdas[0]
    lambda2 = lambdas[1]
    lambda3 = lambdas[2]

    if lambda1.real < 0 and lambda2.real < 0 and lambda3.real < 0:
        print("The dutch roll plus spiral is stable (lambda = " +str(lambda1) + " )")
    elif lambda1.real < 0 or lambda2.real < 0 or lambda3.real < 0:
        print("The dutch roll plus spiral is not stable (lambda = " +str(lambda1) + " )")
    return lambdas


def symmetric_real_eigenvalues_analysis(lambda1, V):
    if isinstance(lambda1, complex):
        print("Given eigenvalue is complex, wrong function used")
    else:
        T_half = np.log(0.5) / lambda1 * chord / V
        tau = -1 / lambda1 * chord / V
        return T_half, tau


def symmetric_complex_eigenvalues_analysis(lambda1, V):
    if lambda1.imag == 0:
        print("Given eigenvalue is real, wrong function used")
    else:
        xi_c = np.real(lambda1)
        eta_c = np.imag(lambda1)
        P = 2 * np.pi / eta_c * chord / V
        T_half = np.log(0.5) / xi_c * chord / V
        return P, T_half


def asymmetric_real_eigenvalues_analysis(lambda1, V):
    if isinstance(lambda1, complex):
        print("Given eigenvalue is complex, wrong function used")
    else:
        T_half = np.log(0.5) / lambda1 * b / V
        return T_half


def asymmetric_complex_eigenvalues_analysis(lambda1, V):
    if isinstance(lambda1, float):
        print("Given eigenvalue is real, wrong function used")
    else:
        xi_c = np.real(lambda1)
        eta_c = np.imag(lambda1)
        P = 2 * np.pi / eta_c * b / V
        T_half = np.log(0.5) / xi_c * b / V
        return P, T_half

