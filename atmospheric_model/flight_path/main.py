import numpy as np
from matplotlib import pyplot as plt
from atmospheric_model import GRAM
import cmath

# Variables
gamma = np.radians(1.5)
C_L_C_D = 1 / np.sin(gamma)
m = 72  # kg
g_u = 8.69  # m/s^2
S = 4.495  # m^2
b = 6  # m
chord = 1  # m
A = b * b / S
e = 1
C_D_0 = np.pi * A * e / (4 * C_L_C_D * C_L_C_D)
C_L_opt = np.sqrt(C_D_0 * np.pi * A * e)
C_D = 2 * C_D_0
p_max = 20 * 10 ** 5  # Pa

# Control and Stability Derivatives
mu_c = 1
mu_b = 1
K_X_2 = 1
K_Y_2 = 1
K_Z_2 = 1
K_XZ = 1
C_X_u = 1
C_X_alpha = 1
C_Y_beta = 1
C_Z_alpha_dot = 1
C_Z_alpha = 1
C_Z_q = 1
C_Z_u = 1
C_Z_0 = 1
C_m_alpha_dot = 1
C_m_q = 1
C_m_alpha = 1
C_m_u = 1
C_z_alpha_dot = 1
C_l_p = 1
C_l_beta = 1
C_l_r = 1
C_n_r = 1
C_n_beta = 1
C_n_p = 1
C_L = C_L_opt


def velocity(rho):
    a = 2 * m * g_u
    b = rho * S * C_L_opt
    V = np.sqrt(a / b)
    return V


def range(h):
    R = C_L_C_D * h
    return R


def time_of_flight(delta_h, rho):
    a = np.sqrt(S / (2 * m * g_u))
    b = C_L_opt ** (3 / 2) / C_D
    t = delta_h * np.sqrt(rho) * a * b
    return t


def eigenvalue_short_period_V_constant():
    A = 2 * mu_c * K_Y_2 * (2 * mu_c - C_Z_alpha_dot)
    B = -2 * mu_c * K_Y_2 * C_Z_alpha - (2 * mu_c + C_Z_q) * C_m_alpha_dot - (2 * mu_c - C_z_alpha_dot) * C_m_q
    C = C_Z_alpha * C_m_q - (2 * mu_c + C_Z_q) * C_m_alpha
    lambda1 = -B / (2 * A) + cmath.sqrt(B * B - 4 * A * C) / (2 * A)
    if lambda1.real < 0:
        print("The short period with constant velocity is stable")
    else:
        print("The short period with constant velocity is not stable")
    return lambda1


def eigenvalue_short_period_V_constant_gamma_constant():
    A = 2 * mu_c * K_Y_2
    B = C_m_alpha_dot + C_m_q
    C = C_m_alpha
    lambda1 = -B / (2 * A) + cmath.sqrt(B * B - 4 * A * C) / (2 * A)
    if lambda1.real < 0:
        print("The short period with constant velocity and constant flight path angle is stable")
    else:
        print("The short period with constant and constant flight path angle velocity is not stable")
    return lambda1


def eigenvalue_phugoid_q_dot_zero_alpha_zero():
    A = -4 * mu_c * mu_c
    B = 2 * mu_c * C_X_u
    C = - C_Z_u * C_Z_0
    lambda1 = -B / (2 * A) + cmath.sqrt(B * B - 4 * A * C) / (2 * A)
    if lambda1.real < 0:
        print("The phugoid with q_dot and angle of attack zero is stable")
    else:
        print("The phugoid with q_dot and angle of attack zero is not stable")
    return lambda1


def eigenvalue_phugoid_q_dot_zero_alpha_dot_zero():
    A = 2 * mu_c * (C_Z_alpha * C_m_q - 2 * mu_c * C_m_alpha)
    B = 2 * mu_c * (C_X_u * C_m_alpha - C_m_u * C_X_alpha) + C_m_q * (C_Z_u * C_X_alpha - C_X_u * C_Z_alpha)
    C = C_Z_0 * (C_m_u * C_Z_alpha - C_Z_u * C_m_alpha)
    lambda1 = -B / (2 * A) + cmath.sqrt(B * B - 4 * A * C) / (2 * A)
    if lambda1.real < 0:
        print("The phugoid with q_dot and alpha_dot zero is stable")
    else:
        print("The phugoid with q_dot and alpha_dot zero is not stable")
    return lambda1


def eigenvalue_heavily_damped_aperiodic_roll():
    A = 4 * mu_b * K_X_2
    lambda1 = C_l_p / A
    if lambda1 < 0:
        print("The aperiodic roll is stable")
    else:
        print("The aperiodic roll is not stable")
    return lambda1


def eigenvalue_dutch_roll_phi_zero():
    A = 8 * mu_b * mu_b * K_Z_2
    B = -2 * mu_b * (C_n_r + 2 * K_Z_2 * C_Y_beta)
    C = 4 * mu_b * C_n_beta + C_Y_beta * C_n_r
    lambda1 = -B / (2 * A) + cmath.sqrt(B * B - 4 * A * C) / (2 * A)
    if lambda1.real < 0:
        print("The dutch roll with phi zero is stable")
    else:
        print("The dutch roll with phi zero is not stable")
    return lambda1


def eigenvalue_dutch_roll_phi_zero_yaw_only():
    A = -2 * mu_b * K_Z_2
    B = 0.5 * C_n_r
    C = -C_n_beta
    lambda1 = -B / (2 * A) + cmath.sqrt(B * B - 4 * A * C) / (2 * A)
    if lambda1.real < 0:
        print("The dutch roll with phi zero and yaw only is stable")
    else:
        print("The dutch roll with phi zero and yaw only is not stable")
    return lambda1


def eigenvalue_aperiodic_spiral():
    A = 2 * C_L * (C_l_beta * C_n_r - C_n_beta * C_l_r)
    B = C_l_p * (C_Y_beta * C_n_r + 4 * mu_b * C_n_beta)
    C = C_n_p * (C_Y_beta * C_l_r + 4 * mu_b * C_l_beta)
    lambda1 = A / (B - C)
    if lambda1 < 0:
        print("The aperiodic spiral is stable")
    else:
        print("The aperiodic spiral is not stable")
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
    lambda4 = lambdas[3]

    if lambda1.real < 0 and lambda2.real < 0 and lambda3.real < 0 and lambda4.real < 0:
        print("The dutch roll plus spiral is stable")
    elif lambda1.real < 0 or lambda2.real < 0 or lambda3.real < 0 or lambda4.real < 0:
        print("The dutch roll plus spiral is not stable")
    return lambdas


def symmetric_real_eigenvalues_analysis(lambda1, V):
    if isinstance(lambda1, complex):
        print("Given eigenvalue is complex, wrong function used")
    else:
        T_half = np.ln(0.5) / lambda1 * chord / V
        tau = -1 / lambda1 * chord / V
        return T_half, tau


def symmetric_complex_eigenvalues_analysis(lambda1, V):
    if isinstance(lambda1, float):
        print("Given eigenvalue is real, wrong function used")
    else:
        xi_c = np.real(lambda1)
        eta_c = np.image(lambda1)
        P = 2 * np.pi / eta_c * chord / V
        T_half = np.ln(0.5) / xi_c * chord / V
        return P, T_half


def asymmetric_real_eigenvalues_analysis(lambda1, V):
    if isinstance(lambda1, complex):
        print("Given eigenvalue is complex, wrong function used")
    else:
        T_half = np.ln(0.5) / lambda1 * b / V
        return T_half


def asymmetric_complex_eigenvalues_analysis(lambda1, V):
    if isinstance(lambda1, float):
        print("Given eigenvalue is real, wrong function used")
    else:
        xi_c = np.real(lambda1)
        eta_c = np.image(lambda1)
        P = 2 * np.pi / eta_c * b / V
        T_half = np.ln(0.5) / xi_c * b / V
        return P, T_half


V = []
h = []
t = []

if __name__ == "__main__":
    # Setting up GRAM
    gram = GRAM.GRAM()
    # gram.runs = 10
    # gram.densityperturbation = 1.0  # Between 0.0 and 2.0, 1.0 = 3 sigma
    gram.run()

    # Finding indeces for the relevant pressures
    index_1 = gram.data.Pressure_Pa[gram.data.Pressure_Pa > 10 ** 4].index
    index_2 = gram.data.Pressure_Pa[gram.data.Pressure_Pa < p_max].index

    # Find common indeces
    index = np.sort(list(set(index_1) & set(index_2)))

    # Find density and heights from indeces
    density = gram.data.Density_kgm3[index]
    height = gram.data.Height_km[index]

    # Find maximum and minimum index to convert to list
    index_min = min(index)
    index_max = max(index)
    points = np.arange(index_min, index_max + 1, 1)

    # Find velocity, height, and flight time per density interval
    for i in points:
        V.append(velocity(density[i]))
        if i < index_max:
            h.append(height[i])
            t.append(time_of_flight((height[i] - height[i + 1]) * 1000, density[i]))

    # Find maximum flight range
    range = range(np.abs(height[index_max]))

    print("The maximum flight range is " + str(range) + " km")
    print("The time of flight is " + str(sum(t) / (60 * 60 * 24)) + " days")

    # Plot altitude vs time
    plt.plot(t, h)
    plt.xlabel("Time [s]")
    plt.ylabel("Altitude [km]")
    plt.title("Altitude vs Time of the Glider")
    plt.show()

    # Find eigenvalues of the symmetric eigenmotions
    lambda_short_1 = eigenvalue_short_period_V_constant()
    lambda_short_2 = eigenvalue_short_period_V_constant_gamma_constant()
    lambda_phugoid_1 = eigenvalue_phugoid_q_dot_zero_alpha_zero()
    lambda_phugoid_2 = eigenvalue_phugoid_q_dot_zero_alpha_dot_zero()
    lambdas_sym = [lambda_short_1, lambda_short_2, lambda_phugoid_1, lambda_phugoid_2]

    # Find eigenvalues of the asymmetric eigenmotions
    lambda_aperiodic_roll = eigenvalue_heavily_damped_aperiodic_roll()
    lambda_dutch_roll_1 = eigenvalue_dutch_roll_phi_zero()
    lambda_dutch_roll_2 = eigenvalue_dutch_roll_phi_zero_yaw_only()
    lambda_aperiodic_spiral = eigenvalue_aperiodic_spiral()
    lambda_dutch_roll_plus_spiral = eigenvalue_dutch_roll_plus_aperiodic_spiral_motion()
    lambdas_asym = [lambda_aperiodic_roll, lambda_dutch_roll_1, lambda_dutch_roll_2, lambda_aperiodic_spiral,
                    lambda_dutch_roll_plus_spiral[0], lambda_dutch_roll_plus_spiral[1],
                    lambda_dutch_roll_plus_spiral[2], lambda_dutch_roll_plus_spiral[3]]

    # Analysis symmetric eigenvalues
    for i in lambdas_sym:
        if isinstance(i, float):
            T_half, tau = symmetric_real_eigenvalues_analysis(i, V[0])
            print("The halving time and time constant of" + str(i) + " are" + str(T_half) + " s" + str(tau) + " s")
        else:
            P, T_half = symmetric_complex_eigenvalues_analysis(i, V[0])
            print("The halving time and period of" + str(i) + " are" + str(T_half) + " s" + str(P) + " s")

    # Analysis asymmetric eigenvalues
    for i in lambdas_asym:
        if isinstance(i, float):
            T_half = asymmetric_real_eigenvalues_analysis(i, V[0])
            print("The halving time" + str(i) + " is" + str(T_half) + " s")
        else:
            P, T_half = asymmetric_complex_eigenvalues_analysis(i, V[0])
            print("The halving time and period of" + str(i) + " are" + str(T_half) + " s" + str(P) + " s")
