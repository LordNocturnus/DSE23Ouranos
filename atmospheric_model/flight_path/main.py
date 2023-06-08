import numpy as np
from matplotlib import pyplot as plt
from atmospheric_model import GRAM
import cmath

# Variables
gamma_glide = np.radians(-1.5)
C_L_C_D = 1 / np.sin(abs(gamma_glide))
m = 72  # kg
g_u = 8.69  # m/s^2
S = 4.495  # m^2
b = 6  # m
chord = 1  # m
A = b * b / S
e = 1
C_D_0 = np.pi * A * e / (4 * C_L_C_D * C_L_C_D)
C_L_opt = np.sqrt(C_D_0 * np.pi * A * e)
C_L_deploy = 0.2
C_D = 2 * C_D_0
C_D_deploy = 0.03
p_max = 20 * 10 ** 5  # Pa
gamma_0_deploy = -np.pi / 2  # rad
V_capsule = 80  # m/s
V_0_deploy = V_capsule  # m/s

# Dimensionless parameters
mu_c = 0.1
mu_b = 0.2
K_X_2 = 0.3
K_Y_2 = 0.4
K_Z_2 = 0.5
K_XZ = 0.6

# Stability Derivatives
C_X_u = 0.7
C_Z_u = 1.3
C_m_u = 1.8

C_X_alpha = 0.8
C_Z_alpha = 1.1
C_m_alpha = 1.7

C_Y_beta = 0.9
C_l_beta = 2.1
C_n_beta = 2.4

C_l_p = 2
C_n_p = 2.5

C_l_r = 2.2
C_n_r = 2.3

# Other Derivatives
C_Z_alpha_dot = 1
C_Z_q = 1.2
C_Z_0 = 1.4
C_m_alpha_dot = 1.5
C_m_q = 1.6
C_z_alpha_dot = 1.9


# Equations for Gliding Flight Calculations
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


# Equations for Glider Deployment

def dV_dt(rho, V, gamma):
    a = -0.5 * rho * V * V * S * C_D / m
    b = -g_u * np.sin(gamma)
    dVdt = a + b
    return dVdt


def dgamma_dt(rho, V, gamma):
    a = 0.5 * rho * V * S * C_L_deploy / m
    b = -g_u * np.cos(gamma) / V
    dgammadt = a + b
    return dgammadt


def v_deploy(V_0, dVdt, dt):
    return V_0 + dVdt * dt


def gamma_deploy(gamma_0, dgammadt, dt):
    return gamma_0 + dgammadt * dt


# Equations for Eigenvalue Calculataions
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
    A = 2 * C_L_opt * (C_l_beta * C_n_r - C_n_beta * C_l_r)
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

    if lambda1.real < 0 and lambda2.real < 0 and lambda3.real < 0:
        print("The dutch roll plus spiral is stable")
    elif lambda1.real < 0 or lambda2.real < 0 or lambda3.real < 0:
        print("The dutch roll plus spiral is not stable")
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


V_flight_path = []
h_flight_path = []
t_flight_path = []

V_deploy = [V_0_deploy]
gamma_deploy_list = [gamma_0_deploy]
h_deploy = []
x_deploy = []
t_deploy = [0]

dt = 0.01
n = 0

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

    ### DEPLOYMENT TRAJECTORY CALCULATIONS ###
    #gamma_deploy_list[-1] < gamma_glide and
    while  t_deploy[-1] < 200:
        n += 1
        rho = gram.data.Density_kgm3[index_min]

        h = V_deploy[-1] * dt * np.sin(gamma_deploy_list[-1])
        x = V_deploy[-1] * dt * np.cos(gamma_deploy_list[-1])

        h_deploy.append(h)
        x_deploy.append(x)

        dVdt = dV_dt(rho, V_deploy[-1], gamma_deploy_list[-1])
        dgammadt = dgamma_dt(rho, V_deploy[-1], gamma_deploy_list[-1])

        V = v_deploy(V_deploy[-1], dVdt, dt)
        gamma = gamma_deploy(gamma_deploy_list[-1], dgammadt, dt)

        V_deploy.append(V)
        gamma_deploy_list.append(gamma)

        t_deploy.append(dt * n)

    plt.subplot(221)
    plt.plot(t_deploy[1:], h_deploy)
    plt.title("Height")
    plt.subplot(222)
    plt.plot(t_deploy[1:], x_deploy)
    plt.title("Horizontal distance")
    plt.subplot(223)
    plt.plot(t_deploy, V_deploy)
    plt.title("Velocity")
    plt.subplot(224)
    plt.plot(t_deploy, gamma_deploy_list)
    plt.title("Gamma")
    plt.show()

    # ### TRAJECTORY CALCULATIONS ###
    # # Find velocity, height, and flight time per density interval
    # for i in points:
    #     V_flight_path.append(velocity(density[i]))
    #     if i < index_max:
    #         h_flight_path.append(height[i])
    #         t_flight_path.append(time_of_flight((height[i] - height[i + 1]) * 1000, density[i]))
    #
    # # Find maximum flight range
    # range = range(np.abs(height[index_max]))
    # total_time = sum(t_flight_path)
    #
    # seconds = total_time
    # seconds_in_day = 60 * 60 * 24
    # seconds_in_hour = 60 * 60
    # seconds_in_minute = 60
    #
    # days = seconds // seconds_in_day
    # hours = (seconds - (days * seconds_in_day)) // seconds_in_hour
    # minutes = (seconds - (days * seconds_in_day) - (hours * seconds_in_hour)) // seconds_in_minute
    #
    # print("The maximum flight range is " + str(range) + " km")
    # print("The time of flight is " + str(days) + " days, " + str(hours) + " hours, and " + str(minutes) + " minutes")
    #
    # # Plot altitude vs time
    # plt.plot(t_flight_path, h_flight_path)
    # plt.xlabel("Time [s]")
    # plt.ylabel("Altitude [km]")
    # plt.title("Altitude vs Time of the Glider")
    # plt.show()
    #
    #
    # ### EIGENVALUE CALCULATIONS ###
    # # Find eigenvalues of the symmetric eigenmotions
    # lambda_short_1 = eigenvalue_short_period_V_constant()
    # lambda_short_2 = eigenvalue_short_period_V_constant_gamma_constant()
    # lambda_phugoid_1 = eigenvalue_phugoid_q_dot_zero_alpha_zero()
    # lambda_phugoid_2 = eigenvalue_phugoid_q_dot_zero_alpha_dot_zero()
    # lambdas_sym = [lambda_short_1, lambda_short_2, lambda_phugoid_1, lambda_phugoid_2]
    #
    # # Find eigenvalues of the asymmetric eigenmotions
    # lambda_aperiodic_roll = eigenvalue_heavily_damped_aperiodic_roll()
    # lambda_dutch_roll_1 = eigenvalue_dutch_roll_phi_zero()
    # lambda_dutch_roll_2 = eigenvalue_dutch_roll_phi_zero_yaw_only()
    # lambda_aperiodic_spiral = eigenvalue_aperiodic_spiral()
    # lambda_dutch_roll_plus_spiral = eigenvalue_dutch_roll_plus_aperiodic_spiral_motion()
    # lambdas_asym = [lambda_aperiodic_roll, lambda_dutch_roll_1, lambda_dutch_roll_2, lambda_aperiodic_spiral,
    #                 lambda_dutch_roll_plus_spiral[0], lambda_dutch_roll_plus_spiral[1],
    #                 lambda_dutch_roll_plus_spiral[2]]
    #
    # # Analysis symmetric eigenvalues
    # for i in np.arange(len(lambdas_sym)):
    #     if lambdas_sym[i].imag == 0:
    #         T_half, tau = symmetric_real_eigenvalues_analysis(lambdas_sym[i].real, V_flight_path[0])
    #         print("The halving time and time constant of " + str(lambdas_sym[i]) + " are " + str(
    #             T_half) + " s and " + str(tau) + " s ")
    #     else:
    #         P, T_half = symmetric_complex_eigenvalues_analysis(lambdas_sym[i], V_flight_path[0])
    #         print("The halving time and period of " + str(lambdas_sym[i]) + " are " + str(T_half) + " s and " + str(
    #             P) + " s ")
    #
    # # Analysis asymmetric eigenvalues
    # for i in np.arange(len(lambdas_asym)):
    #     if lambdas_asym[i].imag == 0:
    #         T_half = asymmetric_real_eigenvalues_analysis(lambdas_asym[i].real, V_flight_path[0])
    #         print("The halving time " + str(lambdas_asym[i]) + " is and " + str(T_half) + " s ")
    #     else:
    #         P, T_half = asymmetric_complex_eigenvalues_analysis(lambdas_asym[i], V_flight_path[0])
    #         print("The halving time and period of " + str(lambdas_asym[i]) + " are " + str(T_half) + " s and " + str(
    #             P) + " s")
