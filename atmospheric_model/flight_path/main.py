import numpy as np

from GlideTrajectory import *
from Eigenvalues import *
from DeploymentTrajectory import *
from atmospheric_model import GRAM
from matplotlib import pyplot as plt

# Lists for Gliding Flight
V_flight_path = []
h_flight_path = []
t_flight_path = []

# Lists for Deployement Trajectory
V_deploy = [V_0_deploy]
gamma_deploy_list = [gamma_0_deploy]
h_deploy = [0]
x_deploy = [0]
s_deploy = [0]
t_deploy = [0]

# Time Interval and Number of Loops
dt = 0.01
n = 0

if __name__ == "__main__":
    # Setting up GRAM
    gram = GRAM.GRAM()
    # gram.runs = 10
    # gram.densityperturbation = 1.0  # Between 0.0 and 2.0, 1.0 = 3 sigma
    gram.run()

    # Finding indexes for the relevant pressures
    index_1 = gram.data.Pressure_Pa[gram.data.Pressure_Pa > 10 ** 4].index
    index_2 = gram.data.Pressure_Pa[gram.data.Pressure_Pa < p_max].index

    # Find common indexes
    index = np.sort(list(set(index_1) & set(index_2)))

    # Find density and heights from indexes
    density = gram.data.Density_kgm3[index]
    height = gram.data.Height_km[index]

    # Find maximum and minimum index to convert to list
    index_min = min(index)
    index_max = max(index)
    points = np.arange(index_min, index_max + 1, 1)

    ### DEPLOYMENT TRAJECTORY CALCULATIONS ###
    # gamma_deploy_list[-1] < gamma_glide and
    while t_deploy[-1] < 100:
        n += 1
        rho = gram.data.Density_kgm3[index_min]

        h = h_deploy[-1] + V_deploy[-1] * dt * np.sin(gamma_deploy_list[-1])
        x = x_deploy[-1] + V_deploy[-1] * dt * np.cos(gamma_deploy_list[-1])
        s = np.sqrt(h * h + x * x)

        h_deploy.append(h)
        x_deploy.append(x)
        s_deploy.append(s)

        dVdt = dV_dt(rho, V_deploy[-1], gamma_deploy_list[-1])
        dgammadt = dgamma_dt(rho, V_deploy[-1], gamma_deploy_list[-1])

        V = v_deploy(V_deploy[-1], dVdt, dt)
        gamma = gamma_deploy(gamma_deploy_list[-1], dgammadt, dt)

        V_deploy.append(V)
        gamma_deploy_list.append(gamma)

        t_deploy.append(dt * n)

    plt.subplot(221)
    plt.plot(t_deploy, h_deploy)
    plt.title("Height")
    plt.subplot(222)
    plt.plot(t_deploy, x_deploy)
    plt.title("Horizontal distance")
    plt.subplot(223)
    plt.plot(t_deploy, V_deploy)
    plt.title("Velocity")
    plt.subplot(224)
    plt.plot(t_deploy, gamma_deploy_list)
    plt.title("Gamma")
    plt.show()

    ### TRAJECTORY CALCULATIONS ###
    # Find velocity, height, and flight time per density interval
    for i in points:
        V_flight_path.append(velocity(density[i]))
        if i < index_max:
            h_flight_path.append(height[i])
            t_flight_path.append(time_of_flight((height[i] - height[i + 1]) * 1000, density[i]))

    # Find maximum flight range
    range = range(np.abs(height[index_max]))
    total_time = sum(t_flight_path)

    seconds = total_time
    seconds_in_day = 60 * 60 * 24
    seconds_in_hour = 60 * 60
    seconds_in_minute = 60

    days = seconds // seconds_in_day
    hours = (seconds - (days * seconds_in_day)) // seconds_in_hour
    minutes = (seconds - (days * seconds_in_day) - (hours * seconds_in_hour)) // seconds_in_minute

    print("The maximum flight range is " + str(range) + " km")
    print("The time of flight is " + str(days) + " days, " + str(hours) + " hours, and " + str(minutes) + " minutes")

    # Plot altitude vs time
    plt.plot(t_flight_path, h_flight_path)
    plt.xlabel("Time [s]")
    plt.ylabel("Altitude [km]")
    plt.title("Altitude vs Time of the Glider")
    plt.show()

    ### EIGENVALUE CALCULATIONS ###
    # Find Eigenvalues of the Symmetric Eigenmotions
    lambda_short_1 = eigenvalue_short_period_V_constant()
    lambda_short_2 = eigenvalue_short_period_V_constant_gamma_constant()
    lambda_phugoid_1 = eigenvalue_phugoid_q_dot_zero_alpha_zero()
    lambda_phugoid_2 = eigenvalue_phugoid_q_dot_zero_alpha_dot_zero()
    lambdas_sym = [lambda_short_1, lambda_short_2, lambda_phugoid_1, lambda_phugoid_2]

    # Find Eigenvalues of the Asymmetric Eigenmotions
    lambda_aperiodic_roll = eigenvalue_heavily_damped_aperiodic_roll()
    lambda_dutch_roll_1 = eigenvalue_dutch_roll_phi_zero()
    lambda_dutch_roll_2 = eigenvalue_dutch_roll_phi_zero_yaw_only()
    lambda_aperiodic_spiral = eigenvalue_aperiodic_spiral()
    lambda_dutch_roll_plus_spiral = eigenvalue_dutch_roll_plus_aperiodic_spiral_motion()
    lambdas_asym = [lambda_aperiodic_roll, lambda_dutch_roll_1, lambda_dutch_roll_2, lambda_aperiodic_spiral,
                    lambda_dutch_roll_plus_spiral[0], lambda_dutch_roll_plus_spiral[1],
                    lambda_dutch_roll_plus_spiral[2]]

    # Analysis Symmetric Eigenvalues
    for i in np.arange(len(lambdas_sym)):
        if lambdas_sym[i].imag == 0:
            T_half, tau = symmetric_real_eigenvalues_analysis(lambdas_sym[i].real, V_flight_path[0])
            print("The halving time and time constant of " + str(lambdas_sym[i]) + " are " + str(
                T_half) + " s and " + str(tau) + " s ")
        else:
            P, T_half = symmetric_complex_eigenvalues_analysis(lambdas_sym[i], V_flight_path[0])
            print("The halving time and period of " + str(lambdas_sym[i]) + " are " + str(T_half) + " s and " + str(
                P) + " s ")

    # Analysis Asymmetric Eigenvalues
    for i in np.arange(len(lambdas_asym)):
        if lambdas_asym[i].imag == 0:
            T_half = asymmetric_real_eigenvalues_analysis(lambdas_asym[i].real, V_flight_path[0])
            print("The halving time " + str(lambdas_asym[i]) + " is and " + str(T_half) + " s ")
        else:
            P, T_half = asymmetric_complex_eigenvalues_analysis(lambdas_asym[i], V_flight_path[0])
            print("The halving time and period of " + str(lambdas_asym[i]) + " are " + str(T_half) + " s and " + str(
                P) + " s")
