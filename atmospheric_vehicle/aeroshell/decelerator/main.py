import os.path

from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
import pandas as pd
import math

from tudatpy.kernel.numerical_simulation import propagation_setup

from atmospheric_model.Entry_model import entry_sim
from atmospheric_model.GRAM import GRAM


def decelerator_sizing(target_time, target_pressure, totalmass, heatshieldmass, capsule_drag_coefficient,
                       parachute_c_ds, shock_load_factors, diameter, alt, lat, lon, speed, flight_path_angle,
                       heading_angle, limit_load, acc=1, steps=5):
    termination_settings = propagation_setup.propagator.dependent_variable_termination(
        dependent_variable_settings=propagation_setup.dependent_variable.altitude("Capsule", "Uranus"),
        limit_value=0.0,
        use_as_lower_limit=True)

    ballistic_dep_vars = entry_sim(totalmass, capsule_drag_coefficient, diameter, alt, lat, lon, speed,
                                   flight_path_angle, heading_angle, [termination_settings], acc=acc)

    print("finished first sim")

    gram = GRAM()
    gram.run()
    target_alt = np.asarray(gram.data.Height_km[gram.data.Pressure_Pa >= target_pressure], dtype=float)[0] * 1000

    termination_settings = propagation_setup.propagator.dependent_variable_termination(
        dependent_variable_settings=propagation_setup.dependent_variable.altitude("Capsule", "Uranus"),
        limit_value=target_alt,
        use_as_lower_limit=True)

    parachute_diameter = diameter

    """plt.plot(dep_vars[:, 0], dep_vars[:, 1])
    plt.show()
    plt.plot(dep_vars[:, 0], dep_vars[:, 2])
    plt.show()
    plt.plot(dep_vars[:, 0], dep_vars[:, 3])
    plt.show()
    plt.plot(dep_vars[:, 0], dep_vars[:, 4])
    plt.show()

    plt.plot(dep_vars[:, 0], dep_vars[:, 11])
    plt.show()
    plt.plot(dep_vars[:, 0], dep_vars[:, 12])
    plt.show()#"""

    for _ in range(0, steps):
        dep_vars = entry_sim(totalmass - heatshieldmass, parachute_c_ds, parachute_diameter,
                             ballistic_dep_vars[-1, -1], ballistic_dep_vars[-1, 2], ballistic_dep_vars[-1, 3],
                             ballistic_dep_vars[-1, 4], ballistic_dep_vars[-1, -3], ballistic_dep_vars[-1, -2],
                             [termination_settings], acc=acc)
        parachute_diameter *= target_time / (dep_vars[-1, 0] - 6000)

    load = ballistic_dep_vars[:, -4][ballistic_dep_vars[:, 4] <= 1000] * np.pi * (parachute_diameter / 2) ** 2 *\
           parachute_c_ds / (totalmass - heatshieldmass) * shock_load_factors

    deploy_alt = ballistic_dep_vars[:, 1][ballistic_dep_vars[:, 4] <= 1000][load <= limit_load]

    termination_settings = propagation_setup.propagator.dependent_variable_termination(
        dependent_variable_settings=propagation_setup.dependent_variable.altitude("Capsule", "Uranus"),
        limit_value=deploy_alt[0],
        use_as_lower_limit=True)

    ballistic_dep_vars = entry_sim(totalmass, capsule_drag_coefficient, diameter, alt, lat, lon, speed,
                                   flight_path_angle, heading_angle, [termination_settings], acc=acc)

    termination_settings = propagation_setup.propagator.dependent_variable_termination(
        dependent_variable_settings=propagation_setup.dependent_variable.altitude("Capsule", "Uranus"),
        limit_value=target_alt,
        use_as_lower_limit=True)

    dep_vars = entry_sim(totalmass - heatshieldmass, parachute_c_ds, parachute_diameter,
                         ballistic_dep_vars[-1, -1], ballistic_dep_vars[-1, 2], ballistic_dep_vars[-1, 3],
                         ballistic_dep_vars[-1, 4], ballistic_dep_vars[-1, -3], ballistic_dep_vars[-1, -2],
                         [termination_settings], ballistic_dep_vars[-1, 0], acc=acc)

    area = 2 * np.pi * (parachute_diameter / 2) ** 2
    chute_weight = area * 35 / 1000

    total_weight = chute_weight

    line_load = np.sqrt(np.sum(np.square(dep_vars[0, 5:8]))) * (totalmass - heatshieldmass)
    line_safe_strength = 325  # N with safety factor of 12
    line_weight_per_meter = 0.013 # kg/m
    line_count = np.ceil(line_load / line_safe_strength)

    distance = parachute_diameter * 1.5
    line_length = np.sqrt((parachute_diameter / 2) ** 2 + distance ** 2) + np.pi / 2 * (parachute_diameter / 2) ** 2
    line_weight = line_length * line_count * line_weight_per_meter

    total_weight += line_weight

    print(total_weight)

    return deploy_alt, total_weight


if __name__ == "__main__":
    decelerator_sizing(4 * 3600, 20*10**5, 500, 250, 1.53, 0.4, 2.0, 4.5, 3.02877105e+07,
                       -6.40748300e-02, -1.63500310e+00 + 2 * np.pi, 1.93919454e+04, np.deg2rad(-30), -2.35413606e+00,
                       100, steps=1)
