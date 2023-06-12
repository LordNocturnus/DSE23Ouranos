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
                       heading_angle, acc=1, steps=5):
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
    print(target_alt)

    termination_settings = propagation_setup.propagator.dependent_variable_termination(
        dependent_variable_settings=propagation_setup.dependent_variable.altitude("Capsule", "Uranus"),
        limit_value=target_alt,
        use_as_lower_limit=True)

    parachute_diameter = np.full_like(parachute_c_ds, diameter, dtype=float)

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
        dep_vars = entry_sim(totalmass - heatshieldmass, parachute_c_ds[-1], parachute_diameter[-1],
                             ballistic_dep_vars[-1, -1], ballistic_dep_vars[-1, 2], ballistic_dep_vars[-1, 3],
                             ballistic_dep_vars[-1, 4], ballistic_dep_vars[-1, -3], ballistic_dep_vars[-1, -2],
                             [termination_settings], acc=acc)
        print(dep_vars[-1, 0] - 6000, target_time)#"""
        parachute_diameter[-1] *= target_time / (dep_vars[-1, 0] - 6000)
        print(parachute_diameter)

    gram.time = dep_vars[:, 0]
    gram.altitudes = dep_vars[:, 1] / 1000
    gram.lat = np.rad2deg(dep_vars[:, 2])
    gram.long = np.rad2deg((dep_vars[:, 3] + 2 * np.pi) % 2 * np.pi)
    gram.run()

    plt.plot(dep_vars[:, 9], dep_vars[:, 1], label="exp")
    plt.plot(gram.data.Density_kgm3, dep_vars[:, 1], label="gram")
    plt.grid()
    plt.legend()
    plt.show()


if __name__ == "__main__":
    decelerator_sizing(2 * 3600, 20*10**5, 500, 250, 1.53, np.asarray([0.5]), np.asarray([1.3]), 4.5, 3.02877105e+07,
                       -6.40748300e-02, -1.63500310e+00 + 2 * np.pi, 1.93919454e+04, np.deg2rad(-30), -2.35413606e+00,
                       steps=1)