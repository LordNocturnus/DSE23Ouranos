import os.path

from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
import pandas as pd
import math

from tudatpy.kernel.numerical_simulation import propagation_setup

from atmospheric_model.Entry_model import entry_sim


def decelerator_sizing(target_time, totalmass, heatshieldmass, capsule_drag_coefficient, diameter, alt, lat, lon, speed,
                       flight_path_angle, heading_angle, acc=1, steps=5):
    termination_altitude_settings = propagation_setup.propagator.dependent_variable_termination(
        dependent_variable_settings=propagation_setup.dependent_variable.altitude("Capsule", "Uranus"),
        limit_value=0.0,
        use_as_lower_limit=True)

    dep_vars = entry_sim(totalmass, capsule_drag_coefficient, diameter, alt, lat, lon, speed, flight_path_angle,
                         heading_angle, [termination_altitude_settings], acc=1)

    plt.plot(dep_vars[:, 0], dep_vars[:, 10])
    plt.grid()
    plt.show()
    plt.plot(dep_vars[:, 1], dep_vars[:, 10])
    plt.grid()
    plt.show()


if __name__ == "__main__":
    decelerator_sizing(2 * 3600, 500, 250, 1.53, 4.5, 3.02877105e+07, -6.40748300e-02, -1.63500310e+00 + 2 * np.pi,
                       1.93919454e+04, np.deg2rad(-30), -2.35413606e+00)