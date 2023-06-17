import os.path

from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
import pandas as pd
import math

from tudatpy.kernel.numerical_simulation import environment_setup, environment, propagation_setup

import atmospheric_model.Entry_model
from atmospheric_model.GRAM import GRAM


def decent_ballistic(totalmass, diameter, alt, lat, lon, speed, flight_path_angle, heading_angle, acc=1):
    drag = atmospheric_model.Entry_model.CapsuleDrag(diameter, 1.125, np.deg2rad(20), lat, lon, acc)

    aero_coefficient_setting = environment_setup.aerodynamic_coefficients.custom_aerodynamic_force_coefficients(
        force_coefficient_function=drag.drag_coefficient,
        reference_area=np.pi * (drag.diameter / 2) ** 2,
        independent_variable_names=[environment.AerodynamicCoefficientsIndependentVariables.mach_number_dependent,
                                    environment.AerodynamicCoefficientsIndependentVariables.altitude_dependent])

    termination_settings = propagation_setup.propagator.dependent_variable_termination(
        dependent_variable_settings=propagation_setup.dependent_variable.mach_number("Capsule", "Uranus"),
        limit_value=2.7,
        use_as_lower_limit=True)

    ballistic_dep_vars = atmospheric_model.Entry_model.entry_sim(totalmass, aero_coefficient_setting, alt, lat, lon,
                                                                 speed, flight_path_angle, heading_angle,
                                                                 [termination_settings], acc=acc)
    return ballistic_dep_vars


def decelerator_decent(target_alt, mass, diameter, c_d, mach, alt, lat, lon, speed,
                       flight_path_angle, heading_angle, acc=1):
    termination_alt_settings = propagation_setup.propagator.dependent_variable_termination(
        dependent_variable_settings=propagation_setup.dependent_variable.altitude("Capsule", "Uranus"),
        limit_value=target_alt,
        use_as_lower_limit=True)  # """

    termination_mach_settings = propagation_setup.propagator.dependent_variable_termination(
        dependent_variable_settings=propagation_setup.dependent_variable.mach_number("Capsule", "Uranus"),
        limit_value=mach,
        use_as_lower_limit=True)#"""

    aero_coefficient_setting = environment_setup.aerodynamic_coefficients.constant(
        np.pi * (diameter / 2) ** 2,
        constant_force_coefficient=[c_d, 0.0, 0.0],
        are_coefficients_in_aerodynamic_frame=True,
        are_coefficients_in_negative_axis_direction=True)#"""

    dep_vars = atmospheric_model.Entry_model.entry_sim(mass, aero_coefficient_setting, alt, lat,
                                                       lon, speed, flight_path_angle, heading_angle,
                                                       [termination_alt_settings, termination_mach_settings], acc=acc)
    return dep_vars


def decelerator_sizing(target_time, target_pressure, totalmass, heatshieldmass, parachute_c_d, shockload, diameter, alt,
                       lat, lon, speed, flight_path_angle, heading_angle, acc=1, steps=5):

    ballistic_dep_vars = decent_ballistic(totalmass, diameter, alt, lat, lon, speed, flight_path_angle, heading_angle,
                                          acc)

    gram = GRAM()
    gram.run()
    target_alt = np.asarray(gram.data.Height_km[gram.data.Pressure_Pa >= target_pressure], dtype=float)[0] * 1000

    parachute_diameter = diameter

    for _ in range(0, steps):

        dep_vars_sonic = decelerator_decent(target_alt, totalmass - heatshieldmass, parachute_diameter,
                                            parachute_c_d[0], 0.9, ballistic_dep_vars[-1, -1],
                                            ballistic_dep_vars[-1, 2], ballistic_dep_vars[-1, 3],
                                            ballistic_dep_vars[-1, 4], ballistic_dep_vars[-1, -3],
                                            ballistic_dep_vars[-1, -2], acc)

        dep_vars = decelerator_decent(target_alt, totalmass - heatshieldmass, parachute_diameter, parachute_c_d[1], 0.0,
                                      dep_vars_sonic[-1, -1], dep_vars_sonic[-1, 2], dep_vars_sonic[-1, 3],
                                      dep_vars_sonic[-1, 4], dep_vars_sonic[-1, -3], dep_vars_sonic[-1, -2], acc)

        parachute_diameter *= target_time / (dep_vars[-1, 0] - dep_vars[dep_vars[:, 1] <= 0.0][0, 0])

    load = ballistic_dep_vars[-1, -4] * np.pi * (parachute_diameter / 2) ** 2 * parachute_c_d[0] / \
           (totalmass - heatshieldmass) * shockload

    area = 2 * np.pi * (parachute_diameter / 2) ** 2
    chute_weight = area * 35 / 1000
    chute_cost = area * 15 * 0.92

    total_weight = chute_weight
    total_cost = chute_cost

    line_safe_strength = 1300 * 4  # N with safety factor of 3
    line_weight_per_meter = 0.013  # kg/m
    line_cost_per_meter = 0.7 #
    line_count = np.ceil(load * (totalmass - heatshieldmass) / line_safe_strength)

    distance = parachute_diameter * 1.5
    line_length = np.sqrt((parachute_diameter / 2) ** 2 + distance ** 2) + np.pi / 2 * (parachute_diameter / 2) ** 2
    line_weight = line_length * line_count * line_weight_per_meter
    line_cost = line_length * line_count * line_cost_per_meter

    total_weight += line_weight
    total_cost += line_cost

    return ballistic_dep_vars[-1, 1] / 1000, total_weight, total_cost, load, dep_vars[dep_vars[:, 1] <= 50000][0, 4]


if __name__ == "__main__":
    print(decelerator_sizing(4 * 3600, 20*10**5, 162.1, 0.0, [0.45, 0.4], 2.0, 3, 3.03327727e+07,
                             5.45941114e-01, -2.33346601e-02, 2.65992642e+04, -5.91036848e-01, -2.96367147e+00, acc=1,
                             steps=1))
