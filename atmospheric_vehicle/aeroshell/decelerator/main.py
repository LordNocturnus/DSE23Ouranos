import os.path

from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
import pandas as pd
import math

from tudatpy.kernel.numerical_simulation import environment_setup, environment, propagation_setup

import atmospheric_model.Entry_model
from atmospheric_model.GRAM import GRAM


def decelerator_sizing(target_time, target_pressure, totalmass, heatshieldmass, parachute_c_d, shockload, diameter, alt,
                       lat, lon, speed, flight_path_angle, heading_angle, acc=1, steps=5):

    drag = atmospheric_model.Entry_model.CapsuleDrag(diameter, 1.125, np.deg2rad(20), lat, lon, acc)

    aero_coefficient_setting = environment_setup.aerodynamic_coefficients.custom_aerodynamic_force_coefficients(
        force_coefficient_function=drag.drag_coefficient,
        reference_area=np.pi * (drag.diameter / 2) ** 2,
        independent_variable_names=[environment.AerodynamicCoefficientsIndependentVariables.mach_number_dependent,
                                    environment.AerodynamicCoefficientsIndependentVariables.altitude_dependent])

    termination_settings = propagation_setup.propagator.dependent_variable_termination(
        dependent_variable_settings=propagation_setup.dependent_variable.mach_number("Capsule", "Uranus"),
        limit_value=2.0,
        use_as_lower_limit=True)

    ballistic_dep_vars = atmospheric_model.Entry_model.entry_sim(totalmass, aero_coefficient_setting, alt, lat, lon,
                                                                 speed, flight_path_angle, heading_angle,
                                                                 [termination_settings], acc=acc)

    print("finished first sim")
    print(ballistic_dep_vars[-1])
    print(ballistic_dep_vars[-2])

    gram = GRAM()
    gram.run()
    target_alt = np.asarray(gram.data.Height_km[gram.data.Pressure_Pa >= target_pressure], dtype=float)[0] * 1000
    print(target_alt)

    termination_altitude_settings = propagation_setup.propagator.dependent_variable_termination(
        dependent_variable_settings=propagation_setup.dependent_variable.altitude("Capsule", "Uranus"),
        limit_value=target_alt,
        use_as_lower_limit=True)

    termination_mach_settings = propagation_setup.propagator.dependent_variable_termination(
        dependent_variable_settings=propagation_setup.dependent_variable.mach_number("Capsule", "Uranus"),
        limit_value=0.9,
        use_as_lower_limit=True)

    parachute_diameter = diameter

    for _ in range(0, steps):

        def drag(var):
            print(var[1])
            if var[0] >= 1.0:
                return 0.45
            return 0.4

        """aero_coefficient_setting = environment_setup.aerodynamic_coefficients.custom_aerodynamic_force_coefficients(
            force_coefficient_function=drag,
            reference_area=np.pi * (parachute_diameter / 2) ** 2,
            independent_variable_names=[environment.AerodynamicCoefficientsIndependentVariables.mach_number_dependent,
                                        environment.AerodynamicCoefficientsIndependentVariables.altitude_dependent])"""

        aero_coefficient_setting = environment_setup.aerodynamic_coefficients.constant(
            np.pi * (parachute_diameter / 2) ** 2,
            constant_force_coefficient=[0.45, 0.0, 0.0],
            are_coefficients_in_aerodynamic_frame=True,
            are_coefficients_in_negative_axis_direction=True
        )

        dep_vars = atmospheric_model.Entry_model.entry_sim(totalmass - heatshieldmass, aero_coefficient_setting,
                                                           ballistic_dep_vars[-2, -1], ballistic_dep_vars[-2, 2],
                                                           ballistic_dep_vars[-2, 3], ballistic_dep_vars[-2, 4],
                                                           ballistic_dep_vars[-2, -3], ballistic_dep_vars[-2, -2],
                                                           [termination_altitude_settings,
                                                            termination_mach_settings], acc=acc)

        print("finished second sim")

        parachute_diameter *= target_time / (dep_vars[-1, 0] - dep_vars[dep_vars[:, 1] <= 0.0][0, 0])
        print(parachute_diameter)

    """load = ballistic_dep_vars[:, -4][ballistic_dep_vars[:, 4] <= 1000] * np.pi * (parachute_diameter / 2) ** 2 *\
           parachute_c_ds / (totalmass - heatshieldmass) * shock_load

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
    chute_cost = area * 15 * 0.92 # $

    total_weight = chute_weight

    line_load = np.sqrt(np.sum(np.square(dep_vars[0, 5:8]))) * (totalmass - heatshieldmass)
    line_safe_strength = 1300 * 4  # N with safety factor of 3
    line_weight_per_meter = 0.013 # kg/m
    line_count = np.ceil(line_load / line_safe_strength)

    distance = parachute_diameter * 1.5
    line_length = np.sqrt((parachute_diameter / 2) ** 2 + distance ** 2) + np.pi / 2 * (parachute_diameter / 2) ** 2
    line_weight = line_length * line_count * line_weight_per_meter

    total_weight += line_weight

    return deploy_alt[0] / 1000, total_weight, max(max(np.sqrt(np.sum(np.square(ballistic_dep_vars[:, 5:8]), axis=1))),
                                                   max(np.sqrt(np.sum(np.square(dep_vars[:, 5:8]), axis=1))))#"""


if __name__ == "__main__":
    print(decelerator_sizing(2 * 3600, 20*10**5, 603, 255, [[1.0, 0.45], [0.0, 0.4]], 2.0, 4.5, 3.03327727e+07,
                             5.45941114e-01, -2.33346601e-02, 2.65992642e+04, -5.91036848e-01, -2.96367147e+00, acc=10,
                             steps=1))
