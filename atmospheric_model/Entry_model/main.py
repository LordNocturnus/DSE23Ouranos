import os.path

from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
import pandas as pd
import math

# Load tudatpy modules
from tudatpy.kernel.interface import spice
from tudatpy.kernel import numerical_simulation
from tudatpy.kernel.numerical_simulation import environment_setup, environment, propagation_setup, propagation
from tudatpy.kernel.astro import element_conversion
from tudatpy.kernel import constants
from tudatpy.util import result2array

from atmospheric_model import GRAM

_path = os.path.dirname(__file__)

spice.load_standard_kernels()
spice.load_kernel(_path + '/../GRAM/GRAM_Suite_1_5/SPICE/spk/satellites/ura116xl.bsp')


# spice.load_kernel(_path+'/Gravity.tpc')


def entry_sim(mass, drag_coefficient, diameter, alt, lat, lon, speed, flight_path_angle, heading_angle,
              termination_settings, simulation_start_epoch=6000.0, max_simulation_time=1 * constants.JULIAN_DAY,
              acc=1):

    # define bodies in simulation
    bodies_to_create = ["Uranus"]

    # create body settings dictionary
    global_frame_origin = "Uranus"
    global_frame_orientation = "J2000"
    body_settings = environment_setup.get_default_body_settings(
        bodies_to_create, global_frame_origin, global_frame_orientation)

    # Modify the gravity field of Uranus to match GRAM
    body_settings.get("Uranus").gravity_field_settings = environment_setup.gravity_field.central_spice()

    print("Setup gravity")

    # Modify the atmosphere to match GRAM
    gram = GRAM.GRAM()
    gram.altitudes = np.arange(7000, -290 - acc, -acc)
    gram.lat = np.full_like(gram.altitudes, np.rad2deg(lat))
    gram.long = np.full_like(gram.altitudes, np.rad2deg((lon + 2 * np.pi) % (2 * np.pi)))
    gram.time = np.zeros_like(gram.altitudes)
    gram.run()

    atmos = {"altitude": np.asarray(gram.data.Height_km) * 1000,
             "density": np.asarray(gram.data.Density_kgm3),
             "pressure": np.asarray(gram.data.Pressure_Pa),
             "temperature": np.asarray(gram.data.Temperature_K),
             "gas_constant": np.asarray(gram.data.SpecificGasConstant_JkgK),
             "specific_heat_ratio": np.asarray(gram.data.SpecificHeatRatio),
             "molar_mass": np.asarray(gram.data.AverageMolecularWeight)}

    atmos = pd.DataFrame(atmos, dtype=float)
    # atmos.to_csv(_path + "/atmos_data.csv", index=False)

    def density(h):
        return np.asarray(gram.data.Density_kgm3[gram.data.Height_km <= h / 1000])[0]

    constant_temperature = np.asarray(gram.data.Density_kgm3[gram.data.Height_km <= 0.0])[0]
    specific_gas_constant = np.asarray(gram.data.Density_kgm3[gram.data.Height_km <= 0.0])[0]
    ratio_of_specific_heats = np.asarray(gram.data.Density_kgm3[gram.data.Height_km <= 0.0])[0]

    body_settings.get("Uranus").atmosphere_settings = environment_setup.atmosphere.custom_constant_temperature(
        density,
        constant_temperature,
        specific_gas_constant,
        ratio_of_specific_heats)

    """body_settings.get("Uranus").atmosphere_settings = environment_setup.atmosphere.exponential(
        40.061 * 1000, 3.788e-01, 76.4, 3456.0, 1.632)

    body_settings.get("Uranus").atmosphere_settings = environment_setup.atmosphere.tabulated(
        _path + "/atmos_data.csv", [environment_setup.atmosphere.AtmosphereDependentVariables.tabulated_density,
                                    environment_setup.atmosphere.AtmosphereDependentVariables.tabulated_pressure,
                                    environment_setup.atmosphere.AtmosphereDependentVariables.tabulated_temperature,
                                    environment_setup.atmosphere.AtmosphereDependentVariables.tabulated_gas_constant,
                                    environment_setup.atmosphere.AtmosphereDependentVariables.tabulated_specific_heat_ratio,
                                    environment_setup.atmosphere.AtmosphereDependentVariables.tabulated_molar_mass])"""

    # Create system of bodies from the body settings
    bodies = environment_setup.create_system_of_bodies(body_settings)

    print("Setup atmosphere")

    bodies.create_empty_body("Capsule")

    # Set mass of vehicle
    bodies.get("Capsule").mass = mass

    # Create aerodynamic coefficient interface settings, and add to vehicle
    aero_coefficient_settings = environment_setup.aerodynamic_coefficients.constant(
        reference_area=np.pi * (diameter / 2) ** 2,
        constant_force_coefficient=[drag_coefficient, 0.0, 0.0])
    environment_setup.add_aerodynamic_coefficient_interface(bodies, "Capsule", aero_coefficient_settings)

    print("Setup Capsule")

    # Define constant angles
    angle_of_attack = np.deg2rad(0.0)
    bank_angle = np.deg2rad(0.0)

    # Define angle function (required for input to rotation settings)

    def angle_function(time):
        return np.asarray([0.0, 0.0, 0.0], dtype=float)

    # Create settings for rotation model
    rotation_model_settings = environment_setup.rotation_model.aerodynamic_angle_based(central_body="Uranus",
                                                                                       base_frame="",
                                                                                       target_frame="Capsule_Fixed",
                                                                                       angle_funcion=angle_function)
    environment_setup.add_rotation_model(bodies, 'Capsule', rotation_model_settings)

    print("Setup rotation")

    # Define bodies that are propagated
    bodies_to_propagate = ["Capsule"]

    # Define central bodies of propagation
    central_bodies = ["Uranus"]

    # Define the accelerations acting on the Capsule (Uranus as a Point Mass, and Uranus's atmosphere)
    accelerations_settings_capsule = dict(Uranus=[propagation_setup.acceleration.point_mass_gravity(),
                                                  propagation_setup.acceleration.aerodynamic()])

    acceleration_settings = {"Capsule": accelerations_settings_capsule}

    acceleration_models = propagation_setup.create_acceleration_models(bodies, acceleration_settings,
                                                                       bodies_to_propagate, central_bodies)

    print("Setup acceleration")

    # Convert the initial state
    initial_uranus_fixed_state = element_conversion.spherical_to_cartesian_elementwise(radial_distance=alt,
                                                                                       latitude=lat,
                                                                                       longitude=lon, speed=speed,
                                                                                       flight_path_angle=flight_path_angle,
                                                                                       heading_angle=heading_angle)

    # Convert the state from the Earth-fixed frame to the inertial frame
    uranus_rotation_model = bodies.get_body("Uranus").rotation_model
    initial_state = environment.transform_to_inertial_orientation(initial_uranus_fixed_state,
                                                                  simulation_start_epoch, uranus_rotation_model)

    print("Setup initial state")

    # Define the list of dependent variables to save during the propagation
    dependent_variables_to_save = [propagation_setup.dependent_variable.altitude("Capsule", "Uranus"),
                                   propagation_setup.dependent_variable.latitude("Capsule", "Uranus"),
                                   propagation_setup.dependent_variable.longitude("Capsule", "Uranus"),
                                   propagation_setup.dependent_variable.airspeed("Capsule", "Uranus"),
                                   propagation_setup.dependent_variable.total_acceleration("Capsule"),
                                   propagation_setup.dependent_variable.mach_number("Capsule", "Uranus"),
                                   propagation_setup.dependent_variable.density("Capsule", "Uranus"),
                                   propagation_setup.dependent_variable.dynamic_pressure("Capsule"),
                                   propagation_setup.dependent_variable.flight_path_angle("Capsule", "Uranus"),
                                   propagation_setup.dependent_variable.heading_angle("Capsule", "Uranus"),
                                   propagation_setup.dependent_variable.relative_distance("Capsule", "Uranus"),]

    print("Setup dependent vars")

    # Define a termination condition to stop after a given time (to avoid an endless skipping re-entry)
    termination_time_settings = propagation_setup.propagator.time_termination(
        simulation_start_epoch + max_simulation_time)
    # Combine the termination settings to stop when one of them is fulfilled
    combined_termination_settings = propagation_setup.propagator.hybrid_termination(
        [termination_time_settings] + termination_settings, fulfill_single_condition=True)

    print("Setup termination")

    # Create numerical integrator settings
    integrator_settings = propagation_setup.integrator.runge_kutta_4(acc)

    # Create the propagation settings
    propagator_settings = propagation_setup.propagator.translational(central_bodies, acceleration_models,
                                                                     bodies_to_propagate, initial_state,
                                                                     simulation_start_epoch, integrator_settings,
                                                                     combined_termination_settings,
                                                                     output_variables=dependent_variables_to_save)

    print("setup propagator")

    # Create the simulation objects and propagate the dynamics
    dynamics_simulator = numerical_simulation.create_dynamics_simulator(bodies, propagator_settings)

    # Extract the resulting simulation dependent variables
    dependent_variables = dynamics_simulator.dependent_variable_history
    # Convert the dependent variables from a dictionary to a numpy array
    return result2array(dependent_variables)


if __name__ == "__main__":
    angle = -45

    termination_altitude_settings = propagation_setup.propagator.dependent_variable_termination(
        dependent_variable_settings=propagation_setup.dependent_variable.altitude("Capsule", "Uranus"),
        limit_value=-150000,
        use_as_lower_limit=True)

    dependent_variables_array = entry_sim(500, 1.53, 4.5, 3.02877105e+07, -6.40748300e-02, -1.63500310e+00 + 2 * np.pi,
                                          1.93919454e+04, np.deg2rad(angle), -2.35413606e+00,
                                          [termination_altitude_settings])

    plt.plot(dependent_variables_array[:, 0], dependent_variables_array[:, 1])
    plt.show()

    gram = GRAM.GRAM()
    gram.time = dependent_variables_array[:, 0]
    gram.altitudes = dependent_variables_array[:, 1] / 1000
    gram.lat = np.rad2deg(dependent_variables_array[:, 2])
    gram.long = np.rad2deg((dependent_variables_array[:, 2] + 2 * np.pi) % (2 * np.pi))
    gram.run()

    plt.plot(np.log(dependent_variables_array[:, -2]), dependent_variables_array[:, 1], label="exp")
    plt.plot(np.log(gram.data.Density_kgm3), dependent_variables_array[:, 1], label="gram")
    plt.grid()
    plt.legend()
    plt.show()

    print("finished")
