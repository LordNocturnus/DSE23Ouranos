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

angle = np.deg2rad(20)
diameter = 4.5

# spice.load_kernel(_path+'/Gravity.tpc')


def entry_sim(mass, aero_coefficient_settings, alt, lat, lon, speed, flight_path_angle, heading_angle,
              termination_settings, simulation_start_epoch=6000.0, max_simulation_time=1 / 2 * constants.JULIAN_DAY,
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
    gram.altitudes = np.arange(-290, 7000 + acc, acc)
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
    atmos.to_csv(_path + "/atmos_data.csv", index=False, header=False, sep=" ")

    body_settings.get("Uranus").atmosphere_settings = environment_setup.atmosphere.tabulated(
        _path + "/atmos_data.csv", [environment_setup.atmosphere.AtmosphereDependentVariables.tabulated_density,
                                    environment_setup.atmosphere.AtmosphereDependentVariables.tabulated_pressure,
                                    environment_setup.atmosphere.AtmosphereDependentVariables.tabulated_temperature,
                                    environment_setup.atmosphere.AtmosphereDependentVariables.tabulated_gas_constant,
                                    environment_setup.atmosphere.AtmosphereDependentVariables.tabulated_specific_heat_ratio,
                                    environment_setup.atmosphere.AtmosphereDependentVariables.tabulated_molar_mass])#"""

    # Create system of bodies from the body settings
    bodies = environment_setup.create_system_of_bodies(body_settings)

    print("Setup atmosphere")

    bodies.create_empty_body("Capsule")

    # Set mass of vehicle
    bodies.get("Capsule").mass = mass

    # Create aerodynamic coefficient interface settings, and add to vehicle

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
    earth_rotation_model = bodies.get_body("Uranus").rotation_model
    initial_state = environment.transform_to_inertial_orientation(initial_uranus_fixed_state,
                                                                  simulation_start_epoch, earth_rotation_model)

    print("Setup initial state")

    # Define the list of dependent variables to save during the propagation
    dependent_variables_to_save = [propagation_setup.dependent_variable.altitude("Capsule", "Uranus"),
                                   propagation_setup.dependent_variable.latitude("Capsule", "Uranus"),
                                   propagation_setup.dependent_variable.longitude("Capsule", "Uranus"),
                                   propagation_setup.dependent_variable.airspeed("Capsule", "Uranus"),
                                   propagation_setup.dependent_variable.total_acceleration_norm("Capsule"),
                                   propagation_setup.dependent_variable.mach_number("Capsule", "Uranus"),
                                   propagation_setup.dependent_variable.density("Capsule", "Uranus"),
                                   propagation_setup.dependent_variable.dynamic_pressure("Capsule"),]

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


class CapsuleDrag:
    def __init__(self, dia, r1, angle, lat, lon, acc=1):
        gram = GRAM.GRAM()
        gram.altitudes = np.arange(7000, -290 - acc, -acc)
        gram.lat = np.full_like(gram.altitudes, np.rad2deg(lat))
        gram.long = np.full_like(gram.altitudes, np.rad2deg((lon + 2 * np.pi) % (2 * np.pi)))
        gram.time = np.zeros_like(gram.altitudes)
        gram.run()
        self.SpecificHeatRatio = sp.interpolate.interp1d(gram.data.Height_km * 1000, gram.data.SpecificHeatRatio)
        self.Pressure_Pa = sp.interpolate.interp1d(gram.data.Height_km * 1000, gram.data.Pressure_Pa)
        self.diameter = dia
        self.area = np.pi * (dia / 2) ** 2
        self.r1 = r1
        self.angle = angle

    def drag_coefficient(self, var): # Mach, specific heat ratio, freestream pressure
        self.mach = var[0]
        self.gamma = self.SpecificHeatRatio(var[1])
        self.p_inf = self.Pressure_Pa(var[1])
        drag = sp.integrate.quad(lambda y: self.p_bar(y) * 2 * np.pi * y, 0.0, self.diameter / 2)[0]
        return [drag / (1 / 2 * self.gamma * self.mach ** 2 * self.area) + 1 / (1 / 2 * self.gamma * self.mach ** 2),
                0.0, 0.0]

    def beta(self, y):
        if y <= np.sin(self.angle) * self.r1:
            return np.sin(y / self.r1)
        else:
            return self.angle

    def p_star_bar(self): # Mach, specific heat ratio, freestream pressure
        return (2 / (self.gamma + 1)) ** (self.gamma / (self.gamma - 1))

    def p_0_stag(self): # Mach, specific heat ratio, freestream pressure
        upper = ((self.gamma + 1) / 2 * self.mach**2) ** (self.gamma / (self.gamma - 1))
        lower = (2*self.gamma / (self.gamma + 1) * self.mach**2 - (self.gamma - 1) /
                 (self.gamma + 1)) ** (1 / (self.gamma - 1))
        return upper / lower

    def p_inf_bar(self):
        return 1 / (1 + self.p_0_stag())

    def s(self, y):
        if y <= np.sin(self.angle) * self.r1:
            return np.arcsin(y / self.r1) * self.r1
        else:
            return self.angle * self.r1 + (y - self.r1 * np.sin(self.angle) / np.cos(self.angle))

    def s_star(self):
        return self.angle * self.r1 + (self.diameter / 2 - self.r1 * np.sin(self.angle) / np.cos(self.angle))

    def p_fd_bar(self, y):
        return 1 - np.exp(-self.Lambda(y)) * (1 - self.p_star_bar()) + 1 / 16 * ((self.s(y) / self.s_star()) ** 2 -
                                                                                 np.exp(-self.Lambda(y)))

    def Lambda(self, y):
        return 5 * np.sqrt(np.log(self.s_star() / self.s(y)))

    def R_N(self, y):
        if y <= np.sin(self.angle) * self.r1:
            return self.r1
        else:
            return y / np.sin(self.angle)

    def R_max(self):
        return self.diameter / (2 * np.sin(self.angle))

    def p_bar(self, y):
        s_ratio = self.s(y) / self.s_star()
        p1 = self.p_inf_bar()
        p2 = (1 - self.p_inf_bar()) * np.cos(self.beta(y)) ** 2
        p3 = (1 - self.p_fd_bar(y)) * ((np.cos(self.beta(y)) ** 2 - self.p_star_bar()) / (1 - self.p_star_bar()))
        p4 = (1 - self.R_N(y) / self.R_max()) * (np.sin(self.beta(y)) ** 2 * (1 - s_ratio) + 1/2 * s_ratio * (
                self.p_fd_bar(y) - 1 + s_ratio * np.sin(self.beta(y)) ** 2 + p3))
        return (p1 + p2 - p3 + p4) * self.p_0_stag()


if __name__ == "__main__":
    drag = CapsuleDrag(4.5, 1.125, np.deg2rad(20), 5.45941114e-01, -2.33346601e-02)
    
    aero_coefficient_setting = environment_setup.aerodynamic_coefficients.custom_aerodynamic_force_coefficients(
        force_coefficient_function=drag.drag_coefficient,
        reference_area=np.pi * (drag.diameter / 2) ** 2,
        independent_variable_names=[environment.AerodynamicCoefficientsIndependentVariables.mach_number_dependent,
                                    environment.AerodynamicCoefficientsIndependentVariables.altitude_dependent])

    termination_altitude_setting = propagation_setup.propagator.dependent_variable_termination(
        dependent_variable_settings=propagation_setup.dependent_variable.altitude("Capsule", "Uranus"),
        limit_value=0.0,
        use_as_lower_limit=True)

    termination_mach_setting = propagation_setup.propagator.dependent_variable_termination(
        dependent_variable_settings=propagation_setup.dependent_variable.mach_number("Capsule", "Uranus"),
        limit_value=2.0,
        use_as_lower_limit=True)

    dependent_variables_array = entry_sim(500, aero_coefficient_setting, 3.03327727e+07, 5.45941114e-01,
                                          -2.33346601e-02, 2.65992642e+04, -5.91036848e-01, -2.96367147e+00,
                                          [termination_altitude_setting, termination_mach_setting], acc=1)

    """plt.plot(dependent_variables_array[:, 0], dependent_variables_array[:, 1])
    plt.grid()
    plt.show()

    plt.plot(dependent_variables_array[:, 1], dependent_variables_array[:, 4])
    plt.grid()
    plt.show()

    plt.plot(dependent_variables_array[:, 1], dependent_variables_array[:, -3])
    plt.grid()
    plt.show()

    plt.plot(dependent_variables_array[:, 1], dependent_variables_array[:, -1])
    plt.grid()
    plt.show()

    c_d = np.zeros_like(dependent_variables_array[:, 1])
    for i,_ in enumerate(c_d):
        c_d[i] = drag.drag_coefficient([dependent_variables_array[i, -3], dependent_variables_array[i, 1]])[0]
    plt.plot(dependent_variables_array[:, -3], c_d)
    plt.grid()
    plt.show()

    print("finished")#"""
