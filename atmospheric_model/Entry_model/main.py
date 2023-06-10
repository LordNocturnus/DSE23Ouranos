import os.path

from matplotlib import pyplot as plt
import numpy as np
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
spice.load_kernel(_path+'/../GRAM/GRAM_Suite_1_5/SPICE/spk/satellites/ura116xl.bsp')
#spice.load_kernel(_path+'/Gravity.tpc')


def entry_sim(mass, drag_coefficient, diameter, alt, lat, lon, speed, flight_path_angle, heading_angle, acc=1):
    # Set simulation start epoch
    simulation_start_epoch = 6000.0

    # Set the maximum simulation time (avoid very long skipping re-entry)
    max_simulation_time = 3 * constants.JULIAN_DAY

    # define bodies in simulation
    bodies_to_create = ["Uranus"]

    # create body settings dictionary
    global_frame_origin = "Uranus"
    global_frame_orientation = "J2000"
    body_settings = environment_setup.get_default_body_settings(
        bodies_to_create, global_frame_origin, global_frame_orientation)

    # Modify the gravity field of Uranus to match GRAM
    body_settings.get("Uranus").gravity_field_settings = environment_setup.gravity_field.central_spice()

    # Modify the atmosphere to match GRAM
    gram = GRAM.GRAM()
    gram.altitudes = np.arange(7000, -290-acc, -acc)
    gram.lat = np.full_like(gram.altitudes, np.rad2deg(lat))
    gram.long = np.full_like(gram.altitudes, np.rad2deg(lon))
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
    atmos.to_csv(_path+"/atmos_data.csv", index=False)


    body_settings.get("Uranus").atmosphere_settings = environment_setup.atmosphere.tabulated(
        _path+"/atmos_data.csv", [environment_setup.atmosphere.AtmosphereDependentVariables.tabulated_density,
                                  environment_setup.atmosphere.AtmosphereDependentVariables.tabulated_pressure,
                                  environment_setup.atmosphere.AtmosphereDependentVariables.tabulated_temperature,
                                  environment_setup.atmosphere.AtmosphereDependentVariables.tabulated_gas_constant,
                                  environment_setup.atmosphere.AtmosphereDependentVariables.tabulated_specific_heat_ratio,
                                  environment_setup.atmosphere.AtmosphereDependentVariables.tabulated_molar_mass])
    # Create system of bodies from the body settings
    bodies = environment_setup.create_system_of_bodies(body_settings)

    bodies.create_empty_body("Capsule")

    # Set mass of vehicle
    bodies.get("Capsule").mass = mass

    # Create aerodynamic coefficient interface settings, and add to vehicle
    aero_coefficient_settings = environment_setup.aerodynamic_coefficients.constant(
        reference_area=np.pi * (diameter / 2) ** 2,
        constant_force_coefficient=[drag_coefficient, 0, 0])
    environment_setup.add_aerodynamic_coefficient_interface(bodies, "Capsule", aero_coefficient_settings)

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


if __name__ == "__main__":
    entry_sim(1550, 1.53, 4.5, 0.0, 0.0)
