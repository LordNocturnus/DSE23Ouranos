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


def entry_sim(mass, drag_coefficient, diameter, lat, lon, acc=1):
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
    gram.lat = np.full_like(gram.altitudes, lat)
    gram.long = np.full_like(gram.altitudes, lon)
    gram.time = np.zeros_like(gram.altitudes)
    gram.run()

    atmos = {"altitude": np.asarray(gram.data.Height_km) * 1000,
             "density": np.asarray(gram.data.Density_kgm3),
             "pressure": np.asarray(gram.data.Pressure_Pa),
             "temperature": np.asarray(gram.data.Temperature_K),
             "gas_constant": np.asarray(gram.data.Height_km),
             "specific_heat_ratio": np.asarray(gram.data.Height_km),
             "molar_mass": np.asarray(gram.data.Height_km),}

    body_settings.get("Uranus").atmosphere_settings = environment_setup.atmosphere.tabulated(
        __file__[:-7]+"atmos_data-csv", ["tabulated_density", "tabulated_pressure", "tabulated_temperature",
                                         "tabulated_gas_constant", "tabulated_specific_heat_ratio",
                                         "tabulated_molar_mass"])
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


if __name__ == "__main__":
    entry_sim(1550, 1.53, 4.5, 0.0, 0.0)
