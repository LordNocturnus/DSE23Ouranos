# Load standard modules
import math
import numpy as np
from matplotlib import pyplot as plt

# Load tudatpy modules
from tudatpy.kernel.interface import spice
from tudatpy.kernel import numerical_simulation
from tudatpy.kernel.numerical_simulation import environment_setup, environment, propagation_setup, propagation
from tudatpy.kernel.astro import element_conversion
from tudatpy.kernel import constants
from tudatpy.util import result2array


def entry_sim(mass, drag_coefficient, diameter):
    # Manually create (empty) settings for a list of bodies, with origin Uranus
    body_settings = environment_setup.BodyListSettings(frame_origin="Uranus")

    # Modify the gravity field of Uranus to match GRAM
    body_settings.get("Uranus").gravity_field_settings = environment_setup.gravity_field.central_spice()
    body_settings.get("Uranus").atmosphere_settings = environment_setup.atmosphere.tabulated(__file__[:-7]+"atmos_data-csv",
                                                                                             ["tabulated_density",
                                                                                              "tabulated_pressure",
                                                                                              "tabulated_temperature",
                                                                                              "tabulated_gas_constant",
                                                                                              "tabulated_specific_heat_ratio",
                                                                                              "tabulated_molar_mass"])
    # Create system of bodies from the body settings
    bodies = environment_setup.create_system_of_bodies(body_settings)

    bodies.create_empty_body("Capsule")

    # Set mass of vehicle
    bodies.get("Capsule").mass = mass

    # Create aerodynamic coefficient interface settings, and add to vehicle
    aero_coefficient_settings = environment_setup.aerodynamic_coefficients.constant(
        reference_area=np.pi * (diameter / 2) ** 2,
        constant_force_coefficient=[drag_coefficient, 0, 0]
    )
    environment_setup.add_aerodynamic_coefficient_interface(
        bodies, "Capsule", aero_coefficient_settings)


if __name__ == "__main__":
    entry_sim(1550, 1.53, 4.5)
