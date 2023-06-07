from matplotlib import pyplot as plt
from tudatpy.kernel.interface import spice
from tudatpy.kernel import numerical_simulation
from tudatpy.kernel.numerical_simulation import environment_setup, environment, propagation_setup, propagation
from tudatpy.kernel.astro import element_conversion
from tudatpy.kernel import constants
from tudatpy.util import result2array
import numpy as np
import scipy as sp
import math

# Create a class for the aerodynamic guidance of the STS, inheriting from 'propagation.AerodynamicGuidance'
class STSAerodynamicGuidance:

    def __init__(self, bodies: environment.SystemOfBodies):

        # Extract the STS and Earth bodies
        self.vehicle = bodies.get_body("STS")
        self.earth = bodies.get_body("Earth")

        # Extract the STS flight conditions, angle calculator, and aerodynamic coefficient interface
        environment_setup.add_flight_conditions( bodies, 'STS', 'Earth')
        self.vehicle_flight_conditions = bodies.get_body("STS").flight_conditions
        self.aerodynamic_angle_calculator = self.vehicle_flight_conditions.aerodynamic_angle_calculator
        self.aerodynamic_coefficient_interface = self.vehicle_flight_conditions.aerodynamic_coefficient_interface

        self.current_time = float("NaN")

    def getAerodynamicAngles(self, current_time: float):
        self.updateGuidance( current_time )
        return np.array([self.angle_of_attack, 0.0, self.bank_angle])

    # Function that is called at each simulation time step to update the ideal bank angle of the vehicle
    def updateGuidance(self, current_time: float):

        if( math.isnan( current_time ) ):
            self.current_time = float("NaN")
        elif( current_time != self.current_time ):
            # Get the (constant) angular velocity of the Earth body
            earth_angular_velocity = np.linalg.norm(self.earth.body_fixed_angular_velocity)
            # Get the distance between the vehicle and the Earth bodies
            earth_distance = np.linalg.norm(self.vehicle.position)
            # Get the (constant) mass of the vehicle body
            body_mass = self.vehicle.mass

            # Extract the current Mach number, airspeed, and air density from the flight conditions
            mach_number = self.vehicle_flight_conditions.mach_number
            airspeed = self.vehicle_flight_conditions.airspeed
            density = self.vehicle_flight_conditions.density

            # Set the current Angle of Attack (AoA). The following line enforces the followings:
            # * the AoA is constant at 40deg when the Mach number is above 12
            # * the AoA is constant at 10deg when the Mach number is below 6
            # * the AoA varies close to linearly when the Mach number is between 12 and 6
            # * a Logistic relation is used so that the transition in AoA between M=12 and M=6 is smoother
            self.angle_of_attack = np.deg2rad(30 / (1 + np.exp(-2*(mach_number-9))) + 10)

            # Update the variables on which the aerodynamic coefficients are based (AoA and Mach)
            current_aerodynamics_independent_variables = [self.angle_of_attack, mach_number]
            # Update the aerodynamic coefficients
            self.aerodynamic_coefficient_interface.update_coefficients(
                current_aerodynamics_independent_variables, current_time)

            # Extract the current force coefficients (in order: C_D, C_S, C_L)
            current_force_coefficients = self.aerodynamic_coefficient_interface.current_force_coefficients
            # Extract the (constant) reference area of the vehicle
            aerodynamic_reference_area = self.aerodynamic_coefficient_interface.reference_area

            # Get the heading, flight path, and latitude angles from the aerodynamic angle calculator
            heading = self.aerodynamic_angle_calculator.get_angle(environment.heading_angle)
            flight_path_angle = self.aerodynamic_angle_calculator.get_angle(environment.flight_path_angle)
            latitude = self.aerodynamic_angle_calculator.get_angle(environment.latitude_angle)

            # Compute the acceleration caused by Lift
            lift_acceleration = 0.5 * density * airspeed ** 2 * aerodynamic_reference_area * current_force_coefficients[2] / body_mass
            # Compute the gravitational acceleration
            downward_gravitational_acceleration = self.earth.gravitational_parameter / (earth_distance ** 2)
            # Compute the centrifugal acceleration
            spacecraft_centrifugal_acceleration = airspeed ** 2 / earth_distance
            # Compute the Coriolis acceleration
            coriolis_acceleration = 2 * earth_angular_velocity * airspeed * np.cos(latitude) * np.sin(heading)
            # Compute the centrifugal acceleration from the Earth
            earth_centrifugal_acceleration = earth_angular_velocity ** 2 * earth_distance * np.cos(latitude) *             (np.cos(latitude) * np.cos(flight_path_angle) + np.sin(flight_path_angle) * np.sin(latitude) * np.cos(heading))

            # Compute the cosine of the ideal bank angle
            cosine_of_bank_angle = ((downward_gravitational_acceleration - spacecraft_centrifugal_acceleration) * np.cos(flight_path_angle) - coriolis_acceleration - earth_centrifugal_acceleration) / lift_acceleration
            # If the cosine lead to a value out of the [-1, 1] range, set to bank angle to 0deg or 180deg
            if (cosine_of_bank_angle < -1):
                self.bank_angle = np.pi
            elif (cosine_of_bank_angle > 1):
                self.bank_angle = 0.0
            else:
                # If the cos is in the correct range, return the computed bank angle
                self.bank_angle = np.arccos(cosine_of_bank_angle)
            self.current_time = current_time


if __name__ == "__main__":

    # Load spice kernels
    spice.load_standard_kernels()

    # Set simulation start epoch
    simulation_start_epoch = 6000.0

    # Set the maximum simulation time (avoid very long skipping re-entry)
    max_simulation_time = 3*constants.JULIAN_DAY


    # ## Environment setup
    #
    # Let’s create the environment for our simulation. This setup covers the creation of (celestial) bodies, vehicle(s), and environment interfaces.
    #

    # Create default body settings for "Earth"
    bodies_to_create = ["Uranus"]

    # Create default body settings for bodies_to_create, with "Earth"/"J2000" as the global frame origin and orientation
    global_frame_origin = "Uranus"
    global_frame_orientation = "J2000"
    body_settings = environment_setup.get_default_body_settings(
        bodies_to_create, global_frame_origin, global_frame_orientation)
    body_settings.get("Uranus").

    # Create system of bodies (in this case only Earth)
    bodies = environment_setup.create_system_of_bodies(body_settings)

    # Create vehicle object and set its constant mass
    bodies.create_empty_body("STS")
    bodies.get_body( "STS" ).set_constant_mass(5.0e3)