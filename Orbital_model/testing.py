# General imports
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple

# Tudat imports
import tudatpy
from tudatpy.kernel.trajectory_design import transfer_trajectory
from tudatpy.kernel import constants
from tudatpy.kernel.numerical_simulation import environment_setup, propagation_setup, propagation
from tudatpy.kernel import numerical_simulation
from tudatpy.kernel.interface import spice
from tudatpy.util import result2array

# Pygmo imports
import pygmo as pg
import os



spice.load_standard_kernels()
path = os.path.dirname(__file__)
spice.load_kernel(path+'/ura111.bsp')
spice.load_kernel(path+'/Gravity.tpc')

print(constants)

print(spice.get_body_cartesian_state_at_epoch("Uranus","Sun","ECLIPJ2000","None",30*constants.JULIAN_YEAR))

print(spice.compute_rotation_matrix_between_frames("ECLIPJ2000","IAU_URANUS",30*constants.JULIAN_YEAR))


'''
simulation_start_epoch = 0

simulation_end_epoch = constants.JULIAN_DAY

bodies = environment_setup.create_simplified_system_of_bodies()

bodies_to_propagate = ["Uranus"]

global_frame_origin = ["Sun"]

# Define accelerations acting on Delfi-C3
#acceleration_settings_Uranus = dict(
#    "Sun":[propagation_setup.acceleration.point_mass_gravity()]
#)

acceleration_settings_Uranus = {"Sun":[propagation_setup.acceleration.point_mass_gravity()]}


acceleration_settings = {"Uranus": acceleration_settings_Uranus}

# Create acceleration models
acceleration_models = propagation_setup.create_acceleration_models(
    bodies, acceleration_settings, bodies_to_propagate, global_frame_origin)


# Create termination settings
termination_settings = propagation_setup.propagator.time_termination(simulation_end_epoch)

# Create numerical integrator settings
fixed_step_size = 10.0
integrator_settings = propagation_setup.integrator.runge_kutta_4(fixed_step_size)

#initial_state="dummy"

initial_state = propagation.get_state_of_bodies(bodies_to_propagate, global_frame_origin, bodies, simulation_start_epoch)

print ("STOP")
# Create propagation settings
propagator_settings = propagation_setup.propagator.translational(
    global_frame_origin,
    acceleration_models,
    bodies_to_propagate,
    initial_state,
    simulation_start_epoch,
    integrator_settings,
    termination_settings)

propagator_settings = propagation_setup.propagator.translational(
    global_frame_origin,
    acceleration_models,
    bodies_to_propagate,
    initial_state,
    simulation_start_epoch,
    integrator_settings,
    termination_settings
)

dynamics_simulator = numerical_simulation.create_dynamics_simulator(
    bodies, propagator_settings)
propagation_results = dynamics_simulator.propagation_results

states = propagation_results.state_history
states_array = result2array(states)

print(states_array[-1])
'''