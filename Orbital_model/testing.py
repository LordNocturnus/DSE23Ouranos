# General imports
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple

# Tudat imports
import tudatpy
from tudatpy.kernel.trajectory_design import transfer_trajectory
from tudatpy.kernel import constants
from tudatpy.kernel.numerical_simulation import environment_setup, propagation_setup, propagation
from tudatpy.kernel import astro
from tudatpy.kernel import numerical_simulation
from tudatpy.kernel.interface import spice
from tudatpy.util import result2array

# Pygmo imports
import pygmo as pg
import os
import scipy.optimize as opt



spice.load_standard_kernels()
path = os.path.dirname(__file__)
spice.load_kernel(path+'/ura111.bsp')
spice.load_kernel(path+'/Gravity.tpc')

#print(constants)

#print(spice.get_body_cartesian_state_at_epoch("Uranus","Sun","ECLIPJ2000","None",30*constants.JULIAN_YEAR))

#print(spice.compute_rotation_matrix_between_frames("ECLIPJ2000","IAU_URANUS",30*constants.JULIAN_YEAR))

dv = 1000

peri1 = 25000000
peri2 = 30000000

gravparam = 5793939212817970.0
v_inf = 4000
r = 1e9

semi_major_1 = -gravparam/(v_inf**2)
ecc1 = 1 -peri1/semi_major_1

arg1 = 0
semi_latus_rectum1 = semi_major_1 * ( ecc1**2 -1)
true_anom1 = np.arccos(semi_latus_rectum1/(ecc1*r)-1/ecc1)
state = np.array([semi_major_1,ecc1,0,arg1,0,true_anom1])

cartesianstate1 = astro.element_conversion.keplerian_to_cartesian(state,gravparam)
cartesianpos1 = cartesianstate1[:3]

semi_major_2 = -gravparam/((v_inf+dv)**2)
ecc2 = 1 - peri2/semi_major_2

def getorbithelper(input,semi_major,ecc,refpos):
    arg = input[0]
    true_anom = input[1]

    refposx = refpos[0]
    refposy = refpos[1]
    state = np.array([semi_major,ecc,0,arg,0,true_anom])
    gravparam = 5793939212817970.0
    cartesianstate = astro.element_conversion.keplerian_to_cartesian(state,gravparam)
    errorx = cartesianstate[0]-refposx
    errory = cartesianstate[1]-refposy
    return errorx,errory

initguess = np.array([0.0,true_anom1])
optimal = opt.fsolve(getorbithelper,initguess,(semi_major_1,ecc1,cartesianpos1))

print(getorbithelper(optimal,semi_major_1,ecc1,cartesianpos1))


print(30362000. + 10.1e6)

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