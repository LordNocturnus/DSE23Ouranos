#generic imports
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple

#Tudat imports
import tudatpy
from tudatpy.kernel.trajectory_design import transfer_trajectory
from tudatpy.kernel import constants
from tudatpy.kernel import numerical_simulation
from tudatpy.kernel.numerical_simulation import environment_setup, propagation_setup
import tudatpy.kernel.astro as astro
from tudatpy.kernel.interface import spice
from tudatpy.util import result2array

#pygmo imports
import scipy.optimize as opt
import pygmo as pg
import os

#func for local testing, ignore
if __name__ == "__main__":
    print("Hello World")


#start of actual program

#loading ephemeris
spice.load_standard_kernels()
print(spice.get_total_count_of_kernels_loaded())
path = os.path.dirname(__file__)
spice.load_kernel(path+'/ura111.bsp')
spice.load_kernel(path+'/Gravity.tpc')

print(spice.get_total_count_of_kernels_loaded())

# Create default body settings for "Uranus"
bodies_to_create = ["Uranus","Titania"]

# Set simulation start and end epochs
simulation_start_epoch = 0.0
simulation_end_epoch = constants.JULIAN_DAY/240

# Create default body settings for bodies_to_create, with "Uranus"/"J2000" as the global frame origin and orientation
global_frame_origin = "Uranus"
global_frame_orientation = "J2000"
body_settings = environment_setup.get_default_body_settings(
    bodies_to_create, global_frame_origin, global_frame_orientation)

# Create system of bodies 
bodies = environment_setup.create_system_of_bodies(body_settings)
central_bodies = ["Uranus"]
#do encounter
bodies.create_empty_body("Orbiter")
bodies.create_empty_body("Capsule")

#bodies_to_propagate = ("Orbiter","Capsule")
bodies_to_propagate_orbiter = ["Orbiter"]
bodies_to_propagate_capsule = ["Capsule"]

# Define accelerations acting on objects
acceleration_settings = dict(
    Uranus=[propagation_setup.acceleration.point_mass_gravity()]
    )


#acceleration_settings = {"Orbiter": acceleration_settings_orbiter,"Capsule": acceleration_settings_capsule}
acceleration_settings_orbiter = {"Orbiter": acceleration_settings}
acceleration_settings_capsule = {"Capsule": acceleration_settings}

# Create acceleration models
acceleration_models_orbiter = propagation_setup.create_acceleration_models(
    bodies, acceleration_settings_orbiter, bodies_to_propagate_orbiter, central_bodies
    )

acceleration_models_capsule = propagation_setup.create_acceleration_models(
    bodies, acceleration_settings_capsule, bodies_to_propagate_capsule, central_bodies
    )

# Set initial conditions for the satellite that will be
# propagated in this simulation. The initial conditions are given in
# Keplerian elements and later on converted to Cartesian elements
uranus_gravitational_parameter = bodies.get("Uranus").gravitational_parameter
initial_state_orbiter = astro.element_conversion.keplerian_to_cartesian_elementwise(
    gravitational_parameter=uranus_gravitational_parameter,
    semi_major_axis=-7500.0e3,
    eccentricity=1.1,
    inclination=np.deg2rad(85.3),
    argument_of_periapsis=np.deg2rad(235.7),
    longitude_of_ascending_node=np.deg2rad(23.4),
    true_anomaly=np.deg2rad(30),)

initial_state_capsule = astro.element_conversion.keplerian_to_cartesian_elementwise(
    gravitational_parameter=uranus_gravitational_parameter,
    semi_major_axis=-7500.0e3,
    eccentricity=1.1,
    inclination=np.deg2rad(4.7),
    argument_of_periapsis=np.deg2rad(235.7),
    longitude_of_ascending_node=np.deg2rad(23.4),
    true_anomaly=np.deg2rad(30),)

# Create termination settings
termination_settings = propagation_setup.propagator.time_termination(simulation_end_epoch)

# Create numerical integrator settings
fixed_step_size = 10.0
integrator_settings = propagation_setup.integrator.runge_kutta_4(fixed_step_size)

# Create propagation settings
propagator_settings_orbiter = propagation_setup.propagator.translational(
    central_bodies,
    acceleration_models_orbiter,
    bodies_to_propagate_orbiter,
    initial_state_orbiter,
    simulation_start_epoch,
    integrator_settings,
    termination_settings
)

propagator_settings_capsule = propagation_setup.propagator.translational(
    central_bodies,
    acceleration_models_capsule,
    bodies_to_propagate_capsule,
    initial_state_capsule,
    simulation_start_epoch,
    integrator_settings,
    termination_settings
)

# Create simulation object and propagate the dynamics
dynamics_simulator_orbiter = numerical_simulation.create_dynamics_simulator(
    bodies, propagator_settings_orbiter
)

dynamics_simulator_capsule = numerical_simulation.create_dynamics_simulator(
    bodies, propagator_settings_capsule
)

# Extract the resulting state history and convert it to an ndarray
states_orbiter = dynamics_simulator_orbiter.state_history
states_orbiter_array = result2array(states_orbiter)

states_capsule = dynamics_simulator_capsule.state_history
states_capsule_array = result2array(states_capsule)

print(
    f"""
Single Earth-Orbiting Satellite Example.
The initial position vector of Orbiter is [km]: \n{
    states_orbiter[simulation_start_epoch][:3] / 1E3}
The initial velocity vector of Orbiter is [km/s]: \n{
    states_orbiter[simulation_start_epoch][3:] / 1E3}
\nAfter {simulation_end_epoch} seconds the position vector of Orbiter is [km]: \n{
    states_orbiter[simulation_end_epoch][:3] / 1E3}
And the velocity vector of the orbiter is [km/s]: \n{
    states_orbiter[simulation_start_epoch][3:] / 1E3}
    """
)
radiusUranus = 25362000

# Define a 3D figure using pyplot
fig = plt.figure(figsize=(6,6), dpi=125)
ax = fig.add_subplot(111, projection='3d')
ax.set_title(f'Trajectories at Uranian encounter')

# Plot the positional state history
ax.plot(states_orbiter_array[:, 1], states_orbiter_array[:, 2], states_orbiter_array[:, 3], label=bodies_to_propagate_orbiter[0], linestyle='-.')
#ax.scatter(0.0, 0.0, 0.0, label="Earth", marker='o', color='blue')

u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x = np.cos(u)*np.sin(v)*radiusUranus
y = np.sin(u)*np.sin(v)*radiusUranus
z = np.cos(v)*radiusUranus
ax.plot_wireframe(x, y, z, label="Uranus", color="blue")

ax.plot(states_capsule_array[:, 1], states_capsule_array[:, 2], states_capsule_array[:, 3], label=bodies_to_propagate_capsule[0], linestyle='-.')
#ax.scatter(0.0, 0.0, 0.0, label="Earth", marker='o', color='red')

# Add the legend and labels, then show the plot
ax.legend()
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_zlabel('z [m]')
plt.show()

#gravity assist program for inspiration

"""

def convert_trajectory_parameters (transfer_trajectory_object: tudatpy.kernel.trajectory_design.transfer_trajectory.TransferTrajectory,
                                   trajectory_parameters: List[float]
                                   ) -> Tuple[ List[float], List[List[float]], List[List[float]] ]:

    # Declare lists of transfer parameters
    node_times = list()
    leg_free_parameters = list()
    node_free_parameters = list()

    # Extract from trajectory parameters the lists with each type of parameters
    departure_time = trajectory_parameters[0]
    times_of_flight_per_leg = trajectory_parameters[1:]

    # Get node times
    # Node time for the intial node: departure time
    node_times.append(departure_time)
    # None times for other nodes: node time of the previous node plus time of flight
    accumulated_time = departure_time
    for i in range(0, transfer_trajectory_object.number_of_nodes - 1):
        accumulated_time += times_of_flight_per_leg[i]
        node_times.append(accumulated_time)

    # Get leg free parameters and node free parameters: one empty list per leg
    for i in range(transfer_trajectory_object.number_of_legs):
        leg_free_parameters.append( [ ] )
    # One empty array for each node
    for i in range(transfer_trajectory_object.number_of_nodes):
        node_free_parameters.append( [ ] )

    return node_times, leg_free_parameters, node_free_parameters


###########################################################################
# CREATE PROBLEM CLASS ####################################################
###########################################################################

class TransferTrajectoryProblem:

    def __init__(self,
                 transfer_trajectory_object: tudatpy.kernel.trajectory_design.transfer_trajectory.TransferTrajectory,
                 departure_date_lb: float, # Lower bound on departure date
                 departure_date_up: float, # Upper bound on departure date
                 legs_tof_lb: np.ndarray, # Lower bounds of each leg's time of flight
                 legs_tof_ub: np.ndarray): # Upper bounds of each leg's time of flight

#        Class constructor.


        self.departure_date_lb = departure_date_lb
        self.departure_date_ub = departure_date_ub
        self.legs_tof_lb = legs_tof_lb
        self.legs_tof_ub = legs_tof_ub

        # Save the transfer trajectory object as a lambda function
        # PyGMO internally pickles its user defined objects and some objects cannot be pickled properly without using lambda functions.
        self.transfer_trajectory_function = lambda: transfer_trajectory_object

    def get_bounds(self) -> tuple:

#        Returns the boundaries of the decision variables.


        # Retrieve transfer trajectory object
        transfer_trajectory_obj = self.transfer_trajectory_function()

        number_of_parameters = self.get_number_of_parameters()

        # Define lists to save lower and upper bounds of design parameters
        lower_bound = list(np.empty(number_of_parameters))
        upper_bound = list(np.empty(number_of_parameters))

        # Define boundaries on departure date
        lower_bound[0] = self.departure_date_lb
        upper_bound[0] = self.departure_date_ub

        # Define boundaries on time of flight between bodies ['Earth', 'Venus', 'Venus', 'Earth', 'Jupiter', 'Saturn']
        for i in range(0, transfer_trajectory_obj.number_of_legs):
            lower_bound[i+1] = self.legs_tof_lb[i]
            upper_bound[i+1] = self.legs_tof_ub[i]

        bounds = (lower_bound, upper_bound)
        return bounds

    def get_number_of_parameters(self):

#        Returns number of parameters that will be optimized


        # Retrieve transfer trajectory object
        transfer_trajectory_obj = self.transfer_trajectory_function()

        # Get number of parameters: it's the number of nodes (time at the first node, and time of flight to reach each subsequent node)
        number_of_parameters = transfer_trajectory_obj.number_of_nodes

        return number_of_parameters

    def fitness(self, trajectory_parameters: List[float]) -> list:

#        Returns delta V of the transfer trajectory object with the given set of trajectory parameters


        # Retrieve transfer trajectory object
        transfer_trajectory = self.transfer_trajectory_function()

        # Convert list of trajectory parameters to appropriate format
        node_times, leg_free_parameters, node_free_parameters = convert_trajectory_parameters(
            transfer_trajectory,
            trajectory_parameters)

        # Evaluate trajectory
        try:
            transfer_trajectory.evaluate(node_times, leg_free_parameters, node_free_parameters)
            delta_v = transfer_trajectory.delta_v

        # If there was some error in the evaluation of the trajectory, use a very large deltaV as penalty
        except:
            delta_v = 1e10

        return [delta_v]
    
###########################################################################
# Define transfer trajectory properties
###########################################################################

# Define the central body
#central_body = "Uranus"

# Define order of bodies (nodes)
#transfer_body_order = ['Earth', 'Venus', 'Venus', 'Earth', 'Jupiter', 'Saturn']

#transfer_body_order = ['Uranus', 'Titania', 'Uranus']


# Define departure orbit
#departure_semi_major_axis = 6371000 + 185000
#departure_eccentricity = 0
departure_semi_major_axis = 311000000
departure_eccentricity = 0.8745980707



# Define insertion orbit
arrival_semi_major_axis = 211000000
arrival_eccentricity = 0.8745980707

# Create simplified system of bodies
bodies = environment_setup.create_simplified_system_of_bodies()

# Define the trajectory settings for both the legs and at the nodes
transfer_leg_settings, transfer_node_settings = transfer_trajectory.mga_settings_unpowered_unperturbed_legs(
    transfer_body_order,
    departure_orbit=(departure_semi_major_axis, departure_eccentricity),
    arrival_orbit=(arrival_semi_major_axis, arrival_eccentricity),minimum_pericenters={'Titania':788400})

# Create the transfer calculation object
transfer_trajectory_object = transfer_trajectory.create_transfer_trajectory(
    bodies,
    transfer_leg_settings,
    transfer_node_settings,
    transfer_body_order,
    central_body)

# Lower and upper bound on departure date
departure_date_lb = 12273 * constants.JULIAN_DAY
departure_date_ub = 12274 * constants.JULIAN_DAY


# List of lower and upper on time of flight for each leg
legs_tof_lb = np.zeros(4)
legs_tof_ub = np.zeros(4)
# Venus first fly-by
legs_tof_lb[0] = 30 * constants.JULIAN_DAY
legs_tof_ub[0] = 6000 * constants.JULIAN_DAY
# Venus second fly-by
legs_tof_lb[1] = 30 * constants.JULIAN_DAY
legs_tof_ub[1] = 6000 * constants.JULIAN_DAY
# Earth fly-by
legs_tof_lb[2] = 30 * constants.JULIAN_DAY
legs_tof_ub[2] = 6000 * constants.JULIAN_DAY
# Jupiter fly-by
legs_tof_lb[3] = 30 * constants.JULIAN_DAY
legs_tof_ub[3] = 6000 * constants.JULIAN_DAY
# Saturn fly-by
#legs_tof_lb[4] = 1000 * constants.JULIAN_DAY
#legs_tof_ub[4] = 6000 * constants.JULIAN_DAY

###########################################################################
# Setup optimization
###########################################################################
# Initialize optimization class
optimizer = TransferTrajectoryProblem(transfer_trajectory_object,
                                        departure_date_lb,
                                        departure_date_ub,
                                        legs_tof_lb,
                                        legs_tof_ub)

# Creation of the pygmo problem object
prob = pg.problem(optimizer)

# To print the problem's information: uncomment the next line
# print(prob)

# Define number of generations per evolution
number_of_generations = 1

# Fix seed
#optimization_seed = 4444
#optimization_seed = np.random.random_integers(0,500)

# Create pygmo algorithm object
#algo = pg.algorithm(pg.de(gen=number_of_generations, seed=optimization_seed, F=0.5))
algo = pg.algorithm(pg.de(gen=number_of_generations, F=0.8))

# To print the algorithm's information: uncomment the next line
# print(algo)

# Set population size
population_size = 20

# Create population
#pop = pg.population(prob, size=population_size, seed=optimization_seed)
pop = pg.population(prob, size=population_size)

###########################################################################
# Run optimization
###########################################################################

# Set number of evolutions
number_of_evolutions = 800
number_of_evolutions = 2000

# Initialize empty containers
individuals_list = []
fitness_list = []

for i in range(number_of_evolutions):

    pop = algo.evolve(pop)

    # individuals save
    individuals_list.append(pop.champion_x)
    fitness_list.append(pop.champion_f)

print('The optimization has finished')

###########################################################################
# Results post-processing
###########################################################################

# Extract the best individual
print('\n########### CHAMPION INDIVIDUAL ###########\n')
print('Total Delta V [m/s]: ', pop.champion_f[0])
best_decision_variables = pop.champion_x/constants.JULIAN_DAY
print('Departure time w.r.t J2000 [days]: ', best_decision_variables[0])
print('Earth-Uranus time of flight [days]: ', best_decision_variables[1])
print('Mars-Venus time of flight [days]: ', best_decision_variables[2])
print('Venus-Earth time of flight [days]: ', best_decision_variables[3])
print('Earth-Uranus time of flight [days]: ', best_decision_variables[4])
#print('Jupiter-Saturn time of flight [days]: ', best_decision_variables[5])

# Plot fitness over generations
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(np.arange(0, number_of_evolutions), np.float_(fitness_list) / 1000, label='Function value: Feval')
# Plot champion
champion_n = np.argmin(np.array(fitness_list))
ax.scatter(champion_n, np.min(fitness_list) / 1000, marker='x', color='r', label='All-time champion', zorder=10)

# Prettify
ax.set_xlim((0, number_of_evolutions))
#ax.set_ylim([4, 25])
ax.grid('major')
ax.set_title('Best individual over generations', fontweight='bold')
ax.set_xlabel('Number of generation')
ax.set_ylabel(r'$\Delta V [km/s]$')
ax.legend(loc='upper right')
plt.tight_layout()
plt.legend()


# Reevaluate the transfer trajectory using the champion design variables
node_times, leg_free_parameters, node_free_parameters = convert_trajectory_parameters(transfer_trajectory_object, pop.champion_x)
transfer_trajectory_object.evaluate(node_times, leg_free_parameters, node_free_parameters)

# Extract the state history
state_history = transfer_trajectory_object.states_along_trajectory(500)
fly_by_states = np.array([state_history[node_times[i]] for i in range(len(node_times))])
state_history = result2array(state_history)
au = 1.5e11

# Plot the state history
fig = plt.figure(figsize=(8,5))
ax = fig.add_subplot(111)
ax.plot(state_history[:, 1] / au, state_history[:, 2] / au)
ax.scatter(fly_by_states[0, 0] / au, fly_by_states[0, 1] / au, color='blue', label='Earth departure')
ax.scatter(fly_by_states[1, 0] / au, fly_by_states[1, 1] / au, color='red', label='Mars fly-by')
ax.scatter(fly_by_states[2, 0] / au, fly_by_states[2, 1] / au, color='yellow', label='Venus fly-by')
ax.scatter(fly_by_states[3, 0] / au, fly_by_states[3, 1] / au, color='blue', label='Earth fly-by')
ax.scatter(fly_by_states[4, 0] / au, fly_by_states[4, 1] / au, color='blue', label='Uranus arrival')
#ax.scatter(fly_by_states[5, 0] / au, fly_by_states[5, 1] / au, color='grey', label='Saturn arrival')
ax.scatter([0], [0], color='orange', label='Sun')
ax.set_xlabel('x wrt Sun [AU]')
ax.set_ylabel('y wrt Sun [AU]')
ax.set_aspect('equal')
ax.legend(bbox_to_anchor=[1, 1])
plt.show()

"""