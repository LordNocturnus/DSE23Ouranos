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

#some placeholder values
duration_of_entry = constants.JULIAN_DAY/24
atmospheric_mission_time = constants.JULIAN_DAY/24*3

#loading ephemeris
spice.load_standard_kernels()
print(spice.get_total_count_of_kernels_loaded())
path = os.path.dirname(__file__)
spice.load_kernel(path+'/ura111.bsp')
spice.load_kernel(path+'/Gravity.tpc')

print(spice.get_total_count_of_kernels_loaded())

# Create default body settings for "Uranus"
bodies_to_create = ["Uranus","Titania"]

radiusUranus = 25362000

# Set simulation start and end epochs
simulation_start_epoch = 0.0
simulation_end_epoch = constants.JULIAN_DAY*15

class orbital_trajectory:

    def __init__(self,initial_state,desired_orbit,duration_of_mission):
        #storing inputs in class
        self.initial_state = initial_state
        self.desired_orbit = desired_orbit

        #setting up basic properties of the orbital simulation as part of the class
        self.global_frame_origin = "Uranus"
        self.global_frame_orientation = "J2000"
        body_settings = environment_setup.get_default_body_settings(
            bodies_to_create, self.global_frame_origin, self.global_frame_orientation)
        self.bodies = environment_setup.create_system_of_bodies(body_settings)
        self.duration_of_mission = duration_of_mission
        self.central_bodies = ["Uranus"]
        self.bodies.create_empty_body("Capsule")
        self.bodies.create_empty_body("Orbiter")
        self.uranus_gravitational_parameter = self.bodies.get("Uranus").gravitational_parameter

    def capsuletrajectory(self,atmosphere_height=5e6,step_size=10):

        bodies_to_propagate_capsule = ["Capsule"]

        # Define accelerations acting on objects
        acceleration_settings = dict(
            Uranus=[propagation_setup.acceleration.point_mass_gravity()]
            )


        #acceleration_settings = {"Orbiter": acceleration_settings_orbiter,"Capsule": acceleration_settings_capsule}
        acceleration_settings_capsule = {"Capsule": acceleration_settings}

        # Create acceleration models
        acceleration_models_capsule = propagation_setup.create_acceleration_models(
            self.bodies, acceleration_settings_capsule, bodies_to_propagate_capsule, self.central_bodies
            )

        # Set initial conditions for the satellite that will be
        # propagated in this simulation. The initial conditions are given in
        # Keplerian elements and later on converted to Cartesian elements
        #uranus_gravitational_parameter = self.bodies.get("Uranus").gravitational_parameter

        
        initial_state_capsule = np.array([1e9,1e9,1e9,0,-380,0])

        print (initial_state_orbiter)

        timenotfound = True
        simulation_end_epoch = 40*constants.JULIAN_DAY

        #doing simulation
        termination_settings = propagation_setup.propagator.time_termination(simulation_end_epoch)
        fixed_step_size = step_size
        integrator_settings = propagation_setup.integrator.runge_kutta_4(fixed_step_size)
        propagator_settings_capsule = propagation_setup.propagator.translational(
            self.central_bodies,
            acceleration_models_capsule,
            bodies_to_propagate_capsule,
            initial_state_capsule,
            simulation_start_epoch,
            integrator_settings,
            termination_settings
        )
        dynamics_simulator_capsule = numerical_simulation.create_dynamics_simulator(
            self.bodies, propagator_settings_capsule
        )

        #postprocessing data

        #getting states from simulator
        states_capsule = dynamics_simulator_capsule.state_history
        states_capsule_array = result2array(states_capsule)

        #detecting when capsule encounters atmosphere
        radius_capsule = np.sqrt( states_capsule_array[:, 1] ** 2 + states_capsule_array[:, 2] ** 2 + states_capsule_array[:, 3] ** 2 )
        altitude_capsule = radius_capsule - radiusUranus - atmosphere_height
        count = 0
        atmospheric_encounter = False
        notfound = True
        while notfound:
            if altitude_capsule[count] < 0:
                atmospheric_encounter = count
                notfound = False
            count += 1
            if count > len(altitude_capsule):
                notfound = False

        if atmospheric_encounter == False:
            print ('atmosphere missed!')
        else: 
            simulation_end_epoch = atmospheric_encounter * step_size

        self.atmospheric_encounter = atmospheric_encounter

        capsule_position = np.array([states_capsule_array[atmospheric_encounter,1],states_capsule_array[atmospheric_encounter,2],states_capsule_array[atmospheric_encounter,3]])

        print ('Atmospheric encounter at',atmospheric_encounter * step_size / constants.JULIAN_DAY,'days')

        capsule_state_at_encounter = states_capsule[atmospheric_encounter*step_size]

        self.states_capsule_array = states_capsule_array

        return astro.element_conversion.cartesian_to_spherical(capsule_state_at_encounter)
    
    def get_capture_delta_v(self):
        initial_state_keplerian = astro.element_conversion.cartesian_to_keplerian(self.initial_state)
        initial_velocity_periapsis = np.sqrt(self.uranus_gravitational_parameter/initial_state_keplerian[0] * (1+initial_state_keplerian[1])/(1-initial_state_keplerian[1]))
        final_velocity_periapsis = np.sqrt(self.uranus_gravitational_parameter/self.desired_orbit[0] * (1+self.desired_orbit[1])/(1-self.desired_orbit[1]))
        return initial_velocity_periapsis - final_velocity_periapsis
    
    def orbiter_initial_trajectory(self,duration_of_entry,capsule_position,step_size=10):

        start_atmospheric_mission = int(self.atmospheric_encounter + duration_of_entry/step_size)

        end_atmsopheric_mission = int(start_atmospheric_mission + atmospheric_mission_time/step_size)

        initial_state_orbiter_keplerian = self.initial_state

        initial_state_captured = np.array([desiredorbit[0],desiredorbit[1],initial_state_orbiter_keplerian[2],initial_state_orbiter_keplerian[3],initial_state_orbiter_keplerian[4],initial_state_orbiter_keplerian[5],])

        delta_mean_anomaly = astro.element_conversion.true_to_mean_anomaly(initial_state_orbiter_keplerian[1],initial_state_orbiter_keplerian[5])

        time_to_periapsis = astro.element_conversion.delta_mean_anomaly_to_elapsed_time(delta_mean_anomaly,self.uranus_gravitational_parameter,initial_state_orbiter_keplerian[0])

        simulation_end_epoch = end_atmsopheric_mission * step_size

        bodies_to_propagate_orbiter = ["Orbiter"]

        # Define accelerations acting on objects
        acceleration_settings = dict(
            Uranus=[propagation_setup.acceleration.point_mass_gravity()]
            )


        #acceleration_settings = {"Orbiter": acceleration_settings_orbiter,"Capsule": acceleration_settings_capsule}
        acceleration_settings_orbiter = {"Orbiter": acceleration_settings}

        # Create acceleration models
        acceleration_models_orbiter = propagation_setup.create_acceleration_models(
            self.bodies, acceleration_settings_orbiter, bodies_to_propagate_orbiter, self.central_bodies
            )

        position_glider = capsule_position[:3]

        #running simulation for orbiter to find periapsis
        termination_settings_before_capture = propagation_setup.propagator.time_termination(time_to_periapsis * step_size)
        termination_settings_after_capture = propagation_setup.propagator.time_termination((end_atmsopheric_mission - time_to_periapsis) * step_size)
        fixed_step_size = step_size
        integrator_settings = propagation_setup.integrator.runge_kutta_4(fixed_step_size)
        propagator_settings_before_capture = propagation_setup.propagator.translational(
            self.central_bodies,
            acceleration_models_orbiter,
            bodies_to_propagate_orbiter,
            initial_state_orbiter,
            simulation_start_epoch,
            integrator_settings,
            termination_settings_before_capture
        )
        propagator_settings_after_capture = propagation_setup.propagator.translational(
            self.central_bodies,
            acceleration_models_orbiter,
            bodies_to_propagate_orbiter,
            initial_state_captured,
            time_to_periapsis,
            integrator_settings,
            termination_settings_after_capture
        )
        dynamics_simulator_before_capture = numerical_simulation.create_dynamics_simulator(
            self.bodies, propagator_settings_before_capture
        )
        dynamics_simulator_after_capture = numerical_simulation.create_dynamics_simulator(
            self.bodies, propagator_settings_after_capture
        )

        states_orbiter_before_capture = dynamics_simulator_before_capture.state_history
        states_orbiter_array_before_capture = result2array(states_orbiter_before_capture)

        states_orbiter_after_capture = dynamics_simulator_after_capture.state_history
        states_orbiter_array_after_capture = result2array(states_orbiter_after_capture)

        states_orbiter_array = np.append(states_orbiter_array_before_capture,states_orbiter_array_after_capture)

        atmospheric_mission_orbiter_states_array = states_orbiter_array[start_atmospheric_mission:]

        print(len(states_orbiter_array),'\n',start_atmospheric_mission)

        atmospheric_mission_orbiter_x = atmospheric_mission_orbiter_states_array[:,1]
        atmospheric_mission_orbiter_y = atmospheric_mission_orbiter_states_array[:,2]
        atmospheric_mission_orbiter_z = atmospheric_mission_orbiter_states_array[:,3]

        telemetry_vector = np.zeros((len(atmospheric_mission_orbiter_x),3))
        telemetry_angle= np.zeros_like(atmospheric_mission_orbiter_x)
        telemetry_distance = np.zeros_like(atmospheric_mission_orbiter_x)

        input_angle = np.zeros_like(atmospheric_mission_orbiter_x)

        for i in range(len(telemetry_vector)):

            telemetry_vector[i,0] = atmospheric_mission_orbiter_x[i] - position_glider[0]
            telemetry_vector[i,1] = atmospheric_mission_orbiter_y[i] - position_glider[1]
            telemetry_vector[i,2] = atmospheric_mission_orbiter_z[i] - position_glider[2]
            telemetry_distance[i] = np.sqrt(telemetry_vector[i,0] ** 2 + telemetry_vector[i,1] ** 2 + telemetry_vector[i,2] ** 2 )
            input_angle[i] = np.inner(telemetry_vector[i],position_glider)/(telemetry_distance[i]*np.linalg.norm(position_glider))
            telemetry_angle[i] = np.arccos(np.inner(telemetry_vector[i],position_glider)/(telemetry_distance[i]*np.linalg.norm(position_glider)))

        return telemetry_distance, telemetry_angle

        

initialstate = np.array([1e9,1e9,1e9,-380,0,0])

desiredorbit = np.array([1,1,1,1,1,1])

time=constants.JULIAN_DAY / 24 * 3

entry_time = constants.JULIAN_DAY / 24

trajectory = orbital_trajectory(initial_state=initialstate,desired_orbit=desiredorbit,duration_of_mission=time)

capsulestate = trajectory.capsuletrajectory(atmosphere_height=5e6,step_size=10)

capture_velocity = trajectory.get_capture_delta_v

capsulestatecartesian = astro.element_conversion.spherical_to_cartesian(capsulestate)

telemetry_distance, telemetry_angle = trajectory.orbiter_initial_trajectory(duration_of_entry=entry_time,capsule_position=capsulestatecartesian,step_size=10)

print('break') 



plt.plot(np.arange(stop=len(telemetry_distance )* 10, start = 0, step = 10), telemetry_distance)
plt.show()

plt.plot(np.arange(stop=len(telemetry_distance) * 10, start = 0, step = 10), telemetry_angle/np.pi * 180)
plt.show()

# Define a 3D figure using pyplot
fig = plt.figure(figsize=(6,6), dpi=125)
ax = fig.add_subplot(111, projection='3d')
ax.set_title(f'Trajectories at Uranian encounter')

# Plot the positional state history
ax.plot(atmospheric_mission_orbiter_x, atmospheric_mission_orbiter_y, atmospheric_mission_orbiter_z, label=bodies_to_propagate_orbiter[0], linestyle='-.')
ax.scatter(position_glider[0],position_glider[1],position_glider[2], label="Glider", marker='x', color='red')

u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x = np.cos(u)*np.sin(v)*radiusUranus
y = np.sin(u)*np.sin(v)*radiusUranus
z = np.cos(v)*radiusUranus
ax.plot_wireframe(x, y, z, label="Uranus", color="blue")

# Add the legend and labels, then show the plot
ax.legend()
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_zlabel('z [m]')
ax.set_aspect('equal', adjustable='box')
plt.show()

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
The initial position vector of Orbiter is [km]: \n{
    states_orbiter[simulation_start_epoch][:3] / 1E3}
The initial velocity vector of Orbiter is [km/s]: \n{
    states_orbiter[simulation_start_epoch][3:] / 1E3}
\nAfter {simulation_end_epoch} seconds the position vector of Orbiter is [km]: \n{
    states_orbiter[simulation_end_epoch][:3] / 1E3}
And the velocity vector of the orbiter is [km/s]: \n{
    states_orbiter[simulation_end_epoch][3:] / 1E3}
    """
)


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
ax.set_aspect('equal', adjustable='box')
plt.show()

#spherical_capsule_data = astro.cartesiantospherical(states_capsule)
radius_capsule = np.sqrt( states_capsule_array[:, 1] ** 2 + states_capsule_array[:, 2] ** 2 + states_capsule_array[:, 3] ** 2 )
altitude_capsule = radius_capsule - radiusUranus

plt.plot(np.arange(stop=len(altitude_capsule) * 10, start = 0, step = 10), altitude_capsule)
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