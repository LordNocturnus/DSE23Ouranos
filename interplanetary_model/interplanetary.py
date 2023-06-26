# General imports
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple

# Tudat imports
import tudatpy
from tudatpy.kernel.trajectory_design import transfer_trajectory
from tudatpy.kernel import constants
from tudatpy.kernel.numerical_simulation import environment_setup
from tudatpy.util import result2array
from tudatpy.kernel.interface import spice
from tudatpy.kernel import astro

# Pygmo imports
import pygmo as pg
import os

def convert_trajectory_parameters (transfer_trajectory_object: tudatpy.kernel.trajectory_design.transfer_trajectory.TransferTrajectory,
                                   trajectory_parameters: List[float], target_periapsis
                                   ) -> Tuple[ List[float], List[List[float]], List[List[float]]]:

    # Declare lists of transfer parameters
    node_times = list()
    leg_free_parameters = list()
    node_free_parameters = []

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
    #node_free_parameters.append([target_periapsis,0,0])

    return node_times, leg_free_parameters, node_free_parameters

###########################################################################
# CREATE PROBLEM CLASS ####################################################
###########################################################################

class TransferTrajectoryProblem:

    def __init__(self,
                 transfer_trajectory_object: tudatpy.kernel.trajectory_design.transfer_trajectory.TransferTrajectory,
                 departure_date_lb: float, # Lower bound on departure date
                 departure_date_ub: float, # Upper bound on departure date
                 legs_tof_lb: np.ndarray, # Lower bounds of each leg's time of flight
                 legs_tof_ub: np.ndarray, # Upper bounds of each leg's time of flight
                 periapsis_height):
        """
        Class constructor.
        """

        self.departure_date_lb = departure_date_lb
        self.departure_date_ub = departure_date_ub
        self.legs_tof_lb = legs_tof_lb
        self.legs_tof_ub = legs_tof_ub
        self.periapsis_height = periapsis_height

        # Save the transfer trajectory object as a lambda function
        # PyGMO internally pickles its user defined objects and some objects cannot be pickled properly without using lambda functions.
        self.transfer_trajectory_function = lambda: transfer_trajectory_object

    def get_bounds(self) -> tuple:
        """
        Returns the boundaries of the decision variables.
        """

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
        """
        Returns number of parameters that will be optimized
        """

        # Retrieve transfer trajectory object
        transfer_trajectory_obj = self.transfer_trajectory_function()

        # Get number of parameters: it's the number of nodes (time at the first node, and time of flight to reach each subsequent node)
        number_of_parameters = transfer_trajectory_obj.number_of_nodes

        return number_of_parameters

    def fitness(self, trajectory_parameters: List[float]) -> list:
        """
        Returns delta V of the transfer trajectory object with the given set of trajectory parameters
        """

        # Retrieve transfer trajectory object
        transfer_trajectory = self.transfer_trajectory_function()

        # Convert list of trajectory parameters to appropriate format
        node_times, leg_free_parameters, node_free_parameters = convert_trajectory_parameters(
            transfer_trajectory,
            trajectory_parameters, self.periapsis_height)

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


class interplanetary_trajectory:
    def __init__(self,transfer_body_order:list[str],departure_orbit:tuple[float,float],target_periapsis:float,launch_interval:tuple[float,float]):

        # Define the central body
        self.central_body = "Sun"

        # Define order of bodies (nodes)
        transfer_body_order = transfer_body_order

        # Define departure orbit
        self.departure_semi_major_axis = departure_orbit[0]
        self.departure_eccentricity = departure_orbit[1]


        arrival_apoapsis = 583000000
        arrival_periapsis = 25380000+10e6

        #arrival_apoapsis = 82664830
        #arrival_periapsis = 82664830

        arrival_semi_major_axis = (arrival_apoapsis + arrival_periapsis) / 1
        arrival_eccentricity = 1-arrival_periapsis/arrival_semi_major_axis

        #verification cassini orbit data
        # Define departure orbit
        #self.departure_semi_major_axis = np.inf
        #self.departure_eccentricity = 0

        # Define insertion orbit
        #arrival_semi_major_axis = 1.0895e8 / 0.02
        #arrival_eccentricity = 0.98


        self.target_periapsis = target_periapsis

        # Create simplified system of bodies
        self.bodies = environment_setup.create_simplified_system_of_bodies()

# Define the trajectory settings for both the legs and at the nodes
        self.transfer_leg_settings, self.transfer_node_settings = transfer_trajectory.mga_settings_unpowered_unperturbed_legs(
            transfer_body_order,
            departure_orbit=(self.departure_semi_major_axis, self.departure_eccentricity),minimum_pericenters = {'Earth': 6678000.0, 'Jupiter': 600000000.0, 'Mars': 3689000.0, 'Mercury': 2740000.0, 'Saturn': 70000000.0, 'Venus': 6351800.0, 'Uranus':target_periapsis},
            arrival_orbit=(arrival_semi_major_axis, arrival_eccentricity))

        # Create the transfer calculation object
        self.transfer_trajectory_object = transfer_trajectory.create_transfer_trajectory(
            self.bodies,
            self.transfer_leg_settings,
            self.transfer_node_settings,
            transfer_body_order,
            self.central_body)
        #transfer_trajectory.print_parameter_definitions(self.transfer_leg_settings,self.transfer_node_settings)

        # Lower and upper bound on departure date
#        self.departure_date_lb = 11000 * constants.JULIAN_DAY
#        self.departure_date_ub = 14000 * constants.JULIAN_DAY
        self.departure_date_lb = launch_interval[0]
        self.departure_date_ub = launch_interval[1]

        # List of lower and upper on time of flight for each leg
        legs_tof_lb = np.zeros(len(transfer_body_order)-1)
        legs_tof_ub = np.zeros(len(transfer_body_order)-1)
        for i in range(len(legs_tof_lb)):
        # Venus first fly-by
            if transfer_body_order[i+1] in ["Jupiter","Saturn"]:
                legs_tof_lb[i] = 400 * constants.JULIAN_DAY
                legs_tof_ub[i] = 2000 * constants.JULIAN_DAY                
            elif (transfer_body_order[i+1] in ['Venus','Earth','Mars'] and transfer_body_order[i] in ['Venus','Earth','Mars']):
                legs_tof_lb[i] = 30 * constants.JULIAN_DAY
                legs_tof_ub[i] = 500 * constants.JULIAN_DAY
            elif transfer_body_order[i+1] in ['Venus','Earth','Mars'] and transfer_body_order[i] in ['Jupiter','Saturn']:
                legs_tof_lb[i] = 100 * constants.JULIAN_DAY
                legs_tof_ub[i] = 1000 * constants.JULIAN_DAY   
            elif transfer_body_order[i+1] == "Uranus":
                legs_tof_lb[i] = 3 * constants.JULIAN_YEAR
                legs_tof_ub[i] = 20 * constants.JULIAN_YEAR

        #V&V cassini bounds
        '''
        # Lower and upper bound on departure date
        self.departure_date_lb = -1000.5 * constants.JULIAN_DAY
        self.departure_date_ub = -0.5 * constants.JULIAN_DAY

        # List of lower and upper on time of flight for each leg
        legs_tof_lb = np.zeros(5)
        legs_tof_ub = np.zeros(5)
        # Venus first fly-by
        legs_tof_lb[0] = 30 * constants.JULIAN_DAY
        legs_tof_ub[0] = 400 * constants.JULIAN_DAY
        # Venus second fly-by
        legs_tof_lb[1] = 100 * constants.JULIAN_DAY
        legs_tof_ub[1] = 470 * constants.JULIAN_DAY
        # Earth fly-by
        legs_tof_lb[2] = 30 * constants.JULIAN_DAY
        legs_tof_ub[2] = 400 * constants.JULIAN_DAY
        # Jupiter fly-by
        legs_tof_lb[3] = 400 * constants.JULIAN_DAY
        legs_tof_ub[3] = 2000 * constants.JULIAN_DAY
        # Saturn fly-by
        legs_tof_lb[4] = 1000 * constants.JULIAN_DAY
        legs_tof_ub[4] = 6000 * constants.JULIAN_DAY
        '''
        self.legs_tof_lb = legs_tof_lb
        self.legs_tof_ub = legs_tof_ub
        

        ###########################################################################
        # Setup optimization
        ###########################################################################
        # Initialize optimization class

    def optimize(self,repeats:float):
        optimizer = TransferTrajectoryProblem(self.transfer_trajectory_object,
                                                self.departure_date_lb,
                                                self.departure_date_ub,
                                                self.legs_tof_lb,
                                                self.legs_tof_ub,
                                                self.target_periapsis)

        # Creation of the pygmo problem object
        prob = pg.problem(optimizer)

        # To print the problem's information: uncomment the next line
        # print(prob)

        # Define number of generations per evolution
        number_of_generations = 1
        bestdv = 1e99
        # Fix seed
        for j in range(repeats):
            #print (j)
            # Fix seed
            optimization_seed = np.random.randint(0,1000)

            optimization_seed = int(optimization_seed)

            # Create pygmo algorithm object
            algo = pg.algorithm(pg.de(gen=number_of_generations, seed=optimization_seed, F=0.55,CR=0.7,variant=4))
            
            # To print the algorithm's information: uncomment the next line
            # print(algo)

            # Set population size
            population_size = 20

            # Create population
            pop = pg.population(prob, size=population_size, seed=optimization_seed)

            ###########################################################################
            # Run optimization
            ###########################################################################

            # Set number of evolutions
            self.number_of_evolutions = 800

            # Initialize empty containers
            individuals_list = []
            self.fitness_list = []

            for i in range(self.number_of_evolutions):

                pop = algo.evolve(pop)

                # individuals save
                individuals_list.append(pop.champion_x)
                self.fitness_list.append(pop.champion_f)
            if pop.champion_f < bestdv:
                bestdv = pop.champion_f
                bestpop = pop

        self.pop = bestpop

    def plot(self):
        ###########################################################################
        # Results post-processing
        ###########################################################################

        # Extract the best individual
        print('\n########### CHAMPION INDIVIDUAL ###########\n')
        print('Total Delta V [m/s]: ', self.pop.champion_f[0])
        best_decision_variables = self.pop.champion_x/constants.JULIAN_YEAR
        flight_times = best_decision_variables[1:]
        flight_time = np.sum(flight_times)
        print('Departure time w.r.t J2000 [days]: ', best_decision_variables[0])
        print('Earth-Venus time of flight [days]: ', flight_time)
        #print('Venus-Earth time of flight [days]: ', best_decision_variables[2])
        #print('Earth-Earth time of flight [days]: ', best_decision_variables[3])
        #print('Earth-Jupiter time of flight [days]: ', best_decision_variables[4])
        #print('Jupiter-Uranus time of flight [days]: ', best_decision_variables[5])
        '''
        # Extract the best individual, V&V
        print('\n########### CHAMPION INDIVIDUAL ###########\n')
        print('Total Delta V [m/s]: ', self.pop.champion_f[0])
        best_decision_variables = self.pop.champion_x/constants.JULIAN_DAY
        print('Departure time w.r.t J2000 [days]: ', best_decision_variables[0])
        print('Earth-Venus time of flight [days]: ', best_decision_variables[1])
        print('Venus-Venus time of flight [days]: ', best_decision_variables[2])
        print('Venus-Earth time of flight [days]: ', best_decision_variables[3])
        print('Earth-Jupiter time of flight [days]: ', best_decision_variables[4])
        print('Jupiter-Saturn time of flight [days]: ', best_decision_variables[5])
        '''
        # Plot fitness over generations
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(np.arange(0, self.number_of_evolutions), np.float_(self.fitness_list) / 1000, label='Function value: Feval')
        # Plot champion
        champion_n = np.argmin(np.array(self.fitness_list))
        ax.scatter(champion_n, np.min(self.fitness_list) / 1000, marker='x', color='r', label='All-time champion', zorder=10)

        # Prettify
        ax.set_xlim((0, self.number_of_evolutions))
        ax.set_ylim([4, 25])
        ax.grid('major')
        ax.set_title('Best individual over generations', fontweight='bold')
        ax.set_xlabel('Number of generation')
        ax.set_ylabel(r'$\Delta V [km/s]$')
        ax.legend(loc='upper right')
        plt.tight_layout()
        plt.legend()
        plt.show()

        # Reevaluate the transfer trajectory using the champion design variables
        node_times, leg_free_parameters, node_free_parameters = convert_trajectory_parameters(self.transfer_trajectory_object, self.pop.champion_x,self.target_periapsis)
        self.transfer_trajectory_object.evaluate(node_times, leg_free_parameters, node_free_parameters)

        # Extract the state history
        state_history = self.transfer_trajectory_object.states_along_trajectory(500)
        fly_by_states = np.array([state_history[node_times[i]] for i in range(len(node_times))])
        self.state_history = result2array(state_history)
        au = 1.5e11

        delta_v_nodes = self.transfer_trajectory_object.delta_v_per_node

        #periapsis_heights = self.transfer_trajectory_object.node_parameters

        print (delta_v_nodes)

        # Plot the state history
        fig = plt.figure(figsize=(8,5))
        ax = fig.add_subplot(111)
        ax.plot(self.state_history[:, 1] / au, self.state_history[:, 2] / au)
        ax.scatter(fly_by_states[0, 0] / au, fly_by_states[0, 1] / au, color='blue', label='Earth departure')
        ax.scatter(fly_by_states[1, 0] / au, fly_by_states[1, 1] / au, color='Red', label='Mars fly-by')
        ax.scatter(fly_by_states[2, 0] / au, fly_by_states[2, 1] / au, color='Brown', label='Jupiter fly-by')
        ax.scatter(fly_by_states[3, 0] / au, fly_by_states[3, 1] / au, color='Blue', label='Uranus arrival')
        #ax.scatter(fly_by_states[4, 0] / au, fly_by_states[4, 1] / au, color='red', label='Jupiter fly-by')
        #ax.scatter(fly_by_states[5, 0] / au, fly_by_states[5, 1] / au, color='grey', label='Saturn arrival')
        ax.scatter([0], [0], color='orange', label='Sun')
        ax.set_xlabel('x wrt Sun [AU]')
        ax.set_ylabel('y wrt Sun [AU]')
        ax.set_aspect('equal')
        ax.legend(bbox_to_anchor=[1, 1])
        plt.show()
    def output(self,time_difference:float):
        spice.load_standard_kernels()
        path = os.path.dirname(__file__)
        spice.load_kernel(path+'/ura111.bsp')
        spice.load_kernel(path+'/Gravity.tpc')
        print("The departure burn will be",self.pop.champion_f[0],'delta v')
        time_of_encounter = self.state_history[-1,0]
        time_of_separation = time_of_encounter-time_difference
        i = 0
        position_Uranus = spice.get_body_cartesian_position_at_epoch("Uranus","Sun","ECLIPJ2000","None",time_of_separation)
        helper_Uranus = spice.get_body_cartesian_position_at_epoch("Uranus","Sun","ECLIPJ2000","None",time_of_separation+1)
        velocity_Uranus = helper_Uranus - position_Uranus
        rotationmatrix = spice.compute_rotation_matrix_between_frames('ECLIPJ2000','IAU_URANUS',time_of_separation)
        position_earth = spice.get_body_cartesian_position_at_epoch("Earth","Sun","ECLIPJ2000","None",time_of_separation)
        distance_earth = np.linalg.norm(position_Uranus - position_earth)
        print ('The distance to earth is',distance_earth)
        print ('The vector from Uranus to Earth is:',(position_earth - position_Uranus)@ rotationmatrix)
        max_distance_earth = distance_earth
        currentdate = time_of_separation
        for i in range(1000):
            position_Uranus = spice.get_body_cartesian_position_at_epoch("Uranus","Sun","ECLIPJ2000","None",currentdate)
            helper_Uranus = spice.get_body_cartesian_position_at_epoch("Uranus","Sun","ECLIPJ2000","None",currentdate+1)
            velocity_Uranus = helper_Uranus - position_Uranus
            rotationmatrix = spice.compute_rotation_matrix_between_frames('ECLIPJ2000','IAU_URANUS',currentdate)
            position_earth = spice.get_body_cartesian_position_at_epoch("Earth","Sun","ECLIPJ2000","None",currentdate)
            distance_earth = np.linalg.norm(position_Uranus - position_earth)
            currentdate += constants.JULIAN_DAY
            if distance_earth > max_distance_earth:
                max_distance_earth = distance_earth

        print ('The maximum distance to Earth during the extended mission is', max_distance_earth)
        notfound = True
        print('state',self.state_history[-1])
        while notfound:
            #if self.state_history[i,0] > time_of_separation:
                cartesianstate = self.state_history[-1,1:]
                notfound = False
                position = cartesianstate[:3]
                velocity = cartesianstate[3:]

                print(position_Uranus,'\n',velocity_Uranus,'\n',velocity)

                rotationvelocitymatrix =   spice.compute_rotation_matrix_derivative_between_frames('ECLIPJ2000','IAU_URANUS',time_of_separation)
                position -= position_Uranus
                velocity -= velocity_Uranus
                position_to_Uranus = position @ rotationmatrix     
                velocity_to_Uranus = velocity @ rotationmatrix #- rotationvelocitymatrix @ np.transpose(rotationmatrix) @ velocity   
                print('The velocity magnitude is',np.linalg.norm(velocity),np.linalg.norm(velocity_to_Uranus))      
                state = np.append(position_to_Uranus,velocity_to_Uranus)
                statekeplerian = astro.element_conversion.cartesian_to_keplerian(state,5793939212817970.0)
                print ('The state of the system in keplerian coordinates is',statekeplerian)
                return state

            #else:
                #i += 1




if __name__ == "__main__":
    print("Hello, world")
    #planets = ['Earth','Venus','Earth','Earth','Jupiter','Uranus']
    #planets = ['Earth','Uranus']
    planets = ['Earth','Mars','Jupiter','Uranus']
    #verification
    #planets = ['Earth', 'Venus', 'Venus', 'Earth', 'Jupiter', 'Saturn']
    earthorbit = (6521000,0)
    periapsis = 30000000
    launching = (30*constants.JULIAN_YEAR,40*constants.JULIAN_YEAR)
    trajectory = interplanetary_trajectory(planets,earthorbit,periapsis,launching)
    trajectory.optimize(25)
    trajectory.plot()
    state = trajectory.output(constants.JULIAN_DAY * 8)
    print(state)