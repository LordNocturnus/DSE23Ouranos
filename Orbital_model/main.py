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



#start of actual program

#some placeholder values
duration_of_entry = constants.JULIAN_DAY/24
atmospheric_mission_time = constants.JULIAN_DAY/24*3

#loading ephemeris
spice.load_standard_kernels()
path = os.path.dirname(__file__)
spice.load_kernel(path+'/ura111.bsp')
spice.load_kernel(path+'/Gravity.tpc')

# Create default body settings for "Uranus"
bodies_to_create = ["Uranus","Titania"]

radiusUranus = 25362000

# Set simulation start and end epochs
simulation_start_epoch = 0.0
simulation_end_epoch = constants.JULIAN_DAY*15

class orbital_trajectory:

    def __init__(self,initial_state):
        #throw error if initial state is incorrect
        if not len(initial_state) == 6:
            raise IndexError("The initial state has the wrong size")
        if np.linalg.norm(initial_state[:3]) < 1e8:
            raise ValueError("The initial state is too close to Uranus")
        if np.linalg.norm(initial_state[3:]) < 10:
            raise ValueError("The initial speed is too low")
        #storing inputs in class
        self.initial_state = initial_state

        #setting up basic properties of the orbital simulation as part of the class
        self.global_frame_origin = "Uranus"
        self.global_frame_orientation = "IAU_URANUS"
        self.radiusUranus = 25362000
        body_settings = environment_setup.get_default_body_settings(
            bodies_to_create, self.global_frame_origin, self.global_frame_orientation)
        self.bodies = environment_setup.create_system_of_bodies(body_settings)
        self.central_bodies = ["Uranus"]
        self.bodies.create_empty_body("Capsule")
        self.bodies.create_empty_body("Orbiter")
        self.uranus_gravitational_parameter = self.bodies.get("Uranus").gravitational_parameter
        print(self.uranus_gravitational_parameter)

    def initial_manoeuvre_capsule(self,desired_periapsis,capsule_retrograde = False):
        if desired_periapsis < 25362000* 0.95:
            raise ValueError("Periapsis location is too low, radisus of Uranus is 25362000 meters")
        radius = np.sqrt(self.initial_state[0]**2+self.initial_state[1]**2+self.initial_state[2]**2)
        initial_state_keplerian = astro.element_conversion.cartesian_to_keplerian(self.initial_state,self.uranus_gravitational_parameter)
        
        periapsis = initial_state_keplerian[0]/(1-initial_state_keplerian[1])

        def helperfunc(a,periapsis,true_anomaly,radius):
            return radius * (1+(1-periapsis/a)*np.cos(true_anomaly))/(1-(1-periapsis/a)**2) - a

        print (initial_state_keplerian)
        final_semi_major_axis = opt.fsolve(helperfunc,initial_state_keplerian[0],(desired_periapsis,initial_state_keplerian[5],radius))
        final_eccentricity = 1-desired_periapsis/final_semi_major_axis

        state = initial_state_keplerian

        state[0] = final_semi_major_axis
        state[1] = final_eccentricity

        self.initial_state_capsule = astro.element_conversion.keplerian_to_cartesian(state,self.uranus_gravitational_parameter)
        if np.inner(self.initial_state_capsule[:3],self.initial_state_capsule[3:]) > 0:
            self.initial_state_capsule[3:] = -self.initial_state_capsule[3:]

        return np.sqrt((self.initial_state[3]-self.initial_state_capsule[3]) ** 2 + (self.initial_state[4]-self.initial_state_capsule[4]) ** 2 + (self.initial_state[5]-self.initial_state_capsule[5]) ** 2 )
    
    def capsule_no_manoeuvre(self):
        self.initial_state_capsule= self.initial_state

    def initial_manoeuvre_orbiter(self,desired_periapsis):
        if desired_periapsis < 25367000:
            raise ValueError("Periapsis location is too low, radisus of Uranus is 25362000 meters. Orbiter would encounter the atmosphere.")
        radius = np.sqrt(self.initial_state[0]**2+self.initial_state[1]**2+self.initial_state[2]**2)
        initial_state_keplerian = astro.element_conversion.cartesian_to_keplerian(self.initial_state,self.uranus_gravitational_parameter)
        
        periapsis = initial_state_keplerian[0]/(1-initial_state_keplerian[1])

        def helperfunc(a,periapsis,true_anomaly,radius):
            return radius * (1+(1-periapsis/a)*np.cos(true_anomaly))/(1-(1-periapsis/a)**2) - a

        
        final_semi_major_axis = opt.fsolve(helperfunc,initial_state_keplerian[0],(desired_periapsis,initial_state_keplerian[5],radius))
        final_eccentricity = 1-desired_periapsis/final_semi_major_axis

        state = initial_state_keplerian

        state[0] = final_semi_major_axis
        state[1] = final_eccentricity

        self.initial_state_orbiter = astro.element_conversion.keplerian_to_cartesian(state,self.uranus_gravitational_parameter)
        if np.inner(self.initial_state_orbiter[:3],self.initial_state_orbiter[3:]) > 0:
            self.initial_state_orbiter[3:] = -self.initial_state_orbiter[3:]

        return np.sqrt((self.initial_state[3]-self.initial_state_orbiter[3]) ** 2 + (self.initial_state[4]-self.initial_state_orbiter[4]) ** 2 + (self.initial_state[5]-self.initial_state_orbiter[5]) ** 2 )

    def capsule_trajectory(self,atmosphere_height=5e6,step_size=10):
        if atmosphere_height < 6000:
            raise ValueError("Atmosphere is too low. Atmosphere height is given in meters.")
        if step_size <= 0:
            raise ValueError("The step size is zero or negative.")
        if not hasattr(self,"initial_state_capsule"):
            raise ValueError("The capsule initial state has not been set. use function .capsule_initial_manoeuvre or .capsule_no_manoeuvre before calling this funciton")
        
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

        initial_state_capsule_keplerian = astro.element_conversion.cartesian_to_keplerian(self.initial_state_capsule,self.uranus_gravitational_parameter)

        delta_mean_anomaly = - astro.element_conversion.true_to_mean_anomaly(initial_state_capsule_keplerian[1],initial_state_capsule_keplerian[5])

        time_to_periapsis = astro.element_conversion.delta_mean_anomaly_to_elapsed_time(delta_mean_anomaly,self.uranus_gravitational_parameter,initial_state_capsule_keplerian[0])

        print (time_to_periapsis/constants.JULIAN_DAY)

        timenotfound = True
        simulation_end_epoch = time_to_periapsis

        #doing simulation
        termination_settings = propagation_setup.propagator.time_termination(simulation_end_epoch)
        fixed_step_size = step_size
        integrator_settings = propagation_setup.integrator.runge_kutta_4(fixed_step_size)
        propagator_settings_capsule = propagation_setup.propagator.translational(
            self.central_bodies,
            acceleration_models_capsule,
            bodies_to_propagate_capsule,
            self.initial_state_capsule,
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
        altitude_capsule = radius_capsule - self.radiusUranus - atmosphere_height
        count = 0
        atmospheric_encounter = False
        notfound = True
        try:

            while notfound:
                if altitude_capsule[count] < 0:
                    atmospheric_encounter = count
                    notfound = False
                count += 1
                if count > len(altitude_capsule):
                    notfound = False
        except:
            altitude = np.min(altitude_capsule)
#            atmospheric_encounter = count -1
            raise ValueError('the periapsis is too high, it is currently',(altitude-25362000)/1000,'Kilometers above 1 bar. The atmospere starts 5000 kilometers above 1 bar')

        if atmospheric_encounter == False:
            raise RuntimeError("The atmosphere was not encountered. Set the periapsis lower.")
        else: 
            simulation_end_epoch = atmospheric_encounter * step_size

        self.atmospheric_encounter = simulation_end_epoch

        print (simulation_end_epoch/constants.JULIAN_DAY)

        capsule_position = np.array([states_capsule_array[atmospheric_encounter,1],states_capsule_array[atmospheric_encounter,2],states_capsule_array[atmospheric_encounter,3]])

        capsule_state_at_encounter = states_capsule[atmospheric_encounter*step_size]

        self.states_capsule_array = states_capsule_array

        return astro.element_conversion.cartesian_to_spherical(capsule_state_at_encounter)
    
    def get_capture_delta_v(self,desired_apoapsis):
        if not hasattr(self,"initial_state_orbiter"):
            raise ValueError("The initial state of the orbiter has not been set. Call .orbiter_initial_manoeuvre before calling this function.")
        initial_state_keplerian = astro.element_conversion.cartesian_to_keplerian(self.initial_state_orbiter,self.uranus_gravitational_parameter)
        initial_velocity_periapsis = np.sqrt(self.uranus_gravitational_parameter/initial_state_keplerian[0] * (1+initial_state_keplerian[1])/(1-initial_state_keplerian[1]))
        periapsis = initial_state_keplerian[0]* (1-initial_state_keplerian[1])

        final_semi_major_axis = (desired_apoapsis + periapsis) / 2
        final_eccentricity = 1 - periapsis/final_semi_major_axis

        captured_state_keplerian = initial_state_keplerian
        captured_state_keplerian[0] = final_semi_major_axis
        captured_state_keplerian[1] = final_eccentricity
        captured_state_keplerian[5] = 0
        captured_orbital_period = 2 * np.pi / np.sqrt(self.uranus_gravitational_parameter) * np.sqrt(final_semi_major_axis ** 3)

        self.captured_state = astro.element_conversion.keplerian_to_cartesian(captured_state_keplerian,self.uranus_gravitational_parameter)

        final_velocity_periapsis = np.sqrt(self.uranus_gravitational_parameter/final_semi_major_axis * (1+final_eccentricity)/(1-final_eccentricity))
        return initial_velocity_periapsis - final_velocity_periapsis, captured_orbital_period
    
    def orbiter_initial_trajectory(self,duration_of_entry,glider_position,step_size=10):
        if not hasattr(self,"initial_state_orbiter"):
            raise ValueError("The initial state of the orbiter has not been set. Call .orbiter_initial_manoeuvre before calling this function.")
        if duration_of_entry <= 0:
            raise ValueError("the entry duration should be positive.")
        if np.linalg.norm(glider_position[:3]) <= 25062000:
            raise ValueError("The capsule posisiton is too low.")
#        if np.linalg.norm(glider_position[:3]) > 30362000:
#            raise ValueError("The capsule posisiton is too high.")
        if abs(glider_position[1]) < 5 or abs(glider_position[2]) < 5:
            raise ValueError("The capsule posisiton is given in the wrong coordinate frame. Cartesian coordinates are expected.")
        if step_size <= 0:
            raise ValueError("The step size should be greater than 0.")

        start_atmospheric_mission = int(self.atmospheric_encounter + duration_of_entry)

        self.start_atmospheric_mission = start_atmospheric_mission

        self.start_atmospheric_mission_index = int(start_atmospheric_mission/step_size)

        end_atmsopheric_mission = int(start_atmospheric_mission + atmospheric_mission_time)

        #print((end_atmsopheric_mission-self.atmospheric_encounter)/step_size)

        self.end_atmospheric_mission_index = int(end_atmsopheric_mission/step_size)

        initial_state_orbiter_keplerian = astro.element_conversion.cartesian_to_keplerian(self.initial_state_orbiter,self.uranus_gravitational_parameter)

        delta_mean_anomaly = - astro.element_conversion.true_to_mean_anomaly(initial_state_orbiter_keplerian[1],initial_state_orbiter_keplerian[5])

        if delta_mean_anomaly < 0:
            delta_mean_anomaly = - delta_mean_anomaly

        time_to_periapsis = astro.element_conversion.delta_mean_anomaly_to_elapsed_time(delta_mean_anomaly,self.uranus_gravitational_parameter,initial_state_orbiter_keplerian[0])

        simulation_end_epoch = end_atmsopheric_mission 

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

        self.position_glider = glider_position[:3]

        #running simulation for orbiter to find periapsis
        termination_settings_before_capture = propagation_setup.propagator.time_termination(time_to_periapsis)
        termination_settings_after_capture = propagation_setup.propagator.time_termination(simulation_end_epoch)
        fixed_step_size = step_size
        integrator_settings = propagation_setup.integrator.runge_kutta_4(fixed_step_size)
        propagator_settings_before_capture = propagation_setup.propagator.translational(
            self.central_bodies,
            acceleration_models_orbiter,
            bodies_to_propagate_orbiter,
            self.initial_state_orbiter,
            simulation_start_epoch,
            integrator_settings,
            termination_settings_before_capture
        )
        propagator_settings_after_capture = propagation_setup.propagator.translational(
            self.central_bodies,
            acceleration_models_orbiter,
            bodies_to_propagate_orbiter,
            self.captured_state,
            time_to_periapsis,
            integrator_settings,
            termination_settings_after_capture
        )
        dynamics_simulator_before_capture = numerical_simulation.create_dynamics_simulator(
            self.bodies, propagator_settings_before_capture
        )

        states_orbiter_before_capture = dynamics_simulator_before_capture.state_history
        states_orbiter_array_before_capture = result2array(states_orbiter_before_capture)

        if simulation_end_epoch > time_to_periapsis:
            dynamics_simulator_after_capture = numerical_simulation.create_dynamics_simulator(
                self.bodies, propagator_settings_after_capture
            )

            states_orbiter_after_capture = dynamics_simulator_after_capture.state_history
            states_orbiter_array_after_capture = result2array(states_orbiter_after_capture)

            states_orbiter_array = np.append(states_orbiter_array_before_capture,states_orbiter_array_after_capture,axis=0)
        else:
            states_orbiter_array = states_orbiter_array_before_capture[:int(simulation_end_epoch/step_size)]

        self.states_orbiter_array = states_orbiter_array

        atmospheric_mission_orbiter_states_array = states_orbiter_array[int(start_atmospheric_mission/step_size):]

        self.atmospheric_mission_orbiter_states_array = atmospheric_mission_orbiter_states_array

        atmospheric_mission_orbiter_x = atmospheric_mission_orbiter_states_array[:,1]
        atmospheric_mission_orbiter_y = atmospheric_mission_orbiter_states_array[:,2]
        atmospheric_mission_orbiter_z = atmospheric_mission_orbiter_states_array[:,3]

        telemetry_vector = np.zeros((len(atmospheric_mission_orbiter_x),3))
        telemetry_angle= np.zeros_like(atmospheric_mission_orbiter_x)
        telemetry_distance = np.zeros_like(atmospheric_mission_orbiter_x)
        

        input_angle = np.zeros_like(atmospheric_mission_orbiter_x)

        for i in range(len(telemetry_vector)):

            telemetry_vector[i,0] = atmospheric_mission_orbiter_x[i] - self.position_glider[0]
            telemetry_vector[i,1] = atmospheric_mission_orbiter_y[i] - self.position_glider[1]
            telemetry_vector[i,2] = atmospheric_mission_orbiter_z[i] - self.position_glider[2]
            telemetry_distance[i] = np.sqrt(telemetry_vector[i,0] ** 2 + telemetry_vector[i,1] ** 2 + telemetry_vector[i,2] ** 2 )
            input_angle[i] = np.inner(telemetry_vector[i],self.position_glider)/(telemetry_distance[i]*np.linalg.norm(self.position_glider))
            telemetry_angle[i] = np.arccos(np.inner(telemetry_vector[i],self.position_glider)/(telemetry_distance[i]*np.linalg.norm(self.position_glider)))

        telemetry_distance_max = np.max(telemetry_distance)

        telemetry_angle_max = np.max(telemetry_angle)

        return telemetry_distance, telemetry_distance_max, telemetry_angle, telemetry_angle_max
    
    def plot_during_atmospheric_mission(self,plot_space_track=True,plot_distance_time=True,plot_angle_time=True):
        if not hasattr(self,"atmospheric_mission_orbiter_states_array"):
            raise AttributeError("The states of the orbiter are missing. Run .orbiter_initial_trajectory before this function")

        # Define a 3D figure using pyplot
        fig = plt.figure(figsize=(6,6), dpi=125)
        ax = fig.add_subplot(111, projection='3d')
        ax.set_title(f'Trajectories at Uranian encounter')

        # Plot the positional state history
        end_index = self.end_atmospheric_mission_index = self.start_atmospheric_mission_index
        ax.plot(self.atmospheric_mission_orbiter_states_array[:end_index,1], self.atmospheric_mission_orbiter_states_array[:end_index,2], self.atmospheric_mission_orbiter_states_array[:end_index,3], label="Orbiter", linestyle='-.')
        ax.scatter(self.position_glider[0],self.position_glider[1],self.position_glider[2], label="Glider", marker='x', color='red')

        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        x = np.cos(u)*np.sin(v)*self.radiusUranus
        y = np.sin(u)*np.sin(v)*self.radiusUranus
        z = np.cos(v)*self.radiusUranus
        ax.plot_wireframe(x, y, z, label="Uranus", color="blue")

        # Add the legend and labels, then show the plot
        ax.legend()
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        ax.set_zlabel('z [m]')
        ax.set_aspect('equal', adjustable='box')
        plt.show()

        plt.plot(np.arange(stop=len(telemetry_distance )* 10, start = 0, step = 10), telemetry_distance)
        plt.show()

        plt.plot(np.arange(stop=len(telemetry_distance) * 10, start = 0, step = 10), telemetry_angle/np.pi * 180)
        plt.show()

    def plot_orbital_trajectory(self):
        if not hasattr(self,"states_capsule_array"):
            raise AttributeError("The states of the capsule are missing. Run .capsule_trajectory before this function")
        if not hasattr(self,"atmospheric_mission_orbiter_states_array"):
            raise AttributeError("The states of the orbiter are missing. Run .orbiter_initial_trajectory before this function")
        # Define a 3D figure using pyplot
        fig = plt.figure(figsize=(6,6), dpi=125)
        ax = fig.add_subplot(111, projection='3d')
        ax.set_title(f'Trajectories at Uranian encounter')

        # Plot the positional state history
        ax.plot(self.states_orbiter_array[:, 1], self.states_orbiter_array[:, 2], self.states_orbiter_array[:, 3], label="Orbiter", linestyle='-.')
        #ax.scatter(0.0, 0.0, 0.0, label="Earth", marker='o', color='blue')

        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        x = np.cos(u)*np.sin(v)*self.radiusUranus
        y = np.sin(u)*np.sin(v)*self.radiusUranus
        z = np.cos(v)*self.radiusUranus
        ax.plot_wireframe(x, y, z, label="Uranus", color="blue")

        ax.plot(self.states_capsule_array[:, 1], self.states_capsule_array[:, 2], self.states_capsule_array[:, 3], label="Capsule", linestyle='-.')
        #ax.scatter(0.0, 0.0, 0.0, label="Earth", marker='o', color='red')

        # Add the legend and labels, then show the plot
        ax.legend()
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        ax.set_zlabel('z [m]')
        ax.set_aspect('equal', adjustable='box')
        plt.show()                

#func for local testing, ignore
if __name__ == "__main__":
    print("Hello World")


    initialstate = np.array([-1.08630339e+10 , 1.24446912e+10 , 7.25305409e+10, -6.57253947e+02 ,7.13997881e+02 , 4.13553122e+03])
    #initialstate = np.array([1e10,1e10,1e10,-4000,-1000,-1000])
    initialstate = np.array([-1.08630339e+10 , 1.24446912e+10 , 7.25305409e+10, 100,-100 , -400])
    
    desiredorbit = np.array([1,1,1,1,1,1])

    time=constants.JULIAN_DAY / 24 * 3

    entry_time = constants.JULIAN_DAY / 24

    trajectory = orbital_trajectory(initial_state=initialstate)

    trajectory.initial_manoeuvre_capsule(25380000)

    capsulestate = trajectory.capsule_trajectory(atmosphere_height=5e6,step_size=10)

    print ('capsule state is',capsulestate)

    capsulestatecartesian = astro.element_conversion.spherical_to_cartesian(capsulestate)

    initialmanoeuvre = trajectory.initial_manoeuvre_orbiter(25380000)

    capturedeltav, orbital_period = trajectory.get_capture_delta_v(58300000000000000)

    print(capturedeltav,orbital_period)

    telemetry_distance, telemetry_distance_max, telemetry_angle, telemetry_angle_max = trajectory.orbiter_initial_trajectory(duration_of_entry=entry_time,glider_position=capsulestatecartesian,step_size=10)

    trajectory.plot_during_atmospheric_mission()

    trajectory.plot_orbital_trajectory()

    print('break') 

