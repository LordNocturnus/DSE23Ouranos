from tudatpy.kernel.numerical_simulation import environment_setup, environment, propagation_setup
import numpy as np
import scipy as sp

import atmospheric_model.Entry_model
from atmospheric_model.GRAM import GRAM


_PICA_density = 0.352 * 100**3 / 1000  # 0.352 - 0.701 g/cm^3
_radius1 = 1.125
_radius2 = 0.126
_bottom_angle = np.deg2rad(20)
_top_angle = np.deg2rad(59.73)

_line_safe_strength = 1300 * 4  # N with safety factor of 3
_line_weight_per_meter = 0.013  # kg/m
_line_cost_per_meter = 0.7  # €/m


class Aeroshell:

    def __init__(self, target_pressure, target_time, glider_mass, diameter, parachute_c_ds, shock_load_factors,
                 adcs_power, pr_power, delta_t, alt, lat, lon, speed, flight_path_angle, heading_angle):
        self.target_pressure = target_pressure
        self.target_time = target_time
        self.glider_mass = glider_mass
        self.heatshield_mass = 0.0
        self.chute_weight = 0.0
        self.diameter = diameter
        self.parachute_c_ds = parachute_c_ds
        self.shock_load_factors = shock_load_factors
        self.adcs_power = adcs_power
        self.pr_power = pr_power
        self.delta_t = delta_t
        self.alt = alt
        self.lat = lat
        self.lon = lon
        self.speed = speed
        self.flight_path_angle = flight_path_angle
        self.heading_angle = heading_angle
        self.gram = GRAM()
        self.gram.run()
        self.target_altitude = np.asarray(self.gram.data.Height_km[self.gram.data.Pressure_Pa >= self.target_pressure],
                                          dtype=float)[0] * 1000
        self.ballistic_dep_vars = np.zeros([1, 12])
        self.sonic_dep_vars = np.zeros([1, 12])
        self.parachute_dep_vars = np.zeros([1, 12])

        self.sonic_alt = alt
        self.sonic_lat = lat
        self.sonic_lon = lon
        self.sonic_speed = speed
        self.sonic_flight_path_angle = flight_path_angle
        self.sonic_heading_angle = heading_angle
        self.parachute_diameter = diameter

        self.parachute_alt = alt
        self.parachute_lat = lat
        self.parachute_lon = lon
        self.parachute_speed = speed
        self.parachute_flight_path_angle = flight_path_angle
        self.parachute_heading_angle = heading_angle

        self.heat_flux_convective = np.zeros([1])
        self.heat_flux_radiative = np.zeros([1])
        self.heat_load = np.zeros([1])
        self.pressure = np.zeros([1])

        self.PICA_thickness = 0.0

        self.load = 0.0
        self.chute_area = 0.0
        self.chute_cost = 0.0
        self.line_weight = 0.0
        self.line_cost = 0.0

    @property
    def mass(self):
        return self.glider_mass + self.heatshield_mass + self.chute_weight + self.line_weight

    def ballistic_entry(self, acc):
        drag = atmospheric_model.Entry_model.CapsuleDrag(self.diameter, 1.125, np.deg2rad(20), self.lat, self.lon, acc)

        capsule_aero_coefficient_setting = \
            environment_setup.aerodynamic_coefficients.custom_aerodynamic_force_coefficients(
                force_coefficient_function=drag.drag_coefficient,
                reference_area=np.pi * (self.diameter / 2) ** 2,
                independent_variable_names=[
                    environment.AerodynamicCoefficientsIndependentVariables.mach_number_dependent,
                    environment.AerodynamicCoefficientsIndependentVariables.altitude_dependent])

        termination_mach_setting = propagation_setup.propagator.dependent_variable_termination(
            dependent_variable_settings=propagation_setup.dependent_variable.mach_number("Capsule", "Uranus"),
            limit_value=2.7,
            use_as_lower_limit=True)
        self.ballistic_dep_vars = atmospheric_model.Entry_model.entry_sim(self.mass, capsule_aero_coefficient_setting,
                                                                          self.alt, self.lat, self.lon, self.speed,
                                                                          self.flight_path_angle, self.heading_angle,
                                                                          [termination_mach_setting], acc=acc)
        self.sonic_alt = self.ballistic_dep_vars[-1, -1]
        self.sonic_lat = self.ballistic_dep_vars[-1, 2]
        self.sonic_lon = self.ballistic_dep_vars[-1, 3]
        self.sonic_speed = self.ballistic_dep_vars[-1, 4]
        self.sonic_flight_path_angle = self.ballistic_dep_vars[-1, -3]
        self.sonic_heading_angle = self.ballistic_dep_vars[-1, -2]

    def sonic_parachute_decent(self, acc):
        termination_alt_settings = propagation_setup.propagator.dependent_variable_termination(
            dependent_variable_settings=propagation_setup.dependent_variable.altitude("Capsule", "Uranus"),
            limit_value=self.target_altitude,
            use_as_lower_limit=True)  # """

        termination_mach_settings = propagation_setup.propagator.dependent_variable_termination(
            dependent_variable_settings=propagation_setup.dependent_variable.mach_number("Capsule", "Uranus"),
            limit_value=0.9,
            use_as_lower_limit=True)  # """

        aero_coefficient_setting = environment_setup.aerodynamic_coefficients.constant(
            np.pi * (self.parachute_diameter / 2) ** 2,
            constant_force_coefficient=[self.parachute_c_ds[0], 0.0, 0.0],
            are_coefficients_in_aerodynamic_frame=True,
            are_coefficients_in_negative_axis_direction=True)  # """

        self.sonic_dep_vars = atmospheric_model.Entry_model.entry_sim(self.mass - self.heatshield_mass,
                                                                      aero_coefficient_setting, self.sonic_alt,
                                                                      self.sonic_lat, self.sonic_lon, self.sonic_speed,
                                                                      self.sonic_flight_path_angle,
                                                                      self.sonic_heading_angle,
                                                                      [termination_alt_settings,
                                                                       termination_mach_settings], acc=acc)

        self.parachute_alt = self.sonic_dep_vars[-1, -1]
        self.parachute_lat = self.sonic_dep_vars[-1, 2]
        self.parachute_lon = self.sonic_dep_vars[-1, 3]
        self.parachute_speed = self.sonic_dep_vars[-1, 4]
        self.parachute_flight_path_angle = self.sonic_dep_vars[-1, -3]
        self.parachute_heading_angle = self.sonic_dep_vars[-1, -2]

    def parachute_decent(self, acc):
        termination_alt_settings = propagation_setup.propagator.dependent_variable_termination(
            dependent_variable_settings=propagation_setup.dependent_variable.altitude("Capsule", "Uranus"),
            limit_value=self.target_altitude,
            use_as_lower_limit=True)  # """

        aero_coefficient_setting = environment_setup.aerodynamic_coefficients.constant(
            np.pi * (self.parachute_diameter / 2) ** 2,
            constant_force_coefficient=[self.parachute_c_ds[1], 0.0, 0.0],
            are_coefficients_in_aerodynamic_frame=True,
            are_coefficients_in_negative_axis_direction=True)  # """

        self.parachute_dep_vars = atmospheric_model.Entry_model.entry_sim(self.mass - self.heatshield_mass,
                                                                          aero_coefficient_setting, self.parachute_alt,
                                                                          self.parachute_lat, self.parachute_lon,
                                                                          self.parachute_speed,
                                                                          self.parachute_flight_path_angle,
                                                                          self.parachute_heading_angle,
                                                                          [termination_alt_settings], acc=acc)

    def simulate_decent(self, acc):
        self.ballistic_entry(acc)
        self.sonic_parachute_decent(acc)
        self.parachute_decent(acc)

    def calculate_entry_heating(self, acc):
        drag = atmospheric_model.Entry_model.CapsuleDrag(self.diameter, 1.125, np.deg2rad(20), self.lat, self.lon, acc)
        self.gram.altitudes = self.ballistic_dep_vars[:, 1] / 1000
        self.gram.time = self.ballistic_dep_vars[:, 0]
        self.gram.lat = np.rad2deg(self.ballistic_dep_vars[:, 2])
        self.gram.long = np.rad2deg((self.ballistic_dep_vars[:, 3] + 2 * np.pi) % (2 * np.pi))
        self.gram.run()

        k = 1 / (np.asarray(self.gram.data.H2mass_pct) / 0.0395 + np.asarray(self.gram.data.Hemass_pct) / 0.0797)
        self.heat_flux_convective = k * self.ballistic_dep_vars[:, 4] ** 3 * (np.asarray(self.gram.data.Density_kgm3) /
                                                                              (np.pi * (self.diameter / 2) ** 2)) ** 0.2
        self.heat_flux_radiative = 9.7632379 ** (-40) * self.diameter ** (-0.17905) * \
                                   np.asarray(self.gram.data.Density_kgm3) ** 1.763827469 * \
                                   self.ballistic_dep_vars[:, 4] ** 10.993852

        q_func = sp.interpolate.interp1d(self.ballistic_dep_vars[:, 0], self.heat_flux_convective +
                                         self.heat_flux_radiative)
        self.heat_load = sp.integrate.quad(lambda x: q_func(x), self.ballistic_dep_vars[0, 0],
                                           self.ballistic_dep_vars[-1, 0])[0]

        drag.gamma = np.asarray(self.gram.data.SpecificHeatRatio)
        drag.mach = self.ballistic_dep_vars[:, 6]
        self.pressure = drag.p_0_stag() * np.asarray(self.gram.data.Pressure_Pa)

    def calculate_heatshield_thickness(self):
        heat_load = self.heat_load / 10000  # convert to J/cm**2
        interface_velocity = self.speed / 1000  # convert to km/s

        self.PICA_thickness = 1.8686 * (heat_load / (interface_velocity ** 2)) ** 0.1879 / 100

    def calculate_parachute(self):
        self.parachute_diameter *= self.target_time / (self.parachute_dep_vars[-1, 0] - self.parachute_dep_vars[
            self.parachute_dep_vars[:, 1] <= 0.0][0, 0])

        self.load = self.ballistic_dep_vars[-1, -4] * np.pi * (self.parachute_diameter / 2) ** 2 * \
                    self.parachute_c_ds[0] / (self.mass - self.heatshield_mass) * self.shock_load_factors

        self.chute_area = 2 * np.pi * (self.parachute_diameter / 2) ** 2
        self.chute_weight = self.chute_area * 35 / 1000
        self.chute_cost = self.chute_area * 15 * 0.92

        line_count = np.ceil(self.load * (self.mass - self.heatshield_mass) / _line_safe_strength)

        distance = self.parachute_diameter * 1.5
        line_length = np.sqrt((self.parachute_diameter / 2) ** 2 + distance ** 2) + np.pi / 2 * (self.parachute_diameter
                                                                                                 / 2) ** 2
        self.line_weight = line_length * line_count * _line_weight_per_meter
        self.line_cost = line_length * line_count * _line_cost_per_meter


if __name__ == "__main__":
    test = Aeroshell(2*10**6, 4 * 3600, 162.1, 3, [0.45, 0.4], 2.0, 0.0, 0.0, 250, 3.03327727e+07, 5.45941114e-01,
                     -2.33346601e-02, 2.65992642e+04, -5.91036848e-01, -2.96367147e+00)
    test.simulate_decent(10)
    test.calculate_parachute()
    print(test.chute_weight + test.line_weight)
