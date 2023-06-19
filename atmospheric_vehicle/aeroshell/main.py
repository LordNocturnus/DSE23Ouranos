from tudatpy.kernel.numerical_simulation import environment_setup, environment, propagation_setup
from matplotlib import pyplot as plt
import numpy as np
import scipy as sp

import atmospheric_model.Entry_model
from atmospheric_model.GRAM import GRAM
from atmospheric_model.Entry_model import CapsuleDrag
from atmospheric_vehicle.aeroshell.heatshield import volume
from atmospheric_vehicle.aeroshell.structure import angle_cone, buckling, volume_truncated, t_pressure


_PICA_density = 0.352 * 100**3 / 1000  # 0.352 - 0.701 g/cm^3

_rho_backshell = 49.7  # https://journals.sagepub.com/doi/pdf/10.1177/0021998313499949 (Q2 selected because strongest)
_sigma_y_backshell = 450 * 10 ** 6  # https://journals.sagepub.com/doi/pdf/10.1177/0021998313499949 (Q2 selected because strongest)
_E_backshell = 130 * 10 ** 9  # https://www.azom.com/article.aspx?ArticleID=6618
_alpha_backshell = 23.6 * 10 ** -6
_back_cost_kg = 40  # Based on research https://www.easycomposites.co.uk/3mm-aluminium-honeycomb
_cost_manu_all = 2063 * 0.92535

_alpha_insulator = 1.2 * 10 ** -7
_sigma_y_insulator = 790 * 10 ** 6
_E_insulator = 173 * 10 ** 9
_rho_insulator = 55

_radius1 = 1.125
_radius2 = 0.126
_bottom_angle = np.deg2rad(20)
_top_angle = np.deg2rad(59.73)
_h_folded_wings = 2
_h_parachute = 0.5
_r_parachute = 0.3
_taper_ratio_bottom = 0.44
_taper_ratio_top = 0.44

_line_safe_strength = 1300 * 4  # N with safety factor of 3
_line_weight_per_meter = 0.013  # kg/m
_line_cost_per_meter = 0.7  # €/m


class Aeroshell:

    def __init__(self, target_pressure, target_time, glider_mass, diameter, parachute_c_ds, shock_load_factors, alt,
                 lat, lon, speed, flight_path_angle, heading_angle):
        self.target_pressure = target_pressure
        self.target_time = target_time
        self.glider_mass = glider_mass
        self.heatshield_mass = 0.0
        self.chute_mass = 4.0
        self.line_mass = 15.0
        self.struct_mass = 0.0
        self.insulator_mass = 0.0
        self.bottom_shell_mass = 0.0
        self.mass_backshell = 0.0
        self.glider_offset_mass = 0.0
        self.diameter = diameter
        self.parachute_c_ds = parachute_c_ds
        self.shock_load_factors = shock_load_factors
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

        self.PICA_thickness = 0.10
        self.insulator_thickness = 0.01
        self.bottom_shell_thickness = 0.01

        self.t_top = 0.01
        self.t_bottom = 0.01
        self.r_top_small = 0.001
        self.r_top_big = self.diameter / 4
        self.r_bottom_small = self.diameter / 4
        self.a_top = np.pi / 2
        self.a_bottom = np.pi / 2

        self.load = 0.0
        self.chute_area = 0.0
        self.chute_cost = 0.0
        self.line_cost = 0.0

        self.cg_y = 0.0
        self.c_p = np.zeros(1)
        self.gram = GRAM()
        self.shield_pressure = CapsuleDrag(self.diameter, _radius1, _bottom_angle, self.lat, self.lon)

    @property
    def mass(self):
        return self.glider_mass + self.heatshield_mass + self.chute_mass + self.line_mass + self.struct_mass + \
               self.insulator_mass + self.bottom_shell_mass + self.mass_backshell + self.glider_offset_mass

    @property
    def cost(self):
        return self.line_cost + self.chute_cost

    @property
    def radius2(self):
        return max(_radius2, (self.PICA_thickness + self.insulator_thickness + self.bottom_shell_thickness) * 1.01)

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
            limit_value=2.5,
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
        self.heat_flux_convective_2 = 2004.2526*1/np.sqrt(self.diameter / 0.6091) * (
            np.asarray(self.gram.data.Density_kgm3) / 1.22522) ** 0.4334341 * (self.ballistic_dep_vars[:, 4] / 3048) ** \
                                      2.9978868 * 1000
        self.heat_flux_radiative = 9.7632379 ** (-40) * self.diameter ** (-0.17905) * \
                                   np.asarray(self.gram.data.Density_kgm3) ** 1.763827469 * \
                                   self.ballistic_dep_vars[:, 4] ** 10.993852 * 1000000

        q_func = sp.interpolate.interp1d(self.ballistic_dep_vars[:, 0], self.heat_flux_convective_2 +
                                         self.heat_flux_radiative)
        self.heat_load = sp.integrate.quad(lambda x: q_func(x), self.ballistic_dep_vars[0, 0],
                                           self.ballistic_dep_vars[-1, 0])[0]

        drag.gamma = np.asarray(self.gram.data.SpecificHeatRatio)
        drag.mach = self.ballistic_dep_vars[:, 6]
        self.pressure = drag.p_0_stag() * np.asarray(self.gram.data.Pressure_Pa)

    def calculate_heatshield_thickness(self):
        heat_load = self.heat_load / 10000  # convert to J/cm**2
        interface_velocity = self.speed / 1000  # convert to km/s

        self.PICA_thickness = 1.8686 * (heat_load / (interface_velocity ** 2)) ** 0.1879 / 100 * 1.117 * 1.5

    def calculate_heatshield_weight(self):
        self.heatshield_mass = volume(self.diameter, _radius1, self.radius2, 0.0, self.PICA_thickness) * _PICA_density

    def calculate_parachute(self):
        self.parachute_diameter *= self.target_time / (self.parachute_dep_vars[-1, 0] - self.parachute_dep_vars[
            self.parachute_dep_vars[:, 1] <= 0.0][0, 0])

        self.chute_area = np.pi * (self.parachute_diameter / 2) ** 2
        self.chute_mass = self.chute_area * 2 * 35 / 1000
        self.chute_cost = self.chute_area * 2 * 15 * 0.92

        self.load = self.ballistic_dep_vars[-1, -4] * self.chute_area * self.parachute_c_ds[0] / \
                    (self.mass - self.heatshield_mass) * self.shock_load_factors

        line_count = np.ceil(self.load * (self.mass - self.heatshield_mass) / _line_safe_strength)

        distance = self.parachute_diameter * 2.5
        line_length = np.sqrt((self.parachute_diameter / 2) ** 2 + distance ** 2) + np.pi / 2 * (self.parachute_diameter
                                                                                                 / 2) ** 2
        self.line_mass = line_length * line_count * _line_weight_per_meter
        self.line_cost = line_length * line_count * _line_cost_per_meter

    def bending_bottom(self):
        """
        Method to calculate the thickness necessary to withstand entry loads
        :param load_entry: Entry loads
        :param l_thermal_shield: Size of the thermal shield
        :param sigma_y: Yield strength of the selected material
        :return: Required minimum thickness to withstand entry loads
        """
        return np.sqrt(max(self.ballistic_dep_vars[:, 5]) * self.mass * self.diameter * 3 / (self.diameter *
                                                                                             _sigma_y_insulator))

    def bending_pressure(self):
        return np.sqrt((3 * max(self.pressure) * self.diameter ** 2 / (2 * self.diameter * _sigma_y_insulator)))

    def calculate_insulator_shell_thickness(self):
        self.insulator_thickness = 0.8 * max(self.bending_bottom(), self.bending_pressure())

    def calculate_insulator_weight(self):
        self.insulator_mass = volume(self.diameter, _radius1, self.radius2, self.PICA_thickness,
                                     self.insulator_thickness) * _rho_insulator

    def calculate_bottom_shell_thickness(self):
        self.bottom_shell_thickness = max(self.bending_bottom(), 1 * 10 ** -3, self.bending_pressure())

    def calculate_bottom_shell_weight(self):
        self.bottom_shell_mass = volume(self.diameter, _radius1, self.radius2, self.PICA_thickness +
                                        self.insulator_thickness, self.bottom_shell_thickness) * _rho_backshell

    def backshell_geometry(self):
        """
        Function that determines all geometrical properties and mass of the backshell. Multiple things are computed,
        during integration, it should be decided what is the most useful return, if one is actually needed. As of now the
        function prints the main geometrical characteristics and mass of the backshell. The backshell is assumed to be
        formed by two truncated cones, a bottom one that has the bottom radius equal to the thermal shield radius and a
        top one that has the top radius equal to the parachute radius (in folded configuration). Taper ratio is selected
        manually to determine the transition from one truncated cone to the other.
        :param peak_load: Peak parachute load experienced during entry (from decelerator analysis)
        :param load_entry: Maximum entry loads experienced by the capsule
        :param p_load: Pressure loads experienced during entry
        :param r_thermal: Thermal shield radius
        :param h_folded_wings: Height of the glider in the folded configuration
        :return: Backshell mass !!!! ADD MASS OF THERMAL SHIELD SUPPORTING STRUCTURE !!!!! (Use volume function from heatshild code)
        """

        # Calculate the big and small radius of the top truncated cone. Limit case is the parachute radius
        self.r_top_small = np.sqrt(max(self.load * (self.mass - self.heatshield_mass) * 1.1 /
                                  _sigma_y_backshell, np.pi * _r_parachute ** 2) / np.pi)
        self.r_top_big = self.r_top_small / (1 - _taper_ratio_top)

        # Calculate the small radius of the bottom truncated cone. Limit case is the radius of the thermal shield
        self.r_bottom_small = (self.diameter / 2) * (1 - _taper_ratio_bottom)

        # Calculate the angle at the base of the truncated cones
        self.a_top = angle_cone(self.r_top_big, self.r_top_small, _h_parachute)
        self.a_bottom = angle_cone((self.diameter / 2), self.r_bottom_small, _h_folded_wings)

        # Calculate thickness based on pressure loads. Use personalised formula
        self.t_top = max(t_pressure(0.1 * 10**5, self.r_top_big, self.a_top, _sigma_y_backshell) * 1.2, 1 * 10 ** -3)
        self.t_bottom = max(t_pressure(0.1 * 10**5, (self.diameter / 2), self.a_bottom, _sigma_y_backshell) * 1.2,
                            1 * 10 ** -3)

        # Calculate thickness based on pressure loads. Use traditional formula for thin walled cylinders
        # t_top = t_hoop(p_load, r_top_big, sigma_y)
        # t_bottom = t_hoop(p_load, (self.diameter / 2), sigma_y)

        # Check buckling
        I_top = 1 / 12 * self.t_top * (_h_parachute / np.cos(self.a_top)) ** 3
        A_top = np.pi * (self.r_top_big + self.r_top_small) * _h_parachute / np.cos(self.a_top)
        buck_top = buckling(_E_backshell, I_top, _h_parachute / np.cos(self.a_top), 0.1 * 10**5 * A_top)
        I_bottom = 1 / 12 * self.t_bottom * _h_folded_wings / np.cos(self.a_bottom)
        A_bottom = np.pi * ((self.diameter / 2) + self.r_bottom_small) * _h_folded_wings / np.cos(self.a_bottom)
        buck_bottom = buckling(_E_backshell, I_bottom, _h_folded_wings / np.cos(self.a_bottom), 0.1 * 10**5 * A_bottom)
        # Calculate the volume of the thin walled structure by subtraction
        volume_top = volume_truncated(self.r_top_big, self.r_top_small, _h_parachute) - \
                     volume_truncated(self.r_top_big - self.t_top, self.r_top_small - self.t_top,
                                      _h_parachute - 2 * self.t_top)
        volume_bottom = volume_truncated((self.diameter / 2), self.r_bottom_small, _h_folded_wings) - \
                        volume_truncated((self.diameter / 2) - self.t_bottom, self.r_bottom_small - self.t_bottom,
                                         _h_folded_wings - 2 * self.t_bottom)
        self.mass_backshell = (volume_top + volume_bottom) * _rho_backshell
        I_backshell = np.pi * ((self.diameter / 2) * 2) ** 4 / 64
        buckling_shield = buckling(_E_backshell, (self.diameter / 2) * 2, I_backshell,
                 self.load * (self.mass - self.heatshield_mass) * np.pi * (self.diameter / 2) ** 2)

        if not buck_bottom or not buck_top or not buckling_shield:
            print(f'Buckling requirements is not satisfied')

    def calculate_cg(self, glider_cg_x):
        self.glider_offset_mass = self.glider_mass * glider_cg_x / self.diameter
        self.r_offset = - _radius1 + _radius1 * (1 - np.cos(_bottom_angle)) + (self.diameter / 2 - self.radius2 *
                                                                               np.cos(np.pi / 2 - _bottom_angle) -
                                                                               _radius1 * np.sin(_bottom_angle)) * \
                        np.tan(_bottom_angle) + self.radius2 * np.sin(np.pi / 2 - _bottom_angle)
        dis_chute = _h_folded_wings + _h_parachute / 2


        bottom = self.chute_mass + self.line_mass + self.glider_mass + self.glider_offset_mass
        top = self.chute_mass * dis_chute + self.line_mass * dis_chute + self.glider_mass * 0.0 + \
              self.glider_offset_mass * 0.0

        self.cg_y = top / bottom

    def c_p_alt(self, mach, alt, p_inf):
        # lower plate
        print(np.sin(self.a_bottom), self.diameter / 4 * _h_folded_wings ** 2, 1 / (3 * np.tan(self.a_bottom)) *
              _h_folded_wings ** 3)
        c_p_top = p_inf * np.sin(self.a_bottom) * (-self.diameter / 4 * _h_folded_wings ** 2 + 1 /
                                                   (3 * np.tan(self.a_bottom)) * _h_folded_wings ** 3)
        c_p_bot = p_inf * np.sin(self.a_bottom) * (-self.diameter / 2 * _h_folded_wings + 1 /
                                                   (2 * np.tan(self.a_bottom)) * _h_folded_wings ** 2)
        # upper plate
        c_p_top += p_inf * np.sin(self.a_bottom) * (-(self.r_top_big / 2 + _h_folded_wings / (2 * np.tan(self.a_top))) *
                                                    ((_h_folded_wings + _h_parachute) ** 2 - _h_folded_wings ** 2) +
                                                    1 / (3 * np.tan(self.a_top)) * ((_h_folded_wings +
                                                                                     _h_parachute) ** 3 -
                                                                                    _h_folded_wings ** 3))
        c_p_bot += p_inf * np.sin(self.a_bottom) * (-(self.r_top_big + _h_folded_wings / np.tan(self.a_top)) *
                                                    _h_parachute + 1 / (2 * np.tan(self.a_top)) *
                                                    ((_h_folded_wings + _h_parachute) ** 2 - _h_folded_wings ** 2))
        t_shield, b_shield = self.shield_pressure.c_p(mach, alt, p_inf)
        print(t_shield, b_shield, t_shield / b_shield)
        print(c_p_top, c_p_bot, c_p_top / c_p_bot)
        return (c_p_top + t_shield) / (c_p_bot + b_shield)

    def calculate_c_p(self):
        self.gram.altitudes = self.ballistic_dep_vars[:, 1] / 1000
        self.gram.time = np.zeros_like(self.ballistic_dep_vars[:, 1])
        self.gram.lat = np.rad2deg(self.ballistic_dep_vars[:, 2])
        self.gram.long = np.rad2deg((self.ballistic_dep_vars[:, 3] + 2 * np.pi) % (2 * np.pi))
        self.gram.run()
        self.c_p = np.zeros(len(self.ballistic_dep_vars))
        for i, v in enumerate(self.ballistic_dep_vars):
            print(v[0])
            self.c_p[i] = self.c_p_alt(v[6], v[1], self.gram.data.Pressure_Pa[i])

    def itterate(self, acc=1, steps=5):
        for i in range(0, steps):
            self.simulate_decent(acc)
            self.calculate_entry_heating(acc)
            self.calculate_heatshield_thickness()
            self.calculate_insulator_shell_thickness()
            self.calculate_bottom_shell_thickness()
            self.calculate_parachute()

            self.calculate_heatshield_weight()
            self.calculate_insulator_weight()
            self.calculate_bottom_shell_weight()
            self.backshell_geometry()

            print(f"Finished Itteration {i+1}")


if __name__ == "__main__":
    test = Aeroshell(2*10**6, 4 * 3600, 162.1, 3, [0.45, 0.4], 1.3, 3.03327727e+07 - 4000000, 5.45941114e-01,
                     -2.33346601e-02, 2.65992642e+04, -5.91036848e-01, -2.96367147e+00)

    test.itterate(1, 1)

    test.calculate_cg(0.36)
    test.calculate_c_p()

    plt.plot(test.ballistic_dep_vars[:, 0], test.c_p)
    plt.show()


    """test.itterate(46, 1, 1)

    plt.plot(test.ballistic_dep_vars[:, 0] - 6000, test.ballistic_dep_vars[:, 1] / 1000, label="Ballistic")
    plt.plot(test.sonic_dep_vars[:, 0] - 6000 + test.ballistic_dep_vars[-1, 0] - 6000,
             test.sonic_dep_vars[:, 1] / 1000, label="Supersonic Parachute")
    plt.plot(test.parachute_dep_vars[test.parachute_dep_vars[:, 1] >= 50000][:, 0] - 6000 +
             test.sonic_dep_vars[-1, 0] - 6000 + test.ballistic_dep_vars[-1, 0] - 6000,
             test.parachute_dep_vars[test.parachute_dep_vars[:, 1] >= 50000][:, 1] / 1000, label="Parachute")
    plt.grid()
    plt.legend()
    plt.xlabel("Time [s]")
    plt.ylabel("Altitude [km]")
    plt.savefig("entry_alt.pdf", dpi='figure', format="pdf", metadata=None,
        bbox_inches=None, pad_inches=0.0,
        facecolor='auto', edgecolor='auto',
        backend=None)
    plt.show()

    plt.semilogy(test.ballistic_dep_vars[:, 4], test.ballistic_dep_vars[:, 1] / 1000, label="Ballistic")
    plt.semilogy(test.sonic_dep_vars[:, 4], test.sonic_dep_vars[:, 1] / 1000, label="Supersonic Parachute")
    plt.semilogy(test.parachute_dep_vars[test.parachute_dep_vars[:, 1] >= 50000][:, 4],
             test.parachute_dep_vars[test.parachute_dep_vars[:, 1] >= 50000][:, 1] / 1000, label="Parachute")
    plt.grid()
    plt.legend()
    plt.ylabel("Logarithmic Altitude [km]")
    plt.xlabel("Velocity [m/s]")
    plt.savefig("entry_vel.pdf", dpi='figure', format="pdf", metadata=None,
                bbox_inches=None, pad_inches=0.0,
                facecolor='auto', edgecolor='auto',
                backend=None)
    plt.show()

    plt.semilogy(test.ballistic_dep_vars[:, 5] / 9.80665, test.ballistic_dep_vars[:, 1] / 1000, label="Ballistic")
    plt.semilogy(test.sonic_dep_vars[:, 5] / 9.80665, test.sonic_dep_vars[:, 1] / 1000, label="Supersonic Parachute")
    plt.semilogy(test.parachute_dep_vars[test.parachute_dep_vars[:, 1] >= 50000][:, 5] / 9.80665,
             test.parachute_dep_vars[test.parachute_dep_vars[:, 1] >= 50000][:, 1] / 1000, label="Parachute")
    plt.grid()
    plt.legend()
    plt.ylabel("Logarithmic Altitude [km]")
    plt.xlabel("Deceleration [g]")
    plt.savefig("entry_dec.pdf", dpi='figure', format="pdf", metadata=None,
                bbox_inches=None, pad_inches=0.0,
                facecolor='auto', edgecolor='auto',
                backend=None)
    plt.show()

    plt.semilogy(test.ballistic_dep_vars[:, 6], test.ballistic_dep_vars[:, 1] / 1000, label="Ballistic")
    plt.semilogy(test.sonic_dep_vars[:, 6], test.sonic_dep_vars[:, 1] / 1000, label="Supersonic Parachute")
    plt.semilogy(test.parachute_dep_vars[test.parachute_dep_vars[:, 1] >= 50000][:, 6],
             test.parachute_dep_vars[test.parachute_dep_vars[:, 1] >= 50000][:, 1] / 1000, label="Parachute")
    plt.grid()
    plt.legend()
    plt.ylabel("Logarithmic Altitude [km]")
    plt.xlabel("Mach number [-]")
    plt.savefig("entry_mach.pdf", dpi='figure', format="pdf", metadata=None,
                bbox_inches=None, pad_inches=0.0,
                facecolor='auto', edgecolor='auto',
                backend=None)
    plt.show()

    #plt.semilogy(test.heat_flux_convective, test.ballistic_dep_vars[:, 1] / 1000, label="Convective heat flux")
    plt.semilogy(test.heat_flux_convective_2 / 1000000, test.ballistic_dep_vars[:, 1] / 1000, label="Convective heat flux")
    plt.semilogy(test.heat_flux_radiative / 1000000, test.ballistic_dep_vars[:, 1] / 1000, label="Radiative heat flux")
    plt.grid()
    plt.legend()
    plt.xlabel("Heat flux [MW/m²]")
    plt.ylabel("Logarithmic Altitude [km]")
    plt.savefig("entry_heatflux.pdf", dpi='figure', format="pdf", metadata=None,
                bbox_inches=None, pad_inches=0.0,
                facecolor='auto', edgecolor='auto',
                backend=None)
    plt.show()

    plt.plot(test.ballistic_dep_vars[:, 0] - 6000, test.ballistic_dep_vars[:, 3], label="Ballistic")
    plt.plot(test.sonic_dep_vars[:, 0] - 6000 + test.ballistic_dep_vars[-1, 0] - 6000,
             test.sonic_dep_vars[:, 3], label="Supersonic Parachute")
    plt.plot(test.parachute_dep_vars[test.parachute_dep_vars[:, 1] >= 50000][:, 0] - 6000 +
             test.sonic_dep_vars[-1, 0] - 6000 + test.ballistic_dep_vars[-1, 0] - 6000,
             test.parachute_dep_vars[test.parachute_dep_vars[:, 1] >= 50000][:, 3], label="Parachute")
    plt.grid()
    plt.legend()
    plt.xlabel("Time [s]")
    plt.ylabel("Longitude [rad]")
    plt.savefig("entry_mach.pdf", dpi='figure', format="pdf", metadata=None,
                bbox_inches=None, pad_inches=0.0,
                facecolor='auto', edgecolor='auto',
                backend=None)
    plt.show()#"""
