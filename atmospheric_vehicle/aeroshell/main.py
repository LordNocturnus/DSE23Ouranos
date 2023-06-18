from tudatpy.kernel.numerical_simulation import environment_setup, environment, propagation_setup
from matplotlib import pyplot as plt
import numpy as np
import scipy as sp

import atmospheric_model.Entry_model
from atmospheric_model.GRAM import GRAM
from atmospheric_vehicle.aeroshell.heatshield import volume

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
        self.chute_weight = 0.0
        self.line_weight = 0.0
        self.struct_mass = 0.0
        self.insulator_weight = 0.0
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

        self.PICA_thickness = 0.0
        self.insulator_thickness = 0.0

        self.load = 0.0
        self.chute_area = 0.0
        self.chute_cost = 0.0
        self.line_cost = 0.0

    @property
    def mass(self):
        return self.glider_mass + self.heatshield_mass + self.chute_weight + self.line_weight + self.struct_mass + \
            self.insulator_weight

    @property
    def cost(self):
        return self.line_cost + self.chute_cost

    @property
    def radius2(self):
        return max(_radius2, (self.PICA_thickness + self.insulator_thickness) * 1.01)

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
        self.chute_weight = self.chute_area * 2 * 35 / 1000
        self.chute_cost = self.chute_area * 2 * 15 * 0.92

        self.load = self.ballistic_dep_vars[-1, -4] * self.chute_area * self.parachute_c_ds[0] / \
                    (self.mass - self.heatshield_mass) * self.shock_load_factors

        line_count = np.ceil(self.load * (self.mass - self.heatshield_mass) / _line_safe_strength)

        distance = self.parachute_diameter * 2.5
        line_length = np.sqrt((self.parachute_diameter / 2) ** 2 + distance ** 2) + np.pi / 2 * (self.parachute_diameter
                                                                                                 / 2) ** 2
        self.line_weight = line_length * line_count * _line_weight_per_meter
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
        self.insulator_weight = volume(self.diameter, _radius1, self.radius2, 0.0, self.PICA_thickness) * _rho_insulator

    def backshell_geometry(self, peak_load, load_entry):
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
        r_top_small = np.sqrt(max(max(max(self.sonic_dep_vars[:, 5]), max(self.parachute_dep_vars[:, 5])) * 1.1 /
                                  _sigma_y_backshell, np.pi * _r_parachute ** 2) / np.pi)
        r_top_big = r_top_small / (1 - _taper_ratio_top)

        # Calculate the small radius of the bottom truncated cone. Limit case is the radius of the thermal shield
        r_bottom_small = (self.diameter / 2) * (1 - _taper_ratio_bottom)

        # Calculate the angle at the base of the truncated cones
        a_top = self.angle_cone(r_top_big, r_top_small, _h_parachute)
        a_bottom = self.angle_cone((self.diameter / 2), r_bottom_small, _h_folded_wings)

        # Calculate thickness based on pressure loads. Use personalised formula
        t_top = max(self.t_pressure(r_top_big, a_top) * 1.2, 1 * 10 ** -3)
        t_bottom = max(self.t_pressure((self.diameter / 2), a_bottom) * 1.2, 1 * 10 ** -3)

        # Calculate thickness based on pressure loads. Use traditional formula for thin walled cylinders
        # t_top = t_hoop(p_load, r_top_big, sigma_y)
        # t_bottom = t_hoop(p_load, (self.diameter / 2), sigma_y)

        # Check buckling
        I_top = 1 / 12 * t_top * (_h_parachute / np.cos(a_top)) ** 3
        A_top = np.pi * (r_top_big + r_top_small) * _h_parachute / np.cos(a_top)
        buck_top = self.buckling(I_top, _h_parachute / np.cos(a_top), p_load * A_top)
        I_bottom = 1 / 12 * t_bottom * _h_folded_wings / np.cos(a_bottom)
        A_bottom = np.pi * ((self.diameter / 2) + r_bottom_small) * _h_folded_wings / np.cos(a_bottom)
        buck_bottom = self.buckling(I_bottom, _h_folded_wings / np.cos(a_bottom), p_load * A_bottom)

        if buck_bottom and buck_top:
            # Calculate the volume of the thin walled structure by subtraction
            volume_top = volume_truncated(r_top_big, r_top_small, _h_parachute) - volume_truncated(r_top_big - t_top,
                                                                                                  r_top_small - t_top,
                                                                                                  _h_parachute - 2 * t_top)
            volume_bottom = volume_truncated((self.diameter / 2), r_bottom_small, _h_folded_wings) - volume_truncated(
                (self.diameter / 2) - t_bottom, r_bottom_small - t_bottom, _h_folded_wings - 2 * t_bottom)

            # Calculate backshell mass
            mass_backshell = (volume_top + volume_bottom) * _rho_backshell
            t_bottom_shell = max(bending_bottom(load_entry, (self.diameter / 2) * 2, _sigma_y_backshell), 1 * 10 ** -3,
                                 bending_pressure(p_load, (self.diameter / 2) * 2, _sigma_y_backshell))
            I_backshell = np.pi * ((self.diameter / 2) * 2) ** 4 / 64
            print(buckling(_E_backshell, (self.diameter / 2) * 2, I_backshell, peak_load * np.pi * (self.diameter / 2) ** 2))
            mass_bottom_shell = volume(t_heatshield, t_bottom_shell, (self.diameter / 2) * 2) * _rho_backshell
            return (mass_backshell + mass_bottom_shell) * 1.2, t_top, t_bottom, t_bottom_shell
        else:
            print(f'Buckling requirements is not satisfied')

    def angle_cone(self, r_big, r_small, height):
        """
        Function that determined the angle between the base and the hypotenuse of the cone
        :param r_big: Bottom radius of the cone
        :param r_small: Top radius of the cone
        :param height: height of the cone
        :return: angle between hypotenuse and base
        """
        if height == 0:
            raise ZeroDivisionError(f'The height has to be a positive value')
        if r_big <= 0 or r_small <= 0 or height < 0:
            raise ValueError(f'The dimensions of the capsule should be positive values')
        return np.arctan((r_big - r_small) / height)

    def t_pressure(self, r_big, a):
        """
        Function to determine the thickness required to withstand pressure loads during entry. Formulas were derived from
        the following website: https://123sanat.com/d/catalogue/ASME-VIII-_8_-div.1-2019.pdf
        :param p: External pressure
        :param r_big: Top radius of cone
        :param a: angle between hypotenuse and base
        :param sigma_y: Yield stress of selected material
        :return: Required minimum thickness
        """
        if a == np.pi / 2:
            raise ZeroDivisionError(f'Capsule side angle equals 90 degrees, change the idealised geometry')
        elif _sigma_y_backshell == 0.6 * abs(self.pressure):
            raise ZeroDivisionError(f'Formula no longer valid due to material properties')
        return abs(self.pressure) * r_big * 2 / (2 * np.cos(a) * (_sigma_y_backshell - 0.6 * abs(self.pressure)))

    def buckling(self, l, I, p_load):
        """
        Function to calculate the buckling limit load
        :param E: Young's Modulus
        :param l: Length of plate
        :param I: Moment of inertia plate
        :param p_load: Pressure loads during entry
        :return: bool if buckling criteria is respected
        """
        if l == 0:
            raise ZeroDivisionError(f'The length of the sheet should be higher than 0')
        elif l < 0 or I <= 0 or _E_backshell <= 0:
            raise ValueError(f'The geometrical and material properties of the sheet should be a positive value')
        return abs(p_load) <= np.pi ** 2 * _E_backshell * I / (l ** 2)

    def itterate(self, struct_mass, acc=1, steps=5):
        self.struct_mass = struct_mass
        for i in range(0, steps):
            self.simulate_decent(acc)
            self.calculate_entry_heating(acc)
            self.calculate_heatshield_thickness()
            self.calculate_insulator_shell_thickness()
            self.calculate_parachute()

            self.calculate_heatshield_weight()
            self.calculate_insulator_weight()

            print(f"Finished Itteration {i+1}")


if __name__ == "__main__":
    test = Aeroshell(2*10**6, 4 * 3600, 162.1, 3, [0.45, 0.4], 1.3, 3.03327727e+07, 5.45941114e-01, -2.33346601e-02,
                     2.65992642e+04, -5.91036848e-01, -2.96367147e+00)

    test.itterate(46, 1, 1)

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
