"""
Script for the orbiter design tool
"""
import Orbiter.propulsion as prop
import Orbiter.structure as str


class Orb:

    def __init__(self):
        self.mass = ...  # Orbiter dry mass
        self.mixture_ratio = 1.65
        self.mass_AV = ...  # Atmospheric vehicle mass (import from AV class)
        self.mass_combined = self.mass + self.mass_AV  # Mass of combined systems
        self.deltaV_transfer = ...  # Combined systems deltaV
        self.deltaV_insertion = ...  # Delta V after splitting up at Uranus
        self.Isp = ...  # Isp of the orbiter thrusters
        self.T = ...  # Orbiter thrust
        self.prop_mass = prop.mass_prop(self.mass, self.deltaV_insertion, self.mass_combined, self.deltaV_transfer, g, self.Isp)
        self.burn_transfer = prop.burntimecombined(self.T, self.mass_combined, self.deltaV_transfer)
        self.burn_insertion = prop.burntimeorbiter(self.T, self.mass, self.deltaV_insertion)
        self.m_fuel = self.prop_mass / (1 + self.mixture_ratio)
        self.m_ox = self.prop_mass - self.m_fuel
        self.prop_properties = [(1 * 10 ** 6, self.m_ox, 1431), (1 * 10 ** 6, self.m_fuel, 874)]
        self.material = [4430, 880 * 10**6, 970 * 10**6, 113.8 * 10**9]
        self.mass_iteration()
        self.dry_mass_final = self.mass_combined + self.m_structure
        self.wet_mass_final = self.dry_mass_final + self.prop_mass

        self.f_lat, self.f_ax = str.natural_frequency(self.l_tanks, self.r_tanks, self.material, self.mass_final, f_ax_min, f_lat_min)



    def mass_iteration(self):
        m_structure = 0.2 * (self.mass_combined)
        diff = 1000
        while diff >= 1:
            self.prop_mass = prop.mass_prop(self.mass + m_structure, self.deltaV_insertion, self.mass_combined, self.deltaV_transfer,
                                            g, self.Isp)
            self.wet_mass = self.mass_combined + m_structure + self.prop_mass
            self.l_tanks, self.r_tanks, self.m_structure = str.final_architecture(self.material, self.prop_properties,
                                                                                  margin, self.wet_mass)
            diff = abs(m_structure - self.m_structure)
            m_structure = self.m_structure


if __name__ == "__main__":
    g = 9.81
    R = 8.314
    margin = 0.2

    acc_axial_tension = 6 * 9.81  # https://www.spacex.com/media/falcon-users-guide-2021-09.pdf
    acc_axial_compr = 2 * 9.81  # https://www.spacex.com/media/falcon-users-guide-2021-09.pdf
    acc_lateral = 2 * 9.81  # Same for compression and tension (https://www.spacex.com/media/falcon-users-guide-2021-09.pdf)
    acc_shock = 1000 * 9.81 # https://www.spacex.com/media/falcon-users-guide-2021-09.pdf
    f_lat_min = 10  # Falcon Heavy user manual
    f_ax_min = 25  # Falcon Heavy user manual

    # Launcher Constraints
    d_fairing = 5.2  # https://www.spacex.com/vehicles/falcon-heavy/
    h_fairing = 13.1  # https://www.spacex.com/vehicles/falcon-heavy/
