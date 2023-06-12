"""
Script for the orbiter design tool
"""
import Orbiter.propulsion as prop


class Orb:

    def __init__(self):
        self.mass = ...  # Orbiter mass
        self.mass_AV = ...  # Atmospheric vehicle mass (import from AV class)
        self.mass_combined = self.mass + self.mass_AV  # Mass of combined systems
        self.deltaV_transfer = ...  # Combined systems deltaV
        self.deltaV_insertion = ...  # Delta V after splitting up at Uranus
        self.Isp = ...  # Isp of the orbiter thrusters
        self.T = ...  # Orbiter thrust
        self.prop_mass = prop.mass_prop(self.mass, self.deltaV_insertion, self.mass_combined, self.deltaV_transfer, g, self.Isp)
        self.burn_transfer = prop.burntimecombined(self.T, self.mass_combined, self.deltaV_transfer)
        self.burn_insertion = prop.burntimeorbiter(self.T, self.mass, self.deltaV_insertion)
        


if __name__ == "__main__":
    g = 9.81
