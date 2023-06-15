import numpy as np

# combined system parameter:
# m_combined = 2898  # total dry mass of both systems in kg
# dV_combined = 1000  # delta V required for the combined system in m/s

# Orbiter system parameters:
# m_orbiter = 2427  # dry mass of only the orbiter in kg
# dV_orbiter = 1000  # delta V required for the orbiter alone in m/s

# Constants
g = 9.81  # gravitational constant in m/s^2

# Engine parameters:
I_sp = 321  # specific impulse in sec
T = 445  # thrust in N
m_flow = 1  # fuel mass flow in kg/s
cost_kg_fuel = 146.54 * 1.1655 * 0.951 / 0.45359237  # FY22 € converted and cost/kg taken from https://www.dla.mil/Portals/104/Documents/Energy/Standard%20Prices/Aerospace%20Prices/E_2017Oct1AerospaceStandardPrices_170913.pdf?ver=2017-09-13-145335-477
cost_kg_ox = 151.93 * 1.1655 * 0.951 / 0.45359237  # FY22 € converted and cost/kg taken from https://www.dla.mil/Portals/104/Documents/Energy/Standard%20Prices/Aerospace%20Prices/E_2017Oct1AerospaceStandardPrices_170913.pdf?ver=2017-09-13-145335-477
m_engine = 4.5
# Calculate required propellant mass:
def mass_prop(m_orbiter, dV_orbiter, m_combined, dV_combined, Isp):
    m_prop_orbiter = m_orbiter * (np.exp(dV_orbiter / (g * Isp)) - 1)
    m_prop_combined = m_combined * (np.exp(dV_combined / (g * Isp)) - 1)
    m_prop = m_prop_orbiter + m_prop_combined
    return m_prop_orbiter, m_prop_combined, m_prop


# Calculate required burn time for combined system manoeuvers in hours:
def burntimecombined(T, m_f, dV_combined, m_orbiter, m_combined, dV_orbiter, Isp):
    a = T / (m_f + mass_prop(m_orbiter, dV_orbiter, m_combined, dV_combined, Isp)[2])
    t_b = dV_combined / a
    t_b = t_b / 60 / 60
    return t_b


# Calculate required burn time for orbiter manoeuvers in hours:
def burntimeorbiter(T, m_orbiter, dV_orbiter, m_combined, dV_combined, Isp):
    a = T / (m_orbiter + mass_prop(m_orbiter, dV_orbiter, m_combined, dV_combined, Isp)[0])
    t_b = dV_orbiter / a
    t_b = t_b / 60 / 60
    return t_b


def total_cost(m_ox, m_fuel, m_engine=m_engine):
    return m_fuel * cost_kg_fuel + m_ox * cost_kg_fuel + 370422.80 # 1.7 * 4.97 * m_engine**0.823 * 1000 * 0.951 * 1.2


if __name__ == "__main__":
    ...
    # print('total dV: ', mass_prop(m_orbiter, dV_orbiter, m_combined, dV_combined, g, I_sp))
    # print(burntimeorbiter(T, m_orbiter, dV_orbiter))
    # print(burntimecombined(T, m_combined, dV_combined))