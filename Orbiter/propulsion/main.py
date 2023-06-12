import numpy as np

# combined system parameter:
m_combined = 2898  # total dry mass of both systems in kg
dV_combined = 1000  # delta V required for the combined system in m/s

# Orbiter system parameters:
m_orbiter = 2427  # dry mass of only the orbiter in kg
dV_orbiter = 1000  # delta V required for the orbiter alone in m/s

# Constants
g = 9.81  # gravitational constant in m/s^2

# Engine parameters:
I_sp = 321  # specific impulse in sec
T = 445  # thrust in N
m_flow = 1  # fuel mass flow in kg/s


# Calculate required propellant mass:
def mass_prop(m_orbiter, dV_orbiter, m_combined, dV_combined, g, I_sp):
    m_prop_orbiter = m_orbiter * (np.exp(dV_orbiter / (g * I_sp)) - 1)
    m_prop_combined = m_combined * (np.exp(dV_combined / (g * I_sp)) - 1)
    m_prop = m_prop_orbiter + m_prop_combined
    return m_prop_orbiter, m_prop_combined, m_prop


# Calculate required burn time for combined system manoeuvers in hours:
def burntimecombined(T, m_f, dV_combined, m_orbiter, m_combined, dV_orbiter, g, I_sp):
    a = T / (m_f + mass_prop(m_orbiter, dV_orbiter, m_combined, dV_combined, g, I_sp)[2])
    t_b = dV_combined / a
    t_b = t_b / 60 / 60
    return t_b


# Calculate required burn time for orbiter manoeuvers in hours:
def burntimeorbiter(T, m_f, dV_orbiter, m_orbiter, m_combined, dV_combined, g, I_sp):
    a = T / (m_f + mass_prop(m_orbiter, dV_orbiter, m_combined, dV_combined, g, I_sp)[0])
    t_b = dV_orbiter / a
    t_b = t_b / 60 / 60
    return t_b


if __name__ == "__main__":
    print('total dV: ', mass_prop(m_orbiter, dV_orbiter, m_combined, dV_combined, g, I_sp))
    print(burntimeorbiter(T, m_flow, dV_orbiter, m_orbiter, m_combined, dV_combined, g, I_sp))
    print(burntimecombined(T, m_flow, dV_combined, m_orbiter, m_combined, dV_orbiter, g, I_sp))
