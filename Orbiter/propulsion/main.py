import numpy as np
g = 9.81 # gravitational constant in m/s^2
I_sp = 321 # specific impulse in sec
m_flow = 2 # fuel mass flow in kg/s
t_b = 4 # burn time in sec
m_f = 2427 # final mass spacecraft without propellant, thus dry mass in kg
dV = 2000 # delta V require in m/s

# Calculate required thrust:
def thrust(m_flow, I_sp, g):
    F_T = m_flow * I_sp * g
    return F_T
# Calculate required propellant mass:
def mass_prop(t_b, g, I_sp):
    m_prop = thrust(m_flow, I_sp, g) * t_b / (g * I_sp)
    m_prop = m_prop / g
    return m_prop

def mass_prop2(m_f, dV, g, I_sp):
    m_prop = m_f * (np.exp(dV / (g * I_sp)) - 1)
    return m_prop



if __name__ == "__main__":

    # print(mass_prop(t_b, g, I_sp))
    print(mass_prop2(m_f, dV, g, I_sp))