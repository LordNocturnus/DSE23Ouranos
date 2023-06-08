from GliderVariables import *


# Equations for Gliding Flight Calculations
def velocity(rho):
    a = 2 * m * g_u
    b = rho * S * C_L_opt
    V = np.sqrt(a / b)
    return V


def range(h):
    R = C_L_C_D * h
    return R


def time_of_flight(delta_h, rho):
    a = np.sqrt(S / (2 * m * g_u))
    b = C_L_opt ** (3 / 2) / C_D
    t = delta_h * np.sqrt(rho) * a * b
    return t
