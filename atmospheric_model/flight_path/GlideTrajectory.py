from GliderVariables import *


# Equations for Gliding Flight Calculations
def velocity(rho, m=m, g=g_u, S=S, C_L=C_L_opt):
    a = 2 * m * g
    b = rho * S * C_L
    V = np.sqrt(a / b)
    return V


def range(h, C_L_C_D=C_L_C_D):
    R = C_L_C_D * h
    return R


def time_of_flight(delta_h, rho, S=S, m=m, g=g_u, C_L=C_L_opt, C_D=C_D_opt):
    a = np.sqrt(S / (2 * m * g))
    b = C_L ** (3 / 2) / C_D
    t = delta_h * np.sqrt(rho) * a * b
    return t
