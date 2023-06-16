from GliderVariables import *


# Equations for Gliding Flight Calculations
def velocity(rho, W_S=W_S, C_L=C_L_opt):
    a = W_S
    b = 2 / rho
    c = 1 / C_L
    V = np.sqrt(a * b * c)
    return V


def range(h, C_L_C_D=C_L_C_D):
    R = C_L_C_D * h
    return R


def time_of_flight(delta_h, rho, W_S = W_S, C_L=C_L_opt, C_D=C_D_opt):
    a = np.sqrt(2/W_S)
    b = C_L ** (3 / 2) / C_D
    t = delta_h * np.sqrt(rho) * a * b
    return t
