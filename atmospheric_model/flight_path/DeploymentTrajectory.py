from GliderVariables import *


# Equations for Glider Deployment

def dV_dt(rho, V, gamma, C_D=C_D_deploy, m=m_glider, g=g_u, S=S):
    a = -0.5 * rho * V * V * S * C_D / m
    b = -g * np.sin(gamma)
    dVdt = a + b
    return dVdt


def dgamma_dt(rho, V, gamma, C_L=C_L_deploy, m=m_glider, g=g_u, S=S):
    a = 0.5 * rho * V * S * C_L / m
    b = -g * np.cos(gamma) / V
    dgammadt = a + b
    return dgammadt


def v_deploy(V_0, dVdt, dt):
    return V_0 + dVdt * dt


def gamma_deploy(gamma_0, dgammadt, dt):
    return gamma_0 + dgammadt * dt
