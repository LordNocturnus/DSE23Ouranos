import numpy as np

# Uranus Variables
g_u = 8.69  # m/s^2
p_max = 20 * 10 ** 5  # Pa

# Glider Variables
m_glider = 269.3  # kg
S = 2.82  # m^2
W_S = 845.85
b = 8.63  # m
chord = 0.3266  # m
A = b * b / S
e = 0.8

# Gliding Variables
# gamma_glide = np.radians(-1.5)
C_L_C_D = 43.2
C_D_0 = 0.0089
C_L_opt = 0.769
C_D_opt = 2 * C_D_0

# Deployment Variables
gamma_0_deploy = np.deg2rad(-70)  # rad
V_capsule = 30  # m/s
V_y0_deploy = V_capsule  # m/s
V_x0_deploy = 0
V_0_deploy = np.sqrt(V_x0_deploy * V_x0_deploy + V_y0_deploy * V_y0_deploy)
C_D_deploy = 0.09
C_L_deploy = 0.9


def calc_glider_parameters(time, altitudedistance, speed, mass, rho, g, cl_0, cl_alpha, cd_0, e, taper):
    cl_cd = speed * time / altitudedistance
    cl = cl_cd * 2 * cd_0
    cd = 2 * cd_0
    alpha = (cl - cl_0) / cl_alpha
    wing_loading = 1 / 2 * rho * speed ** 2 * cl
    surface = mass * g / wing_loading
    AR = cl ** 2 / (cd_0 * np.pi * e)
    b = np.sqrt(surface * AR)
    c_mean = b / AR
    c_r = c_mean * 3 / (2 * (1 + taper + taper ** 2) / (1 + taper))
    c_t = c_r * taper
    return cl_cd, cl, cd, alpha, wing_loading, surface, AR, b, c_mean, c_r, c_t


if __name__ == "__main__":
    print(calc_glider_parameters(12 * 3600, 200000, 200, 269.3, 0.055, 8.7, 0.771, 0.0909, 0.0089, 0.8, 0.4))
