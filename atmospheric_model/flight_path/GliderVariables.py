import numpy as np

# Uranus Variables
g_u = 8.69  # m/s^2
p_max = 20 * 10 ** 5  # Pa

# Glider Variables
m_glider = 274.9  # kg
S = 4.49   # m^2
W_S = 279.3
b = 9.2  # m
chord = 0.9  # m
A = b * b / S
e = 0.8

# Gliding Variables
# gamma_glide = np.radians(-1.5)
C_L_C_D = 76.8
C_D_0 = 0.0083
C_L_opt = 0.45
C_D_opt = 2 * C_D_0

# Deployment Variables
gamma_0_deploy = np.deg2rad(-70)  # rad
V_capsule = 30  # m/s
V_y0_deploy = V_capsule  # m/s
V_x0_deploy = 0
V_0_deploy = np.sqrt(V_x0_deploy * V_x0_deploy + V_y0_deploy * V_y0_deploy)
C_D_deploy = 0.09
C_L_deploy = 0.9
