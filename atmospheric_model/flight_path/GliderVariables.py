import numpy as np

# Uranus Variables
g_u = 8.69  # m/s^2
p_max = 20 * 10 ** 5  # Pa

# Glider Variables
m = 72  # kg
S = 4.495  # m^2
b = 6  # m
chord = 1  # m
A = b * b / S
e = 1

# Gliding Variables
gamma_glide = np.radians(-1.5)
C_L_C_D = 1 / np.sin(abs(gamma_glide))
C_D_0 = np.pi * A * e / (4 * C_L_C_D * C_L_C_D)
C_L_opt = np.sqrt(C_D_0 * np.pi * A * e)
C_D = 2 * C_D_0

# Deployment Variables
gamma_0_deploy = -np.pi / 2  # rad
V_capsule = 80  # m/s
V_0_deploy = V_capsule  # m/s
C_D_deploy = 0.03
C_L_deploy = 0.7498
