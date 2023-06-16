import numpy as np
from atmospheric_model.flight_path.GliderVariables import m_glider, b

# General Paramters
I_x = 625.3693
I_y = 218.2010
I_z = 818.955
J_xz = np.sqrt(I_x * I_x + I_z*I_z)

k_x = np.sqrt(I_x/m_glider)
k_y = np.sqrt(I_y/m_glider)
k_z = np.sqrt(I_z/m_glider)
k_xz = J_xz/m_glider


# Dimensionless parameters
mu_c = 0.1
mu_b = 0.2
K_X_2 = (k_x/b) * (k_x/b)
K_Y_2 = (k_y/b) * (k_y/b)
K_Z_2 =(k_z/b) * (k_z/b)
K_XZ = J_xz/(b*b)

# Stability Derivatives
C_X_u = 0.7
C_Z_u = 1.3
C_m_u = 1.8

C_X_alpha = 0.8
C_Z_alpha = 1.1
C_m_alpha = 1.7

C_Y_beta = 0.9
C_l_beta = 2.1
C_n_beta = 2.4

C_l_p = 2
C_n_p = 2.5

C_l_r = 2.2
C_n_r = 2.3

# Other Derivatives
C_Z_alpha_dot = 1
C_Z_q = 1.2
C_Z_0 = 1.4

C_m_alpha_dot = 1.5
C_m_q = 1.6
