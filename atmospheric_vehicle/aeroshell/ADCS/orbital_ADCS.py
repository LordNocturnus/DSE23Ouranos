import numpy as np
from atmospheric_vehicle.aeroshell import Aeroshell


def adcs_sizing(aeroshell: Aeroshell, theta_dot_spin, max_volume, redundancy=True):
    """
    Input: capsule mass, capsule mmoi, required torque/rotational acceleration
    Output: thruster number, thruster max torque, thruster location
    """
    mass_top_small = np.pi*aeroshell.r_top_small**2*aeroshell.t_top*aeroshell.rho

    mmoi_xx_top_small = 0.5*aeroshell.r_top_small**2*mass_top_small
    mmoi_xx = aeroshell.mmoi[0][0]  # kgm^2
    mmoi_yy = aeroshell.mmoi[1][1]  # kgm^2
    mmoi_zz = aeroshell.mmoi[2][2]  # kgm^2
    mmoi_xy = aeroshell.mmoi[0][1]  # kgm^2
    mmoi_xz = aeroshell.mmoi[0][2]  # kgm^2
    mmoi_yz = aeroshell.mmoi[1][2]  # kgm^2

    vertices = np.array([0, 0],
                        [aeroshell.r_bottom_big, 0],
                        [aeroshell.taper])






# Output: number of thrusters,
# Manouevres: spin stabilisation