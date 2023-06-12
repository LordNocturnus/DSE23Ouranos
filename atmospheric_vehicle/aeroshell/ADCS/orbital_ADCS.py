import numpy as np


def adcs_sizing(mass: float, mmoi: np.matrix, theta_dot, redundancy=True):
    """
    Input: capsule mass, capsule mmoi, required torque/rotational acceleration
    Output: thruster number, thruster max torque, thruster location
    """
    mmoi_xx = mmoi[0][0]  # kgm^2
    mmoi_yy = mmoi[1][1]  # kgm^2
    mmoi_zz = mmoi[2][2]  # kgm^2
    mmoi_xy = mmoi[0][1]  # kgm^2
    mmoi_xz = mmoi[0][2]  # kgm^2
    mmoi_yz = mmoi[1][2]  # kgm^2




    mass_expulsion_torque =


# Output: number of thrusters,
