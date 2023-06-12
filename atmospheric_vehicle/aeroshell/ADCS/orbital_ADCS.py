import matplotlib.pyplot as plt
import numpy as np
# from atmospheric_vehicle.aeroshell import Aeroshell

r_t_small = 0.25
r_t_big = 0.3125
r_b_small = 3.6
r_b_big = 4.5
a_t = 0.12435499454676144
a_b = 0.4636476090008061
theta_d = 1
com = np.array([0, 0.5]).reshape(2, 1)
force_x = 1
force_y = 0
force_z = 0
moment_x = 0
moment_y = 0
moment_z = 0


def adcs_sizing_capsule(r_top_small, r_top_big, r_bottom_small, r_bottom_big, a_top, a_bottom, f_x, f_y, f_z, theta_dot_spin,
                redundancy=True):
    """
    Input: capsule mass, capsule mmoi, required rotational acceleration, required thrust per axis
    Output: thruster number, thruster max torque, thruster location
    """

    # mmoi_xx_top_small = 0.5*aeroshell.r_top_small**2*mass_top_small
    # mmoi_xx = aeroshell.mmoi[0][0]  # kgm^2
    # mmoi_yy = aeroshell.mmoi[1][1]  # kgm^2
    # mmoi_zz = aeroshell.mmoi[2][2]  # kgm^2
    # mmoi_xy = aeroshell.mmoi[0][1]  # kgm^2
    # mmoi_xz = aeroshell.mmoi[0][2]  # kgm^2
    # mmoi_yz = aeroshell.mmoi[1][2]  # kgm^2

    h_bottom = (r_bottom_big-r_bottom_small)/np.tan(a_bottom)
    h_top = (r_top_big-r_t_small)/np.tan(a_top)
    vertices = np.array([[r_bottom_big, 0],
                        [r_bottom_small, h_bottom],
                        [r_top_big, h_bottom]])
    # where can things be placed?


def adcs_sizing_orbiter(r, h, com):

    vertices = np.array([[0,]])
    n = 50
    points = np.array([0, 0]).reshape(2, 1)
    for i in range(int(np.size(vertices)/2)-1):
        x1 = vertices[i][0]
        x2 = vertices[i+1][0]
        y1 = vertices[i][1]
        y2 = vertices[i+1][1]
        x_range = np.linspace(x1, x2, n)
        y_range = np.linspace(y1, y2, n)
        points_new = np.vstack((x_range, y_range))
        points = np.hstack((points, points_new))
    x = points[0][:]
    y = points[1][:]

    print(x, y)

    # print(np.size(points))
    # print(np.size(points[0, :]))
    arm_x = points[0, :] - com[0]*np.ones(np.size(points, 1))
    arm_y = points[1, :] - com[1]*np.ones(np.size(points, 1))
    arms = np.hstack((arm_x, arm_y))

    print(arms)
    plt.scatter(x, y)
    plt.scatter(com[0], com[1], color='red')
    plt.show()


adcs_sizing(r_t_small, r_t_big, r_b_small, r_b_big, a_t, a_b, force_x, force_y, force_z, theta_d)







# Output: number of thrusters,
# Manouevres: spin stabilisation