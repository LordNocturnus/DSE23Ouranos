from mpl_toolkits import mplot3d
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

# cyclindrical coordinates
# x = x
# y = r*cos(theta)
# z = r*sin(theta)


def adcs_sizing_orbiter(r, h, com, isp, thrust, thruster_mass, mmoi, mu_planet, a, phi_orbit, theta_orbit,
                        redundancy=True):
    """
    Function that allows the user to determine the sizing of the ADCS system based characteristics of the orbiter, its
    mission around Uranus, and the expected disturbances encountered in the mission.
    :param r: The radius of the cylinder
    :param h: The height of the cylinder
    :param com: The centre of mass of the spacecraft
    :param isp: The isp of the ADCS thruster considered
    :param thrust: The thrust of the ADCS thruster considered
    :param thruster_mass: Mass of the ADCS thruster considered
    :param mmoi: The matrix containing the mass moments of inertia of the orbiter in kgm^2
    :param mu_planet: The gravitational parameter of the orbited planet in m^3/s^2
    :param a: The semi major axis of the orbit around the planet in meters
    :param phi_orbit: Phase angle of the orbiter in its orbit in radians
    :param theta_orbit: True anomaly of the orbiter in its orbit in radians
    """
    grav_grad = 3/2 * mu_planet/a**3 * np.array([((mmoi[1][1] - mmoi[2][2])*np.sin(2*phi_orbit)),
                                                 ((mmoi[0][0] - mmoi[2][2])*np.sin(2*theta_orbit)),
                                                 0])
    print(grav_grad)







    vertices_2d = np.array([[r, 0],
                            [r, h],
                            [0, h]])
    n = 10  # scales to the power 3
    points = []
    x = vertices_2d[:, 1]  # radii
    # print(x)
    for j in range(len(vertices_2d[:, 0])-1):
        x_range = np.linspace(vertices_2d[j, 1], vertices_2d[j+1, 1], n)  # 50 heights
        r_range = np.linspace(vertices_2d[j, 0], vertices_2d[j+1, 0], n)  # 50 radii
        theta_range = np.linspace(0, 360, n)  # 50 angles
        for x in x_range:
            for radius in r_range:
                for angle in theta_range:
                    x = x
                    y = radius*np.cos(np.radians(angle))
                    z = radius*np.sin(np.radians(angle))
                    points.append([x, y, z])


    # consider possible adcs systems
    # Hydrazine MMH, N2O4,

    vector = []
    for point in points:
        vector.append(point - com)
    print(vector)
    points = np.array(points)
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(points[:, 0], points[:, 1], points[:, 2])
    ax.scatter(com[0], com[1], com[2], color='red', s=100)
    # plt.scatter(com[0], com[1], color='red')
    plt.show()


# adcs_sizing_capsule(r_t_small, r_t_big, r_b_small, r_b_big, a_t, a_b, force_x, force_y, force_z, theta_d)

mmoi = np.array([[1, 0, 0],
                 [0, 1, 0],
                 [0, 1, 0]])
h = 2.13
d = 1
r = d/2
com = np.array([0.5, 0, 0])
isp = 292
adcs_sizing_orbiter(r, h, com)

