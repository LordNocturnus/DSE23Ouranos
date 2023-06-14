import numpy as np



def mmoi(cyl_r, cyl_h, cyl_com, cyl_mass,
         tank_r, tank_com, tank_mass,
         caps_r, caps_h, caps_com, caps_mass, ):
    """"""

    com = (cyl_com*cyl_mass + tank_com*tank_mass + caps_com*caps_mass)/(cyl_mass + tank_mass + caps_mass)
    # cylinder mmoi
    cyl_i_xx = cyl_mass*cyl_r**2
    cyl_i_yy = cyl_mass/12*(6*cyl_r**2+cyl_h**2)
    cyl_i_zz = cyl_i_yy

    cyl_i_xx = cyl_i_xx + cyl_mass*(com[0]-cyl_com[0])**2
    cyl_i_yy = cyl_i_yy + cyl_mass * (com[1] - cyl_com[1]) ** 2
    cyl_i_zz = cyl_i_zz + cyl_mass*(com[2]-cyl_com[2])**2
    print(cyl_i_xx, cyl_i_yy, cyl_i_zz)

    # Tank mmoi
    tank_i_xx, tank_i_yy, tank_i_zz = 2/5*tank_mass*tank_r**2



def grav_grad_torque(mmoi, a, grav_param, phase_angle, true_anomaly):
    """
    :param mmoi: 3x3 mass moment of inertia matrix
    :param a: Semi major axis of the orbited body
    :param grav_param: Gravitational parameter of the orbited planet
    :param phase_angle: Phase angle of the
    :param true_anomaly:
    :return:
    """
    i_xx = mmoi[0][0]
    i_yy = mmoi[1][1]
    i_zz = mmoi[2][2]
    torque_max = 3/2 * grav_param/a**3*np.array([(i_yy-i_zz),
                                           (i_xx-i_zz),
                                           0])

def spin_stab():

    # H = omega * I
    # delta_theta = omega*t = 0.5*alpha*t^2 + omega_0*t
    # omega = alpha*t + omega_0
        # Where alpha is rotational acceleration
        # t is time
        # omega_0 is initial rotational rate

grav_parameter = 5.79394e15
semi_major = 293269000
tank_


if __name__ == "__main__":
    print("Hello World")