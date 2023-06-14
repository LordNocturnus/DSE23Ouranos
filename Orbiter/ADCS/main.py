import numpy as np



def mmoi(cyl_r, cyl_h, cyl_com, cyl_mass, tank_r, tanks_mass, ox_mass, ox_com, fuel_mass, fuel_com,
         caps_r, caps_h, caps_com, caps_mass, tanks_full=True, caps_attached=True):
    """
    :param cyl_r:
    :param cyl_h:
    :param cyl_com:
    :param cyl_mass:
    :param tank_r:
    :param tanks_mass:
    :param ox_mass:
    :param ox_com:
    :param fuel_mass:
    :param fuel_com:
    :param caps_r:
    :param caps_h:
    :param caps_com: Local centre of mass of the capsule, this is later shifted upwards with the height of the cylinder
    :param caps_mass:
    :param tanks_full:
    :param caps_attached:
    :return:
    """
    tank_mass = tanks_mass/2
    caps_com = np.array([(cyl_h + caps_h - caps_com[0]), 0, 0]) - caps_com
    if tanks_full and caps_attached:
        com = (cyl_com * cyl_mass + ox_com * (tank_mass + ox_mass) + fuel_com * (tank_mass + fuel_mass) +
               caps_com * caps_mass)/(cyl_mass + 2 * tank_mass + ox_mass + fuel_mass + caps_mass)
    if tanks_full and not caps_attached:
        com = (cyl_com * cyl_mass + ox_com * (tank_mass + ox_mass) + fuel_com * (tank_mass + fuel_mass))/\
              (cyl_mass + 2 * tank_mass + ox_mass + fuel_mass)
    if caps_attached and not tanks_full:
        com = (cyl_com * cyl_mass + ox_com * tank_mass + fuel_com * tank_mass +
               caps_com * caps_mass) / (cyl_mass + 2 * tank_mass + caps_mass)
    if not caps_attached and not tanks_full:
        com = (cyl_com * cyl_mass + ox_com * tank_mass + fuel_com * tank_mass) / (cyl_mass + 2 * tank_mass)
    print(f"Centre of mass found to be {com}")
    # cylinder mass moment of inertia
    cyl_i_xx = cyl_mass*cyl_r**2
    cyl_i_yy = cyl_mass/12*(6*cyl_r**2+cyl_h**2)
    cyl_i_zz = cyl_i_yy

    cyl_i_xx = cyl_i_xx + cyl_mass * (com[1] - cyl_com[1]) ** 2 + (com[2]-cyl_com[2]) ** 2
    cyl_i_yy = cyl_i_yy + cyl_mass * (com[0] - cyl_com[0]) ** 2 + (com[2]-cyl_com[2]) ** 2
    cyl_i_zz = cyl_i_zz + cyl_mass * (com[0] - cyl_com[0]) ** 2 + (com[1]-cyl_com[1]) ** 2
    cyl_mmoi = np.array([cyl_i_xx, cyl_i_yy, cyl_i_zz])

    if tanks_full:  # Oxidiser tank mass moment of inertia
        ox_tank_i_xx = 2 / 5 * (tank_mass + ox_mass) * tank_r ** 2
        ox_tank_i_yy = ox_tank_i_xx
        ox_tank_i_zz = ox_tank_i_xx
        ox_tank_i_xx = ox_tank_i_xx + (tank_mass + ox_mass) * (com[1] - ox_com[1]) ** 2 + (com[2] - ox_com[2]) ** 2
        ox_tank_i_yy = ox_tank_i_yy + (tank_mass + ox_mass) * (com[0] - ox_com[0]) ** 2 + (com[2] - ox_com[2]) ** 2
        ox_tank_i_zz = ox_tank_i_zz + (tank_mass + ox_mass) * (com[0] - ox_com[0]) ** 2 + (com[1] - ox_com[1]) ** 2
        ox_mmoi = np.array([ox_tank_i_xx, ox_tank_i_yy, ox_tank_i_zz])

        # Fuel moment of inertia
        fuel_tank_i_xx = 2 / 5 * (tank_mass + fuel_mass) * tank_r ** 2
        fuel_tank_i_yy = fuel_tank_i_xx
        fuel_tank_i_zz = fuel_tank_i_xx
        fuel_tank_i_xx = fuel_tank_i_xx + (tank_mass + fuel_mass) * (com[1] - fuel_com[1]) ** 2 + (com[2] - fuel_com[2]) ** 2
        fuel_tank_i_yy = fuel_tank_i_yy + (tank_mass + fuel_mass) * (com[0] - fuel_com[0]) ** 2 + (com[2] - fuel_com[2]) ** 2
        fuel_tank_i_zz = fuel_tank_i_zz + (tank_mass + fuel_mass) * (com[0] - fuel_com[0]) ** 2 + (com[1] - fuel_com[1]) ** 2
        fuel_mmoi = np.array([fuel_tank_i_xx, fuel_tank_i_yy, fuel_tank_i_zz])
    else:
        ox_tank_i_xx = 2 / 3 * tank_mass * tank_r ** 2
        ox_tank_i_yy = ox_tank_i_xx
        ox_tank_i_zz = ox_tank_i_xx
        ox_tank_i_xx = ox_tank_i_xx + tank_mass * (com[1] - ox_com[1]) ** 2 + (com[2] - ox_com[2]) ** 2
        ox_tank_i_yy = ox_tank_i_yy + tank_mass * (com[0] - ox_com[0]) ** 2 + (com[2] - ox_com[2]) ** 2
        ox_tank_i_zz = ox_tank_i_zz + tank_mass * (com[0] - ox_com[0]) ** 2 + (com[1] - ox_com[1]) ** 2
        ox_mmoi = np.array([ox_tank_i_xx, ox_tank_i_yy, ox_tank_i_zz])

        fuel_tank_i_xx, fuel_tank_i_yy, fuel_tank_i_zz = 2 / 3 * tank_mass * tank_r ** 2
        fuel_tank_i_yy = fuel_tank_i_xx
        fuel_tank_i_zz = fuel_tank_i_xx
        fuel_tank_i_xx = fuel_tank_i_xx + tank_mass * (com[1] - fuel_com[1]) ** 2 + (com[2] - fuel_com[2]) ** 2
        fuel_tank_i_yy = fuel_tank_i_yy + tank_mass * (com[0] - fuel_com[0]) ** 2 + (com[2] - fuel_com[2]) ** 2
        fuel_tank_i_zz = fuel_tank_i_zz + tank_mass * (com[0] - fuel_com[0]) ** 2 + (com[1] - fuel_com[1]) ** 2
        fuel_mmoi = np.array([fuel_tank_i_xx, fuel_tank_i_yy, fuel_tank_i_zz])

    caps_i_xx = 0.5 * caps_mass*caps_r ** 2
    caps_i_yy = caps_mass * (3/20 * caps_r ** 2 + 1/10 * caps_h**2)
    caps_i_zz = caps_i_yy

    caps_i_xx = caps_i_xx + caps_mass * ((com[1] - caps_com[1]) ** 2 + (com[2] - caps_com[2])) ** 2
    caps_i_yy = caps_i_yy + caps_mass * ((com[0] - caps_com[0]) ** 2 + (com[2] - caps_com[2])) ** 2
    caps_i_zz = caps_i_zz + caps_mass * ((com[0] - caps_com[0]) ** 2 + (com[1] - caps_com[1])) ** 2

    caps_mmoi = np.array([caps_i_xx, caps_i_yy, caps_i_zz])

    mmoi = cyl_mmoi + ox_mmoi + fuel_mmoi + caps_mmoi
    print(mmoi)


def grav_grad_torque(mmoi, a, grav_param, phase_angle, true_anomaly):
    """
    :param mmoi: a 3-large mass moment of inertia vector
    :param a: Semi major axis of the orbited body
    :param grav_param: Gravitational parameter of the orbited planet
    :param phase_angle: Phase angle of the
    :param true_anomaly:
    :return:
    """
    i_xx = mmoi[0]
    i_yy = mmoi[1]
    i_zz = mmoi[2]
    torque_max = 3/2 * grav_param/a**3*np.array([(i_yy-i_zz),
                                                 (i_xx-i_zz),
                                                 0])

# def spin_stab():

    # H = omega * I
    # delta_theta = omega*t = 0.5*alpha*t^2 + omega_0*t
    # omega = alpha*t + omega_0
        # Where alpha is rotational acceleration
        # t is time
        # omega_0 is initial rotational rate

grav_parameter = 5.79394e15
semi_major = 293269000

cylinder_radius = 0.6
cylinder_height = 3.35
cylinder_com = np.array([cylinder_height/2, 0, 0])
cylinder_mass = 868

tank_radius = 0.5
tanks_mass = 106
oxidiser_mass = 2074
oxidiser_com = np.array([tank_radius/2, 0, 0])
fuelmass = 1257
fuelcom = np.array([tank_radius*1.5, 0, 0])
capsule_r = 2.25
capsule_h = 1.5
capsule_com = np.array([0, 0, 0])  # com of the capsule wrt its own base (thus defined negatively)
capsule_mass = 1230


if __name__ == "__main__":
    mmoi(cylinder_mass, cylinder_height, cylinder_com, cylinder_mass, tank_radius, tanks_mass, oxidiser_mass,
         oxidiser_com, fuelmass, fuelcom, capsule_r, capsule_h, capsule_com, capsule_mass)
    mmoi(cylinder_mass, cylinder_height, cylinder_com, cylinder_mass, tank_radius, tanks_mass, oxidiser_mass,
         oxidiser_com, fuelmass, fuelcom, capsule_r, capsule_h, np.array([0.5, 0, 0]), capsule_mass)
