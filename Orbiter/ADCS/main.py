import numpy as np


def mmoi(cyl_r, cyl_h, cyl_com, cyl_mass, tank_r, tanks_mass, ox_mass, ox_com, fuel_mass, fuel_com,
         caps_r, caps_h, caps_com, caps_mass, tanks_full=True, caps_attached=True, debug=False):
    """
    This function yields the mass moment inertia about the x, y, and z axes. x is aligned longitudinally, y is
    aligned with the antenna on the side of the orbiter, and z is at a right angle with both to produce a right-handed
    coordinate system
    :param cyl_r: Radius of cylindrical spacecraft shell
    :param cyl_h: Height of cylindrical spacecraft shell
    :param cyl_com: Centre of mass of cylindrical spacecraft shell wrt its own base as origin
    :param cyl_mass: Cylindrical spacecraft shell mass (dry mass minus tank mass)
    :param tank_r: Radius of propellant/oxidiser tanks
    :param tanks_mass: Total mass of propellant/oxidiser tanks
    :param ox_mass: Mass of oxidiser
    :param ox_com: Centre of mass of oxidiser (and its tank)
    :param fuel_mass: Mass of fuel
    :param fuel_com: Centre of mass of fuel (and its tank)
    :param caps_r: Base radius of capsule (bottom of cone)
    :param caps_h: Height of capsule from base to top
    :param caps_com: Centre of mass of the capsule wrt its tip (origin at the middle the top)
    :param caps_mass: Mass of the capsule (including glider, of course)
    :param tanks_full: True if the oxidiser and fuel tank are filled with oxidiser and fuel
    :param caps_attached: True if the capsule is not yet separated from the spacecraft
    :param debug: Extra information on the MMOIs if your numbers are off
    :return: Mass moment inertia vector containing I_xx, I_yy, I_zz
    """
    tank_mass = tanks_mass/2
    caps_com = np.array([cyl_h, 0, 0]) + caps_com
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
    # cylinder mass moment of inertia
    cyl_i_xx = cyl_mass * cyl_r ** 2
    cyl_i_yy = cyl_mass / 12 * (6 * cyl_r ** 2 + cyl_h ** 2)
    cyl_i_zz = cyl_i_yy

    cyl_i_xx += cyl_mass * ((com[1] - cyl_com[1]) ** 2 + (com[2]-cyl_com[2]) ** 2)
    cyl_i_yy += cyl_mass * ((com[0] - cyl_com[0]) ** 2 + (com[2]-cyl_com[2]) ** 2)
    cyl_i_zz += cyl_mass * ((com[0] - cyl_com[0]) ** 2 + (com[1]-cyl_com[1]) ** 2)
    cyl_mmoi = np.array([cyl_i_xx, cyl_i_yy, cyl_i_zz])

    if tanks_full:  # Oxidiser tank mass moment of inertia
        ox_tank_i_xx = 2 / 5 * (tank_mass + ox_mass) * tank_r ** 2
        ox_tank_i_yy = ox_tank_i_xx
        ox_tank_i_zz = ox_tank_i_xx
        ox_tank_i_xx += (tank_mass + ox_mass) * ((com[1] - ox_com[1]) ** 2 + (com[2] - ox_com[2]) ** 2)
        ox_tank_i_yy += (tank_mass + ox_mass) * ((com[0] - ox_com[0]) ** 2 + (com[2] - ox_com[2]) ** 2)
        ox_tank_i_zz += (tank_mass + ox_mass) * ((com[0] - ox_com[0]) ** 2 + (com[1] - ox_com[1]) ** 2)
        ox_mmoi = np.array([ox_tank_i_xx, ox_tank_i_yy, ox_tank_i_zz])

        # Fuel moment of inertia
        fuel_tank_i_xx = 2 / 5 * (tank_mass + fuel_mass) * tank_r ** 2
        fuel_tank_i_yy = fuel_tank_i_xx
        fuel_tank_i_zz = fuel_tank_i_xx
        fuel_tank_i_xx += (tank_mass + fuel_mass) * ((com[1] - fuel_com[1]) ** 2 + (com[2] - fuel_com[2]) ** 2)
        fuel_tank_i_yy += (tank_mass + fuel_mass) * ((com[0] - fuel_com[0]) ** 2 + (com[2] - fuel_com[2]) ** 2)
        fuel_tank_i_zz += (tank_mass + fuel_mass) * ((com[0] - fuel_com[0]) ** 2 + (com[1] - fuel_com[1]) ** 2)
        fuel_mmoi = np.array([fuel_tank_i_xx, fuel_tank_i_yy, fuel_tank_i_zz])
    else:
        ox_tank_i_xx = 2 / 3 * tank_mass * tank_r ** 2
        ox_tank_i_yy = ox_tank_i_xx
        ox_tank_i_zz = ox_tank_i_xx
        ox_tank_i_xx += tank_mass * ((com[1] - ox_com[1]) ** 2 + (com[2] - ox_com[2]) ** 2)
        ox_tank_i_yy += tank_mass * ((com[0] - ox_com[0]) ** 2 + (com[2] - ox_com[2]) ** 2)
        ox_tank_i_zz += tank_mass * ((com[0] - ox_com[0]) ** 2 + (com[1] - ox_com[1]) ** 2)
        ox_mmoi = np.array([ox_tank_i_xx, ox_tank_i_yy, ox_tank_i_zz])

        fuel_tank_i_xx = 2 / 3 * tank_mass * tank_r ** 2
        fuel_tank_i_yy = fuel_tank_i_xx
        fuel_tank_i_zz = fuel_tank_i_xx
        fuel_tank_i_xx += tank_mass * ((com[1] - fuel_com[1]) ** 2 + (com[2] - fuel_com[2]) ** 2)
        fuel_tank_i_yy += tank_mass * ((com[0] - fuel_com[0]) ** 2 + (com[2] - fuel_com[2]) ** 2)
        fuel_tank_i_zz += tank_mass * ((com[0] - fuel_com[0]) ** 2 + (com[1] - fuel_com[1]) ** 2)
        fuel_mmoi = np.array([fuel_tank_i_xx, fuel_tank_i_yy, fuel_tank_i_zz])

    if caps_attached:
        caps_i_xx = 0.5 * caps_mass*caps_r ** 2
        caps_i_yy = caps_mass * (3/20 * caps_r ** 2 + 1/10 * caps_h**2)
        caps_i_zz = caps_i_yy

        caps_i_xx += caps_mass * ((com[1] - caps_com[1]) ** 2 + (com[2] - caps_com[2]) ** 2)
        caps_i_yy += caps_mass * ((com[0] - caps_com[0]) ** 2 + (com[2] - caps_com[2]) ** 2)
        caps_i_zz += caps_mass * ((com[0] - caps_com[0]) ** 2 + (com[1] - caps_com[1]) ** 2)

        caps_mmoi = np.array([caps_i_xx, caps_i_yy, caps_i_zz])
    else:
        caps_mmoi = np.zeros(3)
    if debug:
        if caps_attached and tanks_full:
            print(f"---- MMOI with full propellant tanks and attached capsule ----")
        elif caps_attached and not tanks_full:
            print(f"---- How come you ran out of fuel with the capsule still attached? ----")
        elif not caps_attached and tanks_full:
            print(f"---- MMOI after capsule separation with filled propellant tanks ----")
        else:
            print(f"---- MMOI at end of life, with capsule separation and empty fuel tanks ----")
        print(f""
              f"Centre of mass found to be {com}\n"
              f"Cylinder MMOI found to be: {cyl_mmoi}\n"
              f"Oxidiser tank MMOI found to be: {ox_mmoi}\n"
              f"Fuel tank MMOI found to be: {fuel_mmoi}\n"
              f"Capsule MMOI found to be: {caps_mmoi}\n"
              f"Total MMOI found to be: {cyl_mmoi + ox_mmoi + fuel_mmoi + caps_mmoi}\n")
    return cyl_mmoi + ox_mmoi + fuel_mmoi + caps_mmoi


def grav_grad_torque(mmoi, a, grav_param, phase_angle=0, true_anomaly=0):
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
    n2 = grav_param / (a ** 3)  # mean motion
    torque_max = 3 / 2 * n2 * np.array([(i_yy - i_zz), (i_xx - i_zz), 0])
    return torque_max


def mag_torque_max(m, b):
    """
    :param m: Spacecraft residual dipole vector [A m^2]
    :param b: Magnetic field vector in spacecraft coordinates [T]
    :return:
    """
    return m * b * np.ones(3)



# def spin_stab():

    # H = omega * I
    # delta_theta = omega*t = 0.5*alpha*t^2 + omega_0*t
    # omega = alpha*t + omega_0
        # Where alpha is rotational acceleration
        # t is time
        # omega_0 is initial rotational rate


def ad_budget(n_suns, n_gyr, n_sts, n_hors, n_magnet, m_suns, m_gyr, m_sts, m_hors, m_magnet, 
            p_suns, p_gyr, p_sts, p_hors, p_magnet, acc_suns, acc_gyr, acc_sts, acc_hors, acc_magnet):    
    m = n_suns * m_suns + n_gyr * m_gyr + n_sts * m_sts + n_hors * m_hors + n_magnet * m_magnet
    p = n_suns * p_suns + n_gyr * p_gyr + n_sts * p_sts + n_hors * p_hors + n_magnet * p_magnet
    return m, p


def torque_cg_unknown(max_range, thrust):
    thrust_torque = np.array([0, max_range*thrust, max_range*thrust])
    return thrust_torque


def torque_total(max_range, thrust, dipole, mag, moment_oi, a, grav_param):
    torque_point = mag_torque_max(dipole, mag) + grav_grad_torque(moment_oi, a, grav_param)
    torque_burn = torque_cg_unknown(max_range, thrust) + torque_point
    return torque_burn, torque_point


def slew_torque(theta, mmoi: np.array, t):
    """
    :param theta: Slew angle [deg]
    :param mmoi: Mass moment of inertia vector containing I_xx, I_yy, I_zz
    :param t: Time of slew manoeuvre
    :return:
    """
    print(f"Slew torque")
    return 4 * np.radians(theta) * mmoi / t ** 2


def momentum_storage(period, allowable_motion, torque_dist, ):
    ...
magnetic_dipole_uranus = 110e-6

grav_parameter = 5.79394e15
semi_major = 309000000

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
capsule_com = np.array([1, 0, 0])  # com of the capsule wrt its tip
capsule_mass = 1230

mmoi_EOL = mmoi(cylinder_radius, cylinder_height, cylinder_com, cylinder_mass, tank_radius, tanks_mass,
                oxidiser_mass, oxidiser_com, fuelmass, fuelcom, capsule_r, capsule_h, capsule_com, capsule_mass,
                debug=False, tanks_full=False, caps_attached=False)

mmoi_detached_full = mmoi(cylinder_radius, cylinder_height, cylinder_com, cylinder_mass, tank_radius, tanks_mass,
                          oxidiser_mass, oxidiser_com, fuelmass, fuelcom, capsule_r, capsule_h, capsule_com,
                          capsule_mass, debug=False, tanks_full=True, caps_attached=False)

mmoi_launch = mmoi(cylinder_radius, cylinder_height, cylinder_com, cylinder_mass, tank_radius, tanks_mass,
                   oxidiser_mass, oxidiser_com, fuelmass, fuelcom, capsule_r, capsule_h, capsule_com, capsule_mass,
                   debug=True, tanks_full=True, caps_attached=True)
# print(mag_torque_max(magnetic_dipole_uranus, 0.1))


if __name__ == "__main__":
    print(slew_torque(180, mmoi_EOL, 30))
    print(slew_torque(180, mmoi_EOL, 120))

    print(torque_cg_unknown(0.03, 635))
