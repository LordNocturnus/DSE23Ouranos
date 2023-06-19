import numpy as np


rho_uranus_max = 1.66e-9
rho_uranus_min = 6.902e-13
vel_per = 27000
c_d = 2.5
magnetic_dipole_uranus = 110e-6
m_time = 3*365*24*3600
orbital_period = 123*3600

grav_parameter = 5.79394e15
semi_major = 309000000





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
    caps_com_wrt_orbiter = np.array([cyl_h, 0, 0]) + caps_com
    if tanks_full and caps_attached:
        com = (cyl_com * cyl_mass + ox_com * (tank_mass + ox_mass) + fuel_com * (tank_mass + fuel_mass) +
               caps_com_wrt_orbiter * caps_mass)/(cyl_mass + 2 * tank_mass + ox_mass + fuel_mass + caps_mass)
    if tanks_full and not caps_attached:
        com = (cyl_com * cyl_mass + ox_com * (tank_mass + ox_mass) + fuel_com * (tank_mass + fuel_mass)) / \
              (cyl_mass + 2 * tank_mass + ox_mass + fuel_mass)
    if caps_attached and not tanks_full:
        com = (cyl_com * cyl_mass + ox_com * tank_mass + fuel_com * tank_mass +
               caps_com_wrt_orbiter * caps_mass) / (cyl_mass + 2 * tank_mass + caps_mass)
    if not caps_attached and not tanks_full:
        com = (cyl_com * cyl_mass + ox_com * tank_mass + fuel_com * tank_mass) / (cyl_mass + 2 * tank_mass)
    # cylinder mass moment of inertia about base
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
        cyl_area = cyl_h * cyl_r * 2
        cyl_centroid = cyl_h / 2 + caps_h
        caps_area = caps_h * caps_r / 2
        caps_centroid = caps_h * 2 / 3
        area = cyl_area + caps_area
        cop = (caps_centroid * caps_area + cyl_centroid * cyl_area) / area

        caps_i_xx = 3 / 10 * caps_mass * caps_r ** 2
        caps_i_yy = caps_mass * (3 / 20 * caps_r ** 2 + 3 / 80 * caps_h ** 2)
        caps_i_zz = caps_i_yy
        caps_mmoi_local = np.array([caps_i_xx, caps_i_yy, caps_i_zz])
        
        caps_i_xx += caps_mass * ((com[1] - caps_com_wrt_orbiter[1]) ** 2 + (com[2] - caps_com_wrt_orbiter[2]) ** 2)
        caps_i_yy += caps_mass * ((com[0] - caps_com_wrt_orbiter[0]) ** 2 + (com[2] - caps_com_wrt_orbiter[2]) ** 2)
        caps_i_zz += caps_mass * ((com[0] - caps_com_wrt_orbiter[0]) ** 2 + (com[1] - caps_com_wrt_orbiter[1]) ** 2)

        caps_mmoi = np.array([caps_i_xx, caps_i_yy, caps_i_zz])
    else:
        area = cyl_h * cyl_r * 2
        cop = cyl_h/2
        caps_mmoi_local = np.zeros(3)
        caps_mmoi = caps_mmoi_local
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
              f"Total assembly centre of mass at: {com}\n"
              f"Centre of pressure at {cop}\n"
              # f"Cylinder MMOI found to be: {cyl_mmoi}\n"
              # f"Oxidiser tank MMOI found to be: {ox_mmoi}\n"
              # f"Fuel tank MMOI found to be: {fuel_mmoi}\n"
              f"Local capsule MMOI: {caps_mmoi_local}\n"
              f"Total MMOI found to be: {cyl_mmoi + ox_mmoi + fuel_mmoi + caps_mmoi}\n")
    return cyl_mmoi + ox_mmoi + fuel_mmoi + caps_mmoi + caps_mmoi, com, cop, area, caps_mmoi_local


def grav_grad_torque(mmoi, a, grav_param, period):
    """
    :param mmoi: a 3-large mass moment of inertia vector
    :param a: Semi major axis of the orbited body
    :param grav_param: Gravitational parameter of the orbited planet
    :param period: Orbital period
    :return: Maximum torque vector, reaction wheel minimum moment size
    """
    i_xx = mmoi[0]
    i_yy = mmoi[1]
    i_zz = mmoi[2]
    n2 = grav_param / (a ** 3)  # mean motion
    torque_max = 3 / 2 * n2 * np.array([0, abs(i_xx - i_zz), abs(i_xx - i_yy)])
    return torque_max, max(torque_max) * 0.707 / 4 * period


def momentum_storage_mw(torque_dist, theta_a, period):
    return torque_dist / theta_a * period / 4

def mag_torque_max(d, b):
    """
    :param m: Spacecraft residual dipole vector [A m^2]
    :param b: Magnetic field vector in spacecraft coordinates [T]
    :return: Magnetic disturbance torque
    """
    return d * b


def aerodyn_torque(rho, cd, v, mmoi):
    return 0.5 * rho * cd * mmoi[3] * v ** 2 * max(mmoi[2] - mmoi[1])


def torque_s(mmoi, q):
    """
    Solar disturbance torque
    :param q: reflectance factor between 0 and 1
    :param mmoi: mass moment of inertia formula
    :return: solar disturbance torque
    """
    cp_s = mmoi[2]
    cm = mmoi[0][0]
    sun_surface = mmoi[3]
    c = 3e8
    return 4.1 / c * sun_surface * (1 + q) * (cp_s - cm)


def torque_cg_unknown(max_range, thrust):
    return max_range * thrust


def prop_mass(isp, force, r, t):
    return t * force / r / isp / 9.80665

def slew_torque(theta, mmoi: np.array, t):
    """
    :param theta: Slew angle [deg]
    :param mmoi: Mass moment of inertia vector containing I_xx, I_yy, I_zz
    :param t: Time of slew manoeuvre
    :return:
    """
    print(f"Slew torque")
    return 4 * np.radians(theta) * mmoi / t ** 2




cylinder_radius = 0.6
cylinder_height = 1.9
cylinder_mass = 868


# tank_radius = 0.45
# tanks_mass = 19
# oxidiser_mass = 348
# oxidiser_com = np.array([tank_radius, 0, 0])
# fuelmass = 211
# fuelcom = np.array([tank_radius * 3, 0, 0])
# capsule_r = 1.5
# capsule_h = 1.5
# capsule_com = np.array([capsule_h * 3 / 4, 0, 0])  # com of the capsule wrt its tip
# capsule_mass = 579.1


if __name__ == "__main__":
    # mmoi_empty = mmoi(cylinder_radius, cylinder_height, cylinder_com, cylinder_mass, tank_radius, tanks_mass,
    #                   oxidiser_mass, oxidiser_com, fuelmass, fuelcom, capsule_r, capsule_h, capsule_com, capsule_mass,
    #                   debug=True, tanks_full=False, caps_attached=False)
    #
    # mmoi_detached_full = mmoi(cylinder_radius, cylinder_height, cylinder_com, cylinder_mass, tank_radius, tanks_mass,
    #                           oxidiser_mass, oxidiser_com, fuelmass, fuelcom, capsule_r, capsule_h, capsule_com,
    #                           capsule_mass, debug=True, tanks_full=True, caps_attached=False)
    #
    # mmoi_launch = mmoi(cylinder_radius, cylinder_height, cylinder_com, cylinder_mass, tank_radius, tanks_mass,
    #                    oxidiser_mass, oxidiser_com, fuelmass, fuelcom, capsule_r, capsule_h, capsule_com, capsule_mass,
    #                    debug=True, tanks_full=True, caps_attached=True)
    # Report values disturbance torques at Uranus
    # print(grav_grad_torque(mmoi_empty[0], semi_major, grav_parameter, orbital_period))
    # print(mag_torque_max(10, magnetic_dipole_uranus))
    # print(aerodyn_torque(rho_uranus_min, c_d, vel_per, mmoi_empty))
    # print(torque_s(mmoi_empty, 1))
    # print(torque_cg_unknown(0.03, 425))
    # print(mmoi_empty[1])
    print(prop_mass(232, 6.8, 2.94/2, 15000))



    # print(slew_torque(180, mmoi_EOL[0], 30))
    # print(slew_torque(180, mmoi_EOL[0], 120))
    # print(f"CG unknown torque: {torque_cg_unknown(0.03, 400)}"
    #       f"CG unknown momentum: {torque_cg_unknown(0.03, 400) * 2846}")
    # h_rw = rw_sizing(dist_t, orbital_period)
    # h_mw = momentum_storage_mw(dist_t, 0.02, orbital_period)
    # print(h_rw, h_mw)
    # sol_torque = torque_s(mmoi_EOL, 1)
    # print(sol_torque * m_time)

    # print(rcs_fuel(3.4, 218, 4140))
