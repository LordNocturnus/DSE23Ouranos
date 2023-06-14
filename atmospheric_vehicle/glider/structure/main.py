import numpy as np
import matplotlib.pyplot as plt

'WINGS STRUCTURAL ANALYSIS'

def calculate_c_w(c_w_root, c_w_tip, b_w): # estimates different values of chord for different values of spanwise location
    """
    :param c_w_root: root chord
    :param c_w_tip: tip chord
    :param b_w: wing span
    :return: array with values of chord corresponding to different spanwise locations, array of spanwise locations, span increment
    """

    db = 0.1 # [m], span increment used for analysis
    b_range = np.arange(0, b_w / 2 + db, db)  # range of half-span values used for estimation
    c_wing = np.zeros(len(b_range)) # creates array that stores all values of chord corresponding to spanwise location
    for i in range(len(b_range)):
        c_wing[i] = -2 * (c_w_root - c_w_tip) * b_range[i] / b_w + c_w_root
    return c_wing, b_range, db


# In order to calculate the cross sectional properties needed for bending and tensile stresses (hence cross section A, moment of
# inertia I and x location of centroid (section is symmetric wrt z, so no need for  location)) the cross section of the wing has been
# idealised such that it is made by three different hollow parts: a half circle,
# a rectangle and an isosceles triangle. Ask Elia for more details on the geometry

'Wing cross section properties calculations'
#NOTE: for centroids, datum is always taken starting from the leftmost edge of cross section (where half circle begins)
# def calculate_props1(t_t):  # Calculate I, A and x_centroid for section 1 (half circle) as function of thickness t_t
#     Ixx_1 = (np.pi / 16 + 1 / 12) * t**3 * t_t # https://izabellaldgalvan.blogspot.com/2022/09/moment-of-inertia-of-circle.html (for half-circle)
#     Izz_1 = np.pi / 16 * t**3 * t_t # same website as above
#     A_1 = np.pi * t_t * (1 + t) # own work
#     Xcentr_1 = t / 2 * (np.pi - 1) / (np.pi + 1) # https://en.wikipedia.org/wiki/List_of_centroids (for anular sector with alpha=90)
#     return Ixx_1, Izz_1, A_1, Xcentr_1
#
# def calculate_props2(t_t):  # Calculate I, A and x_centroid for section 2 (rectangle) as function of thickness t_t
#     Ixx_2 = 1 / 2 * t**2 * t_t * a_2 + 1 / 6 * t_t * t**3 # https://amaris-has-pacheco.blogspot.com/2022/04/hollow-rectangular-beam-moment-of.html
#     Izz_2 = 1 / 2 * a_2**2 * t_t * t + 1 / 6 * t_t * a_2**3 # same website as above
#     A_2 = 2 * t_t * (t + a_2) # own work
#     Xcentr_2 = (t + a_2) / 2 # simply by looking at geometry
#     return Ixx_2, Izz_2, A_2, Xcentr_2
#
# def calculate_props3(t_t):  # Calculate I, A and x_centroid for section 3 (isosceles triangle) as function of thickness t_t
#     Ixx_3 = 1 / 8 * t**2 * t_t * a_3 + 1 / 24 * t_t * t**3 # https://amesweb.info/section/area-moment-of-inertia-of-isosceles-triangle.aspx (used formula for Iy)
#     Izz_3 = 1 / 6 * a_3**2 *t_t * t + 1 / 18 * a_3**3 * t_t # same website as above (used for Ixx)
#     A_3 = t_t * (t + a_3) # own work
#     Xcentr_3 = t /2 + a_2 + (a_3 / 3 * (t + 2 * a_3)) / (2 * a_3 + 2 * t) # same website as above (did outside triangle-inside triangle analysis)
#     return Ixx_3, Izz_3, A_3, Xcentr_3
#
# def calculate_props_total(Ixx_1, Ixx_2, Ixx_3, Izz_1, Izz_2, Izz_3, A_1, A_2, A_3,Xcentr_1, Xcentr_2, Xcentr_3):
#     Ixx_tot = Ixx_1 + Ixx_2 + Ixx_3
#     A_tot = A_1 + A_2 + A_3
#     Xcentr_tot = (A_1 * Xcentr_1 + A_2 * Xcentr_2 + A_3 * Xcentr_3) / A_tot
#     Izz_tot = (Izz_1 + A_1 * (Xcentr_tot - Xcentr_1) ** 2) + (Izz_2 + A_2 * (Xcentr_tot - Xcentr_2) ** 2) + (Izz_3 + A_3 * (Xcentr_tot - Xcentr_3) ** 2)
#     # parallel axis theorem added to all three parts
#     return Ixx_tot, Izz_tot, A_tot, Xcentr_tot

'Shear stress in xz plane due to torque caused by lift'
def calculate_tau_xz_torque(tau):
    """
    :param tau: critical shear stress
    :return: required thickness to resist load, centroid, torque, enclosed area
    """
    Xcentr_tot = ((np.pi * (1 + t) * t / 2 * (np.pi - 1) / (np.pi + 1) + 2 * (t + a_2) * (t / 2 + a_2 / 2) + (
                t + a_3) * (a_3 / 3 * (t + 2 * a_3) / (2 * a_3 + 2 * t) + t / 2 + a_2)) / (
                              np.pi * (1 + t) + 2 * (t + a_2) + t + a_3))
    L_res = L_distr * ((c_w[0] + c_w[-1]) * b_w / 2) / 2
    Am = t * (np.pi / 4 * t + a_2 + a_3 / 2)
    Ty = np.abs(L_res * ( Xcentr_tot - t / 2 - a_2 / 2))
    t_t = Ty / (2 * Am * tau)
    return t_t, Xcentr_tot, Ty, Am


'Normal stress in xz plane due to moment in the x direction caused by lift'
def calculate_sigma_bend_xz(sigma):
    """
    :param sigma: critical normal stress
    :return: required thickness to resist load, moment, moment of inertia wrt x
    """
    Mx = -2 / 3 * (c_w_root - c_w_tip) / b_w * L_distr * b_range**3 + c_w_root / 2 * L_distr * b_range**2
    t_t = (Mx * t / 2) / ((1 / 2 * t**2 * a_2 + 1 / 6 * t**3 + (np.pi / 16 + 1 / 12) * t**3 + 1 / 8 * t**2 * a_3
                           + 1 / 24 * t**3) * sigma)
    Ixx = t_t * (1 / 2 * t**2 * a_2 + 1 / 6 * t**3 + (np.pi / 16 + 1 / 12) * t**3 + 1 / 8 * t**2 * a_3 + 1 / 24 * t**3)
    return t_t, Mx, Ixx


'Shear stress in xz plane due to lift shear force'
def calculate_tau_xz_shear(t_t): # gives already q_max for every cross section
    """
    :param t_t: required thickness to resist load (to be iterated upon)
    :return: maximum shear flow in cross section
    """
    q_max = np.zeros(len(c_w)-1)
    for i in range(len(c_w)-1):
        S_shear = (c_w[i] + c_w[i+1]) * db / 2
        Vy = L_distr * S_shear
        A = np.array([[2 * (np.sqrt(2) + 1) / (t[i + 1] * G * t_t), -2 / (t[i + 1] * G * t_t), 0, -1],
                      [-1 / (2 * G * t_t * a_2[i + 1]), (a_2[i + 1] + t[i + 1]) / (G * t_t * a_2[i + 1] * t[i + 1]), -1 / (2 * G * t_t * a_2[i + 1]), -1],
                      [0, -1 / (a_3[i + 1] * G * t_t), (1 / (a_3[i + 1] * G * t_t) + np.sqrt(t[i + 1]**2 / 4 + a_3[i + 1]**2) / (t[i + 1] * a_3[i + 1] * G * t_t)), -1],
                      [t[i + 1]**2 / 2, 2 * a_2[i + 1] * t[i + 1], t[i + 1] * a_3[i + 1], 0]])
        s = np.array([[-Vy / (t[i + 1]**2 * G * t_t)],
                       [-Vy / (2 * G * t_t * a_2[i + 1] * t[i + 1])],
                       [-Vy / (2 * t[i + 1] * a_3[i + 1] * G * t_t)],
                       [-Vy * a_2[i + 1]]])
        solution = np.linalg.solve(A, s)
        qs0_1 = solution[0]
        qs0_2 = solution[1]
        qs0_3 = solution[2]
        q26 = Vy / (2 * t[i + 1]) + qs0_1 - qs0_2
        q35 = Vy / (2 * t[i + 1]) + qs0_2 - qs0_3
        if np.sum(q26) > np.sum(q35):
            q_max[i] = q26
        else:
            q_max[i] = q35
    return q_max


def thickness_tau_xz_shear(tau_yield): # calculates thickness necessary to resist shear caused by shear force in xz plane
    """
    :param tau_yield: critical shear stress
    :return: True if shear stress is survived, False if not, together with corresponding thickness
    """
    t_t_shear = np.zeros(len(c_w) - 1)
    safe_thick = 1.3
    for i in range(len(c_w) - 1):
        t_t0 = 0.000001
        q_maxs = calculate_tau_xz_shear(t_t0)
        q_max = q_maxs[i]
        tau = q_max / t_t0
        while tau > tau_yield:
            t_t0 += 0.000001
            q_max = calculate_tau_xz_shear(t_t0)[i]
            tau = q_max / t_t0
        t_t_shear[i] = t_t0
    if 0.5 * 10**-3 > safe_thick * max(t_t_shear):
        return True
    else:
        print(f'False, required thickness is: {max(t_t_shear)}')


'Stress due to thermal differences in atmosphere'
def calculate_sigma_thermal(T_space, T_01, T_1, T_20, sigma):
    """
    :param T_space: temperature in space
    :param T_01: temperature at 0.1 bar
    :param T_1: temperature at 1 bar
    :param T_20: temperature at 20 bar
    :param sigma: critical normal stress
    :return: True if load is survived, False if not, together with value of stress
    """
    temperatures = np.array([[T_space], [T_01], [T_1], [T_20]]) # array of most important temperature values encountered
    delta_T = np.zeros(len(temperatures) - 1)
    for i in range(len(temperatures) - 1):
        delta_T[i] = temperatures[i + 1] - temperatures[i]
    max_delta = max(delta_T)
    sigma_thermal = alpha * max_delta * E
    if sigma_thermal < sigma:
        return True
    else:
        print(f'False, stress due to thermal expansion is: {sigma_thermal / 10 ** 6} MPa')


'Normal stress due to buckling'
def calculate_sigma_xz_buckling(E_spar):  # calculates thickness of I-spars required to resist buckling
    """
    :param E_spar: elastic modulus of spars (assumed to be the same material as the wing skin for simplicity)
    :return: thickness required to resist load
    """
    t_I1 = np.zeros(len(c_w) - 1)
    t_I2 = np.zeros(len(c_w) - 1)
    for i in range(len(c_w)-1):
        S_shear = (c_w[i] + c_w[i+1]) * db / 2
        L = L_distr * S_shear
        t_I1[i] = (t[i] + np.sqrt(t[i] ** 2 - 96 * L / (np.pi ** 2 * E_spar))) / 4
        t_I2[i] = (t[i] - np.sqrt(t[i] ** 2 - 96 * L / (np.pi ** 2 * E_spar))) / 4
    if t_I1.any() > 1 * 10 ** -2:
        return t_I2
    elif t_I2.any() > 1 * 10 ** -2:
        return t_I1


'Shear stress due to lift shear force in yz plane'
def calculate_tau_yz_shear(tau): # calculates thickness necessary to resist shear caused by shear force in xz plane
    """
    :param tau: critical shear stress
    :return: thickness required to survive load, surface area of wing
    """
    S = (c_w[0] + c_w[-1]) * (b_w / 2) / 2
    L = -L_distr * S
    qs0 = L / (2 * t[-1]) * (1 / 2 - (t[-1] / (6 * (b_w + t[-1] / 3))) + b_w / (2 * (b_w + t[-1] / 3)))
    qmax23 = L / (4 * (b_w + t[-1] / 3)) - (L * b_w) / (t[-1] * (b_w + t[-1] / 3)) + qs0
    q12 = -L * (b_w) / (t[-1] * (b_w + t[-1] / 3)) + qs0
    qmax = max(qmax23, q12)
    t_t = qmax / tau
    return t_t, S


'Normal stress due to differences in pressure'
def calculate_hoop_stress_wing(sigma): # calculates hoop stress by assuming wing shape to be a truncated cone
    """
    :param sigma: critical normal stress
    :return: required thickness to survive load
    """
    delta_p = (20 - 1) * 101325  # [Pa], difference in pressure between maximum outside pressure and pressure inside wing
    # (wing keeps same pressure inside that it had when wing was "closed", hence Earth sea level atmospheric pressure)
    a = np.arctan((c_w[0] - c_w[-1]) / (2 * b_w))
    t_t = delta_p * c_w[0] * 2 / (2 * np.cos(a) * (sigma - 0.6 * delta_p))
    return t_t


'Mass of the wings based on thickness and density'
def calculate_mass_wings(t_t, rho): # calculates mass of wings according to maximum thickness needed to sustain loads
    """
    :param t_t: thickness of wings (given from loads/stresses analysis)
    :param rho: density of material used
    :return: mass of the wings
    """
    # find volume given by tip cross section
    Area1 = ((t[-1] / 2 + 1) * np.pi + 2 * (a_2[-1] + t[-1]) + (t[-1] + 2 * np.sqrt((t[-1] / 2) ** 2 + a_2[-1] ** 2))) * t_t
    Volume1 = Area1 * b_w
    # find volume given by root cross section
    Area2 = ((t[0] / 2 + 1) * np.pi + 2 * (a_2[0] + t[0]) + (t[0] + 2 * np.sqrt((t[0] / 2) ** 2 + a_2[0] ** 2))) * t_t
    Volume2 = Area2 * b_w
    # average volume
    Volume = (Volume1 + Volume2) / 2
    mass_wings = Volume * rho
    return mass_wings




'FUSELAGE STRUCTURAL ANALYSIS'

def thickness_torque(f_tail, l_tail, r_fuselage, tau_y):
    """
    Function to calculate the minimum thickness required to withstand the torque loads generated by the tail
    :param f_tail: Aerodynamics load generated by the tail
    :param l_tail: Tail span
    :param r_fuselage: Selected radius of the fuselage
    :param tau_y: Yield shear stress of the selected material
    :return: Minimum required thickness for the fuselage based on torque
    """
    return f_tail * (l_tail / 2 + r_fuselage) / (2 * np.pi * r_fuselage ** 2 * tau_y)


def thickness_hoop(delta_p, r_fuselage, sigma_y):
    """
    Function to calculate the minimum thickness required to withstand the pressure loads from the atmosphere
    :param delta_p: Atmospheric pressure encountered during operations at Uranus (worst case scenario 10 bar)
    :param r_fuselage: Selected radius of the fuselage
    :param sigma_y: Yield stress of the selected material
    :return: Minimum required thickness for the fuselage based on outside pressure
    """
    return delta_p * r_fuselage / sigma_y


def buckling_check(t_fuselage, r_fuselage, delta_p, E, l_fuselage):
    """
    Function to check if the critical buckling strength of the fuselage (based on previously calculated parameters) is
    big enough to withstand the atmospheric loads
    :param t_fuselage: Fuselage thickness
    :param r_fuselage: Fuselage radius
    :param delta_p: Atmospheric pressure encountered during operations at Uranus (worst case scenario 10 bar)
    :param E: Young modulus of the selected material
    :param l_fuselage: Fuselage length
    :return: Print summary statement indicating if check has been passed or not
    """
    sigma_c = (9 * (t_fuselage / r_fuselage) ** 1.6 + 0.16 * (t_fuselage / l_fuselage) ** 1.3) * E
    if sigma_c >= delta_p:
        print(f'Buckling check passed')
    else:
        print(f'The idealised structure would fail in buckling')


def fuselage_architecture(f_tail, l_tail, r_fuselage, l_fuselage, tau_y, delta_p, sigma_y):
    """
    Function to calculate the final architecture of the fuselage (thickness and mass) based on load applied. Buckling
    check is performed
    :param f_tail: Aerodynamics load generated by the tail
    :param l_tail: Tail span
    :param r_fuselage: Fuselage radius
    :param l_fuselage: Fuselage length
    :param tau_y: Yield shear stress of the selected material
    :param delta_p: Atmospheric pressure encountered during operations at Uranus (worst case scenario 10 bar)
    :param sigma_y: Yield stress of the selected material
    :return: mass and thickness of the fuselage
    """
    t_torque = thickness_torque(tail_load, l_tail, r_fuselage, tau_y)
    t_hoop = thickness_hoop(delta_p, r_fuselage, sigma_y)
    t_fuselage = round(max(t_hoop, t_torque, 0.8 * 10**-3), 4)
    buckling_check(t_fuselage, r_fuselage, p_atmos, E_metal, l_fuselage)
    mass_fuselage = 2 * np.pi * r_fuselage * l_fuselage * t_fuselage * density_metal
    return mass_fuselage, t_fuselage


if __name__ == "__main__":

    'WINGS ANALYSIS'


    'Material properties (Ti-6Al-4V, Titanium-Aluminum alloy)'  # source: Mechanics of Materials textbook (besides tau, look in notes)
    E = 120 * 10 ** 9  # [Pa]
    G = 44 * 10 ** 9  # [Pa]
    sigma_yield = 924 * 10 ** 6 * 1.1  # [Pa], normal yield stress, used as failure criterion with margin of safety
    tau_yield = 508.2 * 10 ** 6 * 1.1  # [Pa], shear yield stress, used as failure criterion with margin of safety
    rho = 4430  # [kg/m^3]
    alpha = 9.1 * 10 ** -6  # [1/K], coefficient of thermal expansion of material selected

    'Wing properties'
    c_w_root = 0.785  # [m], root chord, given from Luigi
    c_w_tip = 0.713  # [m], tip chord, given from Luigi
    b_w = 6  # [m], wing span, given from Luigi
    L_distr = 224.482  # [N/m^2], found from MSc Thesis (given from Luigi)

    c_w = calculate_c_w(c_w_root, c_w_tip, b_w)[0]  # array of chord values that will be used for later estimations
    b_range = calculate_c_w(c_w_root, c_w_tip, b_w)[1]
    db = calculate_c_w(c_w_root, c_w_tip, b_w)[2]

    'Atmospheric/physical properties'
    T_space = 2.73  # [K], average temperature in space. Source: https://phys.libretexts.org/Bookshelves/University_Physics/Book%3A_University_Physics_(OpenStax)/Book%3A_University_Physics_II_-_Thermodynamics_Electricity_and_Magnetism_(OpenStax)/01%3A_Temperature_and_Heat/1.04%3A_Thermal_Expansion
    T_01 = 53  # [K], temperature at 0.1 bar inside Uranus's atmosphere
    T_1 = 76  # [K], temperature at 1 bar inside Uranus's atmosphere
    T_20 = 193  # [K], temperature at 20 bar inside Uranus's atmosphere

    'Simplified dimensions used for wing cross section properties calculations'
    A = 0.28  # percentage value of segment a_2 with respect to the chord (found from geometrical sketch of airfoil)
    B = 0.08  # percentage value of segment t with respect to the chord (same as above)
    C = 0.68  # percentage value of segment a_3 with respect to the chord (same as above)
    a_2 = c_w * A  # [m], value used for the simplified airfoil cross section: will be found by analysing the actual geometry of the
    # airfoil and simplifying the central section into a hollow rectangle, set as dummy for now
    t = c_w * B  # [m], maximum thickness of airfoil. Set as dummy for now
    a_3 = c_w * C  # [m], same reasoning as for a_2, set as dummy for now

    'Safety factors'  # to be applied as safety margin
    safe_thick = 1.3  # applied to thicknesses (if calculated)
    safe_load = 1.1  # applied to loads (if calculated)

    'Results of calculations'
    t_t_sigma_bend_xz = safe_thick * calculate_sigma_bend_xz(sigma_yield)[0]
    t_t_sigma_bend_xz = max(np.delete(t_t_sigma_bend_xz, 0))
    print('The minimum thickness required for the wings to sustain normal stress loads due to bending in the xz plane'
          ' applied in the x direction is: ', t_t_sigma_bend_xz * 10 ** 3, ' mm')

    t_t_torque_xz = max(safe_thick * calculate_tau_xz_torque(tau_yield)[0])
    print('The minimum thickness required for the wings to sustain shear loads due to torques in the xz plane is: ',
          t_t_torque_xz * 10 ** 3, ' mm')

    print(f'Shear stress due to shear force check passed: {thickness_tau_xz_shear(tau_yield)}')

    print(f'Normal stress due to temperature changes check passed: {calculate_sigma_thermal(T_space, T_01, T_1, T_20, sigma_yield)}')

    t_t_buckling_xz = max(calculate_sigma_xz_buckling(E)) * safe_thick
    print(f'The minimum thickness required to resist buckling in spars is: {t_t_buckling_xz * 10 ** 3} mm')

    t_t_shear_yz = safe_thick * calculate_tau_yz_shear(tau_yield)[0]
    print('The minimum thickness required for the wings to sustain shear loads due to shear force in the yz plane is: ',
          t_t_shear_yz * 10 ** 3, ' mm')

    t_t_hoop = safe_thick * calculate_hoop_stress_wing(sigma_yield)
    print('The minimum thickness required for the wings to sustain hoop stress is: ',
          t_t_hoop * 10 ** 3, ' mm')

    max_t_t = max(t_t_sigma_bend_xz, t_t_torque_xz, t_t_buckling_xz, t_t_shear_yz, t_t_hoop)
    print('The minimum thickness value needed for the wings structure to survive loads is: ', max_t_t * 10 ** 3, 'mm')

    mass_wings = calculate_mass_wings(max_t_t, rho)
    print('The mass of the wings is: ', mass_wings, 'kg')


    'FUSELAGE ANALYSIS'

    # Chose composite material based on glider on earth https://avpay.aero/company/aeros/product/new-aeros-ac-21-motor-glider-for-sale/
    # Metal properties taken from wing design. Metal thickness found in https://www.sciencedirect.com/topics/materials-science/ti-6al-4v
    # CFRP properties: https://pdf.sciencedirectassets.com/305927/1-s2.0-S2214785321X00399/1-s2.0-S2214785320357618/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEEIaCXVzLWVhc3QtMSJGMEQCIAzM0u1KVyUv1rRKTmBOPN2OD3MD4awU%2F6V%2Fi7t1e2v3AiApbysS31YUijV5u0TFSbYDR2EpEQOkQriedGRWygrR%2Fyq7BQiL%2F%2F%2F%2F%2F%2F%2F%2F%2F%2F8BEAUaDDA1OTAwMzU0Njg2NSIM1BjKg1IEVW4b5wX5Ko8F7K1j5pX0AQ%2Fgx%2FIu%2FjMN66VC8%2Btd2yb%2FsGGriIzQwM4gvevRQfJiYND7qGlQ7gjg7ymH3q597xDU2nOdi729e6os18c5%2BrIoKmKoMb2mGFAuRnhOmTJ7I%2Fsgoi8NRsrg77I4WiGYCFBWq4TjGfE31y19L0A%2BucWQuxKJL2LHtJWWhgjOIDM1JBjfDUAKz15VM3Luni8SCT4Xql4nttVydLzyadRH4ThQBKj7nSvtRnZsFohPlwENCr6glYQF7RaV3mwRsSnxG%2BdFKXMdAh6Emn6hjEML9SEvtGMwTOQkPVmSP2PUNSJkR0%2FoaSNH4WJhBmp1xbKRmLaDnye8qhCQmVzvz7aneKg1jttbbgMd3Q%2FFBXty%2Fm3giBWzFG206aZhgXgpHzZH0B0CARQv%2BIGdBE63dZpGwrPdLNN7ZRRKzy91Q1LLbiKMNmlEpp3QiOle8rECEY5M969gpIDITa7fkhwPzOePGveGrBzokZ%2F2n1kfEB7RceabjF5I5GsW86RLlmOXB4e14a6f8XBhdKPB6mSWEUDVn9BBbrMWjATKJVYok0zR6xvdx%2F1sxBpEb%2B2eu6UvHmBd1g1qppeXUJu7YGOwOhbjE%2FHwp6Szf%2B0Hnix0%2B0y5ave6G6r%2FQP2HHSxpTEt85vHZVUMS5fBB0Z0GuAZpIFixKTvhQyAjvSdFrQWd0NK6j9mnJb6PjCtZBOcMnqsRqMdg5n6U9uPXIoS%2B6lvOnUzxR%2BPyS%2FnTFR4dpDhzgeTGWthairhBk4cV3AMaHjXJ1lfeP75X4ZUENqeJCdO0DJZhKiWelgyXtAve9EWhL0lflemZpMhppX99hdl8wvHYsG2bhnsvPATwURft7jHgZS4IbdX5DYlS6%2FNndzCO54ukBjqyAaxM%2B2w4A18IYQZ8q873yp148wx8euJMOuMW8AxjdeLnC3F9atotltoC9qe1Y0xbO1UJNLU%2B4yqjBw7yIqqltLX3NhCyaj4W3qeNHetgmpCgtcvafphbP3gxdtgBORxwlBeAaNgy7OxTt4%2F9q6xq5vFrI4F5OGth%2FFq%2FDeSo0xYF%2Fd9e65mG6TBCSJU9bW%2FiIvCx3A5GSKluAQ5mcN%2FEstX0KjeOzgM1Ss4%2FaM3GmA%2BeNiw%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20230609T100622Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYWSCMYTFE%2F20230609%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=c8ae11c9a34cff7baf36b0d2dfe67b59b1c24c10a3ecbb73df3399f5eaf1fb80&hash=8fdd270f536e1f29521b0d59b6485cedc89b1a3468c4ab92f667594d8f47a3f0&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S2214785320357618&tid=spdf-f50fac15-cc70-429e-affa-a548257d29b0&sid=d9b148212649a84cae6a0426731e1314c615gxrqb&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=0b0a520b505a59030252&rr=7d488fc07a67b93f&cc=nl
    sigma_y_CFRP = 750 * 10 ** 6  # Yield strength CFRP file:///C:/Users/mzamb/AppData/Local/Microsoft/Windows/INetCache/IE/771OFKP1/Prepreg_Technology_[1].pdf
    E_CFRP = 37 * 10 ** 9  # Young Modulus CFRP https://matweb.com/search/datasheettext.aspx?matguid=39e40851fc164b6c9bda29d798bf3726
    tau_y_CFRP = 86.7 * 10 ** 6  # Shear yield strength CFRP
    density_CFRP = 2100  # Density CFRP
    sigma_y_metal = 455 * 10 ** 6  # Yield strength metal
    tau_y_metal = 172 * 10 ** 6  # Shear yield strength metal
    E_metal = 73 * 10 ** 9  # Young modulus metal
    density_metal = 2770  # Density metal
    p_atmos = 1000000 * 1.5  # Atmospheric pressure with safety factor (ADSEE reader launcher)
    l_battery = 0.3  # Battery length
    l_obc = 0.1  # OBC length https://blog.satsearch.co/2021-01-21-spotlight-how-to-choose-a-satellite-on-board-computer-obc
    tail_load = 224.48 * 0.2 * 1.1  # Aerodynamic loads on the tail
    l_tail = 1.3  # Tail span
    r_fuselage = 0.56 / 2  # Fuselage radius
    l_fuselage = 0.49 + l_obc + l_battery + 0.785  # Fuselage length (0.49 comes from payload sizing)

    m_fuselage, t_fuselage = fuselage_architecture(tail_load, l_tail, r_fuselage, l_fuselage, tau_y_metal, p_atmos,
                                                   sigma_y_metal)
    print('Fuselage mass: ', m_fuselage, 'kg')
    print('Fuselage thickness: ', t_fuselage * 10**3, 'mm')
    print('Total structure mass: ', m_fuselage + mass_wings, 'kg')
