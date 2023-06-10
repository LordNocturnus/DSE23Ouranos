import numpy as np
import matplotlib.pyplot as plt



def calculate_c_w(c_w_root, c_w_tip, b_w): # estimates different values of chord for different values of spanwise location
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
def calculate_tau_xz(tau):
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
    Mx = -2 / 3 * (c_w_root - c_w_tip) / b_w * L_distr * b_range**3 + c_w_root / 2 * L_distr * b_range**2
    t_t = (Mx * t / 2) / ((1 / 2 * t**2 * a_2 + 1 / 6 * t**3 + (np.pi / 16 + 1 / 12) * t**3 + 1 / 8 * t**2 * a_3
                           + 1 / 24 * t**3) * sigma)
    Ixx = t_t * (1 / 2 * t**2 * a_2 + 1 / 6 * t**3 + (np.pi / 16 + 1 / 12) * t**3 + 1 / 8 * t**2 * a_3 + 1 / 24 * t**3)
    return t_t, Mx, Ixx


def calculate_tau_xz_shear(t_t): # gives already q_max for every cross section
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


def thickness_tau_xz_shear(tau_yield):
    t_t_shear = np.zeros(len(c_w) - 1)
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
    return 1 * 10**-3 > max(t_t_shear)

if __name__ == "__main__":
    'Material properties (Ti-6Al-4V, Titanium-Aluminum alloy)'  # source: Mechanics of Materials textbook
    E = 120 * 10 ** 9  # [Pa]
    G = 44 * 10 ** 9  # [Pa]
    sigma_yield = 924 * 10 ** 6  # [Pa], normal yield stress, used as failure criterion
    tau_yield = 508.2 * 10 ** 6  # [Pa], shear yield stress, used as failure criterion
    rho = 4430  # [kg/m^3]

    'Wing properties'
    c_w_root = 0.785  # [m], root chord, given from Luigi
    c_w_tip = 0.713  # [m], tip chord, given from Luigi
    b_w = 6  # [m], wing span, given from Luigi
    L_distr = 224.482  # [N/m^2], found from MSc Thesis (given from Luigi)

    c_w = calculate_c_w(c_w_root, c_w_tip, b_w)[0]  # array of chord values that will be used for later estimations
    b_range = calculate_c_w(c_w_root, c_w_tip, b_w)[1]
    db = calculate_c_w(c_w_root, c_w_tip, b_w)[2]

    'Simplified dimensions used for wing cross section properties calculations'
    A = 0.35  # percentage value of segment a_2 with respect to the chord
    B = 0.5  # percentage value of segment t with respect to the chord: once the airfoil is selected, can be found on airfoiltools.com
    C = 0.4  # percentage value of segment a_3 with respect to the chord
    a_2 = c_w * A  # [m], value used for the simplified airfoil cross section: will be found by analysing the actual geometry of the
    # airfoil and simplifying the central section into a hollow rectangle, set as dummy for now
    t = c_w * B  # [m], maximum thickness of airfoil. Set as dummy for now
    a_3 = c_w * C  # [m], same reasoning as for a_2, set as dummy for now

    'Safety factor(s)'  # to be applied as safety margin
    safe_thick = 1.3  # applied to thicknesses calculated for all load cases HAS TO BE DELETED, WILL BE USED FACTOR FOR LOADS RATHER THAN THICKNESSES
    safe_load = 1.1  # applied to loads as safety factor REMEMBER TO APPLY THIS AND NOT ONE ABOVE
    t_t_sigma_bend_xz = safe_thick * calculate_sigma_bend_xz(sigma_yield)[0]
    t_t_sigma_bend_xz = max(np.delete(t_t_sigma_bend_xz, 0))
    print('The minimum thickness required for the wings to sustain normal stress loads due to bending in the xz plane'
          ' applied in the x direction is: ', t_t_sigma_bend_xz * 10 ** 3, ' mm')

    t_t_tau_xz = max(safe_thick * calculate_tau_xz(tau_yield)[0])
    print('The minimum thickness required for the wings to sustain shear loads due to torques in the xz plane is: ',
          t_t_tau_xz * 10 ** 3, ' mm')

    print(f'Shear stress due to shear force check passed: {thickness_tau_xz_shear(tau_yield)}')


