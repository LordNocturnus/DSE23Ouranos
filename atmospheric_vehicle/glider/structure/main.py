import numpy as np
import matplotlib.pyplot as plt
'Masses of different subsystems'
mass_PW = 1  # [kg], mass of power subsystem, set as dummy for now
mass_DH = 1  # [kg], mass of data handling subsystem, set as dummy for now
mass_TTC = 1  # [kg], mass of communication subsystem, set as dummy for now
mass_TM = 1  # [kg], mass of thermal subsystem, set as dummy for now
mass_ADCS = 1  # [kg], mass of ADCS subsystem, set as dummy for now
mass_PLD = 1  # [kg], mass of payload subsystem, set as dummy for now

'Volumes of different subsystems'
vol_PW = 1  # [m^3], volume of power subsystem, set as dummy for now
vol_DH = 1  # [m^3], volume of data handling subsystem, set as dummy for now
vol_TTC = 1  # [m^3], volume of communication subsystem, set as dummy for now
vol_TM = 1  # [m^3], volume of thermal subsystem, set as dummy for now
vol_ADCS = 1  # [m^3], volume of ADCS subsystem, set as dummy for now
vol_PLD = 1  # [m^3], volume of payload subsystem, set as dummy for now

'Wing properties'
c_w_root = 0.785 # [m], root chord, given from Luigi
c_w_tip = 0.713 # [m], tip chord, given from Luigi
b_w = 2 # [m], wing span, given from Luigi
L_distr = 224.482 # [N/m^2], found from MSc Thesis (given from Luigi)

def calculate_c_w(c_w_root, c_w_tip, b_w): # estimates different values of chord for different values of spanwise location
    b_range = np.arange(0, b_w / 2 + 0.1, 0.1)  # range of half-span values used for estimation
    c_wing = np.zeros(len(b_range)) # creates array that stores all values of chord corresponding to spanwise location
    for i in range(len(b_range)):
        c_wing[i] = -2 * (c_w_root - c_w_tip) * b_range[i] / b_w + c_w_root
    return c_wing

c_w = calculate_c_w(c_w_root, c_w_tip, b_w) # array of chord values that will be used for later estimations
#print(c_w)

'Simplified dimensions used for wing cross section properties calculations'
a_2 = c_w * 0.1 # [m], value used for the simplified airfoil cross section: will be found by analysing the actual geometry of the
# airfoil and simplifying the central section into a hollow rectangle, set as dummy for now. Shall be given as percentage of chord
a_3 = c_w * 0.2 # [m], same reasoning as above, set as dummy for now. Shall be given as percentage of chord
t = c_w * 0.05 # [m], maximum thickness of airfoil: once the airfoil is selected, can be found on airfoiltools.com as
# percentage of the chord. Set as dummy for now

# In order to calculate the cross sectional properties needed for bending and tensile stresses (hence cross section A, moment of
# inertia I and x location of centroid (section is symmetric wrt z, so no need for  location)) the cross section of the wing has been
# idealised such that it is made by three different hollow parts: a half circle,
# a rectangle and an isosceles triangle. Ask Elia for more details on the geometry

'Wing cross section properties calculations'
#NOTE: for centroids, datum is always taken starting from the leftmost edge of cross section (where half circle begins)
def calculate_props1(t_t):  # Calculate I, A and x_centroid for section 1 (half circle) as function of thickness t_t
    Ixx_1 = (np.pi / 16 + 1 / 12) * t**3 * t_t # https://izabellaldgalvan.blogspot.com/2022/09/moment-of-inertia-of-circle.html (for half-circle)
    Izz_1 = np.pi / 16 * t**3 * t_t # same website as above
    A_1 = np.pi * t_t * (1 + t) # own work
    Xcentr_1 = t / 2 * (np.pi - 1) / (np.pi + 1) # https://en.wikipedia.org/wiki/List_of_centroids (for anular sector with alpha=90)
    return Ixx_1, Izz_1, A_1, Xcentr_1

def calculate_props2(t_t):  # Calculate I, A and x_centroid for section 2 (rectangle) as function of thickness t_t
    Ixx_2 = 1 / 2 * t**2 * t_t * a_2 + 1 / 6 * t_t * t**3 # https://amaris-has-pacheco.blogspot.com/2022/04/hollow-rectangular-beam-moment-of.html
    Izz_2 = 1 / 2 * a_2**2 * t_t * t + 1 / 6 * t_t * a_2**3 # same website as above
    A_2 = 2 * t_t * (t + a_2) # own work
    Xcentr_2 = (t + a_2) / 2 # simply by looking at geometry
    return Ixx_2, Izz_2, A_2, Xcentr_2

def calculate_props3(t_t):  # Calculate I, A and x_centroid for section 3 (isosceles triangle) as function of thickness t_t
    Ixx_3 = 1 / 8 * t**2 * t_t * a_3 + 1 / 24 * t_t * t**3 # https://amesweb.info/section/area-moment-of-inertia-of-isosceles-triangle.aspx (used formula for Iy)
    Izz_3 = 1 / 6 * a_3**2 *t_t * t + 1 / 18 * a_3**3 * t_t # same website as above (used for Ixx)
    A_3 = t_t * (t + a_3) # own work
    Xcentr_3 = t /2 + a_2 + (a_3 / 3 * (t + 2 * a_3)) / (2 * a_3 + 2 * t) # same website as above (did outside triangle-inside triangle analysis)
    return Ixx_3, Izz_3, A_3, Xcentr_3

def calculate_props_total(Ixx_1, Ixx_2, Ixx_3, Izz_1, Izz_2, Izz_3, A_1, A_2, A_3,Xcentr_1, Xcentr_2, Xcentr_3):
    Ixx_tot = Ixx_1 + Ixx_2 + Ixx_3
    A_tot = A_1 + A_2 + A_3
    Xcentr_tot = (A_1 * Xcentr_1 + A_2 * Xcentr_2 + A_3 * Xcentr_3) / A_tot
    Izz_tot = (Izz_1 + A_1 * (Xcentr_tot - Xcentr_1)**2) + (Izz_2 + A_2 * (Xcentr_tot - Xcentr_2)**2) + (Izz_3 + A_3 * (Xcentr_tot - Xcentr_3)**2)
    # parallel axis theorem added to all three parts
    return Ixx_tot, Izz_tot, A_tot, Xcentr_tot
















if __name__ == "__main__":
    print("Hello World")