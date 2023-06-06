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
c_w_root = 1 # [m], root chord, set as dummy for now
c_w_tip = 0.5 # [m], tip chord, set as dummy for now
b_w = 2 # [m], wing span, set as dummy for now


def calculate_c_w(c_w_root, c_w_tip, b_w): # estimates different values of chord for different values of spanwise location
    b_range = np.arange(0, b_w / 2 + 0.1, 0.1)  # range of half-span values used for estimation
    c_wing = np.zeros(len(b_range)) # creates array that stores all values of chord corresponding to spanwise location
    for i in range(len(b_range)):
        c_wing[i] = -2 * (c_w_root - c_w_tip) * b_range[i] / b_w + c_w_root
    return c_wing

c_w = calculate_c_w(c_w_root, c_w_tip, b_w) # array of chord values that will be used for later estimations
#print(c_w)

a_2 = c_w * 0.1 # [m], value used for the simplified airfoil cross section: will be found by analysing the actual geometry of the airfoil and
        # simplifying the central section into a hollow rectangle, set as dummy for now. Shall be given as percentage of chord
a_3 = c_w * 0.2 # [m], same reasoning as above, set as dummy for now. Shall be given as percentage of chord

'In order to calculate the cross sectional properties needed for bending and tensile stresses (hence cross section A and moment of'
'inertia I) the cross section of the wing has been idealised such that it is made by three different hollow parts: a half circle, '
'a rectangle and an isosceles triangle. Ask Elia for more details on the geometry'
def calculate_Ixx1(t_t):  # Calculate I for section 1 (half circle) as function of thickness t_t
    t = 1 # [m], maximum thickness of airfoil: once the airfoil is selected, can be found on airfoiltools.com as function of the chord
    Ixx_1 = np.pi / 16 * t**3 * t_t
    return Ixx_1

def calculate_Ixx2(t_t, a_2):

if __name__ == "__main__":
    print("Hello World")