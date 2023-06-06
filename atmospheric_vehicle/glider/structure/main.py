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

def calculate_A_wing(c_w): # estimates area of wing section as function of chord
    'Analysis performed for airfoil FX 62-K-153/20 (same one used in MSc Thesis), if changes then remember to change properties!'
    K_A = 0.6 # constant value valid for approx. all commonly used airfoils, found on MIT document
    t_w = 0.153 * c_w # max thickness, value found on airfoiltools.com
    A_wing = K_A * t_w * c_w
    return A_wing

A_w = calculate_A_wing(c_w)
#print(A_w)

def calculate_I_wing(c_w): # estimates moment of inertia of wing section as function of chord
    'Analysis performed for airfoil FX 62-K-153/20 (same one used in MSc Thesis), if changes then remember to change properties!'
    K_I = 0.036 # constant value valid for approx. all commonly used airfoils, found on MIT document
    t_w = 0.153 * c_w  # max thickness, value found on airfoiltools.com
    h_w = 0.041 * c_w # maximum camber, value found on airfoiltools.com
    I_wing = K_I * c_w * t_w * (t_w**2 + h_w**2)
    return I_wing

I_w = calculate_I_wing(c_w)
print(I_w)

if __name__ == "__main__":
    print("Hello World")