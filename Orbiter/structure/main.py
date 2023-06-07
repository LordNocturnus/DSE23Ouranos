# Tool to determine the sizing of the orbiter structure
import numpy as np


def axial_loads(axial, sigma_y):
    """
    Function to determine the base area of the orbiter based on launch loads
    :param axial: Axial peak loads during launch
    :param sigma_y: Yield stress of the selected material
    :return: Required minimum area to withstand the loads
    """
    return np.sqrt(axial / (sigma_y * np.pi))

def lateral_loads(lateral, I, sigma_y, r):
    """
    Function to calculate the minimum required length of the orbiter based on launch loads.
    Worst case scenario assumed, with lateral loads fully applied at the tip of the structure.
    :param lateral: Lateral loads during launch
    :param I: Moment of inertia of the structure
    :param sigma_y: Yield strength of the selected material
    :param r: Radius of the orbiter structure
    :return: Minimum required length to witshtand lateral loads
    """
    return 6 * lateral / (sigma_y * r)

def t_hoop_sphere(p, r, sigma_y):
    """
    Function to calculate the thickness of the orbiter due to hoop stress in sphere. Need to be used for
    spherical caps in case are present or for the entire structure if spherical tanks are implemented.
    :param p: Pressure required to store propellant
    :param r: Radius of the propellant tank (therefore orbiter)
    :param sigma_y: Yield strength of the selected material
    :return: Required minimum thickness to withstand pressure loads
    """
    return p * r / (2 * sigma_y)

def t_hoop_cylind(p, r, sigma_y):
    """
    Function to calculate the thickness of the orbiter due to hoop stress in cylinder. To be used
    in case the propellant tanks need to have a cylindrical geometry with spherical caps.
    :param p: Pressure required to store propellant
    :param r: Radius of the propellant tank (therefore orbiter)
    :param sigma_y: Yield strength of the selected material
    :return: Required minimum thickness to withstand pressure loads
    """
    return p * r / sigma_y


def volume_tank(molar, R, T, p, m_total):
    n = molar * m_total
    return n * R * T / p


if __name__ == '__main__':
    # Constants
    g_earth = 9.81
    R = 8.314

    # Mass Budget
    m_subsystems = 472  # From preliminary design of concept 2
    m_structure = ...

    # Propellant Properties
    prop_pressure = ...
    prop_temp = ...
    prop_molar = ...

    # Material Properties
    material = {'Ti-6Al-4V': []}

    # Launch Loads
    acc_axial_tension = 6 * 9.81  # https://www.spacex.com/media/falcon-users-guide-2021-09.pdf
    acc_axial_compr = 2 * 9.81  # https://www.spacex.com/media/falcon-users-guide-2021-09.pdf
    acc_lateral = 2 * 9.81  # Same for compression and tension (https://www.spacex.com/media/falcon-users-guide-2021-09.pdf)
    acc_shock = 1000 * 9.81 # https://www.spacex.com/media/falcon-users-guide-2021-09.pdf


    # Launcher Constraints
    d_fairing = 5.2  # https://www.spacex.com/vehicles/falcon-heavy/
    h_fairing = 13.1  # https://www.spacex.com/vehicles/falcon-heavy/

    # Aeroshell Properties
    l_aeroshell = ...
    d_aeroshell = ...



