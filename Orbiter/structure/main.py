# Tool to determine the sizing of the orbiter structure
import numpy as np
import matplotlib.pyplot as plt

# Constants

# Launch Loads
acc_axial_tension = 6 * 9.81  # https://www.spacex.com/media/falcon-users-guide-2021-09.pdf
acc_axial_compr = 2 * 9.81  # https://www.spacex.com/media/falcon-users-guide-2021-09.pdf
acc_lateral = 2 * 9.81  # Same for compression and tension (https://www.spacex.com/media/falcon-users-guide-2021-09.pdf)
acc_shock = 1000 * 9.81  # https://www.spacex.com/media/falcon-users-guide-2021-09.pdf
f_lat_min = 10  # Falcon Heavy user manual
f_ax_min = 25  # Falcon Heavy user manual

# Launcher Constraints
d_fairing = 5.2  # https://www.spacex.com/vehicles/falcon-heavy/
h_fairing = 13.1  # https://www.spacex.com/vehicles/falcon-heavy/


def axial_loads(axial, sigma_y, r):
    """
    Function to determine the base area of the orbiter based on launch loads
    :param axial: Axial peak loads during launch
    :param sigma_y: Yield stress of the selected material
    :return: Required minimum radius to withstand the loads
    """
    return sigma_y * 0.9 >= axial / (np.pi * r**2)

def lateral_loads(lateral, sigma_y, r, l):
    """
    Function to calculate the minimum required length of the orbiter based on launch loads.
    Worst case scenario assumed, with lateral loads fully applied at the tip of the structure.
    :param lateral: Lateral loads during launch
    :param sigma_y: Yield strength of the selected material
    :param r: Radius of the orbiter structure
    :return: Minimum required length to withstand lateral loads
    """
    return sigma_y * 0.9 >= 6 * lateral / (r * l)

def t_hoop_sphere(p, r, sigma_y, margin):
    """
    Function to calculate the thickness of the orbiter due to hoop stress in sphere. Need to be used for
    spherical caps in case are present or for the entire structure if spherical tanks are implemented.
    :param p: Pressure required to store propellant
    :param r: Radius of the propellant tank (therefore orbiter)
    :param sigma_y: Yield strength of the selected material
    :return: Required minimum thickness to withstand pressure loads
    """
    return p * r / (2 * sigma_y) * (1 + margin)

def t_hoop_cylind(p, r, sigma_y, margin):
    """
    Function to calculate the thickness of the orbiter due to hoop stress in cylinder. To be used
    in case the propellant tanks need to have a cylindrical geometry with spherical caps.
    :param p: Pressure required to store propellant
    :param r: Radius of the propellant tank (therefore orbiter)
    :param sigma_y: Yield strength of the selected material
    :return: Required minimum thickness to withstand pressure loads
    """
    return p * r / sigma_y * (1 + margin)


def geometry_mass(material, propellant, margin, plot=True):
    """
    Function to calculate all the geometrical property of the tank and its mass. Tank is
    assumed to be cylindrical with two semi spherical caps
    :param material: string with the name of the material selected for the tank
    :param propellant: list of properties of the propellant
    :param margin: safety margin for the thickness calculations
    :param plot: bool to determine if plotting is necessary or not
    :return: total structure mass, length of the tank, radius of the tank, thickness of spherical
    and cylindrical sections
    """

    # Calculate volume to store in the cylindrical part of the tank
    v_tot = propellant[1] / propellant[2]

    # Generate range of possible value for radius to use for calculations
    r = np.arange(0.1, 1.5, 0.05)

    # Calculate all geometrical properties based on geometry and propellant pressurisation
    l = (v_tot - 4/3 * np.pi * r**3) / (np.pi * r**2)
    r = r[l >= 0]
    l = l[l >= 0]
    t_caps = t_hoop_sphere(propellant[0], r, material[2], margin)
    t_cylind = t_hoop_cylind(propellant[0], r, material[2], margin)

    # Calculate final material volume and mass
    v_material = 2 * np.pi * r * t_cylind * l + 4 * np.pi * r**2 * t_caps
    m_structure = v_material * material[0]

    # Possibility to plot radius vs mass or radius vs length
    if plot:
        plt.subplot(131)
        plt.plot(r, m_structure)
        plt.subplot(132)
        plt.plot(r, l)
        plt.subplot(133)
        plt.plot(r, m_structure)
        plt.show()

    return m_structure[-1], l[-1], r[-1], max(t_caps[-1], 1 * 10**(-3)), max(t_cylind[-1], 1 * 10**(-3))

def final_architecture(material, propellant, margin, m_subsystems):
    """
    Function that defines the final architecture of the tanks based on the selected material, subsystem
    mass and propellant selected. Main function to run. Stress checks are performed.
    :param material: string with the name of the material selected for the tank
    :param propellant: list of properties of the propellant
    :param margin: safety margin for the thickness calculations
    :param m_subsystems: mass of the subsystems of the orbiter
    :return: prints the final properties of the orbiter structure
    """
    m_tank_o, l_o, r_o, t_caps_o, t_cylind_o = geometry_mass(material, propellant[0], margin, False)
    m_tank_f, l_f, r_f, t_caps_f, t_cylind_f = geometry_mass(material, propellant[1], margin, False)
    m_tot_structure = m_tank_o + m_tank_f
    l_tot = l_o + l_f + r_f * 2 + r_o * 2
    m_total = m_subsystems + m_tot_structure

    axial_check_tens = axial_loads(acc_axial_tension * m_total, material[1], min(r_f, r_o))
    axial_check_compr = axial_loads(acc_axial_compr * m_total, material[1], min(r_f, r_o))
    axial_shock = axial_loads(acc_shock * m_total, material[2], min(r_f, r_o))
    if not axial_shock or not axial_check_compr or not axial_check_tens:
        print(f'DANGER! Stress checks not passed:\n'
              f'Axial Tension --> {axial_check_tens}\n'
              f'Axial Compression --> {axial_check_compr}\n'
              f'Axial Shock --> {axial_shock}')
    # else:
    #     print(f'Stress checks passed')
    #     print(f'----------------------------------\n'
    #           f'Total Mass: {m_tot_structure}\n'
    #           f'Total Length: {l_tot}\n'
    #           f'Radius Oxidiser: {r_o}\n'
    #           f'Length oxidiser: {l_o}\n'
    #           f'Radius Fuel: {r_f}\n'
    #           f'Length Fuel: {l_f}\n'
    #           f'Thickness Oxidiser: {t_cylind_o, t_caps_o}\n'
    #           f'Thicknesses Fuel: {t_cylind_f, t_caps_f}\n'
    #           f'----------------------------------')
    return l_tot, max(r_f, r_o), m_tot_structure

def natural_frequency(l_tot, r, material, m_tot):
    I = 1/12 * 2 * r * l_tot ** 3
    f_lat = 0.276 * np.sqrt((material[3] * I) / (m_tot * l_tot**3 + 0.236 * m_tot * l_tot**3))
    f_ax = 0.160 * np.sqrt((np.pi * r**2 * material[3]) / (m_tot * l_tot + 0.333 * m_tot * l_tot))
    if f_ax > f_ax_min and f_lat > f_lat_min:
        print(f'Frequency checks passed')

    else:
        print(f'DANGER! Resonance could occur')
    return f_lat, f_ax






if __name__ == '__main__':
    ...
    # l_structure, r_structure, m_structure = final_architecture(material, propellant, margin, m_subsystems)
    # f_lat, f_ax = natural_frequency(l_structure, r_structure, material, m_subsystems + m_structure, f_ax_min, f_lat_min)



