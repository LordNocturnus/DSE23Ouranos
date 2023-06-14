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

cost_kg = 45 * 0.92535  # https://www.navstarsteel.com/titanium-grade-5-sheet.html
manuf_cost = 132868.896  # Based on cost found in source and previous cst/kg of alluminum (check thermal) https://www.osti.gov/servlets/purl/62626

def axial_loads(axial, sigma_y, r):
    """
    Function to determine the base area of the orbiter based on launch loads
    :param axial: Axial peak loads during launch
    :param sigma_y: Yield stress of the selected material
    :return: Required minimum radius to withstand the loads
    """
    if r > 0 and sigma_y >= 0:
        return sigma_y * 0.9 >= abs(axial) / (np.pi * r**2)
    elif r < 0:
        raise ValueError(f'The radius needs to be a positive value')
    elif r == 0:
        raise ZeroDivisionError(f'The raidus cannot be 0')
    else:
        raise ValueError(f'The yield strength needs to be a positive value')

def lateral_loads(lateral, sigma_y, l):
    """
    Function to calculate the minimum required length of the orbiter based on launch loads.
    Worst case scenario assumed, with lateral loads fully applied at the tip of the structure.
    :param lateral: Lateral loads during launch
    :param sigma_y: Yield strength of the selected material
    :param r: Radius of the orbiter structure
    :return: Minimum required length to withstand lateral loads
    """
    if l > 0 and sigma_y >= 0:
        return sigma_y * 0.9 >= 6 * abs(lateral) / (l * l)
    elif l < 0:
        raise ValueError(f'The length of the orbiter has to be a positive value')
    else:
        raise ValueError(f'The yield strength of the material needs to be a positive value')

def t_hoop_sphere(p, r, sigma_y, margin):
    """
    Function to calculate the thickness of the orbiter due to hoop stress in sphere. Need to be used for
    spherical caps in case are present or for the entire structure if spherical tanks are implemented.
    :param p: Pressure required to store propellant
    :param r: Radius of the propellant tank (therefore orbiter)
    :param sigma_y: Yield strength of the selected material
    :return: Required minimum thickness to withstand pressure loads
    """
    if isinstance(r, int):
        if r >= 0 and sigma_y > 0 and margin >= 0:
            return abs(p) * r / (2 * sigma_y) * (1 + margin)
        elif r < 0:
            raise ValueError(f'The radius needs to have a positive value')
        elif sigma_y < 0:
            raise ValueError(f'The yield strength of the material needs to be a positive value')
        elif sigma_y == 0:
            raise ZeroDivisionError(f'The yield strength of a material cannot be 0')
        else:
            raise ValueError(f'The margin of safety needs to be a positive value')
    else:
        if any(r) >= 0 and sigma_y > 0 and margin >= 0:
            return abs(p) * r / (2 * sigma_y) * (1 + margin)
        elif any(r) < 0:
            raise ValueError(f'The radius needs to have a positive value')
        elif sigma_y < 0:
            raise ValueError(f'The yield strength of the material needs to be a positive value')
        elif sigma_y == 0:
            raise ZeroDivisionError(f'The yield strength of a material cannot be 0')
        else:
            raise ValueError(f'The margin of safety needs to be a positive value')
def t_hoop_cylind(p, r, sigma_y, margin):
    """
    Function to calculate the thickness of the orbiter due to hoop stress in cylinder. To be used
    in case the propellant tanks need to have a cylindrical geometry with spherical caps.
    :param p: Pressure required to store propellant
    :param r: Radius of the propellant tank (therefore orbiter)
    :param sigma_y: Yield strength of the selected material
    :return: Required minimum thickness to withstand pressure loads
    """
    if isinstance(r, int):
        if r >= 0 and sigma_y > 0 and margin >= 0:
            return abs(p) * r / sigma_y * (1 + margin)
        elif r < 0:
            raise ValueError(f'The orbiter radius needs to be a positive value')
        elif sigma_y < 0:
            raise ValueError(f'The yield strength of the material needs to be a positive value')
        elif sigma_y == 0:
            raise ZeroDivisionError(f'The yield strength of a material cannot be 0')
        else:
            raise ValueError(f'The margin of safety needs to be a positive value')
    else:
        if any(r) >= 0 and sigma_y > 0 and margin >= 0:
            return abs(p) * r / sigma_y * (1 + margin)
        elif any(r) < 0:
            raise ValueError(f'The orbiter radius needs to be a positive value')
        elif sigma_y < 0:
            raise ValueError(f'The yield strength of the material needs to be a positive value')
        elif sigma_y == 0:
            raise ZeroDivisionError(f'The yield strength of a material cannot be 0')
        else:
            raise ValueError(f'The margin of safety needs to be a positive value')


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
    if isinstance(material, list) and isinstance(propellant, tuple) and margin >= 0:
        # Calculate volume to store in the tanks
        v_prop = propellant[1] / propellant[2]
        v_tot = v_prop + pressurising_gas(propellant, v_prop)

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
            plt.subplot(121)
            plt.plot(r, m_structure)
            plt.subplot(122)
            plt.plot(r, l)
            plt.show()

        return m_structure[-1], l[-1], r[-1], max(t_caps[-1], 1 * 10**(-3)), max(t_cylind[-1], 1 * 10**(-3)), v_tot
    elif margin < 0:
        raise ValueError(f'The margin of safety needs to be a positive value')
    else:
        raise TypeError(f'The variable "propellant" should be a tuple and the variable "material" should be a list')

def minimum_thickness(m_AV, sigma_y, r):
    if m_AV >= 0 and r > 0 and sigma_y > 0:
        return acc_axial_tension * m_AV / (sigma_y * 2 * np.pi * r)
    elif r < 0:
        raise ValueError(f'The orbiter radius needs to be a positive value')
    elif sigma_y < 0:
        raise ValueError(f'The yield strength of the material needs to be a positive value')
    elif m_AV < 0:
        raise ValueError(f'The mass of the atmospheric vehicle cannot be negative')
    else:
        raise ZeroDivisionError(f'The yield strength and radius cannot be 0 valued')


def pressurising_gas(propellant, v_o):
    """
    Method to calculate the quantity of pressurising gas necessary in the propellant tanks.
    Formula taken from https://propulsion-skrishnan.com/pdf/N2O4-MMH%20Upper%20Stage%20Thruster.pdf
    :param propellant: Propellant properties, either the fuel or oxidiser
    :param v_o: Volume of the propellant, either fuel or oxidiser
    :return: Total volume of the pressurising gas
    """
    if isinstance(propellant, tuple):
        if propellant[0] > 0:
            return (0.9 * 10**6 / propellant[0]) * v_o / (1 - 0.9 * 10**6 / propellant[0])
        elif propellant[0] == 0:
            raise ZeroDivisionError('The propellant should be stored at a pressure higher than 0 Pa')
        else:
            raise ValueError(f'The propellant should be stored at a positive valued pressure')
    else:
        raise TypeError('The variable "propellant" should be a list with the properties of the propellant')


def final_architecture(material, propellant, margin, m_AV):
    """
    Function that defines the final architecture of the tanks based on the selected material, subsystem
    mass and propellant selected. Main function to run. Stress checks are performed.
    :param material: string with the name of the material selected for the tank
    :param propellant: list of properties of the propellant
    :param margin: safety margin for the thickness calculations
    :param m_subsystems: mass of the subsystems of the orbiter
    :return: prints the final properties of the orbiter structure
    """
    if isinstance(material, list) and isinstance(propellant, list):
        m_tank_o, l_o, r_o, t_caps_o, t_cylind_o, v_o = geometry_mass(material, propellant[0], margin, False)
        m_tank_f, l_f, r_f, t_caps_f, t_cylind_f, v_f = geometry_mass(material, propellant[1], margin, False)
        t_min_o = minimum_thickness(m_AV, material[1], r_o)
        t_min_f = minimum_thickness(m_AV, material[1], r_f)
        t_caps_o, t_cylind_o = max(t_caps_o, t_min_o), max(t_cylind_o, t_min_o)
        t_caps_f, t_cylind_f = max(t_caps_f, t_min_f), max(t_cylind_f, t_min_f)
        m_tanks = m_tank_o + m_tank_f
        m_tot_structure = ((r_o + r_f) * max(t_cylind_f, t_cylind_o) * 2 + 4 * r_o**2 * t_cylind_o +
                                                 4 * r_f**2 * t_cylind_f + 4 * max(r_f, r_o)**2 * 2 * t_cylind_f)\
                                                 * material[0] # Calculating the total mass of all metal surfaces on the orbiter
        l_tot = l_o + l_f + r_f * 2 + r_o * 2

        axial_check_tens = axial_loads(acc_axial_tension * m_AV, material[1], min(r_f, r_o))
        axial_check_compr = axial_loads(acc_axial_compr * m_AV, material[1], min(r_f, r_o))
        lateral_check = lateral_loads(acc_lateral * m_AV, material[1], l_tot)
        axial_shock = axial_loads(acc_shock * m_AV, material[2], min(r_f, r_o))
        if not axial_shock or not axial_check_compr or not axial_check_tens or not lateral_check:
            print(f'DANGER! Stress checks not passed:\n'
                  f'Axial Tension --> {axial_check_tens}\n'
                  f'Axial Compression --> {axial_check_compr}\n'
                  f'Axial Shock --> {axial_shock}')
        else:
            print(f'Stress checks passed')
        return l_tot, max(r_f, r_o), m_tot_structure, t_cylind_o, t_cylind_f, t_caps_o, t_caps_f, m_tanks
    else:
        raise TypeError(f'The variable "propellant" should be a list[tuple, ...] and the variable "material" should be a list[float, ...]')

def natural_frequency(l_tot, r, t, material, m_orb, m_AV):
    """
    Mathod to calculate the axial and lateral frequency of the idealised structure
    :param l_tot: Total length of the orbiter
    :param r: Radius of the orbiter
    :param material: Material list with material properties
    :param m_tot: Total mass of the orbiter
    :return: Lateral and axial frequency if they comply with the launcher constraints
    """
    if m_orb == 0 and m_AV == 0:
        raise ZeroDivisionError(f'One of the two masses should have a value different than 0')
    elif m_orb < 0 or m_AV < 0:
        raise ValueError(f'The mass cannot have a negative value')
    elif r <= 0 or t <= 0:
        raise ValueError(f'The dimensions of the orbiter need to have positive values')
    I = np.pi * t * (2 * r) ** 3 / 8
    f_lat = 0.276 * np.sqrt((material[3] * I) / (m_AV * l_tot**3 + 0.236 * m_orb * l_tot**3))
    f_ax = 0.160 * np.sqrt((np.pi * r**2 * material[3]) / (m_AV * l_tot + 0.333 * m_orb * l_tot))
    if f_ax > f_ax_min and f_lat > f_lat_min:
        print(f'Frequency checks passed')

    else:
        print(f'DANGER! Resonance could occur')
    return f_lat, f_ax

def total_cost(mass_tot):
    """
    Method to calcualte the total cost of the structure. Price taken from:
    :param mass_tot: Total orbiter mass
    :return: Total cost
    """
    return mass_tot * cost_kg + manuf_cost






if __name__ == '__main__':
    ...
    # Material Properties: https://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MTP641


