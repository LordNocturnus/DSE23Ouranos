"""Tool to implement the thermal management of the orbiter.. All geometrical aspects should be
derived from the orbiter structure tool.
"""
import numpy as np

# Constants
boltzman = 5.67 * 10 ** (-8)

# Planet list [Sun distance, radius, albedo factor, radiating temperature, *closest approach during gravity assist]
# Radius and distance to sun https://www.jpl.nasa.gov/edu/pdfs/scaless_reference.pdf
# Albedo and Temperature  ADSEE reader --> ECSS-S standards
planets_list = {'Uranus': [2872500000, 51118 / 2, 0.51, 58.2, 20],
                'Mars': [227900000, 6792 / 2, 0.15, 210.1, 200],
                'Jupiter': [778600000, 142984 / 2, 0.52, 109.5, 30011]}

# Inputs
alpha = 0.09  # Absorptivity (Aluminized Kapton foil from SMAD or ADSEE I reader)
epsilon = 0.8  # Emissivity (Aluminized Kapton foil from SMAD or ADSEE I reader)
T_operational = 283.15  # Operational temperature of payload instruments in K

# RTG properties
l_rtg = 1.14
w_rtg = 0.422
A_rtg = np.pi * w_rtg * l_rtg
p_rtg_tot = 4500
A_single_l = 0.05 * w_rtg * np.pi
# mass_louvres = 0.001 * np.pi * n_rtg * l_rtg * w_rtg * 2700  # Assumed to be an alluminum plate https://ntrs.nasa.gov/api/citations/20190028943/downloads/20190028943.pdf

# Alluminum properties https://material-properties.org/aluminium-thermal-properties-melting-point-thermal-conductivity-expansion/
cost_kg_l = 2.2 * 0.951  # https://markets.businessinsider.com/commodities/aluminum-price
cost_kg_shield = 14 * 0.951  # Rogers DK, Plucinski J. Low-cost carbon foams for thermal protection and reinforcement applications
                               # https://pdf.sciencedirectassets.com/271508/1-s2.0-S0008622309X00059/1-s2.0-S0008622309000591/mainext.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEKf%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJHMEUCIQDxP5dUHkIZpVDWdoLmplrWPGNnWOInPf%2F%2BRfNdalzoNwIgBjgGdocVkdQk30a%2FvmZsPmCCD0xrptYd2QSciD%2BG3fMquwUI7%2F%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FARAFGgwwNTkwMDM1NDY4NjUiDOnXkI7dqY7H4rUtXiqPBRlrcIFgnKjUPRnMvYAH47d3nwI2xOqp7YVGonEkT6n2Dh6Z%2ByStkcrgKtFNJ93lmfovgtuoGaIPS3lsOWXyG30IkQySAqQF71T0lfLFxsl2M6414gBrQV5AMdZpbPY4tUZ%2FTuly3d7Zry%2B%2FjegSNGxntTX%2BCszhkpA8ccHccOK77lq3JyNiSxqBNUIEhayZR9Nq%2F1idYa3Q4pAmWSvlOPUGNwWt%2B8uAVbc7FysokvLila6CD1riiMB%2FWOBb7GfQ3i0H%2BkQRJm65opIZVhkcuBGMWXaA62Vnnidy9xJhy5Vcvwy0qn%2F05Vov9FXzY%2BfAXfbznGYE2wT%2FI2M8HXbA6dnwexWgtSyALwNWdY2pYIdFpFPI4LvmXN6T88OO3twBIQCkWUjnka38iigw5oRbpTVP70MDu9GthR7ah%2BSxAtKiqKWy3gRRAnYPNuElWKqS1ZgJj492GFdzkl65j4KezLNDmLG7WmP9jFpJaCPWh9hXgsy65MFTDqTM4mqVjLXcbYHHjFS2ftNhGY4%2BIrMqKiAoBRXew58E6vDIUqwgeqDvezeZy0hJukoZxXqIKQKgHh%2By%2F2NMfhrIxZykNtvFHxSRccFdI5ijiLYFu5XmkfJ70%2FT7G%2B8U3Sm5hXHSJQFpx1abGu1OPzjf0Ubf7bXNKaKNCRd7lBEwUe58LKYXgRT%2FU2Oc9WgH2YuQLIl3UbyrIvZ%2FFtrJIQ133BvCqC5CBrw%2Ba8EwSX5SnNcQ8JLPxNcHS103J3%2BISePizRfIcCv7ovltmp0psDb7zqDRuI6tA67iFRrbK20Xv%2FSTmiwLSMf1EBQOvKCcbn3GPaWU4MY8V7xzBUg5rPp4elRu%2B0adXuLI67or9ecJkXnX4Eypzmgw2fOhpAY6sQH1S0c0oZ1LtPQ%2FK6Eha3mV8hc0KhyPjH97h29Q3ZbKXDjg%2BfZwCLXkC7WLYa8U%2BsccIggXruJfn%2FWI8SWjnnfkV8O7txCowOjPLvENJAP%2Fp4dLCBeaqquriTw8xjTwDL4o6vTBZ3oJEKOJenFT8yNiFgWFcFhv3oPWby4hqX7h%2Fi1PKtCBKbWL9EybOi%2B5lTxE3UIhNN4R3t3pb5c0yzkUnNQ5nJHfo%2FKkG3PSeBZa5uY%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20230613T151007Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYZFBG2VPP%2F20230613%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=0a8dd66e91273ba35e7cea5c66912228ace3c89b0262df90f02bb6a1faee5de6&hash=0a7caaa0557d0cb4cafecdad718772c429de7fe7fe83f676d4fa142050549394&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0008622309000591&tid=spdf-adc79775-744e-41d2-9e9f-c1e607edc507&sid=c010b063594b854fba0adef864f4dd9f4f48gxrqb&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=0b0a520b520055575307&rr=7d6b422e7fa30e14&cc=nl
# Constants
boltzman = 5.67 * 10 ** (-8)

def solar_intensity(r, planet, planet_list=planets_list):
    """
    Function to calculate the solar radiation on the spacecraft. Power of sun and distance of Uranus
    are taken, respectively, from SMAD and online research.
    :param r: Orbiter distance from the center of Uranus in km
    :param planet: Chosen planet where to calculate the solar intensity
    :return: Solar intensity J at Uranus
    """
    J_sun = 3.856 * 10**26 / (4 * np.pi * (planet_list[planet][0] * 1000 + r * 1000)**2)
    return J_sun


def albedo(J, r, planet, planet_list=planets_list):
    """
    Function to calculate the albedo flux received by the spacecraft in orbit. The albedo factor and radius of
    Uranus are derived, respectively, from ECSS-E-ST-10-04C and online research
    :param J: Solar intensity
    :param r: Radius of Uranus in km
    :param planet: Chosen planet where to calculate the albedo radiations
    :return: Albedo flux on orbiter
    """
    F = (planet_list[planet][1] / (r + planet_list[planet][1]))**2
    J_albedo = planet_list[planet][2] * J * F
    return J_albedo


def power_absorbed(r, A_rec, alpha, epsilon, planet, planet_list=planets_list, solar=True):
    """
    Function to compute the total power received by the orbiter. It is assumed the worst case scenario,
    therefore all radiations are present (solar, albedo, IR)
    :param r: Orbital distance (closest to planet) in km
    :param A_rec: Area of spacecraft receiving the radiations (assumed to be rectangular)
    :param alpha: Absorptivity of orbiter
    :param epsilon: Emissivity of orbiter
    :param planet: Chosen planet where to calculate the absorbed power
    :param solar: bool to determine if solar radiations should be taken into account or not
    :return: Total received heat
    """
    F = (planet_list[planet][1] / (r + planet_list[planet][1])) ** 2
    J_s = solar_intensity(r, planet, planet_list)
    J_a = albedo(J_s, r, planet, planet_list)
    J_ir = 5.67 * 10**(-8) * planet_list[planet][3]**4 * F  # ECSS

    P_s = alpha * J_s * A_rec
    P_a = alpha * J_a * A_rec
    P_ir = epsilon * J_ir * A_rec
    if solar:
        P_abs = P_s + P_a + P_ir
        return P_abs
    else:
        P_abs = P_a + P_ir
        return P_abs


def power_emitted(A_emit, epsilon, T_op):
    """
    Function computed the heat leaving the orbit during operations
    :param A_emit: Orbiter area that emits radiations
    :param epsilon: Emissivity of the orbiter material
    :param T_op: Operational temperature of the payload on board in K
    :return: Total orbiter emitted power
    """
    P_emit = epsilon * A_emit * (T_op**4) * 5.67 * 10**(-8)
    return P_emit


def power_dissipated(power_em, power_abs):
    """
    Function to calculate the required dissipated power in order to achieve thermal balance. A positive
    value means that the power has to be released while a negative power means that the power needs to
    be kept in the system.
    :param power_em: Orbiter emitted power
    :param power_abs: Orbiter absorbed power
    :return: Orbiter dissipated power
    """
    P_diss = power_abs - power_em
    return P_diss

def distance_rtg(n_rtg, p_rtg, p_diss, A_rec, alpha):
    """
    Function to calculate the required distance the RTG has to stay away from the orbiter subsystems.
    This is based on the required dissipated power of the orbiter
    :param n_rtg: Number of RTGs on board the orbiter
    :param p_rtg: Total power of each rtg
    :param p_diss: Dissipated power required to achieve thermal balance
    :param A_rec: Orbiter receiving area
    :param alpha: absorptivity
    :return: Value for the S/C to RTG distance
    """
    d_rtg = np.sqrt(A_rec * alpha * p_rtg * n_rtg / (4 * np.pi * (-1 * p_diss)))
    return d_rtg

def louvres_area(p_diss, A_rec, alpha, d_rtg_uranus, A_rtg, p_rtg_tot, n_rtg, A_single_l):
    """
    Function that determines the number of louvres that should be closed in order to produce the required
    dissipated power with the RTGs
    :param p_diss: Dissipated power of specific mission phase
    :param A_rec: Orbiter receiving area
    :param alpha: Absorptivity
    :param d_rtg_uranus: Distance S/C - RTG based on conditions at Uranus
    :param A_rtg: Surface area RTG
    :param p_rtg_tot: Total power of one RTG
    :param n_rtg: Number of RTGs on board
    :param A_single_l: Surface area of a single louvre blind
    :return: Percentage of surface are that should be covered by louvers
    """
    J_rtg_venus = - p_diss / (A_rec * alpha)
    p_rtg_venus = J_rtg_venus * 4 * np.pi * d_rtg_uranus**2
    a_tot_lv = A_rtg * p_rtg_venus / p_rtg_tot
    n_lv = np.ceil((A_rtg * n_rtg - a_tot_lv) / A_single_l)
    return n_lv


def power_phases(A_rec, A_emit, n_rtg, planet_list=planets_list, alpha=alpha, epsilon=epsilon, p_rtg_tot=p_rtg_tot, A_single_l=A_single_l, T_operational=T_operational, A_rtg=A_rtg):
    """
    Function that determines the rtg distance at Uranus and the number of closed louvres cells for the
    different mission phases. RTG distance is only computed for Uranus because this is the driving
    environment since we will perform mission there.
    :param planet_list: List of all planets characteristics [distance, radius, albedo factor, radiating temperature]
    :param r: Orbit radius at Uranus
    :param A_rec: Orbiter receiving area
    :param A_emit: Orbiter emitting area
    :param alpha: Absorptivity
    :param epsilon: Emissivity
    :param n_rtg: Number of RTGs
    :param p_rtg_tot: Total power of single RTG
    :param A_single_l: Surface area of a single louvre blind
    :return: Distance S/C - RTG and louvres percentage area for different phases.
    """
    areas = []
    d_rtg = 0
    for planet in planet_list:
        r = planet_list[planet][4]
        solar = False if planet == 'Venus' else True
        power_abs = power_absorbed(r, A_rec, alpha, epsilon, planet, planet_list, solar=solar)
        power_em = power_emitted(A_emit, epsilon, T_operational)
        power_diss = power_dissipated(power_em, power_abs)
        if planet == 'Uranus':
            d_rtg = distance_rtg(n_rtg, p_rtg_tot, power_diss, A_rec, alpha)
        else:
            areas.append((f'{planet}', louvres_area(power_diss, A_rec, alpha, d_rtg, A_rtg, p_rtg_tot, n_rtg, A_single_l)))
            #m_heater = - power_diss / (6 / 0.0001) * 2
    return d_rtg, areas

def total_cost(m_louvres, m_shield):
    return m_louvres * cost_kg_l + m_shield * cost_kg_shield + 1.70 * 50.6 * (m_louvres + m_shield)**0.707 * 1000 * 0.951 * 1.61


if __name__ == "__main__":
    # d_rtg, areas = power_phases(planets_list, r_orbit, A_rec, A_emit, alpha, epsilon, n_rtg, p_rtg_tot, A_single_l)
    # print(areas)
    power_diss = power_dissipated()

