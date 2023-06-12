"""
Script for the orbiter design tool
"""
import numpy as np
import Orbiter.propulsion as prop
import Orbiter.structure as strt
import Orbiter.comms as comms
import Orbiter.power as pwr
import Orbiter.thermal as thm


class Orb:

    def __init__(self, P_comms, d_antenna=5):
        # Payload
        self.T_operational = 283.15  # Operational temperature of payload instruments in K
        self.m_payload = ...
        
        # Data Handling

        # Comms
        self.d_antenna = d_antenna
        self.DR = ...
        self.f_dl = ...
        self.wavelength_dl = c / (self.f_dl * 10 ** 9)
        self.f_ul = self.f_dl * TurnAroundRatio
        self.wavelength_ul = c / (self.f_ul * 10 ** 9)
        self.EbN0_dl = comms.downlink(self.P_comms, L_l, L_r, L_a, self.DR, Tnoisedown, k)
        self.EBN0_ul = comms.uplink(self.f_ul, P_gs, L_l, L_r, L_a, DR_ul, Tnoiseup, k)
        self.m_comms = ...

        # Power
        self.t_mission = ...  # Mission Timeline Tool
        self.P_comms = P_comms
        self.P_prop = ...
        self.P_adcs = ...
        self.P_dh = ...
        self.P_payload = ...
        self.P_thermal = 0
        self.P_pw = ...
        self.P_req = self.P_comms + self.P_pw + self.P_dh + self.P_adcs + self.P_payload + self.P_thermal + self.P_prop
        self.n_rtg = pwr.numberRTG(self.P_req)
        self.m_power = pwr.massRTG(mass_RTG)
        self.cost_rtg = pwr.costRTG(costRTG1)

        # Structure and Prop
        self.mass = ...  # Orbiter dry mass
        self.mixture_ratio = 1.65
        self.mass_AV = ...  # Atmospheric vehicle mass (import from AV class)
        self.mass_combined = self.mass + self.mass_AV  # Mass of combined systems
        self.deltaV_transfer = ...  # Combined systems deltaV
        self.deltaV_insertion = ...  # Delta V after splitting up at Uranus
        self.Isp = ...  # Isp of the orbiter thrusters
        self.T = ...  # Orbiter thrust
        self.prop_mass = prop.mass_prop(self.mass, self.deltaV_insertion, self.mass_combined, self.deltaV_transfer, g, self.Isp)
        self.burn_transfer = prop.burntimecombined(self.T, self.mass_combined, self.deltaV_transfer)
        self.burn_insertion = prop.burntimeorbiter(self.T, self.mass, self.deltaV_insertion)
        self.m_fuel = self.prop_mass / (1 + self.mixture_ratio)
        self.m_ox = self.prop_mass - self.m_fuel
        self.prop_properties = [(1 * 10 ** 6, self.m_ox, 1431), (1 * 10 ** 6, self.m_fuel, 874)]
        self.material = [4430, 880 * 10**6, 970 * 10**6, 113.8 * 10**9]
        self.mass_iteration()
        self.dry_mass_final = self.mass_combined + self.m_structure
        self.wet_mass_final = self.dry_mass_final + self.prop_mass
        self.f_lat, self.f_ax = strt.natural_frequency(self.l_tanks, self.r_tanks, self.material, self.mass_final, f_ax_min, f_lat_min)


        # Thermal
        self.A_rec = np.pi * self.l_tanks * self.r_tanks
        self.A_emit = 2 * np.pi * self.r_tanks * self.l_tanks + np.pi * self.r_tanks ** 2
        self.m_louvres = 0.001 * np.pi * self.n_rtg * l_rtg * w_rtg * 2700
        self.d_rtg, self.n_l_closed = thm.power_phases(planets_list, r_orbit, self.A_rec, self.A_emit, alpha, epsilon, self.n_rtg, p_rtg_tot, A_single_l)
        self.m_thermal = self.m_louvres + ...


    def mass_iteration(self):
        m_structure = 0.2 * (self.mass_combined)
        diff = 1000
        while diff >= 1:
            self.prop_mass = prop.mass_prop(self.mass + m_structure, self.deltaV_insertion, self.mass_combined, self.deltaV_transfer,
                                            g, self.Isp)
            self.wet_mass = self.mass_combined + m_structure + self.prop_mass
            self.l_tanks, self.r_tanks, self.m_structure = strt.final_architecture(self.material, self.prop_properties,
                                                                                  margin, self.wet_mass)
            diff = abs(m_structure - self.m_structure)
            m_structure = self.m_structure

    
    
    
    
if __name__ == "__main__":
    g = 9.81
    R = 8.314
    margin = 0.2
    boltzman = 5.67 * 10 ** (-8)  # Boltzamn constant for thermal
    k = 1.38 * 10 ** (-23)  # boltzmann constant for comms

    # --- STRUCTURE ---

    acc_axial_tension = 6 * 9.81  # https://www.spacex.com/media/falcon-users-guide-2021-09.pdf
    acc_axial_compr = 2 * 9.81  # https://www.spacex.com/media/falcon-users-guide-2021-09.pdf
    acc_lateral = 2 * 9.81  # Same for compression and tension (https://www.spacex.com/media/falcon-users-guide-2021-09.pdf)
    acc_shock = 1000 * 9.81 # https://www.spacex.com/media/falcon-users-guide-2021-09.pdf
    f_lat_min = 10  # Falcon Heavy user manual
    f_ax_min = 25  # Falcon Heavy user manual

    # Launcher Constraints
    d_fairing = 5.2  # https://www.spacex.com/vehicles/falcon-heavy/
    h_fairing = 13.1  # https://www.spacex.com/vehicles/falcon-heavy/

    # --- COMMS ---
    c = 300000000  # speed of light in m/s
    earthRadius = 6371000.  # radius Earth in m
    L_a = -0.5  # atmospheric attenuation in dB
    AU = 149597870691  # one AU in m
    d_EarthSC = 20.8  # max distance between SC and Earth in AU
    Tnoisedown = 424  # Noise temperature in K
    Tnoiseup = 763  # Noise temperature in K
    TurnAroundRatio = 3599 / 3344

    eta_antenna = 0.55
    L_l = 0.9
    L_r = 0.75
    PointingAccuracy = 0.0572958
    P_gs = 800
    d_gs = 70
    DR_ul = 25000

    # --- POWER ---
    P_0 = 300  # Begin of life power one GPHS-RTG in W
    tau1 = 87.7  # half life fuel in years
    mass_RTG = 55.9  # mass of one GPHS-RTG in kg
    costRTG1 = 145699633.36  # cost of one GPHS-RTG in FY$2022, This is the highest value. It could be around 130 million as well

    # --- THERMAL ---
    planets_list = {'Uranus': [2872500000, 51118 / 2, 0.51, 58.2],
                    'Venus': [108200000, 12104 / 2, 0.65, 227, 200],
                    'Earth': [149600000, 12756 / 2, 0.44, 255, 200],
                    'Mars': [227900000, 6792 / 2, 0.15, 210.1, 200],
                    'Jupiter': [778600000, 142984 / 2, 0.52, 109.5, 200]}
    r_orbit = 200  # From Mission Timeline Tool
    alpha = 0.09  # Absorptivity (Aluminized Kapton foil from SMAD or ADSEE I reader)
    epsilon = 0.8  # Emissivity (Aluminized Kapton foil from SMAD or ADSEE I reader)
    l_rtg = 1.14
    w_rtg = 0.422
    A_rtg = np.pi * w_rtg * l_rtg
    p_rtg_tot = 4500
    A_single_l = 0.05 * w_rtg * np.pi