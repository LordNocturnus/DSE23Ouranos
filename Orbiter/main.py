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

    def __init__(self):
        # Payload
        self.T_operational = 283.15  # Operational temperature of payload instruments in K
        self.m_payload = 26.77
        
        # Data Handling
        self.m_dh = 0

        # Power
        self.t_mission = 25  # Mission Timeline Tool
        self.P_comms = 60
        self.P_prop = 46
        self.P_adcs = 37.5
        self.P_dh = 37.5
        self.P_payload = 38.2
        self.P_thermal = 0
        self.P_pw = 30

        # Comms
        self.f_dl = 32
        self.f_ul = self.f_dl * TurnAroundRatio
        self.DR = 8000
        self.EbN0_dl = comms.downlink(self.P_comms, L_l, L_r, L_a, self.DR, Tnoisedown, k)
        self.EBN0_ul = comms.uplink(self.f_ul, P_gs, L_l, L_r, L_a, DR_ul, Tnoiseup, k)
        self.m_comms = 80

        # Structure and Prop
        self.mass = 2427  # Orbiter dry mass
        self.mixture_ratio = 1.65
        self.mass_AV = 471  # Atmospheric vehicle mass (import from AV class)
        self.mass_combined = self.mass + self.mass_AV  # Mass of combined systems
        self.deltaV_transfer = 1000  # Combined systems deltaV
        self.deltaV_insertion = 1000  # Delta V after splitting up at Uranus
        self.Isp = 321  # Isp of the orbiter thrusters
        self.T = 445  # Orbiter thrust
        self.material = [4430, 880 * 10**6, 970 * 10**6, 113.8 * 10**9]

        # ADCS
        self.m_adcs = 0

        # Iteration
        self.iteration()
        self.total_dry_mass = self.mass_combined
        self.total_wet_mass = self.total_dry_mass + self.prop_mass
        self.burn_transfer = prop.burntimecombined(self.T, self.orbiter_mass, self.deltaV_insertion, self.total_dry_mass, self.deltaV_transfer)
        self.burn_insertion = prop.burntimeorbiter(self.T, self.mass, self.deltaV_insertion)
        self.f_lat, self.f_ax = strt.natural_frequency(self.l_tanks, self.r_tanks, self.material, self.wet_mass_final, f_ax_min, f_lat_min)




    def mass_prop(self, m_dry):
        self.prop_mass = prop.mass_prop(m_dry, self.deltaV_insertion, self.mass_combined, self.deltaV_transfer,
                                        g, self.Isp)[2]
        self.m_fuel = self.prop_mass / (1 + self.mixture_ratio)
        self.m_ox = self.prop_mass - self.m_fuel
        self.prop_properties = [(1 * 10 ** 6, self.m_ox, 1431), (1 * 10 ** 6, self.m_fuel, 874)]
        self.wet_mass = m_dry + self.prop_mass
        self.l_tanks, self.r_tanks, self.m_structure = strt.final_architecture(self.material, self.prop_properties,
                                                                              margin, self.wet_mass)

    def power(self):
        self.n_rtg = pwr.numberRTG(self.P_req)
        self.m_power = pwr.massRTG(mass_RTG)
        self.cost_rtg = pwr.costRTG(costRTG1)

    def thermal(self):
        self.A_rec = np.pi * self.l_tanks * self.r_tanks
        self.A_emit = 2 * np.pi * self.r_tanks * self.l_tanks + np.pi * self.r_tanks ** 2
        self.m_louvres = 0.001 * np.pi * self.n_rtg * l_rtg * w_rtg * 2700
        self.d_rtg, self.n_l_closed = thm.power_phases(planets_list, r_orbit, self.A_rec, self.A_emit, alpha, epsilon,
                                                       self.n_rtg, p_rtg_tot, A_single_l)
        self.m_thermal = self.m_louvres + 50


    def iteration(self):
        diff = 1000
        while diff > 1:
            self.mass_prop(self.mass_combined)
            self.P_req = 500
            self.power()
            self.thermal()
            self.orbiter_mass = self.m_structure + self.m_power + self.m_thermal + self.m_payload + self.m_dh + self.m_comms + self.m_adcs
            new_mass_combined = self.mass_AV + self.orbiter_mass
            self.P_req = self.P_comms + self.P_pw + self.P_dh + self.P_adcs + self.P_payload + self.P_thermal + self.P_prop
            diff = abs(new_mass_combined - self.mass_combined)
            print(diff)
            self.mass_combined = new_mass_combined

if __name__ == "__main__":
    g = 9.81
    R = 8.314
    margin = 0.2
    boltzman = 5.67 * 10 ** (-8)  # Boltzamn constant for thermal
    k = 1.38 * 10 ** (-23)  # boltzmann constant for comms

    # --- STRUCTURE ---

    # acc_axial_tension = 6 * 9.81  # https://www.spacex.com/media/falcon-users-guide-2021-09.pdf
    # acc_axial_compr = 2 * 9.81  # https://www.spacex.com/media/falcon-users-guide-2021-09.pdf
    # acc_lateral = 2 * 9.81  # Same for compression and tension (https://www.spacex.com/media/falcon-users-guide-2021-09.pdf)
    # acc_shock = 1000 * 9.81  # https://www.spacex.com/media/falcon-users-guide-2021-09.pdf
    # f_lat_min = 10  # Falcon Heavy user manual
    # f_ax_min = 25  # Falcon Heavy user manual

    # Launcher Constraints
    # d_fairing = 5.2  # https://www.spacex.com/vehicles/falcon-heavy/
    # h_fairing = 13.1  # https://www.spacex.com/vehicles/falcon-heavy/

    # --- COMMS ---
    L_a = -0.5  # atmospheric attenuation in dB
    Tnoisedown = 424  # Noise temperature in K
    Tnoiseup = 763  # Noise temperature in K
    TurnAroundRatio = 3599 / 3344
    L_l = 0.9
    L_r = 0.75
    P_gs = 800
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

    orbiter = Orb()
