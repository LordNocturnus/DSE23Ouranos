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
        self.m_dh = 20

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
        self.m_comms = 120

        # Structure and Prop
        self.mass = 2427  # Orbiter dry mass
        self.mixture_ratio = 1.65
        self.mass_AV = 500  # Atmospheric vehicle mass (import from AV class)
        self.mass_combined = self.mass + self.mass_AV  # Mass of combined systems
        self.deltaV_transfer = 1000  # Combined systems deltaV
        self.deltaV_insertion = 1000  # Delta V after splitting up at Uranus
        self.Isp = 321  # Isp of the orbiter thrusters
        self.T = 445  # Orbiter thrust
        self.material = [4430, 880 * 10**6, 970 * 10**6, 113.8 * 10**9]

        # ADCS
        self.m_adcs = 30

        # Iteration
        self.iteration()
        self.total_dry_mass = self.mass_combined
        self.total_wet_mass = self.total_dry_mass + self.prop_mass
        self.burn_transfer = prop.burntimecombined(self.T, self.mass, self.deltaV_insertion, self.total_dry_mass, self.deltaV_transfer)
        self.burn_insertion = prop.burntimeorbiter(self.T, self.mass, self.deltaV_insertion, self.total_dry_mass, self.deltaV_transfer)
        self.f_lat, self.f_ax = strt.natural_frequency(self.l_tanks, self.r_tanks, self.material, self.total_wet_mass)



    def mass_prop(self, m_dry):
        self.prop_mass = prop.mass_prop(m_dry, self.deltaV_insertion, self.mass_combined, self.deltaV_transfer,
                                        g, self.Isp)[2]
        self.m_fuel = self.prop_mass / (1 + self.mixture_ratio)
        self.m_ox = self.prop_mass - self.m_fuel
        self.prop_properties = [(3 * 10 ** 6, self.m_ox, 1431), (3 * 10 ** 6, self.m_fuel, 874)]
        self.wet_mass = m_dry + self.prop_mass
        self.l_tanks, self.r_tanks, self.m_structure = strt.final_architecture(self.material, self.prop_properties,
                                                                              margin, self.wet_mass, self.mass_AV)

    def power(self):
        self.n_rtg = pwr.numberRTG(self.P_req, self.t_mission)
        self.m_power = 1.2 * (pwr.massRTG(mass_RTG, self.P_req, self.t_mission) + 25)  # 25 kg is an estimate for PDU and regulators
        self.cost_rtg = pwr.costRTG(costRTG1, self.P_req, self.t_mission)

    def thermal(self):
        self.A_rec = np.pi * self.l_tanks * self.r_tanks
        self.A_emit = 2 * np.pi * self.r_tanks * self.l_tanks + np.pi * self.r_tanks ** 2
        self.m_louvres = 0.001 * np.pi * self.n_rtg * l_rtg * w_rtg * 2700
        self.d_rtg, self.n_l_closed = thm.power_phases(self.A_rec, self.A_emit, self.n_rtg)
        self.m_thermal = self.m_louvres + np.pi * self.r_tanks**2 * 0.1143 * 400  # 0.1143 thickness of shield https://science.nasa.gov/technology/technology-highlights/heat-shield-protect-mission-to-sun
                                                                                  # 400 is density of carbon phoam https://www.cfoam.com/wp-content/uploads/Carbon-Foams-amp16111p029-3.pdf


    def iteration(self):
        diff = 1000
        while diff > 1 * 10**-3:
            self.mass_prop(self.mass)
            self.P_req = 500
            self.power()
            self.thermal()
            new_orbiter_mass = 1.1 * (self.m_structure + self.m_power + self.m_thermal + self.m_payload + self.m_dh + self.m_comms + self.m_adcs)  # 1.1 is 10% margin based on smad page 316
            self.P_req = self.P_comms + self.P_pw + self.P_dh + self.P_adcs + self.P_payload + self.P_thermal + self.P_prop
            diff = abs(new_orbiter_mass - self.mass)
            self.mass = new_orbiter_mass
        self.mass_combined = self.mass + self.mass_AV

    def __str__(self):
        return f'Orbiter Dry Mass: {self.mass}\n' \
               f'Total Dry Mass: {self.total_dry_mass}\n' \
               f'Orbiter Wet Mass: {self.total_wet_mass}\n' \
               f'Propellant Mass: {self.prop_mass}\n' \
               f'Atmospheric Vehicle Mass: {self.mass_AV}\n' \
               f'Structure Mass: {self.m_structure}\n' \
               f'Comms Mass: {self.m_comms}\n' \
               f'Payload Mass: {self.m_payload}\n' \
               f'DH Mass: {self.m_dh}\n' \
               f'ADCS Mass: {self.m_adcs}\n' \
               f'Power Mass: {self.m_power}\n' \
               f'Thermal Mass: {self.m_thermal}\n' \
               f'Length Tanks: {self.l_tanks}\n' \
               f'Radius Tanks: {self.r_tanks}'

if __name__ == "__main__":
    g = 9.81
    R = 8.314
    margin = 0.2
    boltzman = 5.67 * 10 ** (-8)  # Boltzamn constant for thermal
    k = 1.38 * 10 ** (-23)  # boltzmann constant for comms

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
    l_rtg = 1.14
    w_rtg = 0.422

    orbiter = Orb()
    print(orbiter.r_tanks)
    print(str(orbiter))
