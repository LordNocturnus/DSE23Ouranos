"""
Script for the orbiter design tool
"""
import numpy as np
import Orbiter.propulsion as prop
import Orbiter.structure as strt
import Orbiter.comms as comms
import Orbiter.power as pwr
import Orbiter.thermal as thm
import Orbiter.ADCS as adcs

g = 9.80665
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

# --- ADCS --- #
capsule_radius = 1.5
capsule_height = 1.5
capsule_com = np.array([capsule_height * 3 / 4, 0, 0])
capsule_mass = 548.0

rho_uranus_max = 1.66e-9
rho_uranus_min = 6.902e-13
vel_per = 27000
c_d = 2.5
m_time = 3*365*24*3600

magnetic_dipole_uranus = 110e-6
orbital_period = 123*3600
grav_parameter = 5.79394e15
semi_major = 309000000
class Orb:

    def __init__(self, optimisation=True):
        # Payload
        self.T_operational = 283.15  # Operational temperature of payload instruments in K
        self.m_payload = 75.3
        
        # Data Handling
        self.m_dh = 1

        # Power
        self.t_mission = 23.6  # Mission Timeline Tool
        self.t_payload = 150  # Payload lifetime in months
        self.P_comms = 120
        self.P_prop = 35
        self.P_adcs = 150
        self.P_dh = 46
        self.P_payload = 147.04
        self.P_thermal = 0
        self.P_pw = 25

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
        self.mass_AV = 1000  # Atmospheric vehicle mass (import from AV class)
        self.mass_combined = self.mass + self.mass_AV  # Mass of combined systems
        self.deltaV_transfer = 170  # Combined systems deltaV
        self.deltaV_insertion = 2800 + 500 + 100  # Delta V after splitting up at Uranus, Moon discovery and ADCS
        self.Isp = 321  # Isp of the orbiter thrusters
        self.T = 425  # Orbiter thrust
        self.m_engine = 4.3  # Main engine mass in kg
        self.material = [4430, 880 * 10**6, 970 * 10**6, 113.8 * 10**9]

        # COST  page 297 SMAD
        self.launch_cost = 150000000
        self.cost_dh = 41400000  # Arnaud
        self.cost_comms = comms.total_cost(self.m_comms)
        self.cost_ADCS = ...
        self.cost_payload = 0  # (328 * self.m_payload**0.426 * self.P_payload**0.414 * self.t_payload**0.375 * 1000) * 1.34 * 0.951 * 1.39
            #1000000 + 31000000  # Magnetometer (https://gi.copernicus.org/preprints/gi-2017-53/gi-2017-53-AC1-supplement.pdf,
                                                # Camera (https://www.bhphotovideo.com/explora/photography/features/cameras-on-37-interplanetary-spacecraft)
        self.cost_ground_station = 2.01 * self.cost_dh  # SMAD
        self.cost_ops = 5980638 * self.t_mission + self.cost_ground_station  # SMAD


        # Iteration
        if optimisation:
            self.iteration()
            self.total_dry_mass = self.mass_combined
            self.wet_mass = self.total_dry_mass + self.prop_mass
            self.burn_transfer = prop.burntimecombined(self.T, self.mass, self.deltaV_transfer, self.total_dry_mass, self.deltaV_insertion, self.Isp)
            self.burn_insertion = prop.burntimeorbiter(self.T, self.mass, self.deltaV_insertion, self.total_dry_mass, self.deltaV_transfer, self.Isp)
            self.f_lat, self.f_ax = strt.natural_frequency(self.l_tanks, self.r_tanks, max(self.t_cy_o, self.t_cy_f), self.material, self.mass, self.mass_AV)
            self.cost_orbiter = self.cost_str + self.cost_thermal + self.cost_rtg + self.cost_comms + self.cost_prop + \
                                self.cost_dh + self.cost_payload + self.cost_test_assembly
            self.total_cost = (self.cost_orbiter + self.cost_test_assembly + self.cost_ops + self.launch_cost)\
                               * 1.2 / (10**6)  # Nasa Green Book

    def mass_prop(self, m_dry):
        self.prop_mass = prop.mass_prop(m_dry, self.deltaV_insertion, self.mass_combined, self.deltaV_transfer,
                                        self.Isp)[2] * 1.25
        self.m_fuel = self.prop_mass / (1 + self.mixture_ratio)
        self.m_ox = self.prop_mass - self.m_fuel
        self.prop_properties = [(3 * 10 ** 6, self.m_ox, 1431), (3 * 10 ** 6, self.m_fuel, 874)]
        self.wet_mass = m_dry + self.prop_mass
        self.l_tanks, self.r_tanks, self.m_structure, self.t_cy_o, self.t_cy_f, self.t_caps_o, self.t_caps_f, self.m_tanks = strt.final_architecture(self.material, self.prop_properties,
                                                                              margin, self.mass_AV)
        self.m_propulsion = self.m_tanks + self.prop_mass
        self.mainengine_burntime = self.prop_mass * self.Isp * 9.80665 / self.T
    def ADCS(self):
        cylinder_com = np.array([self.l_tanks / 2, 0, 0])
        self.mmoi = adcs.mmoi(self.r_tanks, self.l_tanks, np.array([self.l_tanks / 2, 0, 0]),
                              self.wet_mass - self.m_propulsion, self.r_tanks, self.m_tanks, self.m_ox,
                              np.array([self.r_tanks, 0, 0]), self.m_fuel, np.array([3 * self.r_tanks, 0, 0]),
                              capsule_radius, capsule_height, capsule_com, capsule_mass, tanks_full=False,
                              caps_attached=False, debug=False)
        self.reactionwheel_H = adcs.grav_grad_torque(self.mmoi, semi_major, grav_parameter, orbital_period,
                                                     debug=True) + \
                               adcs.mag_torque_max(magnetic_dipole_uranus, 1, orbital_period, debug=True)
        self.angular_momentum = abs(adcs.torque_s(self.mmoi, 1, debug=True)) * m_time + \
                                adcs.aerodyn_torque(rho_uranus_min, c_d, vel_per, self.mmoi, debug=True) * 220 * 15 * 60
        self.m_adcs_fuel = adcs.prop_mass(229, 4.33, self.l_tanks / 2, self.mainengine_burntime) + \
                           adcs.prop_mass(229, 4.33, self.l_tanks / 2, self.angular_momentum / (4.33 * self.r_tanks))
        self.m_adcs = self.m_adcs_fuel * 1.1 + 51.7
        self.mmoi_capsule = adcs.mmoi(self.r_tanks, self.l_tanks, np.array([self.l_tanks / 2, 0, 0]),
                                      self.mass - self.m_propulsion, self.r_tanks, self.m_tanks, self.m_ox,
                                      np.array([self.r_tanks, 0, 0]), self.m_fuel, np.array([3 * self.r_tanks, 0, 0]),
                                      capsule_radius, capsule_height, capsule_com, capsule_mass, tanks_full=True,
                                      caps_attached=True)[4]

    def power(self):
        self.n_rtg = pwr.numberRTG(self.P_req, self.t_mission)[1]
        self.m_power = 1.2 * (pwr.massRTG(self.P_req, self.t_mission) + 25)  # 25 kg is an estimate for PDU and regulators
        self.cost_rtg = pwr.costRTG(self.P_req, self.t_mission)

    def thermal(self):
        self.A_rec = self.l_tanks * self.r_tanks * 2
        self.A_emit = 2 * np.pi * self.r_tanks * self.l_tanks + np.pi * self.r_tanks ** 2
        self.d_rtg, self.m_radiator, self.m_louvres = thm.power_phases(self.A_rec, self.A_emit, self.n_rtg, T_operational=self.T_operational)
        self.m_kapton = self.A_emit * 0.001 * 1.55 * 1000
        self.m_thermal = self.m_louvres + self.m_radiator + self.m_kapton

    def iteration(self):
        diff = 1000
        while diff > 1 * 10 ** -3:
            self.mass_prop(self.mass)
            self.P_req = self.P_comms + self.P_pw + self.P_dh + self.P_adcs + self.P_payload + self.P_thermal + self.P_prop
            self.ADCS()
            self.power()
            self.thermal()
            new_orbiter_mass = self.m_structure + self.m_tanks + self.m_power + self.m_thermal + self.m_payload + self.m_dh + self.m_comms + self.m_adcs + self.m_engine
            diff = abs(new_orbiter_mass - self.mass)
            self.mass = new_orbiter_mass
        self.mass *= 1.25  # Nasa Green Book
        self.mass_combined = self.mass + self.mass_AV
        self.cost_prop = prop.total_cost(self.m_ox, self.m_fuel)
        self.cost_str = strt.total_cost(self.m_structure + self.m_tanks)
        self.cost_thermal = thm.total_cost(self.m_louvres, self.m_thermal)
        self.cost_test_assembly = 10.4 * 1000 * 1.7 * 0.951 * self.mass

    def __str__(self):
        return f'Total dry mass: {self.total_dry_mass}\n' \
               f'Total wet mass: {self.wet_mass}\n' \
               f'Total cost: {self.total_cost}\n' \
               f'Total Power: {self.P_req}\n' \
               f'Burn time transfer: {self.burn_transfer}\n' \
               f'Burn orbit insertion: {self.burn_insertion}'

    def mass_breakdwon(self):
        print(f'Orbiter Dry Mass: {self.mass}\n'
               f'Total Dry Mass: {self.total_dry_mass}\n'
               f'Orbiter Wet Mass: {self.wet_mass}\n'
               f'Propellant Mass: {self.prop_mass}\n'
               f'Atmospheric Vehicle Mass: {self.mass_AV}\n'
               f'Structure Mass: {self.m_structure}\n'
               f'Comms Mass: {self.m_comms}\n'
               f'Payload Mass: {self.m_payload}\n'
               f'DH Mass: {self.m_dh}\n'
               f'ADCS Mass: {self.m_adcs}\n'
               f'Power Mass: {self.m_power}\n'
               f'Thermal Mass: {self.m_thermal}\n'
               f'Tanks Mass: {self.m_tanks}')

    def cost_breakdown(self):
        print(f'Total cost: {self.total_cost} M€\n'
              f'Orbiter Cost: {self.cost_orbiter / 10**6} M€ \n'
              f'Launch cost: {self.launch_cost / 10**6} M€\n'
              f'Operational cost: {self.cost_ops / 10**6} M€\n'
              f'Ground station cost {self.cost_ground_station / 10**6} M€\n'
              f'Assembly cost: {self.cost_test_assembly / 10**6} M€\n'
              f'Payload cost: {self.cost_payload / 10**6} M€\n'
              f'Data Handling cost: {self.cost_dh / 10**6} M€\n'
              f'Propulsion cost: {self.cost_prop / 10**6} M€\n'
              f'Communication cost: {self.cost_comms / 10**6} M€\n'
              f'Power cost: {self.cost_rtg / 10**6} M€\n'
              f'Thermal cost: {self.cost_thermal / 10**6} M€\n'
              f'Structure cost: {self.cost_str / 10**6} M€\n')


if __name__ == "__main__":
    orbiter = Orb()
    print(orbiter.mass_breakdwon())
    print(orbiter.m_adcs_fuel)
