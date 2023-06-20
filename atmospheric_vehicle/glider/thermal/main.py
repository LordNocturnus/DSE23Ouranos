import numpy as np
import matplotlib.pyplot as plt

area_structural_crosssection_payload = 0.0001  # [m^2]
area_structural_crosssection_subsys = 0.0001  # [m^2]
area_payload = 0.579224  # [m^2]
area_subsystem = 0.565  # [m^2]
temp_vacuum = 3
temp_low = 53
temp_high = 193
temp_payload = 300  # temp at which the payload is held 253 min, 276.5 med, 300 max
temp_subsys = 293  # temp at which the subsystems are held 273 min, 283 med, 293 max
t_skin_body = 0.0009  # [m]
eps = 0.03  # found in SMAD page 422 of second edition
alpha_s = 0.775  # found in SMAD page 422 of second edition
k_mat = 7.1
k_atmosphere = 0.00052
l_structure = 0.19
l_structure_ss = 0.19
l_atmosphere_payload = 0.138
l_atmosphere_subsystems = 0.18
alpha_IR = eps
power = 255.6 * 0
efficiency = 0.99

def thermal_equilibrium(area_structure, area_structure_ss, area_payload, area_subsystems,
                        temp_atm, temp_payload, temp_subsystems, t_skin, epsilon,
                        alpha_s, k, k_atm, l_str, l_str_ss, l_atm_pl, l_atm_ss, alpha_IR, power, eta=0.85, j_s=3.69,
                        temp_IR=58.2, debug=False):
    """
    Function that determines the thermal equilibrium temperature of the glider
    :param area_structure: Cross sectional area of payload mounting structure [m^2]
    :param area_structure_ss: Cross sectional area of subsystem mounting structure [m^2]
    :param area_body_payload: Exposed payload surface area [m^2]
    :param area_payload: Area of the volume containing the payload [m^2]
    :param area_subsystems: Area of the volume containing the glider subsystems [m^2]
    :param temp_atm: Atmospheric temperature
    :param temp_payload: Temperature of the payload [K]
    :param temp_subsystems: Temperature of the subsystems [K]
    :param t_skin: Skin thickness fuselage [m]
    :param epsilon: Body emissivity [-] 0.4 for Ti Al
    :param alpha_s: Solar absorptivity [-]
    :param k: Coefficient of thermal conductivity [W/mK]
    :param k_atm: Coefficient of thermal conductivity atmosphere [W/mK]
    :param l_str: Length of the suspension rod [m]
    :param alpha_IR: IR absorptivity [-]
    :param power: Power usage of the glider [W]
    :param eta: Power usage efficiency [-]
    :param j_s: Solar intensity [W/m^2]
    :param temp_IR: Effective radiating temperature of planet [K]
    """

    sigma = 5.670374e-8  # Stefan-Boltzmann constant
    temp_out = temp_atm

    # Payload calculations
    q_dot_radiation_pl = epsilon * sigma * area_payload * temp_payload**4  # Heat loss due to black body radiation
    q_dot_cond_str_pl = k * area_structure * (temp_payload-temp_out) / l_str  # Heat loss due to conduction through str
    q_dot_cond_atm_pl = k_atm * area_payload * (temp_payload - temp_out) / l_atm_pl  # Heat loss through air
    q_dot_solar_pl = 0 * alpha_s * j_s * area_payload/2  # Solar heating
    q_dot_IR_pl = 0 * alpha_IR * sigma * temp_IR**4 * area_payload/2  # infrared radiation heating
    q_ineff_pl = 0 * power * (1 - eta)  # Internal heating from inefficiency (payload specific)
    if temp_out == 3:
        q_dot_cond_atm_pl = 0
    q_req_pl = (q_dot_radiation_pl + q_dot_cond_str_pl + q_dot_cond_atm_pl - q_dot_solar_pl - q_dot_IR_pl - q_ineff_pl)
    # Subsystem calculations
    q_dot_radiation_ss = epsilon * sigma * area_subsystems * temp_subsystems ** 4  # Heat loss from black body radiation
    q_dot_cond_str_ss = k * area_structure_ss * (temp_subsystems - temp_out) / l_str_ss  # Heat loss through structure
    q_dot_cond_atm_ss = k_atm * area_subsystems * (temp_subsystems - temp_out) / l_atm_ss  # Heat loss through air
    q_dot_solar_ss = 0 * alpha_s * j_s * area_subsystems / 2  # Solar heating
    q_dot_IR_ss = 0 * alpha_IR * sigma * temp_IR ** 4 * area_subsystems / 2  # infrared radiation heating
    q_ineff_ss = power * (1 - eta)  # Internal heating from inefficiency (subsystems specific)
    if temp_out == 3:
        q_dot_cond_atm_ss = 0
    q_req_ss = (q_dot_radiation_ss + q_dot_cond_str_ss + q_dot_cond_atm_ss -
                q_dot_solar_ss - q_dot_IR_ss - q_ineff_ss)
    q_req = q_req_ss + q_req_pl
    q_req_max = q_req * 1.1
    if debug:
        print(f"---- Heat Flux Terms ----\n"
              f"q_radiation_pl: {round(q_dot_radiation_pl, 3)} [W]     q_radiation_ss: {round(q_dot_radiation_ss, 3)} [W]\n"
              f"q_cond_str_pl: {round(q_dot_cond_str_pl, 3)} [W]       q_cond_str_ss: {round(q_dot_cond_str_ss, 3)} [W]\n"
              f"q_dot_cond_atm_pl: {round(q_dot_cond_atm_pl, 5)}       q_dot_cond_atm_ss: {round(q_dot_cond_atm_ss, 5)}\n"
              f"q_solar_pl: {round(q_dot_solar_pl, 3)} [W]             q_solar_ss: {round(q_dot_solar_ss, 3)} [W]\n"
              f"q_IR_pl: {round(q_dot_IR_pl, 3)} [W]                   q_IR_ss: {round(q_dot_IR_ss, 3)} [W]\n"
              f"q_ineff_pl: {round(q_ineff_pl, 3)} [W]                 q_ineff_ss: {round(q_ineff_ss, 3)} [W]\n"
              f"Heating required: {q_req} [W]")

    # heater mass estimation
    mass_high = q_req_max/2e4 * 2  # [kg]
    mass_low = q_req_max/6e4 * 2  # [kg]
    # print(f"Mass of the heater is approx {round(mass_low, 4)}-{round(mass_high, 4)} [kg]\n")
    return q_req
    # battery heating requirement


if __name__ == "__main__":
    thermal_equilibrium(area_structural_crosssection_payload, area_structural_crosssection_subsys, area_payload,
                        area_subsystem, temp_high, temp_payload, temp_subsys,
                        t_skin_body, eps, alpha_s, k_mat, k_atmosphere, l_structure, l_structure_ss,
                        l_atmosphere_payload, l_atmosphere_subsystems, alpha_IR, power, efficiency, debug=True)

temperatures_pl = np.linspace(253, 300, 20)
temperatures_ss = np.linspace(273, 293, 20)
temperatures_air = np.linspace(temp_low, temp_high, 20)
# y = []
# for t in temperatures_ss:
#     y.append(thermal_equilibrium(area_structural_crosssection_payload, area_structural_crosssection_subsys,
#                                  area_payload, area_subsystem, temp_low, t, temp_subsys, t_skin_body, eps, alpha_s,
#                                  k_mat, k_atmosphere, l_structure, l_structure_ss, l_atmosphere_payload,
#                                  l_atmosphere_subsystems, alpha_IR, power, efficiency))

# y = []
# for air_temperature in temperatures_air:
#     y.append(thermal_equilibrium(area_structural_crosssection_payload, area_structural_crosssection_subsys,
#                                  area_payload, area_subsystem, air_temperature, temp_payload, temp_subsys,
#                                  t_skin_body, eps, alpha_s, k_mat, k_atmosphere, l_structure, l_structure_ss,
#                                  l_atmosphere_payload, l_atmosphere_subsystems, alpha_IR, power, efficiency))
# plt.plot(y, temperatures_ss)
# plt.show()