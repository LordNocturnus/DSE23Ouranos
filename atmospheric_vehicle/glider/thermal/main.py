def thermal_equilibrium(area_structure, area_structure_ss, area_body_payload, area_payload, area_subsystems,
                        temp_atm_min, temp_atm_max, temp_payload, temp_subsystems, t_skin, epsilon,
                        alpha_s, k, l_str, l_str_ss, alpha_IR, power, eta=0.85, j_s=3.69, temp_IR=58.2):
    """
    Function that determines the thermal equilibrium temperature of the glider
    :param area_structure: Cross sectional area of payload mounting structure [m^2]
    :param area_structure_ss: Cross sectional area of subsystem mounting structure [m^2]
    :param area_body_payload: Exposed payload surface area [m^2]
    :param area_payload: Area of the volume containing the payload [m^2]
    :param area_subsystems: Area of the volume containing the glider subsystems [m^2]
    :param temp_atm_min: Minimum temperature of the atmosphere [K]
    :param temp_atm_max: Maximum temperature of the atmosphere [K]
    :param temp_payload: Temperature of the payload [K]
    :param temp_subsystems: Temperature of the subsystems [K]
    :param t_skin: Skin thickness fuselage [m]
    :param epsilon: Body emissivity [-] 0.4 for Ti Al
    :param alpha_s: Solar absorptivity [-]
    :param k: Coefficient of thermal conductivity [W/mK]
    :param l_str: Length of the suspension rod [m]
    :param alpha_IR: IR absorptivity [-]
    :param power: Power usage of the glider [W]
    :param eta: Power usage efficiency [-]
    :param j_s: Solar intensity [W/m^2]
    :param temp_IR: Effective radiating temperature of planet [K]
    """

    sigma = 5.670374e-8  # Stefan-Boltzmann constant
    temp_out = temp_atm_min
    
    # Payload calculations
    q_dot_radiation = epsilon * sigma * area_payload * temp_payload**4  # Heat loss due to black body radiation
    q_dot_cond_str = k * area_structure * (temp_payload-temp_out) / l_str  # Heat loss due to conduction through mount

    q_dot_solar = alpha_s * j_s * area_payload/2  # Solar heating
    q_dot_IR = alpha_IR * sigma * temp_IR**4 * area_payload/2  # infrared radiation heating
    q_ineff_pl = power * (1 - eta)  # Internal heating from inefficiency (payload specific)


    q_req_pl = (q_dot_radiation + q_dot_cond_str - q_dot_solar - q_dot_IR - q_ineff_pl)



    # Subsystem calculations
    q_dot_radiation_ss = epsilon * sigma * area_subsystems * temp_subsystems ** 4  # Heat loss from black body radiation
    q_dot_cond_str_ss = k * area_structure_ss * (temp_subsystems - temp_out) / l_str_ss  # Heat loss through structure

    q_dot_solar_ss = alpha_s * j_s * area_subsystems / 2  # Solar heating
    q_dot_IR_ss = alpha_IR * sigma * temp_IR ** 4 * area_subsystems / 2  # infrared radiation heating
    q_ineff_ss = power * (1 - eta)  # Internal heating from inefficiency (subsystems specific)

    q_req_ss = (q_dot_radiation_ss + q_dot_cond_str_ss - q_dot_solar_ss - q_dot_IR_ss - q_ineff_ss)
    q_req = q_req_ss + q_req_pl
    q_req_max = q_req*1.1
    print(f"Required power: {q_req_max} [W]")

    print(f"---- Heat Flux Terms ----\n"
          f"q_radiation_pl: {round(q_dot_radiation, 3)} [W]     q_radiation_ss: {round(q_dot_radiation_ss, 3)} [W]\n"
          f"q_cond_pl: {round(q_dot_cond_str, 3)} [W]\n         q_cond_ss: {round(q_dot_cond_str_ss, 3)} [W]\n"
          f"q_solar_pl: {round(q_dot_solar, 3)} [W]\n           q_solar_ss: {round(q_dot_solar_ss, 3)} [W]\n"
          f"q_IR_pl: {round(q_dot_IR, 3)} [W]\n                 q_IR_ss: {round(q_dot_IR_ss, 3)} [W]\n"
          f"q_ineff_pl: {round(q_ineff_pl, 3)} [W]\n            q_ineff_ss: {round(q_ineff_ss, 3)} [W]\n"
          f"Max heating required: {q_req_max} [W]")

    # heater mass estimation
    mass_high = q_req_max/2e4 * 2  # [kg]
    mass_low = q_req_max/6e4 * 2  # [kg]
    print(f"Mass of the heater is approx {round(mass_low, 4)}-{round(mass_high, 4)} [kg]\n")

    # battery heating requirement


if __name__ == "__main__":
    area_s = 0.00001  # [m^2]
    area_s_ss = 0.0001  # [m^2]
    area_body_payload = 0.1  # [m^2]
    area_payload = 0.8
    area_subsystem = 1
    temp_low = 53
    temp_high = 193
    temp_payload = 245
    temp_subsys = 273
    t_skin_body = 0.0009  # [m]
    eps = 0.88  # found in SMAD page 422 of second edition
    alpha_s = 0.775  # found in SMAD page 422 of second edition
    k_mat = 7.1
    l_structure = 0.19
    l_structure_ss = 0.19
    alpha_IR = eps
    power = 134.3
    efficiency = 0.75

    thermal_equilibrium(area_s, area_s_ss, area_body_payload, area_payload, area_subsystem,
                        temp_low, temp_high, temp_payload, temp_subsys,
                        t_skin_body, 0.129, 0.448, k_mat, l_structure, l_structure_ss, 0.129, power, efficiency)
    thermal_equilibrium(area_s, area_s_ss, area_body_payload, area_payload, area_subsystem,
                        temp_low, temp_high, temp_payload, temp_subsys,
                        t_skin_body, 0.472, 0.766, k_mat, l_structure, l_structure_ss, 0.472, power, efficiency)
