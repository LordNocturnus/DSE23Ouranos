def thermal_equilibrium(area_structure, area_body_payload, area_payload, area_subsystems,
                        temp_atm_min, temp_atm_max, temp_payload, temp_subsystems,
                        t_skin, epsilon, alpha_s, k, l_str, alpha_IR, power, eta=0.85, j_s=3.69, temp_IR=58.2):
    """
    Function that determines the thermal equilibrium temperature of the glider
    :param area_structure: Cross sectional area of the structure [m^2]
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
    q_dot_radiation = epsilon * sigma * area_payload * temp_payload**4
    q_dot_cond_str = k * area_structure * (temp_payload-temp_out) / l_str

    q_dot_solar = alpha_s * j_s * area_payload/2
    q_dot_IR = alpha_IR * sigma * temp_IR**4 * area_payload/2

    # Internal heating from inefficiency
    q_ineff = power * (1 - eta)

    # internal temperature
    print(f"q_solar: {q_dot_solar} [W]\n"
          f"q_IR: {q_dot_IR} [W]\n"
          f"q_radiation: {q_dot_radiation} [W]\n"
          f"q_cond: {q_dot_cond_str} [W]\n"
          f"q_ineff: {q_ineff} [W]\n")
    q_req = (q_dot_radiation + q_dot_cond_str - q_dot_solar - q_dot_IR - q_ineff)
    q_req_max = q_req*1.1
    print(f"Required power: {q_req_max} [W]")


    # mass estimation
    mass_high = q_req_max/2e4 * 2  # [kg]
    mass_low = q_req_max/6e4 * 2  # [kg]
    print(f"Mass of the heater is approx {round(mass_low, 4)}-{round(mass_high, 4)} [kg]\n")

    # battery heating requirement







if __name__ == "__main__":
    area_s = 0.00001  # [m^2]
    area_body_payload = 0.1  # [m^2]
    area_payload = 0.8
    area_subsystem = 1
    temp_low = 53
    temp_high = 193
    temp_payload = 245
    temp_subsys = 273
    eps = 0.88  # found in SMAD page 422 of second edition
    alpha_s = 0.775  # found in SMAD page 422 of second edition
    k_mat = 7.1
    t_skin_body = 0.0009  # [m]
    l_structure = 0.19
    alpha_IR = eps
    power = 255
    efficiency = 0.75

    thermal_equilibrium(area_s, area_body_payload, area_payload, area_subsystem,
                        temp_low, temp_high, temp_payload, temp_subsys,
                        t_skin_body, 0.129, 0.448, k_mat, l_structure, 0.129, power, efficiency)
    thermal_equilibrium(area_s, area_body_payload, area_payload, area_subsystem,
                        temp_low, temp_high, temp_payload,  temp_subsys,
                        t_skin_body, 0.472, 0.766, k_mat, l_structure, 0.472, power, efficiency)