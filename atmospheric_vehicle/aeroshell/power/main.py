dod = 1
n_bat = 0.9
n_cab = 0.9
vol_energy = 1653.65e3
spec_energy = 711.3  # https://www.greencarcongress.com/2023/04/20230429-li.html
margin = 0.5
days = 3

adcs_p_req = 20
pr_p_req = 20


def f_tot_p_req(adcs_p_req, pr_p_req):
    tot_p_req = adcs_p_req + pr_p_req
    return tot_p_req


def f_t(days):
    t = days * 24
    return t


def f_m_bat(tot_p_req, t, dod, n_bat, n_cab, spec_energy):
    m_bat = tot_p_req * t / dod / n_bat / n_cab / spec_energy
    return m_bat


def f_v_bat(tot_p_req, t, dod, n_bat, n_cab, vol_energy):
    v_bat = tot_p_req * t / dod / n_bat / n_cab / vol_energy
    return v_bat


def power_sizing(adcs_p_req, pr_p_req, flight_time):
    tot_p_req = f_tot_p_req(adcs_p_req, pr_p_req)
    t = f_t(flight_time / (24 * 3600))
    m_bat = f_m_bat(tot_p_req, t, dod, n_bat, n_cab, spec_energy)
    m_bat = m_bat * (1 + margin)
    v_bat = f_v_bat(tot_p_req, t, dod, n_bat, n_cab, spec_energy)
    v_bat = v_bat * (1 + margin)
    return m_bat, v_bat


if __name__ == "__main__":
    tot_p_req = f_tot_p_req(adcs_p_req, pr_p_req)
    t = f_t(days)
    m_bat = f_m_bat(tot_p_req, t, dod, n_bat, n_cab, spec_energy)
    m_bat = m_bat * (1 + margin)
    v_bat = f_v_bat(tot_p_req, t, dod, n_bat, n_cab, spec_energy)
    v_bat = v_bat * (1 + margin)
    print("Required energy [kWh]:", tot_p_req * t / 1000)
    print("Battery mass [kg]:", m_bat)
    print("Battery mass [kg]:", v_bat)
