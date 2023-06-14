import math as m
dod = 1
n_bat = 0.9
n_cab = 0.9
vol_energy = 1000e3
spec_energy = 710 #https://www.greencarcongress.com/2023/04/20230429-li.html
margin = 0.2
days=3
nom_watt = 0.5
capacity = 19*2

dh_p_req = 0
ttc_p_req = 20
adcs_p_req = 0
pl_p_req = 96.1
tm_p_req = 30


def f_tot_p_req(dh_p_req, ttc_p_req,
                adcs_p_req, pl_p_req, tm_p_req):
    tot_p_req = dh_p_req + ttc_p_req + adcs_p_req + pl_p_req + tm_p_req
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

def f_n_bat(tot_p_req, t, dod, n_bat, n_cab, capacity):
    n_bat = m.ceil(tot_p_req * t / dod / n_bat / n_cab / capacity)
    return n_bat

if __name__ == "__main__":
    tot_p_req = f_tot_p_req(dh_p_req, ttc_p_req, adcs_p_req, pl_p_req, tm_p_req)
    t = f_t(days)
    print("Required energy [kWh]:", tot_p_req * t / dod / n_bat / n_cab/1000)
    n_bat = f_n_bat(tot_p_req, t, dod, n_bat, n_cab, capacity)
    print("Battery number [n]:", n_bat)
    print("Battery volume [m3]:", ((33e-3)/2)**2*m.pi*59e-3*n_bat)
    print("Battery mass [kg]:", n_bat * 72e-3)

