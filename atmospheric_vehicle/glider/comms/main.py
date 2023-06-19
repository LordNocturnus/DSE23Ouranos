import numpy as np

n_tx = 0.933  # -
n_rx = 0.55  # -


def unit_to_db(unit):
    return 10 * np.log10(unit)


def db_to_unit(db):
    return 10 ** (db / 10)


def l_fs(freq, r):
    wavelength = 3e8 / freq
    l_fs = (20 * np.log10(4 * np.pi * r / wavelength))
    return l_fs


def s(r):
    s = np.pi * r * r
    return s


def p_rx(p_tx, g_tx, l_tx, l_fs, l_m, g_rx, l_rx):  # input in dB
    p_rx = p_tx + g_tx - l_tx - l_fs - l_m + g_rx - l_rx
    return p_rx


def e_b(p_rx, bit_rate):
    e_b = p_rx / bit_rate * (2351/100)
    return e_b


def n_0(t_s):  # first estimation, update with https://deepspace.jpl.nasa.gov/dsndocs/810-005/105/105E.pdf
    return 1.380e-23 * t_s


def energy_per_bit_criteria(e_b, n_0):
    print(unit_to_db(e_b / n_0))
    print("Energy per bit criteria:", (e_b / n_0 >= 10))
    return e_b/n_0, e_b / n_0 >= 3


if __name__ == "__main__":
    g_tx_ao =1  # dB
    print("g_tx:", g_tx_ao)
    l_tx = unit_to_db(1 / n_tx)  # dB
    print("l_tx:", -l_tx)
    freq_ao = 330e6
    #freq_ao = 400e6  # Hz
    max_r_ao = 2.8e12
    #max_r_ao = 252401691  # m
    l_fs_ao = l_fs(freq_ao, max_r_ao)  # dB
    print("l_fs:", -l_fs_ao)
    r_rx_ao = 2.4  # m
    s_rx_ao = s(250)+30*s(45/2)+39*s(12.5)+2*s(50)+s(76/2)
    #s_rx_ao = s(r_rx_ao)  # m^2
    g_rx_ao = unit_to_db(s_rx_ao)  # dB
    print("g_rx:", g_rx_ao)
    l_rx = unit_to_db(1 / n_rx)  # dB
    print("l_rx:", -l_rx)
    l_m =0 # dB
    p_tx_ao = 55 # W
    p_tx_ao = unit_to_db(p_tx_ao)  # dB
    print("p_tx:", p_tx_ao)
    p_rx_ao = p_rx(p_tx_ao, g_tx_ao, l_tx, l_fs_ao, l_m, g_rx_ao, l_rx)
    bit_rate_ao = 100
    #bit_rate_ao = 2351#*(1.481*1.1+1)  # bps
    print("p_rx:", db_to_unit(p_rx_ao))
    e_b_ao = e_b(db_to_unit(p_rx_ao), bit_rate_ao)
    print("e_b:", e_b_ao)
    print("e_b [dB]:", unit_to_db(e_b_ao))
    n_0 = n_0(103.34)
    print("n_0:", n_0)
    print("n_0 [dB]:", unit_to_db(n_0))
    energy_per_bit_criteria(e_b_ao, n_0)
    connected =  [(0, 48900), (89690, 118030), (210580, 238890), (279800, 321630), (326270, 349940)]
    nonconnected = [(48900, 89690), (118030, 210580), (238890, 279800), (321630, 326270)]
    #for i in range(len(connected)):
     #   print(connected[i][1]-connected[i][0])
    #print(".")
    #for i in range(len(nonconnected)):
       # print(nonconnected[i][1]-nonconnected[i][0])
    #for i in range(len(nonconnected)):
     #   print(((nonconnected[i][1]-nonconnected[i][0])/(connected[i+1][1]-connected[i+1][0])))