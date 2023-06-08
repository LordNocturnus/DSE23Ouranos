import numpy as np
n_tx = 0.18 #-
n_rx = 0.6 #-

def unit_to_db(unit):
    return 20*np.log10(unit)

def db_to_unit(db):
    return 10**(db/20)

def l_fs(freq, r):
    wavelength = 3e8 / freq
    l_fs = (20*np.log10(4*np.pi*r/wavelength))
    return l_fs

def s(r):
    s = 2*np.pi*r*r
    return s

def p_rx(p_tx, g_tx, l_tx, l_fs, l_m, g_rx, l_rx): #input in dB
    p_rx = p_tx + g_tx - l_tx - l_fs - l_m + g_rx - l_rx
    return p_rx

def e_b(p_rx, bit_rate):
    e_b = p_rx/bit_rate
    return e_b

def n_0(t_s): #first estimation, update with https://deepspace.jpl.nasa.gov/dsndocs/810-005/105/105E.pdf
    return 1.380e-23*t_s

def energy_per_bit_criteria(e_b, n_0):
    return print("Energy per bit criteria:", (e_b / n_0 >= 10))


if __name__ == "__main__":
    g_tx_ao = unit_to_db(2) #dB
    l_tx = unit_to_db(1/n_tx) #dB
    freq_ao = 3e9 #Hz
    max_r_ao = 10000 #km
    l_fs_ao = l_fs(freq_ao, max_r_ao*1000) #dB
    r_rx_ao = 0.01 #m
    s_rx_ao = s(r_rx_ao) #m^2
    g_rx_ao = unit_to_db(s_rx_ao) #dB
    l_rx = unit_to_db(1/n_rx) #dB
    l_m = 0
    p_tx_ao = 1000 #W
    p_tx_ao = unit_to_db(p_tx_ao) #dB
    p_rx_ao = p_rx(p_tx_ao, g_tx_ao, l_tx, l_fs_ao, l_m, g_rx_ao, l_rx)
    bit_rate_ao = 2455 #bps
    print(db_to_unit(p_rx_ao))
    e_b_ao = e_b(db_to_unit(p_rx_ao), bit_rate_ao)
    print(e_b_ao)
    n_o = n_0(1000)
    print(n_o)
    energy_per_bit_criteria(e_b_ao, n_o)