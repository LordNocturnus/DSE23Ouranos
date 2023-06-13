import numpy as np

n_tx = 0.18  # -
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
    e_b = p_rx / bit_rate
    return e_b


def n_0(t_s):  # first estimation, update with https://deepspace.jpl.nasa.gov/dsndocs/810-005/105/105E.pdf
    return 1.380e-23 * t_s


def energy_per_bit_criteria(e_b, n_0):
    print(unit_to_db(e_b / n_0))
    print("Energy per bit criteria:", (e_b / n_0 >= 10))
    return e_b/n_0, e_b / n_0 >= 10


if __name__ == "__main__":
    g_tx_ao = unit_to_db(2)  # dB
    print("g_tx:", g_tx_ao)
    l_tx = unit_to_db(1 / n_tx)  # dB
    print("l_tx:", l_tx)
    freq_ao = 400e6  # Hz
    max_r_ao = 10000  # km
    l_fs_ao = l_fs(freq_ao, max_r_ao * 1000)  # dB
    print("l_fs:", l_fs_ao)
    r_rx_ao = 0.5  # m
    s_rx_ao = s(r_rx_ao)  # m^2
    g_rx_ao = unit_to_db(s_rx_ao)  # dB
    print("g_rx:", g_rx_ao)
    l_rx = unit_to_db(1 / n_rx)  # dB
    print("l_rx:", l_rx)
    l_m = unit_to_db(1 / 0.9)
    p_tx_ao = 20  # W
    p_tx_ao = unit_to_db(p_tx_ao)  # dB
    print("p_tx:", p_tx_ao)
    p_rx_ao = p_rx(p_tx_ao, g_tx_ao, l_tx, l_fs_ao, l_m, g_rx_ao, l_rx)
    bit_rate_ao = 2455  # bps
    print("p_rx:", db_to_unit(p_rx_ao))
    e_b_ao = e_b(db_to_unit(p_rx_ao), bit_rate_ao)
    print("e_b:", e_b_ao)
    print("e_b [dB]:", unit_to_db(e_b_ao))
    n_0 = n_0(103.15)
    print("n_0:", n_0)
    print("n_0 [dB]:", unit_to_db(n_0))
    energy_per_bit_criteria(e_b_ao, n_0)
