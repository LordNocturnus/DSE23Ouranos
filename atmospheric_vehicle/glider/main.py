def unit_to_db(unit):
    return 20*np.log10(unit)

def db_to_unit(db):
    return 10**(db/20)

def l_fs(freq, r):
    wavelength = 3e8 / freq
    l_fs = (20*np.log10(4*pi*r/wavelength))
    return lfs

def p_tx(min_p_rx, g_tx, l_tx, l_fs, l_m, g_rx, l_rx): #input in dB
    p_tx = min_p_rx - g_tx + l_tx + l_fs() + l_m - g_rx + l_rx
    return p_tx

def e_b(min_p_rx, bit_rate):
    return min_p_rx/bit_rate

def n_0(t_s): #first estimation, update with https://deepspace.jpl.nasa.gov/dsndocs/810-005/105/105E.pdf
    return 1.380e-23*t_s

def energy_per_bit_criteria(min_p_rx, bit_rate, t_s):
    e_b = e_b(min_p_rx, bit_rate)
    n_0 = n_0(t_s)
    return print("Energy per bit criteria:"+(e_b / n_0 >= 10))


if __name__ == "__main__":
    print("Hello World")