import matplotlib as plt




def stability_cg(C_l_alpha_h, C_l_alpha_ah, de_dalpha, l_h, chord, V_h, V, x_ac_bar, SM):
    S_h_S_stab = []

    S_h_S_stab = (1/((C_l_alpha_h/C_l_alpha_ah)*(1-de_dalpha)*(l_h/chord)*((V_h/V)**2))) - ((x_ac_bar - SM)/((C_l_alpha_h/C_l_alpha_ah)*(1-de_dalpha)*(l_h/chord)*((V_h/V)**2)))

    return S_h_S_stab

def controllability_cg(x_bar_cg, C_L_h, C_L_Ah, l_h, chord, V_h, V, C_m_ac, x_bar_ac):
    S_h_S_cont = []

    S_h_S_cont = (x_bar_cg/((C_L_h/C_L_Ah)*(l_h/chord)*((V_h/V)**2)))+(((C_m_ac/C_L_Ah)-x_bar_ac)/(((C_L_h/C_L_Ah)*(l_h/chord)*((V_h/V)**2))))

    return S_h_S_cont

def trim(C_m_ac, C_L_h, x_ac, chord):
    C_L_Ah = []
    d_x_cg = []
    for i in range(0, 91, 1):
    numbers = [round(x, 2) for x in range(0, 91, 1)]
    d_x_cg = [x / 100 for x in numbers]


        C_L_Ah += (-C_m_ac + C_L_h) * (chord / (d_x_cg[i] - x_ac))


    plt(d_x_cg,C_L_Ah, 'r=+')

    return C_L_Ah
trim(0.1,0.9, 5, 10)
###################################################################

if __name__ == "__main__":
    print("Hello World")
