#Estimating Wing Geometry
#Parameters (wing assumed to be trapezoidal)
b = 6 #[m]
W_S_avg = 224.482 #[N/m^2]
m = 112 #[kg]
g_ur = 9.01 #[m/s^2]
c_r =0.785 #[m]
Re_low_bar = 3.21*(10**6)
Re_high_bar = 5.41*(10**6)
Viscosity_low_bar = 2.29*(10**(-5)) #[m^2/s]
Viscosity_high_bar = 2.51*(10**(-6)) #[m^2/s]
a_low_bar = 560.9 #[m/s]
a_high_bar = 992.2 #[m/s]
Temp_low_bar = 54 #[k]
Temp_high_bar = 193 #[k]
alpha_low_range = -8.5 #[deg]
alpha_high_range = 9.5 #[deg]
V_glider_low_bar = 80 #[m/s]
V_glider_high_bar= 0 #[m/s]
rho_low_bar = 0.05 #[kg/m^3]
rho_high_bar = 3.012 #[kg/m^3]

#determine Mach Number at pressure levels considered:
M_low_bar = V_glider_low_bar/a_low_bar
M_high_bar = V_glider_high_bar/a_high_bar

#Weight, dynamic pressure determination:
W = m*g_ur
p_dynamic_low_bar = 0.5*rho_low_bar*(V_glider_low_bar**2)
p_dynamic_high_bar = 0.5*rho_high_bar*(V_glider_high_bar**2)

#determine Wing area and Tip Chord for Wing:
S_avg = (1/W_S_avg)*(W) #[m^2]
c_t = ((S_avg/2)-((b/4)*c_r))*(4/b) #[m]
c_avg = (S_avg/2)/(b/2)


#Stability and Control derivatives + Parameters:
x_cg = 0
x_n_fixed = 0
C_m_delta_e = 0
C_m_alpha = 0
C_m_0 = 0
C_N_alpha = 0
alpha_0 = 0 #[deg]
alpha = []

#Stability:
def stability_cg(C_l_alpha_h, C_l_alpha_ah, de_dalpha, l_h, chord, V_h, V, x_ac_bar, SM):
    S_h_S_stab = []

    S_h_S_stab = (1/((C_l_alpha_h/C_l_alpha_ah)*(1-de_dalpha)*(l_h/chord)*((V_h/V)**2))) - ((x_ac_bar - SM)/((C_l_alpha_h/C_l_alpha_ah)*(1-de_dalpha)*(l_h/chord)*((V_h/V)**2)))

    return S_h_S_stab

#Control:
def controllability_cg(x_bar_cg, C_L_h, C_L_Ah, l_h, chord, V_h, V, C_m_ac, x_bar_ac):
    S_h_S_cont = []

    S_h_S_cont = (x_bar_cg/((C_L_h/C_L_Ah)*(l_h/chord)*((V_h/V)**2)))+(((C_m_ac/C_L_Ah)-x_bar_ac)/(((C_L_h/C_L_Ah)*(l_h/chord)*((V_h/V)**2))))

    return S_h_S_cont

#Elevator trim:
delta_e = -(1/C_m_delta_e)*(C_m_0 + (W/(p_dynamic_low_bar*S_avg))*((x_cg-x_n_fixed)/c_avg))
delta_e = -(1/C_m_delta_e)*(C_m_0 + (W/(p_dynamic_high_bar*S_avg))*((x_cg-x_n_fixed)/c_avg))
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

if __name__ == "__main__":
    print("Hello World")