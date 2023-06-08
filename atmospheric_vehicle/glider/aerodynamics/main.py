#Estimating Wing Geometry
#Parameters (wing assumed to be trapezoidal)
b_w = 6 #[m]
b_h = 1.3025 #[m]
W_S_avg = 224.482 #[N/m^2]
m = 112 #[kg]
g_ur = 9.01 #[m/s^2]
c_r_w = 0.785 #[m]
c_r_h = 1.5 #[m]
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
S_avg_h = 1.357 #[m]

#determine Mach Number at pressure levels considered:
M_low_bar = V_glider_low_bar/a_low_bar
M_high_bar = V_glider_high_bar/a_high_bar

#Weight, dynamic pressure determination:
W = m*g_ur
p_dynamic_low_bar = 0.5*rho_low_bar*(V_glider_low_bar**2)
p_dynamic_high_bar = 0.5*rho_high_bar*(V_glider_high_bar**2)

#determine Wing area and Tip Chord for Wing:
S_avg_w = (1/W_S_avg)*(W) #[m^2]
c_t_w = ((S_avg_w/2)-((b_w/4)*c_r_w))*(4/b_w) #[m]
c_avg_w = (S_avg_w/2)/(b_w/2)
A_w = (b_w**2)/S_avg_w

c_t_h = ((S_avg_h/2)-((b_h/4)*c_r_h))*(4/b_h)
c_avg_h = (S_avg_h/2)/(b_h/2)
A_h = (b_h**2)/S_avg_h


#Stability and Control derivatives + Parameters:
x_n_fixed = 0
C_m_delta_e = 0
C_m_alpha = 0
C_m_0 = 0
C_N_alpha = 0
alpha_0 = 0 #[deg]
alpha = []
C_l_alpha_ah = 4.02
C_l_alpha_h = 4.02
de_da = []
l_h = 2.765 #[m]
l_w = 1.5 #[m]

V_h_over_V = 1
SM = 0.15 #[%]
x_ac = l_w - 0.75*c_avg_w
x_cg = l_w - (c_avg_w - 0.333)


#determining de/da:
def downwash_vs_alpha(A):
    de_da = 4/(A+2)
    return de_da

#Stability:
def stability_cg(C_l_alpha_h, C_l_alpha_ah, de_dalpha, l_h, chord, V_h_over_V, x_ac_bar, SM, x_cg):
    S_h_S_stab = (x_cg/((C_l_alpha_h/C_l_alpha_ah)*(1-de_dalpha)*(l_h/chord)*(V_h_over_V**2))) - ((x_ac_bar - SM)/((C_l_alpha_h/C_l_alpha_ah)*(1-de_dalpha)*(l_h/chord)*((V_h_over_V)**2)))

    return S_h_S_stab

#Control:
def controllability_cg(x_bar_cg, C_L_h, C_L_Ah, l_h, chord, V_h, V, C_m_ac, x_bar_ac):
    S_h_S_cont = []

    S_h_S_cont = (x_bar_cg/((C_L_h/C_L_Ah)*(l_h/chord)*((V_h/V)**2)))+(((C_m_ac/C_L_Ah)-x_bar_ac)/(((C_L_h/C_L_Ah)*(l_h/chord)*((V_h/V)**2))))

    return S_h_S_cont

#Elevator trim:
'''''
delta_e = -(1/C_m_delta_e)*(C_m_0 + (W/(p_dynamic_low_bar*S_avg))*((x_cg-x_n_fixed)/c_avg))
delta_e = -(1/C_m_delta_e)*(C_m_0 + (W/(p_dynamic_high_bar*S_avg))*((x_cg-x_n_fixed)/c_avg))
'''''
print(stability_cg(C_l_alpha_h, C_l_alpha_ah, downwash_vs_alpha(A_w), l_h, c_avg_w, V_h_over_V, x_ac, SM, x_cg))

if __name__ == "__main__":
    print("Hello World")