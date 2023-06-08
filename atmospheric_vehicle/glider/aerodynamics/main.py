#Estimating Wing Geometry
#Parameters (wing assumed to be trapezoidal)

class aerodynamicproperties():
    def __init__(self):
        self.b_w = 6 #[m]
        self.b_h = 1.3025 #[m]
        self.W_S_avg = 224.482 #[N/m^2]
        self.m = 112 #[kg]
        self.g_ur = 9.01 #[m/s^2]
        self.c_r_w = 0.785 #[m]
        self.c_r_h = 1.5 #[m]
        self.Re_low_bar = 3.21*(10**6)
        self.Re_high_bar = 5.41*(10**6)
        self.Viscosity_low_bar = 2.29*(10**(-5)) #[m^2/s]
        self.Viscosity_high_bar = 2.51*(10**(-6)) #[m^2/s]
        self.a_low_bar = 560.9 #[m/s]
        self.a_high_bar = 992.2 #[m/s]
        self.Temp_low_bar = 54 #[k]
        self.Temp_high_bar = 193 #[k]
        self.alpha_low_range = -8.5 #[deg]
        self.alpha_high_range = 9.5 #[deg]
        self.V_glider_low_bar = 80 #[m/s]
        self.V_glider_high_bar= 0 #[m/s]
        self.rho_low_bar = 0.05 #[kg/m^3]
        self.rho_high_bar = 3.012 #[kg/m^3]
        self.S_avg_h = 1.357 #[m]

        #determine Mach Number at pressure levels considered:
        self.M_low_bar = self.V_glider_low_bar/self.a_low_bar
        self.M_high_bar = self.V_glider_high_bar/self.a_high_bar

        #Weight, dynamic pressure determination:
        W = self.m*self.g_ur
        self.p_dynamic_low_bar = 0.5*self.rho_low_bar*(self.V_glider_low_bar**2)
        self.p_dynamic_high_bar = 0.5*self.rho_high_bar*(self.V_glider_high_bar**2)

        #determine Wing area and Tip Chord for Wing:
        self.S_avg_w = (1/self.W_S_avg)*(W) #[m^2]
        self.c_t_w = ((self.S_avg_w/2)-((self.b_w/4)*self.c_r_w))*(4/self.b_w) #[m]
        self.c_avg_w = (self.S_avg_w/2)/(self.b_w/2)
        self.A_w = (self.b_w**2)/self.S_avg_w

        self.c_t_h = ((self.S_avg_h/2)-((self.b_h/4)*self.c_r_h))*(4/self.b_h)
        self.c_avg_h = (self.S_avg_h/2)/(self.b_h/2)
        self.A_h = (self.b_h**2)/self.S_avg_h


        #Stability and Control derivatives + Parameters:
        self.x_n_fixed = 0
        self.C_m_delta_e = 0
        self.C_m_alpha = 0
        self.C_m_0 = 0
        self.C_N_alpha = 0
        self.alpha_0 = 0 #[deg]
        self.alpha = []
        self.C_l_alpha_ah = 4.02
        self.C_l_alpha_h = 4.02
        self.de_da = []
        self.l_h = 2.765 #[m]
        self.l_w = 1.5 #[m]

        self.V_h_over_V = 1
        self.SM = 0.15 #[%]
        self.x_ac = self.l_w - 0.75*self.c_avg_w
        self.x_cg = self.l_w - (self.c_avg_w - 0.333)
        self.de_da = 4/(self.A_w+2)    



    #Stability:
    def stability_cg(self):
        S_h_S_stab = (self.x_cg/((self.C_l_alpha_h/self.C_l_alpha_ah)*
                                 (1-self.de_da)*(self.l_h/self.c_avg_w)*(self.V_h_over_V**2)))-( 
                                 (self.x_ac - self.SM)/((self.C_l_alpha_h/self.C_l_alpha_ah)*
                                 (1-self.de_da)*(self.l_h/self.c_avg_w)*((self.V_h_over_V)**2)))

        return S_h_S_stab

    #Control:
    def controllability_cg(self,x_bar_cg, C_L_h, C_L_Ah, l_h, chord, V_h, V, C_m_ac, x_bar_ac):
        S_h_S_cont = []

        S_h_S_cont = (x_bar_cg/((C_L_h/C_L_Ah)*(l_h/chord)*((V_h/V)**2)))+(((C_m_ac/C_L_Ah)-x_bar_ac)/(((C_L_h/C_L_Ah)*(l_h/chord)*((V_h/V)**2))))

        return S_h_S_cont

#Elevator trim:
'''''
delta_e = -(1/C_m_delta_e)*(C_m_0 + (W/(p_dynamic_low_bar*S_avg))*((x_cg-x_n_fixed)/c_avg))
delta_e = -(1/C_m_delta_e)*(C_m_0 + (W/(p_dynamic_high_bar*S_avg))*((x_cg-x_n_fixed)/c_avg))
'''''


if __name__ == "__main__":
    print("Hello World")
    glider = aerodynamicproperties()

    #print(glider.stability_cg(C_l_alpha_h, C_l_alpha_ah, downwash_vs_alpha(A_w), l_h, c_avg_w, V_h_over_V, x_ac, SM, x_cg))
    print(glider.stability_cg())