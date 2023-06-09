#Library Import
import numpy as np

#Estimating Wing Geometry
#Parameters (wing assumed to be trapezoidal)


class aerodynamicproperties():
    def __init__(self):
        self.b_w = 6 #[m]
        self.b_h = 1.3025 #[m]
        self.W_S_avg = 224.482 #[N/m^2]
        self.m = 112 #[kg]
        self.oswald = 0.8
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
        self.C_L_over_C_D = 38.202 # value obtained from Noah's atmospheric flight model
        self.C_L_max_range = 0.45 # value obtained from literature
        self.C_m_0 = 0
        self.C_m_ac = self.C_L_max_range/self.C_L_over_C_D

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
        self.C_m_alpha = -0.009
        self.C_m_0 = 0
        self.C_N_w_alpha = 0
        self.C_N_h_alpha = 0
        self.alpha_0 = 0 #[deg]
        self.alpha = []
        self.C_l_alpha_ah = 4.02
        self.C_l_alpha_h = 4.02
        self.l_h = 2.765 #[m]
        self.l_w = 1.5 #[m]
        self.C_L_h = 0.45
        self.C_L_ah = 0.45
        self.C_X_0 = 0
        self.C_Z_0 = 0
        self.d_C_D_over_d_V = 0
        self.d_C_L_over_d_V = 0
        self.C_m_u = 0
        self.delta_C_L = 0
        self.C_Z_0 = 0 # C_Z at 0 alpha
        self.C_l_p = 0 # read from table if necessary
        self.C_Y_v_alpha = 0
        self.dtheta_dbeta = 0
        self.sweep = 0 #[deg]
        self.beta = 0 #[rad]
        self.dL_dV = 0
        self.L = 0 #[N]
        self.N = 0 #[N]
        self.p = 0 #[rad/s or deg/s] roll rate
        self.C_Y_r_v = 0
        self.r = 0 ##[rad/s or deg/s] yaw rate
        self.C_L = 0
        self.Y_v = 0 #[N]
        self.N = 0 #[N]


        self.V_h_over_V = 1
        self.SM = 0.15 #[%]
        self.x_ac = self.l_w - 0.75*self.c_avg_w
        self.x_cg = self.l_w - (self.c_avg_w - 0.333)
        self.de_da = 4/(self.A_w+2)
        self.z_cg = 0  # [m]
    def mass_properties(self,mass,x_cg,I_xx,I_yy,I_zz,J_xz):
        self.m = mass
        self.x_cg = x_cg
        self.K_x_2 = np.sqrt(I_xx/mass)
        self.K_y_2 = np.sqrt(I_yy/mass)
        self.K_z_2 = np.sqrt(I_zz/mass)
        self.K_x_z = J_xz / mass
        self.W = self.mass * 9.01 #[m/s^2]

    def mainwing(self,wingspan,wing_loading,taper_ratio,x_w):

        self.b_w = wingspan
        self.S_w = wing_loading * self.m
        self.c_avg_w = self.S_w/wingspan
        self.c_t_w = self.c_avg_w /(taper_ratio+1)  
        self.c_r_w = taper_ratio * self.c_t_w      
        self.V = 0 #[m^3]
        self.x_w = x_w

    def horizontal_tail(self,wingspan,surface_area,taper_ratio,l_h,z_h):
        self.b_h = wingspan
        self.S_h = surface_area
        self.c_avg_h = surface_area/wingspan
        self.c_t_h = self.c_avg_h /(taper_ratio+1)  
        self.c_r_h = taper_ratio * self.c_t_h     
        self.l_h = l_h
        self.z_h = z_h

    def vertical_tail(self, l_v, x_v, S_v, V_v, z_v):
        self.l_v = l_v #[m]
        self.x_v = x_v  # [m]
        self.z_v = z_v  # [m]
        self.S_v = S_v #[m^2]
        self.V_v = V_v #[m^3]

    def atmospheric_properties(self,speed_of_sound,reynolds_number,viscosity,temperature,density,velocity):
        self.Re_low_bar = reynolds_number[0]
        self.Re_high_bar = reynolds_number[-1]
        self.Viscosity_low_bar = viscosity[0]
        self.Viscosity_high_bar = viscosity[-1]
        self.a_low_bar = speed_of_sound[0]
        self.a_high_bar = speed_of_sound[-1]
        self.Temp_low_bar = temperature[0]
        self.Temp_high_bar = temperature[-1]
        self.V_glider_low_bar = velocity[0]
        self.V_glider_high_bar= velocity[-1]
        self.rho_low_bar = density[0]
        self.rho_high_bar = density[-1]       
        self.p_dynamic_low_bar = 0.5*self.rho_low_bar*(self.V_glider_low_bar**2)
        self.p_dynamic_high_bar = 0.5*self.rho_high_bar*(self.V_glider_high_bar**2)

    #Stability:
    def stability_cg(self):
        S_h_S_stab = (self.x_cg/((self.C_l_alpha_h/self.C_l_alpha_ah)*
                                 (1-self.de_da)*(self.l_h/self.c_avg_w)*(self.V_h_over_V**2)))-( 
                                 (self.x_ac - self.SM)/((self.C_l_alpha_h/self.C_l_alpha_ah)*
                                 (1-self.de_da)*(self.l_h/self.c_avg_w)*((self.V_h_over_V)**2)))

        return S_h_S_stab

    #Control:
    def controllability_cg(self):

        S_h_S_cont = ((self.x_cg/self.c_avg_w)/
                      ((self.C_L_h/self.C_L_ah)*(self.l_h/self.c_avg_w)*
                       (self.V_h_over_V**2)))+(((self.C_m_ac/self.C_L_ah)-(self.x_ac/self.c_avg_w))/(((self.C_L_h/self.C_L_ah)*(self.l_h/self.c_avg_w)*(self.V_h_over_V**2))))

        return S_h_S_cont

    #Stability Derivatives with respect to airspeed
    def stab_derivatives_speed(self):
        C_X_u_lowbar = 2*self.C_X_0 - self.d_C_D_over_d_V*self.V_glider_low_bar
        C_Z_u_lowbar = 2 * self.C_Z_0 - self.d_C_L_over_d_V * self.V_glider_low_bar

        C_X_u_highbar = 2 * self.C_X_0 - self.d_C_D_over_d_V * self.V_glider_high_bar
        C_Z_u_highbar = 2 * self.C_Z_0 - self.d_C_L_over_d_V * self.V_glider_high_bar

        return C_X_u_lowbar, C_X_u_highbar, C_Z_u_lowbar, C_Z_u_highbar

    # Stability Derivatives with respect to angle of attack
    def stab_derivatives_alpha(self,S_h_S):
        C_X_alpha = self.C_L_ah*(1-((2*self.C_l_alpha_ah)/(np.pi*self.A_w*self.oswald)))
        C_Z_alpha = -self.C_N_w_alpha - self.C_N_h_alpha*(1-self.de_da)*(self.V_h_over_V**2)*(S_h_S)
        C_m_alpha = self.C_N_w_alpha*((self.x_cg-self.l_w)/self.c_avg_w) - self.C_N_h_alpha*(1-self.de_da)*(self.V_h_over_V**2)*(S_h_S)*(self.l_h/self.c_avg_w)

        return C_X_alpha, C_Z_alpha, C_m_alpha

    # Stability Derivatives with respect to sideslip angle
    def stab_derivatives_beta(self):
        C_Y_beta  = self.C_Y_v_alpha*(1-self.dtheta_dbeta)*(self.V_h_over_V**2)*(self.S_v/self.S_w)
        C_l_beta_lowbar = (2*self.dL_dV)/(self.rho_low_bar*(self.V_glider_low_bar**2)*self.S_w*np.sin(2*self.sweep)*self.beta) + C_Y_beta*(((self.z_v - self.z_cg)/self.b_w)*np.cos(self.alpha_0) - ((self.x_v - self.x_cg)/self.b_w)*np.sin(self.alpha_0))
        C_n_beta = self.C_Y_v_alpha*(1-self.dtheta_dbeta)*((self.V_v/self.V)**2)*((self.S_v*self.l_v)/(self.S_w*self.b_w))

        C_l_beta_highbar = (2 * self.dL_dV) / (self.rho_high_bar * (self.V_glider_high_bar ** 2) * self.S_w * np.sin(
            2 * self.sweep) * self.beta) + C_Y_beta * (((self.z_v - self.z_cg) / self.b_w) * np.cos(self.alpha_0) - (
                    (self.x_v - self.x_cg) / self.b_w) * np.sin(self.alpha_0))

        return C_Y_beta, C_l_beta_lowbar, C_l_beta_highbar, C_n_beta

    # Stability Derivatives with respect to roll rate
    def stab_derivatives_rollrate(self):
        C_n_p_lowbar = self.N/(((self.p*self.b_w)/2*self.V_glider_low_bar)*0.5*self.rho_low_bar*(self.V_glider_low_bar**2)*self.S_w*self.b_w)
        C_n_p_highbar = self.N/(((self.p*self.b_w)/2*self.V_glider_high_bar)*0.5*self.rho_high_bar*(self.V_glider_high_bar**2)*self.S_w*self.b_w)
        #value should be between -1 & -0.1

        return C_n_p_lowbar, C_n_p_highbar

    # Stability Derivatives with respect to yaw rate
    def stab_derivatives_yawrate(self):
        C_Y_r_v = 2*self.C_Y_v_alpha*((self.V_v/self.V)**2)*((self.S_v*self.l_v)/(self.S_w*self.b_w))
        C_l_r_lowbar = self.L/(((self.r*self.b_w)/(2*self.V_glider_low_bar))*0.5*self.rho_low_bar*(self.V_glider_low_bar**2)*self.S_w*self.b_w) + self.C_Y_r_v*(((self.z_v-self.z_cg)/self.b_w*np.cos(self.alpha_0))-((self.x_v-self.x_cg)/self.b_w*np.sin(self.alpha_0)))
        C_l_r_highbar = self.L/(((self.r*self.b_w)/(2*self.V_glider_high_bar))*0.5*self.rho_high_bar*(self.V_glider_high_bar**2)*self.S_w*self.b_w) + self.C_Y_r_v*(((self.z_v-self.z_cg)/self.b_w*np.cos(self.alpha_0))-((self.x_v-self.x_cg)/self.b_w*np.sin(self.alpha_0)))
        C_n_r = -C_Y_r_v*(self.l_v/self.b_w)

        return C_Y_r_v, C_l_r_lowbar, C_l_r_highbar, C_n_r

    #Other Derivatives
    def stab_derivatives_others(self):
        C_Z_alpha_dot = -self.C_N_h_alpha*(self.V_h_over_V**2)*self.de_da*((self.S_h*self.l_h)/(self.S_w*self.c_avg_w))
        C_Z_q = -2*self.C_N_h_alpha*(self.V_h_over_V**2)*((self.S_h*self.l_h)/(self.S_w*self.c_avg_w))
        C_m_alpha_dot = -self.C_N_h_alpha*(self.V_h_over_V**2)*self.de_da*((self.S_h*self.l_h**2)/(self.S_w*self.c_avg_w**2))
        C_m_q = 1.1*C_Z_q*0.5*(self.l_h/self.c_avg_w)

        return C_Z_alpha_dot, C_Z_q, C_m_alpha_dot, C_m_q


#Elevator trim:
'''''
delta_e = -(1/C_m_delta_e)*(C_m_0 + (W/(p_dynamic_low_bar*S_avg))*((x_cg-x_n_fixed)/c_avg))
delta_e = -(1/C_m_delta_e)*(C_m_0 + (W/(p_dynamic_high_bar*S_avg))*((x_cg-x_n_fixed)/c_avg))
'''''


if __name__ == "__main__":
    print("Hello World")
    glider = aerodynamicproperties()

    #print(glider.stability_cg(C_l_alpha_h, C_l_alpha_ah, downwash_vs_alpha(A_w), l_h, c_avg_w, V_h_over_V, x_ac, SM, x_cg))
    '''
    print(glider.stability_cg())
    print(glider.controllability_cg())
    print(glider.stability_cg()*glider.S_avg_w)
    print(glider.controllability_cg() * glider.S_avg_w)
    print(glider.S_avg_w)
    '''