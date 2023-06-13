#Library imports
import numpy as np
from numpy.linalg import eig
import matplotlib.pyplot as plt

class Reentry():
    def __init__(self):
        self.m = 4976 #[kg]
        self.Ixx = 2 #[kg/m^2]
        self.Iyy = 2 #[kg/m^2]
        self.Izz = 2 #[kg/m^2]
        self.Ixz = 1 #[kg/m^2]
        self.S_ref = 12 #[m^2] needs to be
        self.b_ref = 3.9  #[m] assumed to be the diameter of the heatshield
        self.c_ref = 0.75*self.b_ref #[m] taken to be the height of the capsule
        self.V_0 =  11000 #[m/s]
        self.rho_0 = 1.225 #[kg/m^3]
        self.a_0 = 340.2 #[m/s]
        self.M_0 = self.V_0/self.a_0
        self.q_bar_0 = 0.5*self.rho_0*self.V_0**2
        self.C_D0 = 0.72
        self.D_0 = self.q_bar_0 * self.S_ref * self.C_D0
        self.g_0 = 9.80665 #[m/s^2]
        self.gamma_0 = np.deg2rad(-9.536)  # [rad]
        self.mu_0 = np.deg2rad(1) #rad
        self.R_0 = 220000  # [m]
        self.gamma_prime_0 = 0
        self.theta_0 = np.deg2rad(1) #[rad]
        self.dCD_dCM = 1
        self.dCD_dM = 1
        self.dCD_dalpha = 1
        self.dCL_dCM = 1
        self.C_L0 = 0.16
        self.L_0 = self.q_bar_0*self.S_ref*self.C_L0
        self.dCL_dAlpha = 1
        self.dCS_dBeta = 1
        self.C_L = 1
        self.dCl_dBeta = 1
        self.dCm_dAlpha = 1
        self.dCm_dM = 0.5
        self.dCn_dBeta = 1
        self.dCL_dM = 1
        self.alpha_0 = np.deg2rad(-24.5) #[rad]
        self.dCl_dDeltaA = 0
        self.dCn_dDeltar = 0
        self.dCn_dDeltaA = 0
        self.dCm_dDeltaE = 0
        self.p_0 = ((self.g_0 / self.V_0) * np.cos(self.gamma_0) * np.sin(self.mu_0)) * np.sin(self.alpha_0) + (((self.L_0) / (self.m * self.V_0)) * np.tan(self.gamma_0) * np.sin(self.mu_0)) * np.cos(self.alpha_0)
        self.q_0 = (self.L_0 / (self.m * self.V_0)) - (self.g_0 / self.V_0) * np.cos(self.gamma_0) * np.cos(self.mu_0)
        self.r_0 = -((self.g_0 / self.V_0) * np.cos(self.gamma_0) * np.sin(self.mu_0)) * np.cos(self.alpha_0) + ((self.L_0 / (self.m * self.V_0)) * np.tan(self.gamma_0) * np.sin(self.mu_0)) * np.sin(self.alpha_0)
        self.A_matrix_entries()
        self.B_matrix_entries()
        self.x_matrix_entries()
        self.u_matrix_entries()
        self.state_space_matrices()
        self.eigen_estimation()


        #Individual components of State Space Matrix:
    def A_matrix_entries(self):
        self.a_VV = (-1/(self.m*self.V_0))*(self.M_0*self.dCD_dM*self.q_bar_0*self.S_ref +2*self.D_0)
        self.a_Vgamma = -self.g_0*np.cos(self.gamma_0)
        self.a_VR = 2*(self.g_0/self.R_0)*np.sin(self.gamma_0)
        self.a_Vp = 0
        self.a_Vq = 0
        self.a_Vr = 0
        self.a_Valpha = (-1/self.m)*(self.dCD_dalpha)*self.q_bar_0*self.S_ref
        self.a_Vbeta = 0
        self.a_Vtheta = 0
        #Creating array
        self.a_V = np.array([self.a_VV, self.a_Vgamma, self.a_VR, self.a_Vp, self.a_Vq, self.a_Vr, self.a_Valpha, self.a_Vbeta, self.a_Vtheta])

        self.a_gammaV = (1/self.V_0)*(-self.gamma_prime_0 + ((2*self.V_0)/self.R_0)*np.cos(self.gamma_0)) + (np.cos(self.theta_0)/(self.m*self.V_0**2))*(self.M_0*self.dCL_dM*self.q_bar_0*self.S_ref + 2*self.L_0)
        self.a_gammagamma = -((self.V_0/self.R_0)-(self.g_0/self.V_0))*np.sin(self.gamma_0)
        self.a_gammaR = (((2*self.g_0)/self.V_0) - (self.V_0/self.R_0))*(np.cos(self.gamma_0)/self.R_0)
        self.a_gammap = 0
        self.a_gammaq = 0
        self.a_gammar = 0
        self.a_gammaAlpha = (np.cos(self.theta_0)/(self.m*self.V_0))*self.dCL_dAlpha*self.q_bar_0*self.S_ref
        self.a_gammabeta = -(np.sin(self.theta_0)/(self.m*self.V_0))*self.dCS_dBeta*self.q_bar_0*self.S_ref
        self.a_gammatheta = -(self.L_0/(self.m*self.V_0))*np.sin(self.theta_0)
        # Creating array
        self.a_gamma = np.array([self.a_gammaV, self.a_gammagamma, self.a_gammaR, self.a_gammap, self.a_gammaq, self.a_gammar, self.a_gammaAlpha, self.a_gammabeta, self.a_gammatheta])


        self.a_RV = np.sin(self.gamma_0)
        self.a_Rgamma = self.V_0*np.cos(self.gamma_0)
        self.a_RR = 0
        self.a_Rp = 0
        self.a_Rq = 0
        self.a_Rr = 0
        self.a_Ralpha = 0
        self.a_Rbeta = 0
        self.a_Rtheta = 0
        # Creating array
        self.a_R = np.array([self.a_RV, self.a_Rgamma, self.a_RR, self.a_Rp, self.a_Rq, self.a_Rr, self.a_Ralpha, self.a_Rbeta, self.a_Rtheta])


        self.a_pV = 0
        self.a_pgamma = 0
        self.a_pR = 0
        self.a_pp = ((self.Ixx-self.Iyy+self.Izz)/(self.Ixx*self.Izz - self.Ixz*self.Ixz))*self.q_0
        self.a_pq = ((self.Ixx-self.Iyy+self.Izz)/(self.Ixx*self.Izz - self.Ixz*self.Ixz))*self.p_0 + (((self.Iyy-self.Izz)*self.Izz - self.Ixz**2)/(self.Ixx*self.Izz - self.Ixz*self.Ixz))*self.r_0
        self.a_pr = (((self.Iyy-self.Izz)*self.Izz - self.Ixz**2)/(self.Ixx*self.Izz - self.Ixz*self.Ixz))*self.q_0
        self.a_palpha = 0
        self.a_pbeta = (1/self.Ixx)*self.dCl_dBeta*self.q_bar_0*self.S_ref*self.b_ref
        self.a_ptheta = 0
        # Creating array
        self.a_p = np.array([self.a_pV, self.a_pgamma, self.a_pR, self.a_pp, self.a_pq, self.a_pr, self.a_palpha, self.a_pbeta, self.a_ptheta])


        self.a_qV = (self.M_0/(self.Iyy*self.V_0))*self.dCm_dM*self.q_bar_0*self.S_ref*self.c_ref
        self.a_qgamma = 0
        self.a_qR = 0
        self.a_qp = -2*(self.Ixz/self.Iyy)*self.p_0 + ((self.Izz-self.Ixx)/self.Iyy)*self.r_0
        self.a_qq = 0
        self.a_qr = ((self.Izz-self.Ixx)/self.Iyy)*self.p_0 +2*(self.Ixz/self.Iyy)*self.r_0
        self.a_qalpha = (1/self.Iyy)*self.dCm_dAlpha*self.q_bar_0*self.S_ref*self.c_ref
        self.a_qBeta = 0
        self.a_qtheta = 0
        # Creating array
        self.a_q = np.array([self.a_qV, self.a_qgamma, self.a_qR, self.a_qp, self.a_qq, self.a_qr, self.a_qalpha, self.a_qBeta, self.a_qtheta])


        self.a_rV = 0
        self.a_rgamma = 0
        self.a_rR = 0
        self.a_rp = (((self.Ixx-self.Iyy)*self.Ixx - self.Ixz**2)/(self.Ixx*self.Izz - self.Ixz*self.Ixz))*self.q_0
        self.a_rq = (((self.Ixx-self.Iyy)*self.Ixx - self.Ixz**2)/(self.Ixx*self.Izz - self.Ixz*self.Ixz))*self.p_0 + (((-self.Ixx + self.Iyy - self.Izz)*self.Ixz)/(self.Ixx*self.Izz - self.Ixz*self.Ixz))*self.r_0
        print(self.a_rq)
        self.a_rr = (((-self.Ixx + self.Iyy - self.Izz)*self.Ixz)/(self.Ixx*self.Izz - self.Ixz*self.Ixz))*self.q_0
        self.a_ralpha = 0
        self.a_rbeta = (1/self.Izz)*self.dCn_dBeta*self.q_bar_0*self.S_ref*self.b_ref
        self.a_rtheta = 0
        # Creating array
        self.a_r = np.array([self.a_rV, self.a_rgamma, self.a_rR, self.a_rp, self.a_rq, self.a_rr, self.a_ralpha, self.a_rbeta, self.a_rtheta])


        self.a_alphaV = -(self.g_0/(self.V_0**2))*np.cos(self.gamma_0)*np.cos(self.theta_0) - (1/(self.m*self.V_0**2))*(self.M_0*self.dCL_dM + self.C_L)*self.q_bar_0*self.S_ref
        self.a_alphagamma = -(self.g_0/self.V_0)*np.sin(self.gamma_0)*np.cos(self.theta_0)
        self.a_alphaR = -(2*self.g_0/(self.R_0*self.V_0))*np.cos(self.gamma_0)*np.cos(self.theta_0)
        self.a_alphap = 0
        self.a_alphaq = 1
        self.a_alphar = 0
        self.a_alphaalpha = (-1/(self.m*self.V_0))*self.dCL_dAlpha*self.q_bar_0*self.S_ref
        self.a_alphabeta = 0
        self.a_alphatheta = (-self.g_0/self.V_0)*np.cos(self.gamma_0)*np.sin(self.theta_0)
        # Creating array
        self.a_alpha = np.array([self.a_alphaV, self.a_alphagamma, self.a_alphaR, self.a_alphap, self.a_alphaq, self.a_alphar, self.a_alphaalpha, self.a_alphabeta, self.a_alphatheta])


        self.a_betaV = (self.g_0/(self.V_0**2))*np.cos(self.gamma_0)*np.sin(self.theta_0)
        self.a_betagamma = (self.g_0/self.V_0)*np.sin(self.gamma_0)*np.sin(self.theta_0)
        self.a_betaR = (2*self.g_0/(self.R_0*self.V_0))*np.cos(self.gamma_0)*np.sin(self.theta_0)
        self.a_betap = np.sin(self.alpha_0)
        self.a_betaq = 0
        self.a_betar = -np.cos(self.alpha_0)
        self.a_betaalpha = 0
        self.a_betabeta = -(1/(self.m*self.V_0))*self.dCS_dBeta*self.q_bar_0*self.S_ref
        self.a_betatheta = -(self.g_0/self.V_0)*np.cos(self.gamma_0)*np.cos(self.theta_0)
        # Creating array
        self.a_beta = np.array([self.a_betaV, self.a_betagamma, self.a_betaR, self.a_betap, self.a_betaq, self.a_betar, self.a_betaalpha, self.a_betabeta, self.a_betatheta])


        self.a_thetaV = ((np.tan(self.gamma_0)*np.sin(self.theta_0))/(self.m*self.V_0**2))*(self.M_0*self.dCL_dM + self.C_L)*self.q_bar_0*self.S_ref
        self.a_thetagamma = (self.L_0/(self.m*self.V_0))*np.sin(self.theta_0)
        self.a_thetaR = 0
        self.a_thetap = -np.cos(self.alpha_0)
        self.a_thetaq = 0
        self.a_thetar = -np.sin(self.alpha_0)
        self.a_thetaalpha = ((np.tan(self.gamma_0)*np.sin(self.theta_0))/(self.m*self.V_0))*self.dCL_dAlpha*self.q_bar_0*self.S_ref
        self.a_thetabeta = ((np.tan(self.gamma_0)*np.cos(self.theta_0))/(self.m*self.V_0))*self.dCS_dBeta*self.q_bar_0*self.S_ref - (self.L_0/(self.m*self.V_0)) + (self.g_0/self.V_0)*np.cos(self.gamma_0)*np.cos(self.theta_0)
        self.a_thetatheta = np.tan(self.gamma_0)*np.cos(self.theta_0)*(self.L_0/(self.m*self.V_0))
        # Creating array
        self.a_theta = np.array([self.a_thetaV, self.a_thetagamma, self.a_thetaR, self.a_thetap, self.a_thetaq, self.a_thetar, self.a_thetaalpha, self.a_thetabeta, self.a_thetatheta])
        self.A = np.transpose(np.array([self.a_V, self.a_gamma, self.a_R, self.a_p, self.a_q, self.a_r, self.a_alpha, self.a_beta, self.a_theta]))
    def x_matrix_entries(self):
        self.Delta_V = 0.5
        self.Delta_gamma = 0.5
        self.Delta_R = 50
        self.Delta_p = 0
        self.Delta_q = 0
        self.Delta_r = 0
        self.Delta_alpha = 0
        self.Delta_beta = 0
        self.Delta_theta = 0.5
        self.x = np.array([self.Delta_V, self.Delta_gamma, self.Delta_R, self.Delta_p, self.Delta_q, self.Delta_r, self.Delta_alpha, self.Delta_beta, self.Delta_theta])

    def B_matrix_entries(self):
        self.b_Ve = 0
        self.b_Va = 0
        self.b_Vr = 0
        self.b_Vx = 0
        self.b_Vy = 0
        self.b_Vz = 0
        # Creating array
        self.b_V = np.array([self.b_Ve, self.b_Va, self.b_Vr, self.b_Vx, self.b_Vy, self.b_Vz])

        self.b_gammae = 0
        self.b_gammaa = 0
        self.b_gammar = 0
        self.b_gammax = 0
        self.b_gammay = 0
        self.b_gammaz = 0
        # Creating array
        self.b_gamma = np.array([self.b_gammae, self.b_gammaa, self.b_gammar, self.b_gammax, self.b_gammay, self.b_gammaz])


        self.b_Re = 0
        self.b_Ra = 0
        self.b_Rr = 0
        self.b_Rx = 0
        self.b_Ry = 0
        self.b_Rz = 0
        # Creating array
        self.b_R = np.array([self.b_Re, self.b_Ra, self.b_Rr, self.b_Rx, self.b_Ry, self.b_Rz])


        self.b_pa = (1/self.Ixx)*self.dCl_dDeltaA*self.q_bar_0*self.S_ref*self.b_ref
        self.b_pe = 0
        self.b_pr = 0
        self.b_px = self.Izz/(self.Ixx*self.Izz - self.Ixz*self.Ixz)
        self.b_py = 0
        self.b_pz = self.Ixz/(self.Ixx*self.Izz - self.Ixz*self.Ixz)
        # Creating array
        self.b_p = np.array([self.b_pe, self.b_pa, self.b_pr, self.b_px, self.b_py, self.b_pz])


        self.b_qe = (1/self.Iyy)*self.dCm_dDeltaE*self.q_bar_0*self.S_ref*self.c_ref
        self.b_qa = 0
        self.b_qr = 0
        self.b_qx = 0
        self.b_qy = 1/self.Iyy
        self.b_qz = 0
        # Creating array
        self.b_q = np.array([self.b_qe, self.b_qa, self.b_qr, self.b_qx, self.b_qy, self.b_qz])


        self.b_re = 0
        self.b_ra = (1/self.Izz)*self.dCn_dDeltaA*self.q_bar_0*self.S_ref*self.b_ref
        self.b_rr = (1/self.Izz)*self.dCn_dDeltar*self.q_bar_0*self.S_ref*self.b_ref
        self.b_rx = self.Ixz/(self.Ixx*self.Izz - self.Ixz*self.Ixz)
        self.b_ry = 0
        self.b_rz = self.Ixx/(self.Ixx*self.Izz - self.Ixz*self.Ixz)
        # Creating array
        self.b_r = np.array([self.b_re, self.b_ra, self.b_rr, self.b_rx, self.b_ry, self.b_rz])


        self.b_alphae = 0
        self.b_alphaa = 0
        self.b_alphar = 0
        self.b_alphax = 0
        self.b_alphay = 0
        self.b_alphaz = 0
        # Creating array
        self.b_alpha = np.array([self.b_alphae, self.b_alphaa, self.b_alphar, self.b_alphax, self.b_alphay, self.b_alphaz])


        self.b_betae = 0
        self.b_betaa = 0
        self.b_betar = 0
        self.b_betax = 0
        self.b_betay = 0
        self.b_betaz = 0
        # Creating array
        self.b_beta = np.array([self.b_betae, self.b_betaa, self.b_betar, self.b_betax, self.b_betay, self.b_betaz])


        self.b_thetae = 0
        self.b_thetaa = 0
        self.b_thetar = 0
        self.b_thetax = 0
        self.b_thetay = 0
        self.b_thetaz = 0
        # Creating array
        self.b_theta = np.array([self.b_thetae, self.b_thetaa, self.b_thetar, self.b_thetax, self.b_thetay, self.b_thetaz])
        self.B = np.array([self.b_V, self.b_gamma, self.b_R, self.b_p, self.b_q, self.b_r, self.b_alpha, self.b_beta, self.b_theta])

    def u_matrix_entries(self):
        self.Delta_deltae = 0
        self.Delta_deltaa = 0
        self.Delta_deltar = 0
        self.Delta_M_Tx = 1000
        self.Delta_M_Ty = 1000
        self.Delta_M_Tz = 1000
        self.u = np.array([self.Delta_deltae, self.Delta_deltaa, self.Delta_deltar, self.Delta_M_Tx, self.Delta_M_Ty, self.Delta_M_Tz])


    def state_space_matrices(self):
        self.x_prime = np.matmul(self.A, self.x) + np.matmul(self.B, self.u)

    def eigen_estimation(self):
        self.eigen_value, self.eigen_vector = eig(self.A)
        self.Period = []
        self.T_half = []
        self.damping_ratio = []
        self.nat_freq = []

        for i in range(len(self.eigen_value)):
            if self.eigen_value[i].imag == 0:
                self.T_half.append(np.log(0.5) / (self.eigen_value[i].real))
                #print(self.T_half)
                self.damping_ratio.append(-1.00)
                #print(self.damping_ratio)
                self.nat_freq.append(self.eigen_value[i].real)
                #print(self.nat_freq)
                self.Period.append(0.00)
            else:
                self.damping_ratio.append(-self.eigen_value[i].real / (np.sqrt(self.eigen_value[i].real ** 2 + self.eigen_value[i].imag ** 2)))
                #print(self.damping_ratio)
                self.nat_freq.append(np.sqrt(self.eigen_value[i].real ** 2 + self.eigen_value[i].imag ** 2))
                #print(self.nat_freq)
                self.Period.append(2*np.pi/(self.eigen_value[i].imag))
                #print(self.Period)
                self.T_half.append(np.log(0.5)/(self.eigen_value[i].real))
                #print(self.T_half)
    def eigenvalue_plot(self):
        x = self.eigen_value.real
        y= self.eigen_value.imag
        plt.plot(x,y,"ro")
        plt.grid(True)
        plt.title("Root Locus Plot")
        plt.xlabel("Real")
        plt.ylabel("Imaginary")
        plt.show()



if __name__ == "__main__":
    print("Hello World")
    test = Reentry()
    #test.eigenvalue_plot()

    #test.eigen_estimation()
    #print(test.Period)
    #print(test.eigen_value)
    # print(test.A.shape)
    # print(test.B.shape)
    #print(test.a_VV)
    # print(test.x.shape)
    # print(test.u.shape)
    #print(test.Period)
    #print(test.T_half)
    #print(test.damping_ratio)
    #print(test.nat_freq)
    #print(test.eigen_vector)

