#Library imports
import numpy as np



class Reentry:
    def __init__(self):
        #individual components of State Space Matrix:
    def A_matrix_entries(self):
        self.a_vv = (-1/(self.m*self.V_0))*(self.M_0*dCD_dCM*self.q_bar_0*self.S_ref +2*self.D_0)
        self.a_vgamma = -self.g_o*np.cos(self.gamma_0)
        self.a_V_R = 2*(self.g_0/self.R_0)*np.sin(self.gamma_0)
        self.a_Vp = 0
        self.a_Vq = 0
        self.a_Vr = 0
        self.a_valpha = (-1/self.m)*(dCD_dalpha)*self.q_bar_0*self.S_ref
        self.a_Vbeta = 0
        self.a_Vtheta = 0
        self.a_gammaV = (1/self.V_0)*(-self.gamma_prime_0 + 2*(self.V_0/self.R_0)*np.cos(self.gamma_0)) + (self.theta_0/(self.m*(self.V_0**2)))*(self.M_0*self.dCL_dCM*self.q_bar_0*self.S_ref + 2*self.L_0)
        self.a_gammagamma = -((self.V_0/self.R_0)-(self.g_0/self.V_0))*np.sin(self.gamma_0)
        self.a_gammaR = (((2*self.g_0)/self.V_0) - (self.V_0/self.R_0))*(np.cos(self.gamma_0)/self.R_0)
        self.a_gammaAlpha = (np.cos(self.theta_0)/(self.m*self.V_0))*self.dCL_dAlpha*self.q_bar_0*self.S_ref
        self.a_gammabeta = -(np.sin(self.theta_0)/(self.m*self.V_0))*self.dCS_dBeta*self.q_bar_0*self.S_ref
        self.a_gammatheta = -(self.L_0/(self.m*self.V_0))*np.sin(self.theta_0)
        self.a_gammap = self.a_gammaq = self.a_gammar = 0
        self.a_RV = np.sin(self.gamma_0)
        self.a_Rgamma = self.V_0*np.cos(self.gamma_0)
        self.a_RR =self.a_Rp =self.a_Rq = self.a_Rr =self.a_Ralpha = self.a_Rbeta = self.a_Rtheta = 0
        self.a_pbeta = (1/self.Ixx)*self.dCl_dBeta*self.q_bar_0*self.S_ref*self.b_ref
        self.a_pV = self.a_pgamma = self.a_pR = self.a_pp = self.a_pq = self.a_pr = self.a_palpha = self.a_ptheta = 0
        self.a_qV = (self.M_0/(self.Iyy*self.V_0))*self.dCm_dM*self.q_bar_0*self.S_ref*self.c_ref
        self.a_qalpha = (1/self.Iyy)*self.dCm_dAlpha*self.q_bar_0*self.S_ref*self.c_ref
        self.a_qgamma= self.a_qR =self.a_qp = self.a_qq = self.a_qr = self.a_qBeta = self.a_qtheta = 0
        self.a_rbeta = (1/self.Izz)*self.dCn_dBeta*self.q_bar_0*self.S_ref*self.b_ref
        self.a_rV = self.a_rgamma= self.a_rR = self.a_rp = self.a_rq =self.a_rr = self.a_ralpha = self.a_rtheta = 0
        self.a_alphaV = -(self.g_0/(self.V_0**2))*np.cos(self.gamma_0)*np.cos(self.theta_0) - (1/(self.m*self.V_0))*(self.M_0*self.dCL_dM + self.C_L)*self.q_bar_0*self.S_ref
        self.a_alphagamma = -(self.g_0/self.V_0)*np.sin(self.gamma_0)*np.cos(self.theta_0)
        self.a_alphaR -(2*self.g_0/(self.R_0*self.V_0))*np.cos(self.gamma_0)*np.cos(self.theta_0)
        self.a_alphaq = 1
        self.a_alphaalpha = (-1/(self.m*self.V_0))*self.dCL_dalpha*self.q_bar_0*self.S_ref
        self.a_alphatheta = (-self.g_0/self.V_0)*np.cos(self.gamma_0)*np.sin(self.theta_0)
        self.a_alphap = self.a_alphar = self.a_alphabeta = 0
        self.a_betaV = (self.g_0/(self.V_0**2))*np.cos(self.gamma_0)*np.sin(self.theta_0)
        self.a_betagamma = (self.g_0/self.V_0)*np.sin(self.gamma_0)*np.sin(self.theta_0)
        self.a_betaR = (2*self.g_0/(self.R_0*self.V_0))*np.cos(self.gamma_0)*np.sin(self.theta_0)
        self.a_betap = np.sin(self.alpha_0)
        self.betar = -np.cos(alpha_0)
        self.betabeta = -(1/(self.m*self.V_0))*self.dCS_dBeta*self.q_bar_0*self.S_ref
        self.betatheta = -(self.g_0/self.V_0)*np.cos(self.gamma_0)*np.cos(self.theta_0)
        self.a_betaq = self.a_betaalpha = 0
        self.a_thetaV = ((np.tan(self.gamma_0)*np.sin(self.theta_0))/(self.m*self.V_0**2))*(self.M_0*self.dCL_dM + self.C_L)*self.q_bar_0*self.S_ref
        self.a_thetagamma = (self.L_0/(self.m*self.V_0))*np.sin(self.theta_0)
        self.a_thetap = -np.cos(self.alpha_0)
        self.a_thetar = -np.sin(self.alpha_0)
        self.a_thetaalpha = ((np.tan(self.gamma_0)*np.sin(self.theta_0))/(self.m*self.V_0))*self.dCL_dAlpha*self.q_bar_0*self.S_ref
        self.a_thetabeta = ((np.tan(self.gamma_0)*np.cos(self.theta_0))/(self.m*self.V_0))*self.dCS_dBeta*self.q_bar_0*self.S_ref - (self.L_0/(self.m*self.V_0)) + (self.g_0/self.V_0)*np.cos(self.gamma_0)*np.cos(theta_0)
        self.a_thetatheta = np.tan(self.gamma_0)*np.cos(self.theta_0)*(self.L_0/(self.m*self.V_0))
        self.a_thetaR = self.a_thetaq = 0

    def B_matrix_entries(self):
        self.b_Ve = self.b_Va =self.b_Vr = self.b_Vx = self.b_Vy =self.b_Vz = 0
        self.b_gammae = self.b_gammaa = self.b_gammar = self.b_gammax = self.b_gammay = self.b_gammaz = 0
        self.b_Re = self.b_Ra = self.b_Rr = self.b_Rx = self.b_Ry = self.b_Rz = 0
        self.b_pa = (1/self.Ixx)*self.dCl_dDeltaA*self.q_bar_0*self.S_ref*self.b_ref
        self.b_px = 1/self.Ixx
        self.b_pe = self.b_pr = self.b_py = self.b_pz = 0
        self.b_qe = (1/self.Iyy)*self.dCm_dDeltaE*self.q_bar_0*self.S_ref*self.c_ref
        self.b_qy = 1/self.Iyy
        self.b_qa = self.b_qr = self.b_qx = self.b_qz = 0
        self.b_ra = (1/self.Izz)*self.dCn_dDeltaA*self.q_bar_0*self.S_ref*self.b_ref
        self.b_rr = (1/self.Izz)*self.dCn_dDeltar*self.q_bar_0*self.S_ref*self.b_ref
        self.b_rz = 1/self.Izz
        self.b_re = self.b_rx = self.b_ry = 0
        self.b_alphae = self.b_alphaa = self.b_alphar = self.b_alphax = self.b_alphay = self.b_alphaz = 0
        self.b_etae = self.b_betaa = self.b_betar = self.b_betax = self.b_betay = self.b_betaz = 0
        self.b_thetae = self.b_thetaa = self.b_thetar = self.b_thetax = self.b_thetay = self.b_thetaz = 0





if __name__ == "__main__":
    print("Hello World")