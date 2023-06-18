import numpy as np
import matplotlib.pyplot as plt

# falcon_heavy = {
#     "wet_mass_b": 0,  # booster wet mass
#     "str_mass_b": 0,  # booster structural mass
#     "isp_b": 0,  # booster isp
#     "b_number": 2,  # number of boosters
#     "wet_mass_1": 433300,  # first stage wet mass
#     "str_mass_1": 22200,  # first stage structural mass
#     "isp_1": 288,  # first stage isp
#     "wet_mass_2": 111500,  # second stage wet mass
#     "str_mass_2": 4000,  # second stage structural mass
#     "isp_2": 348,  # second stage isp
#     "payload_mass": 3500,  # payload mass
#     "shell_mass": 1700,  # fairing mass
#     "deltaV": 8400  # transfer delta V
# }


class Launcher:
    def __init__(self, deltaV, shell, payload,
                 isp_2_vac, str_mass_2, wet_mass_2,
                 isp_1_vac, isp_1_sea, str_mass_1, wet_mass_1, tb_1,
                 booster_number=None, isp_boost_vac=None, isp_boost_sea=None,
                 str_mass_boost=None, wet_mass_boost=None, tb_boost=None):
        self.dv = deltaV
        self.m_shell = shell
        self.m_pl = payload

        self.isp_2_vac = isp_2_vac
        self.ve_2_vac = self.isp_2_vac * 9.80665
        self.m_str_2 = str_mass_2
        self.m_wet_2 = wet_mass_2
        self.m_prop_2 = self.m_wet_2 - self.m_str_2

        self.isp_1_vac = isp_1_vac
        self.isp_1_sea = isp_1_sea
        self.ve_1_vac = self.isp_1_vac * 9.80665
        self.ve_1_sea = self.isp_1_sea * 9.80665

        self.m_str_1 = str_mass_1
        self.m_wet_1 = wet_mass_1
        self.m_prop_1 = self.m_wet_1 - self.m_str_1
        self.tb_1 = tb_1

        if booster_number:
            self.nr_b = booster_number
            self.isp_b_sea = isp_boost_sea
            self.isp_b_vac = isp_boost_vac
            self.ve_b_sea = self.isp_b_sea * 9.80665
            self.ve_b_vac = self.isp_b_vac * 9.80665
            self.m_str_b = str_mass_boost
            self.m_wet_b = wet_mass_boost
            self.m_prop_b = self.m_wet_b - self.m_str_b
            self.tb_b = tb_boost

    def dv_for_payload(self, dv_goal, m_pl_correct, name):
        """
        For a given launcher, using payload to LEO, calculate payload for input delta V.
        dv_goal: Delta V beyond LEO
        NOTE: The calculator only works downwards, i.e., the final destination has to be closer by than the dV being
        calculated for. It answers the question: "if I decrease my payload, how much further can I go?"
        """
        # how much fuel does stage 2 use for a transfer?

        m_prop_trans = (np.exp(self.dv / self.ve_2_vac) - 1) * (self.m_pl + self.m_str_2)
        m_prop_leo = self.m_prop_2 - m_prop_trans  # how much fuel did the circularisation burn take?
        dv_leo = self.ve_2_vac*np.log((self.m_prop_2 + self.m_pl + self.m_shell +
                                       self.m_str_2)/(m_prop_trans + self.m_pl + self.m_shell + self.m_str_2))
        # print(f"Delta V of the second stage to LEO: {round(dv_leo)} [m/s]")
        m_stage_2 = self.m_pl + self.m_wet_2 + self.m_shell

        if self.nr_b:  # if we have a booster
            tb_core = self.tb_1 - self.tb_b
            m_prop_core = tb_core/self.tb_1 * self.m_prop_1
            dv_core = self.ve_1_vac * np.log((m_stage_2 + self.m_str_1 + m_prop_core)/(m_stage_2 + self.m_str_1))
            m_dot_1 = self.m_prop_1/self.tb_1
            m_dot_b = self.nr_b * self.m_prop_b/self.tb_b

            # We average the isp by mass flow
            ve_1b_sea = self.ve_1_sea * m_dot_1/(m_dot_1 + m_dot_b) + self.ve_b_sea * m_dot_b/(m_dot_1 + m_dot_b)
            ve_1b_vac = self.ve_1_vac * m_dot_1 / (m_dot_1 + m_dot_b) + self.ve_b_vac * m_dot_b / (m_dot_1 + m_dot_b)
            ve_avg = 0.2 * ve_1b_sea + 0.8 * ve_1b_vac

            dv_1b = ve_avg * np.log((m_stage_2 + self.m_str_1 + self.nr_b * self.m_str_b + self.m_prop_1 +
                                     self.nr_b * self.m_prop_b)/(m_stage_2 + self.m_str_1 + self.nr_b * self.m_str_b +
                                                                 m_prop_core))

        delta_v = dv_1b + dv_core + dv_leo
        print(f"Delta V basic: {delta_v}")
        
        


        m_pl_new = [self.m_pl]
        dv_new = [delta_v]
        x = True
        while x and m_pl_new[-1] > 0:  # stage 2 to LEO
            m_pl_new.append(m_pl_new[-1] - 10)
            m_1b_wet = (m_pl_new[-1] + self.m_shell + self.m_str_2 + self.m_prop_2 +
                       self.m_str_1 + self.m_prop_1 + self.nr_b * (self.m_prop_b + self.m_str_b))
            m_1b_dry = (m_pl_new[-1] + self.m_shell + self.m_str_2 + self.m_prop_2 +
                        self.m_str_1 + self.nr_b * self.m_str_b + (self.m_prop_1 - self.tb_b * m_dot_1))
            dv_1b_new = ve_avg * np.log(m_1b_wet/m_1b_dry)

            m_core_wet = (m_pl_new[-1] + self.m_shell + self.m_str_2 + self.m_prop_2 + self.m_str_1 +
                         (self.m_prop_1 - self.tb_b * m_dot_1))
            m_core_dry = (m_pl_new[-1] + self.m_shell + self.m_str_2 + self.m_prop_2 + self.m_str_1)
            dv_core_new = self.ve_1_vac * np.log(m_core_wet/m_core_dry)

            m_2_wet = (m_pl_new[-1] + self.m_shell + self.m_str_2 + self.m_prop_2)
            m_2_dry = (m_pl_new[-1] + self.m_shell + self.m_str_2)
            dv_2_new = self.ve_2_vac * np.log(m_2_wet/m_2_dry)

            dv_new.append(dv_1b_new + dv_core_new + dv_2_new)

            if m_pl_new[-1] == self.m_pl-10:
                print(dv_new, delta_v)
            if dv_new[-1] > delta_v + dv_goal:
                x = False
                print(f"Payload of {round(m_pl_new[-1])} [kg] for {round(dv_new[-1]-delta_v)} [m/s] with {name}")

        v_inf = []
        dv = []
        v_leo = 7800
        mu_earth = 3.986e14
        for i in dv_new:
            dv.append(i-delta_v)
            # print(2*mu_earth/200000)
            # print((v_leo + dv[-1])**2 - 2*mu_earth/66000000)
            # print(v_leo + dv[-1])
            # print(np.sqrt((v_leo + dv[-1])**2 - 2*mu_earth/200000))
            v_inf.append((v_leo + dv[-1])**2 - 2*mu_earth/6600000)
            # print(v_inf[-1])

        fig, ax1 = plt.subplots()
        ax1.plot(m_pl_new, dv, color='lightblue', label="Delta V vs payload mass")
        ax1.set_xlabel("Payload mass [kg]")
        ax1.set_ylabel("Delta V budget from LEO [m/s]")
        # ax1.legend()

        ax2 = ax1.twinx()
        ax2.plot(m_pl_new, v_inf, color='blue', label="C3 vs payload mass")
        ax2.set_ylabel("C3 m2/s2")
        plt.title("Payload mass vs Delta V, C3")
        # plt.legend()
        plt.figlegend(loc='right', bbox_to_anchor=(0.8, 0.8))
        plt.show()


# falcon_heavy_pluto = Launcher(8200, 1700, 3500,  # general data
#                               348, 4000, 111500,  # second stage data
#                               311, 282, 22200, 433100, 187,  # first stage data
#                               2, 311, 282, 22200, 433100, 154)  # booster data

falcon_heavy_mars = Launcher(3600, 1700, 16800,  # general data
                             348, 4000, 111500,  # second stage data
                             311, 282, 22200, 433100, 187,  # first stage data
                             2, 311, 282, 22200, 433100, 154)  # booster data

falcon_heavy_leo = Launcher(0, 1700, 63800,  # general data
                            348, 4000, 111500,  # second stage data
                            311, 282, 22200, 433100, 187,  # first stage data
                            2, 311, 282, 22200, 433100, 154)  # booster data

falcon_heavy_gto = Launcher(2300, 1700, 26700,  # general data
                            348, 4000, 111500,  # second stage data
                            311, 282, 22200, 433100, 187,  # first stage data
                            2, 311, 282, 22200, 433100, 154)  # booster data

# vulcan_centaur_LEO = Launcher(0, 1700, 27200,
#                               454, 1, 1,
#                               340, 315, 1, 1, 1,
#                               6, 280, 280,
#                               5182, 47852, 87)



# falcon_heavy_pluto.dv_for_payload(5000)
falcon_heavy_mars.dv_for_payload(7000, 16800, "FH_Mars")
# falcon_heavy_gto.dv_for_payload(7000, 16800, "FH_GTO")
falcon_heavy_leo.dv_for_payload(7900, 16800, "FH_LEO")

