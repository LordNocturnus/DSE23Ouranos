import numpy as np

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

    def dv_for_payload(self, dv_goal):
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
        print(f"Delta V of the second stage to LEO: {round(dv_leo)} [m/s]")
        m_stage_2 = self.m_pl + self.m_wet_2 + self.m_shell

        m_pl_attempt = self.m_pl
        x = True
        while x and m_pl_attempt > 0:  # stage 2 to LEO

            m_pl_attempt -= 10
            # Prop to LEO for lower payload
            m_prop_leo_attempt = (1 - 1 / np.exp(dv_leo / self.ve_2_vac)) * \
                                 (self.m_prop_2 + self.m_str_2 + self.m_shell + m_pl_attempt)
            m_prop_trans_attempt = self.m_prop_2 - m_prop_leo_attempt  # Prop available beyond LEO
            dv_trans_attempt = self.ve_2_vac * np.log((m_pl_attempt + self.m_str_2 + m_prop_trans_attempt) /
                                                      (m_pl_attempt + self.m_str_2))  # dV from this propellant
            if dv_trans_attempt >= dv_goal:
                x = False
                print(f"Delta V for transfer of {round(m_pl_attempt)} [kg] payload: {round(dv_trans_attempt)} [m/s]\n")
                      # f"Fuel for LEO manoeuvre: {round(m_prop_leo_attempt)} of {round(self.m_prop_2)} [kg]\n")
        """
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
                                     self.nr_b * self.m_prop_b)/(m_stage_2 + self.m_str_1 + self.nr_b * self.m_str_b))
            
        """


falcon_heavy_pluto = Launcher(8200, 1700, 3500,  # general data
                              348, 4000, 111500,  # second stage data
                              311, 282, 22200, 433100, 187,  # first stage data
                              2, 311, 282, 22200, 433100, 154)  # booster data

falcon_heavy_mars = Launcher(2900, 1700, 16800,  # general data
                             348, 4000, 111500,  # second stage data
                             311, 282, 22200, 433100, 187,  # first stage data
                             2, 311, 282, 22200, 433100, 154)  # booster data

falcon_heavy_leo = Launcher(0, 1700, 63800,  # general data
                            348, 4000, 111500,  # second stage data
                            311, 282, 22200, 433100, 187,  # first stage data
                            2, 311, 282, 22200, 433100, 154)  # booster data

falcon_heavy_gto = Launcher(2000, 1700, 26700,  # general data
                            348, 4000, 111500,  # second stage data
                            311, 282, 22200, 433100, 187,  # first stage data
                            2, 311, 282, 22200, 433100, 154)  # booster data



# falcon_heavy_pluto.dv_for_payload(5000)
falcon_heavy_mars.dv_for_payload(4000)
falcon_heavy_leo.dv_for_payload(4000)
falcon_heavy_gto.dv_for_payload(4000)

