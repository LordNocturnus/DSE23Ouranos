import unittest
from atmospheric_vehicle.aeroshell.ADCS import Reentry

class MyTestCase(unittest.TestCase):
    def test_A_matrix(self):
        a_VV_true = -548.7583462347663
        a_VV_false = self.ADCS.a_VV
        self.assertAlmostEqual(a_VV_true, a_VV_false)
        a_Vgamma_true = -9.671138811
        a_Vgamma_false = self.ADCS.a_Vgamma
        self.assertAlmostEqual(a_Vgamma_true, a_Vgamma_false)
        a_VR_true = -0.000014769
        a_VR_false = self.ADCS.a_VR
        self.assertAlmostEqual(a_VR_true, a_VR_false)
        a_Valpha_true = -178727.8939
        a_Valpha_false = self.ADCS.a_Valpha
        self.assertAlmostEqual(a_Valpha_true, a_Valpha_false, 4)

        a_gammaV_true = 0.048241747
        a_gammaV_false = self.ADCS.a_gammaV
        self.assertAlmostEqual(a_gammaV_true, a_gammaV_false,4)
        a_gammagamma_true = 0.008135669
        a_gammagamma_false = self.ADCS.a_gammagamma
        self.assertAlmostEqual(a_gammagamma_true, a_gammagamma_false)
        a_gammaR_true = -0.000000216
        a_gammaR_false = self.ADCS.a_gammaR
        self.assertAlmostEqual(a_gammaR_true, a_gammaR_false)
        a_gammaalpha_true = 16.24799035
        a_gammaalpha_false = self.ADCS.a_gammaAlpha
        self.assertAlmostEqual(a_gammaalpha_true, a_gammaalpha_false,2)
        a_gammabeta_true = -0.283566531
        a_gammabeta_false = self.ADCS.a_gammabeta
        self.assertAlmostEqual(a_gammabeta_true, a_gammabeta_false)
        a_gammatheta_true = -0.045370645
        a_gammatheta_false = self.ADCS.a_gammatheta
        self.assertAlmostEqual(a_gammatheta_true, a_gammatheta_false)

        a_RV_true = -0.165667275
        a_RV_false = self.ADCS.a_RV
        self.assertAlmostEqual(a_RV_true, a_RV_false)
        a_Rgamma_true = 10847.99875
        a_Rgamma_false = self.ADCS.a_Rgamma
        self.assertAlmostEqual(a_Rgamma_true, a_Rgamma_false,5)

        a_pbeta_true = 1734232500
        a_pbeta_false = self.ADCS.a_pbeta
        self.assertAlmostEqual(a_pbeta_true, a_pbeta_false,4)
        a_pp_true = 1.732532931
        a_pp_false = self.ADCS.a_pp
        self.assertAlmostEqual(a_pp_true, a_pp_false)
        a_pq_true = -0.005676816
        a_pq_false = self.ADCS.a_pq
        self.assertAlmostEqual(a_pq_true, a_pq_false)
        a_pr_true = -0.866266465
        a_pr_false = self.ADCS.a_pr
        self.assertAlmostEqual(a_pr_true, a_pr_false)

        a_qV_true = 1911631.944
        a_qV_false = self.ADCS.a_qV
        self.assertAlmostEqual(a_qV_true, a_qV_false,3)
        a_qalpha_true = 1300674375
        a_qalpha_false = self.ADCS.a_qalpha
        self.assertAlmostEqual(a_qalpha_true, a_qalpha_false)
        a_qp_true = 0.006941861
        a_qp_false = self.ADCS.a_qp
        self.assertAlmostEqual(a_qp_true, a_qp_false)
        a_qr_true = 0.003146726
        a_qr_false = self.ADCS.a_qr
        self.assertAlmostEqual(a_qr_true, a_qr_false)

        a_rp_true = -0.866266465
        a_rp_false = self.ADCS.a_rp
        self.assertAlmostEqual(a_rp_true, a_rp_false)
        a_rq_true = 2.161363353*(10**(-4))
        a_rq_false = self.ADCS.a_rq
        self.assertAlmostEqual(a_rq_true, a_rq_false)
        a_rr_true = -1.732532931
        a_rr_false = self.ADCS.a_rr
        self.assertAlmostEqual(a_rr_true, a_rr_false)
        a_rbeta_true = 1734232500
        a_rbeta_false = self.ADCS.a_rbeta
        self.assertAlmostEqual(a_rbeta_true, a_rbeta_false)

        a_alphaV_true = -4.92372827*(10**(-2))
        a_alphaV_false = self.ADCS.a_alphaV
        self.assertAlmostEqual(a_alphaV_true, a_alphaV_false)
        a_alphagamma_true = 0.000147672
        a_alphagamma_false = self.ADCS.a_alphagamma
        self.assertAlmostEqual(a_alphagamma_true, a_alphagamma_false)
        a_alphaR_true = -7.99145938*(10**(-9))
        a_alphaR_false = self.ADCS.a_alphaR
        self.assertAlmostEqual(a_alphaR_true, a_alphaR_false)
        a_alphaalpha_true = -1.624799035*(10**1)
        a_alphaalpha_false = self.ADCS.a_alphaalpha
        self.assertAlmostEqual(a_alphaalpha_true, a_alphaalpha_false)
        a_alphatheta_true = -1.53440587*(10**(-5))
        a_alphatheta_false = self.ADCS.a_alphatheta
        self.assertAlmostEqual(a_alphatheta_true, a_alphatheta_false)

        a_betaV_true = 1.394914421*(10**(-9))
        a_betaV_false = self.ADCS.a_betaV
        self.assertAlmostEqual(a_betaV_true, a_betaV_false)
        a_betagamma_true = -2.57762679 * (10 ** (-6))
        a_betagamma_false = self.ADCS.a_betagamma
        self.assertAlmostEqual(a_betagamma_true, a_betagamma_false)
        a_betaR_true = 1.39491442 * (10 ** (-10))
        a_betaR_false = self.ADCS.a_betaR
        self.assertAlmostEqual(a_betaR_true, a_betaR_false)
        a_betap_true = -4.146932427 * (10 ** (-1))
        a_betap_false = self.ADCS.a_betap
        self.assertAlmostEqual(a_betap_true, a_betap_false)
        a_betar_true = -9.09961271 * (10 ** (-1))
        a_betar_false = self.ADCS.a_betar
        self.assertAlmostEqual(a_betar_true, a_betar_false)
        a_betabeta_true = -1.624799035 * (10 ** (1))
        a_betabeta_false = self.ADCS.a_betabeta
        self.assertAlmostEqual(a_betabeta_true, a_betabeta_false)
        a_betatheta_true = -8.79060532 * (10 ** (-4))
        a_betatheta_false = self.ADCS.a_betatheta
        self.assertAlmostEqual(a_betatheta_true, a_betatheta_false)

        a_thetaV_true = -1.44353885 * (10 ** (-4))
        a_thetaV_false = self.ADCS.a_thetaV
        self.assertAlmostEqual(a_thetaV_true, a_thetaV_false)
        a_thetagamma_true = 4.537064503 * (10 ** (-2))
        a_thetagamma_false = self.ADCS.a_thetagamma
        self.assertAlmostEqual(a_thetagamma_true, a_thetagamma_false)
        a_thetap_true = -9.09961271 * (10 ** (-1))
        a_thetap_false = self.ADCS.a_thetap
        self.assertAlmostEqual(a_thetap_true, a_thetap_false)
        a_thetar_true = 4.146932427 * (10 ** (-1))
        a_thetar_false = self.ADCS.a_thetar
        self.assertAlmostEqual(a_thetar_true, a_thetar_false)
        a_thetaalpha_true = -4.7635942 * (10 ** (-2))
        a_thetaalpha_false = self.ADCS.a_thetaalpha
        self.assertAlmostEqual(a_thetaalpha_true, a_thetaalpha_false)
        a_thetabeta_true = -5.327860686
        a_thetabeta_false = self.ADCS.a_thetabeta
        self.assertAlmostEqual(a_thetabeta_true, a_thetabeta_false)
        a_thetatheta_true = -4.36649806 * (10 ** (-1))
        a_thetatheta_false = self.ADCS.a_thetatheta
        self.assertAlmostEqual(a_thetatheta_true, a_thetatheta_false)

    def test_B_matrix(self):
        b_px_true = 2/3
        b_px_false = self.ADCS.b_px
        self.assertAlmostEqual(b_px_true, b_px_false)
        b_pz_true = 1 / 3
        b_pz_false = self.ADCS.b_pz
        self.assertAlmostEqual(b_pz_true, b_pz_false)

        b_qy_true = 1 / 2
        b_qy_false = self.ADCS.b_qy
        self.assertAlmostEqual(b_qy_true, b_qy_false)


        b_rx_true = 1 / 3
        b_rx_false = self.ADCS.b_rx
        self.assertAlmostEqual(b_rx_true, b_rx_false)
        b_rz_true = 2 / 3
        b_rz_false = self.ADCS.b_rz
        self.assertAlmostEqual(b_rz_true, b_rz_false)



    def setUp(self) -> None:
        self.ADCS = Reentry()

if __name__ == '__main__':
    unittest.main()
