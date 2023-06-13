from unittest import TestCase
from Orbiter.comms.main import *


class Test(TestCase):
    def test_space_loss(self, wavelength=590 * 10 ** -9, d=15.6, AU=1375683.9):
        space_loss_true = -293.2000221
        space_loss_numerical = SpaceLoss(wavelength, d, AU)
        self.assertAlmostEqual(space_loss_numerical, space_loss_true)

    def test_gain(self, D=5.3, wavelength=590 * 10 ** -9, eta=0.76):
        gain_true = 147.8196105
        gain_numerical = gain(D, wavelength, eta)
        self.assertAlmostEqual(gain_numerical, gain_true)

    def test_halfpowerbeamwidth(self, D=5.6, f=5.084745763 * 10 ** 14):
        half_power_beam_width_true = 7.79245283 * 10 ** -15
        half_power_beam_width_numerical = halfpowerbeamwidth(D, f)
        self.assertAlmostEqual(half_power_beam_width_numerical, half_power_beam_width_true)

    def test_pointing_loss(self, alpha=7.79245283 * 10 ** -15, pointingAccuracy=2 * 10 ** -14):
        poiting_loss_true = -79.04836166
        poiting_loss_numerical = pointingLoss(alpha, pointingAccuracy)
        self.assertAlmostEqual(poiting_loss_numerical, poiting_loss_true)

    def test_downlink(self, P=150, L_l=0.6, L_r=0.8, L_a=-0.1, DR=7656, Tnoise=200, k=2.6 * 10 ** -23, d_ant=2.5,
                      d_ground=10.6,
                      wavelength=590 * 10 ** -9, eta_ant=0.8, PointAccSC=2 * 10 ** -14, freq=5.084745763 * 10 ** 14,
                      d=15.6, AU=1375683.9):
        E_bN_0_true = 167.1436941
        E_bN_0_numerical = downlink(P, L_l, L_r, L_a, DR, Tnoise, k, d_ant, d_ground, wavelength, eta_ant, PointAccSC,
                                    freq, d, AU)
        self.assertAlmostEqual(E_bN_0_numerical, E_bN_0_true, 5)

    def test_uplink(self,f_gs=5.472488039*10**14, P_gs=250, L_l=0.6, L_r=0.8, L_a=-0.1, uplinkDR=9628, Tnoise=200,
                    k=2.6*10**-23, wavelength=5.48*10**-7, eta_ant=0.8, d_ground=10.6,d_ant=2.5, d=15.6, AU=1375683.9,
                    PointAccSC=2*10**-14):
        E_bN_0_true = 186.392728
        E_bN_0_numerical = uplink(f_gs, P_gs, L_l, L_r, L_a, uplinkDR, Tnoise, k, wavelength, eta_ant,
                                  d_ground,d_ant, d, AU, PointAccSC)
        self.assertAlmostEqual(E_bN_0_numerical, E_bN_0_true,6)
