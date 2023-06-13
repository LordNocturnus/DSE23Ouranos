from unittest import TestCase
from Orbiter.comms.main import *


class Test(TestCase):
    def test_space_loss(self, wavelength=590 * 10 ** -9, d=15.6, AU=1375683.9):
        space_loss_true = -293.2000221
        space_loss_numerical = SpaceLoss(wavelength, d, AU)
        self.assertAlmostEquals(space_loss_numerical, space_loss_true)

    def test_gain(self, D=5.3, wavelength=590 * 10 ** -9, eta=0.76):
        gain_true = 147.8196105
        gain_numerical = gain(D, wavelength, eta)
        self.assertAlmostEquals(gain_numerical, gain_true)

    def test_halfpowerbeamwidth(self, D=5.6, f=5.084745763 * 10 ** 14):
        half_power_beam_width_true = 7.79245283 * 10 ** -15
        half_power_beam_width_numerical = halfpowerbeamwidth(D, f)
        self.assertAlmostEquals(half_power_beam_width_numerical, half_power_beam_width_true)

    def test_pointing_loss(self, alpha=7.79245283 * 10 ** -15, pointingAccuracy=2 * 10 ** -14):
        poiting_loss_true = -79.04836166
        poiting_loss_numerical = pointingLoss(alpha, pointingAccuracy)
        self.assertAlmostEquals(poiting_loss_numerical, poiting_loss_true)

    # def test_downlink(self):
    #
    #
    # def test_uplink(self):
    #
