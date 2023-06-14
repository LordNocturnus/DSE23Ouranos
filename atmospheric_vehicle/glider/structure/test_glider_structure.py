from unittest import TestCase
from atmospheric_vehicle.glider.structure.main import *

c_w_root = 0.785  # [m], root chord, given from Luigi
c_w_tip = 0.713  # [m], tip chord, given from Luigi
b_w = 6  # [m], wing span, given from Luigi
c_w = calculate_c_w(c_w_root, c_w_tip, b_w)[0]  # array of chord values that will be used for later estimations
b_range = calculate_c_w(c_w_root, c_w_tip, b_w)[1]
db = calculate_c_w(c_w_root, c_w_tip, b_w)[2]

class GliderStructureTests(TestCase):
    def test_calculate_c_w1(self):
        self.assertEqual(len(c_w), len(b_range))
