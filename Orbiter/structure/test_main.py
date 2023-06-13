from unittest import *
from unittest import TestCase

from Orbiter.structure import *


class Test(TestCase):

    def test_axial_loads(self):
        self.assertFalse(axial_loads(1000, 10, 1))
        self.assertTrue(axial_loads(0, 1000, 1))

    @expectedFailure
    def test_axial_loads_failures(self):
        self.assertTrue(axial_loads(1000, 100, 0))

    @expectedFailure
    def test_axial_loads_failures2(self):
        self.assertTrue(axial_loads(1000, -1, 1))

    @expectedFailure
    def test_axial_loads_failures3(self):
        self.assertTrue(axial_loads(-200, 100, 1))

    def test_lateral_loads(self):
        self.assertTrue(lateral_loads(100, 1000, 1))
        self.assertTrue(lateral_loads(-100, 1000, 1))

    @expectedFailure
    def test_lateral_loads_failure(self):
        self.assertFalse(lateral_loads(100, -100, 1))

    @expectedFailure
    def test_lateral_loads_failure1(self):
        self.assertFalse(lateral_loads(100, 100, -1))

    @expectedFailure
    def test_lateral_loads_failure3(self):
        self.assertTrue(lateral_loads(100, -100, -1))

    def test_t_hoop_sphere(self):
        self.assertEqual(0.6, round(t_hoop_sphere(100, 4, 1000, 2), 4))
