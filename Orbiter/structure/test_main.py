from unittest import *
from unittest import TestCase

from Orbiter.structure import *


class Test_Units(TestCase):

    def test_axial_loads1(self):
        self.assertFalse(axial_loads(1000, 10, 1))
        self.assertTrue(axial_loads(0, 1000, 1))
        self.assertTrue(axial_loads(200, 100, 1))

    def test_axial_loads2(self):
        with self.assertRaises(ZeroDivisionError):
            a = axial_loads(1000, 100, 0)

    def test_axial_loads3(self):
        with self.assertRaises(ValueError):
            a = axial_loads(1000, -1, 1)

    def test_lateral_loads1(self):
        self.assertTrue(lateral_loads(100, 1000, 1))
        self.assertTrue(lateral_loads(-100, 1000, 1))

    def test_lateral_loads2(self):
        with self.assertRaises(ValueError):
            a = lateral_loads(100, -100, 1)

    def test_lateral_loads3(self):
        with self.assertRaises(ValueError):
            a = lateral_loads(100, 100, -1)

    def test_t_hoop_sphere1(self):
        self.assertEqual(0.6, round(t_hoop_sphere(100, 4, 1000, 2), 4))
        self.assertEqual(0.6, round(t_hoop_sphere(-100, 4, 1000, 2), 4))
        self.assertEqual(0.2, round(t_hoop_sphere(100, 4, 1000, 0), 4))
        self.assertEqual(0, t_hoop_sphere(0, 1, 1, 2))
        self.assertEqual(0, t_hoop_sphere(100, 0, 1, 2))

    def test_t_hoop_sphere2(self):
        with self.assertRaises(ValueError):
            a = t_hoop_sphere(100, -4, 1000, 2)

    def test_t_hoop_sphere3(self):
        with self.assertRaises(ValueError):
            a = t_hoop_sphere(100, 4, -1000, 2)

    def test_t_hoop_sphere4(self):
        with self.assertRaises(ValueError):
            a = t_hoop_sphere(100, 4, 1000, -2)

    def test_t_hoop_sphere5(self):
        with self.assertRaises(ZeroDivisionError):
            a = t_hoop_sphere(100, 4, 0, 2)

    def test_t_hoop_cylind1(self):
        self.assertEqual(1.2, round(t_hoop_cylind(100, 4, 1000, 2), 4))
        self.assertEqual(1.2, round(t_hoop_cylind(-100, 4, 1000, 2), 4))
        self.assertEqual(0, t_hoop_cylind(0, 4, 1000, 2))
        self.assertEqual(0, t_hoop_cylind(100, 0, 1000, 2))
        self.assertEqual(0.4, t_hoop_cylind(100, 4, 1000, 0))

    def test_t_hoop_cylind2(self):
        with self.assertRaises(ValueError):
            a = t_hoop_cylind(100, -4, 1000, 1)

    def test_t_hoop_cylind3(self):
        with self.assertRaises(ValueError):
            a = t_hoop_cylind(100, 4, 1000, -1)

    def test_t_hoop_cylind4(self):
        with self.assertRaises(ZeroDivisionError):
            a = t_hoop_cylind(100, 4, 0, 2)

    def test_minimum_thickness1(self):
        self.assertEqual(936.785995, round(minimum_thickness(10000, 100, 1), 6))
        self.assertEqual(0, minimum_thickness(0, 100, 1))

    def test_minimum_thickness2(self):
        with self.assertRaises(ZeroDivisionError):
            a = minimum_thickness(10000, 0, 1)

    def test_minimum_thickness3(self):
        with self.assertRaises(ZeroDivisionError):
            a = minimum_thickness(10000, 100, 0)

    def test_minimum_thickness3(self):
        with self.assertRaises(ValueError):
            a = minimum_thickness(-1, 100, 1)

    def test_pressurising_gas(self):
        propellant = [100 * 10**6]
        self.assertEqual(0.090817, round(pressurising_gas(propellant, 10), 6))
