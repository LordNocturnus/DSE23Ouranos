from unittest import TestCase
from atmospheric_vehicle.aeroshell.structure.main import *


class Test(TestCase):
    def test_angle_cone1(self):
        self.assertEqual(0.2915, round(angle_cone(5, 2, 10), 4))
        self.assertEqual(0, angle_cone(5, 5, 10))

    def test_angle_cone2(self):
        with self.assertRaises(ZeroDivisionError):
            a = angle_cone(1, 1, 0)

    def test_angle_cone3(self):
        with self.assertRaises(ValueError):
            a = angle_cone(0, 1, 1)

    def test_angle_cone4(self):
        with self.assertRaises(ValueError):
            a = angle_cone(1, 0, 1)

    def test_angle_cone5(self):
        with self.assertRaises(ValueError):
            a = angle_cone(-1, 1, 1)

    def test_angle_cone6(self):
        with self.assertRaises(ValueError):
            a = angle_cone(1, 1, -1)

    def test_t_pressure1(self):
        self.assertEqual(0.7522, round(t_pressure(10 * 10 ** 6, 5, np.pi / 4, 100 * 10 ** 6), 4))
        self.assertEqual(0, t_pressure(10, 0, np.pi / 4, 1))
        self.assertEqual(0.7522, round(t_pressure(-10 * 10 ** 6, 5, np.pi / 4, 100 * 10 ** 6), 4))

    def test_t_pressure2(self):
        with self.assertRaises(ZeroDivisionError):
            a = t_pressure(10, 5, np.pi / 2, 1)

    def test_t_pressure3(self):
        with self.assertRaises(ZeroDivisionError):
            a = t_pressure(10, 5, np.pi / 4, 6)

    def test_volume_truncated(self):
        self.assertEqual(51.31, round(volume_truncated(5, 3, 1), 2))

    def test_volume_truncated2(self):
        with self.assertRaises(ValueError):
            a = volume_truncated(0, 1, 1)

    def test_volume_truncated3(self):
        with self.assertRaises(ValueError):
            a = volume_truncated(1, 0, 1)

    def test_volume_truncated4(self):
        with self.assertRaises(ValueError):
            a = volume_truncated(1, 1, -1)

    def test_buckling1(self):
        self.assertTrue(buckling(100, 1, 10, 1000))
        self.assertTrue(buckling(100, 1, 10, -1000))
        self.assertTrue(buckling(1, 1, 1, np.pi ** 2))

    def test_buckling2(self):
        with self.assertRaises(ZeroDivisionError):
            a = buckling(1, 0, 1, 1)

    def test_buckling3(self):
        with self.assertRaises(ValueError):
            a = buckling(-1, 1, 1, 1)

    def test_buckling4(self):
        with self.assertRaises(ValueError):
            a = buckling(1, 1, -1, 1)

    def test_buckling5(self):
        with self.assertRaises(ValueError):
            a = buckling(1, -1, 1, 1)

