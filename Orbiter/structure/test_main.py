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

    def test_pressurising_gas1(self):
        propellant = (100 * 10 ** 6, 1)
        self.assertEqual(0.090817, round(pressurising_gas(propellant, 10), 6))

    def test_pressurising_gas2(self):
        propellant = 1
        with self.assertRaises(TypeError):
            a = pressurising_gas(propellant, 10)

    def test_pressurising_gas3(self):
        propellant = (-100 * 10 ** 6, 1)
        with self.assertRaises(ValueError):
            a = pressurising_gas(propellant, 10)

    def test_pressurising_gas4(self):
        propellant = (0, 1)
        with self.assertRaises(ZeroDivisionError):
            a = pressurising_gas(propellant, 10)

    def test_natural_frequency1(self):
        material = [1, 1, 1, 9 * 10 ** 9]
        self.assertEqual(48.725, round(natural_frequency(3, 1, 0.001, material, 100, 10)[0], 3))
        self.assertEqual(2360.541, round(natural_frequency(3, 1, 0.001, material, 100, 10)[1], 3))

    def test_natural_frequency2(self):
        material = [1, 1, 1, 9 * 10 ** 9]
        with self.assertRaises(ZeroDivisionError):
            a = natural_frequency(3, 1, 0.001, material, 0, 0)[0]

    def test_natural_frequency3(self):
        material = [1, 1, 1, 9 * 10 ** 9]
        with self.assertRaises(ValueError):
            a = natural_frequency(3, 1, 0.001, material, -100, 1)[0]

    def test_natural_frequency4(self):
        material = [1, 1, 1, 9 * 10 ** 9]
        with self.assertRaises(ValueError):
            a = natural_frequency(3, 1, 0.001, material, 1, -100)[0]

    def test_natural_frequency5(self):
        material = [1, 1, 1, 9 * 10 ** 9]
        with self.assertRaises(ValueError):
            a = natural_frequency(3, -1, 0.001, material, 1, 100)[0]

    def test_natural_frequency6(self):
        material = [1, 1, 1, 9 * 10 ** 9]
        with self.assertRaises(ValueError):
            a = natural_frequency(3, 1, -0.001, material, 1, 100)[0]


class Test_Integration(TestCase):
    def test_geometry_mass(self):
        material = [4430, 880 * 10 ** 6, 970 * 10 ** 6, 113.8 * 10 ** 9]
        propellant = (3 * 10 ** 6, 200, 1431)
        m, l, r, t_caps, t_cyl, v_tot = geometry_mass(material, propellant, 0.2, False)
        self.assertEqual(np.pi * r ** 2 * l + 4 / 3 * np.pi * r ** 3, v_tot)
        self.assertGreater(t_caps, t_cyl / 2)
        self.assertGreaterEqual(t_caps, 1 * 10**-3)
        self.assertGreaterEqual(t_cyl, 1 * 10**-3)

    def test_geometry_mass2(self):
        material = [4430, 880 * 10 ** 6, 970 * 10 ** 6, 113.8 * 10 ** 9]
        with self.assertRaises(TypeError):
            a, args = geometry_mass(material, 1, 0.2, False)

    def test_geometry_mass3(self):
        propellant = (3 * 10 ** 6, 200, 1431)
        with self.assertRaises(TypeError):
            a, args = geometry_mass(1, propellant, 0.2, False)

    def test_geometry_mass4(self):
        propellant = (3 * 10 ** 6, 200, 1431)
        material = [4430, 880 * 10 ** 6, 970 * 10 ** 6, 113.8 * 10 ** 9]
        with self.assertRaises(ValueError):
            a, args = geometry_mass(material, propellant, -1, False)

    def test_final_architecture1(self):
        propellant = [(3 * 10 ** 6, 200, 1431), (3 * 10 ** 6, 200, 1431)]
        material = [4430, 880 * 10 ** 6, 970 * 10 ** 6, 113.8 * 10 ** 9]
        l_tot, r, m, t_cy_o, t_cy_f, t_cap_o, t_cap_f = final_architecture(material, propellant, 0.2, 100)
        self.assertEqual(30, round(m, 0))
        self.assertGreaterEqual(t_cy_o, t_cap_o)
        self.assertGreaterEqual(t_cy_f, t_cap_f)
        self.assertEqual(0.052, round((l_tot - r * 4) / 2, 3))

    def test_final_architecture2(self):
        with self.assertRaises(TypeError):
            propellant = [1, 2, 3, 4]
            material = [4430, 880 * 10 ** 6, 970 * 10 ** 6, 113.8 * 10 ** 9]
            a, args = final_architecture(material, propellant, 0.2, 100)

    def test_final_architecture3(self):
        with self.assertRaises(TypeError):
            propellant = [(3 * 10 ** 6, 200, 1431), (3 * 10 ** 6, 200, 1431)]
            material = 1
            a, args = final_architecture(material, propellant, 0.2, 100)





