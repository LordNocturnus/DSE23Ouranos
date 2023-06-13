from unittest import TestCase
from Orbiter.main import *

orbiter = Orb(False)


class TestOrb(TestCase):
    def test_mass_prop(self):
        orbiter.mass_prop(100)
        self.assertEqual(100, orbiter.wet_mass - orbiter.prop_mass)
        self.assertEqual(orbiter.m_ox, orbiter.prop_mass - orbiter.m_fuel)
        self.assertGreaterEqual(orbiter.t_cy_f, orbiter.t_caps_f)
        self.assertGreaterEqual(orbiter.t_cy_o, orbiter.t_caps_o)

    def test_mass_prop2(self):
        orbiter.mixture_ratio = 1
        orbiter.mass_prop(100)
        self.assertEqual(orbiter.m_ox, orbiter.m_fuel)

    def test_mass_prop3(self):
        with self.assertRaises(ValueError):
            orbiter.mass_AV = -3
            orbiter.mass_prop(100)
            orbiter.mass_AV = 200

    def test_power1(self):
        orbiter.P_req = 500
        orbiter.power()
        self.assertEqual((55.9 * orbiter.n_rtg + 25) * 1.2, orbiter.m_power)

    def test_power2(self):
        orbiter.P_req = 500
        orbiter.power()
        self.assertEqual(145699633.36 * orbiter.n_rtg, orbiter.cost_rtg)

    def test_thermal(self):
        orbiter.thermal()
        self.assertEqual(orbiter.r_tanks * orbiter.l_tanks, orbiter.A_rec)
        self.assertEqual(4, len(orbiter.n_l_closed))
        self.assertGreaterEqual(orbiter.n_l_closed[0], orbiter.n_l_closed[-1])
        self.assertEqual(4.08 * orbiter.n_rtg, round(orbiter.m_louvres, 2))

