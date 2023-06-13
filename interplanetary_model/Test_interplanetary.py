from unittest import TestCase
from interplanetary_model.arrival_departure import *
from interplanetary_model.interplanetary import *





class Test(TestCase):
    def test_earth_departure(self):
        a_numerical,e_numerical = earth_departure(1e6,8269.761087)
        a_compared = 1e9
        e_compared = 1.0
        self.assertGreater(a_compared,a_numerical)
        self.assertAlmostEqual(e_numerical,e_compared)
