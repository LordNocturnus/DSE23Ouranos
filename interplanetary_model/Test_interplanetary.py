from unittest import TestCase
from interplanetary_model.arrival_departure import *
#from interplanetary_model.interplanetary import *





class Test(TestCase):
    def test_earth_departure_extremecase(self):
        a_numerical,e_numerical = earth_departure(1e6,8269.761087)
        a_compared = 1e9
        e_compared = 1.0
        self.assertGreater(a_compared,a_numerical)
        self.assertGreater(e_numerical,e_compared)
    def test_earth_departure_value(self):
        a_numerical,e_numerical = earth_departure(1e6,20035.0307)
        a_analytical = -993024.4145
        e_compared = 1.0
        self.assertLess(a_analytical-a_numerical,0.0001)
        self.assertGreater(e_numerical,e_compared)        
