import os.path

from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
import pandas as pd
import math

# Load tudatpy modules
from tudatpy.kernel.interface import spice
from tudatpy.kernel import numerical_simulation
from tudatpy.kernel.numerical_simulation import environment_setup, environment, propagation_setup, propagation
from tudatpy.kernel.astro import element_conversion
from tudatpy.kernel import constants
from tudatpy.util import result2array

from atmospheric_model import GRAM

_path = os.path.dirname(__file__)

spice.load_standard_kernels()
spice.load_kernel(_path + '/../../../atmospheric_model/GRAM/GRAM_Suite_1_5/SPICE/spk/satellites/ura116xl.bsp')


def decelerator_sizing(target_time, totalmass, heatshieldmass, capsule_drag_coefficient, diameter, alt, lat, lon, speed,
                       flight_path_angle, heading_angle, acc=1, steps=5):
    pass

if __name__ == "__main__":
    decelerator_sizing(2 * 3600, 500, 250, 1.53, 4.5, 3.02877105e+07, -6.40748300e-02, -1.63500310e+00 + 2 * np.pi,
                       1.93919454e+04, np.deg2rad(-30), -2.35413606e+00)