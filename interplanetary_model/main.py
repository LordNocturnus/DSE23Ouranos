from astropy import units as u
from matplotlib import pyplot as plt
import numpy as np
from typing import List, Tuple

# Tudat imports
import tudatpy
from tudatpy.kernel.trajectory_design import transfer_trajectory
from tudatpy.kernel import constants
from tudatpy.kernel.numerical_simulation import environment_setup
from tudatpy.util import result2array

# Pygmo imports
import pygmo as pg


if __name__ == "__main__":
    print("Hello, world")
