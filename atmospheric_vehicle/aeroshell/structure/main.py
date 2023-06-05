import numpy as np

# Atmospheric Properties
atm_pressure = ...


# Glider Parameters
d_glider = ...
h_glider = ...

# Launcher Parameters
d_launcher = ...


# Capsule Properties
cd_capsule = ...
d_capsule = ...


# Material Properties
rho_caps = ...


def drag(rho, v_entry, d, cd):
    return 1/2 * rho * v_entry**2 * np.pi * d**2 * cd / 4


def lift(rho, v_entry, d, cl):
    return 1/2 * rho * v_entry**2 * np.pi * d**2 * cl / 4

def axial_stress(L, D, m):
    ...




if __name__ == "__main__":
    print("Hello World")