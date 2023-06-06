import numpy as np
from matplotlib import pyplot as plt

# Variables
S = 12  # m^2
rho = 0.42  # kg/m^3
C_D = 0.05
C_L = 0.7
m = 72  # kg
V_0 = 23.6  # m/s
gamma = np.radians(-15)  # rad
chi = 0  # rad
V_x_wind = 250  # m/s
g_u = 8.69  # m/s^2


# Functions
def velocity(t):
    z = rho * S * C_D / (2 * m)
    a = z * V_0 - g_u * np.sin(gamma)
    b = 2 * z * V_0 * V_0
    V = a / b * np.exp(-b * t) + (2 * z * V_0 * V_0 * V_0 - a) / b
    return V  # m/s


def range(t):
    z = rho * S * C_D / (2 * m)
    a = z * V_0 - g_u * np.sin(gamma)
    b = 2 * z * V_0 * V_0
    x = np.cos(gamma) * np.cos(chi) * (-a / (b * b) * np.exp(-b * t) + (2 * z * V_0 * V_0 * V_0 - a) / b) + V_x_wind * t
    return x  # m


def height(t):
    z = rho * S * C_D / (2 * m)
    a = z * V_0 - g_u * np.sin(gamma)
    b = 2 * z * V_0 * V_0
    h = (-a / (b * b) * np.exp(-b * t) + (2 * z * V_0 * V_0 * V_0 - a) / b) * np.sin(gamma)
    return h  # m


t = np.linspace(0, 10, 100)

if __name__ == "__main__":
    V = velocity(t)
    x = range(t)
    h = height(t)
    plt.plot(t, V, label="Velocity")
    plt.plot(t, x, label="Range")
    plt.plot(t, h, label="Depth")
    plt.legend()
    plt.show()
