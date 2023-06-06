import numpy as np
from matplotlib import pyplot as plt

# Variables
S = 1  # m^2
rho = 1  # kg/m^3
C_D = 1
C_L = 1
m = 1  # kg
V_0 = 1  # m/s
g = 1  # m/s^2
gamma = 1  # rad
chi = 1  # rad
V_x_wind = 1  # m/s
g_u = 1  # m/s^2


# Functions
def velocity(t):
    z = rho * S * C_D / m
    a = V_0 * (1 + 2 / z)
    b = 1 / z / V_0
    c = g_u * np.sin(gamma) / (z * V_0)
    d = np.exp(-z * V_0 * V_0 * t)
    k = (g_u * np.sin(gamma) + 2 * V_0 * V_0 * V_0 - V_0) / (2 * V_0 * V_0)
    V = (a - b + c) * d + k
    return V


def range(t):
    z = rho * S * C_D / m
    a = np.cos(gamma) * np.cos(chi)
    b = V_0 * (1 + 2 / z)
    c = -1 / (z * V_0 * V_0)
    d = -1 / (z * V_0)
    e = g_u * np.sin(gamma) / (z * V_0)
    f = np.exp(-z * V_0 * V_0 * t)
    g = (g_u * np.sin(gamma) + 2 * V_0 * V_0 * V_0 - V_0) / (2 * V_0 * V_0)
    x = a * ((b * c + d + e) * f + g * t) + V_x_wind * t
    return x


def height(t):
    z = rho * S * C_D / m
    a = -np.sin(gamma) / (z * V_0 * V_0)
    b = V_0 * (1 + 2 / z)
    c = -1 / (z * V_0)
    d = g_u * np.sin(gamma) / (z * V_0)
    e = np.exp(-z * V_0 * V_0 * t)
    h = a * (b + c + d) * e
    return h


t = np.linspace(0, 10, 100)

if __name__ == "__main__":
    # print("Hello World")
    V = velocity(t)
    x = range(t)
    h = height(t)

    plt.plot(t, V, label="Velocity")
    plt.plot(t, x, label="Range")
    plt.plot(t, h, label="Depth")
    plt.legend()
    plt.show()
