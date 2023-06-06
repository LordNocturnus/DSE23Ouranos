import numpy as np
from matplotlib import pyplot as plt

# Variables
S = 12  # m^2
rho = 0.42  # kg/m^3
C_D = 0.05
C_L = 0.7
m = 72  # kg
V_0 = 23.6  # m/s
chi = 0  # rad
V_x_wind = 250  # m/s
g_u = 8.69  # m/s^2


# Functions
def velocity(t, gamma):
    z = rho * S * C_D / (2 * m)
    a = z * V_0 - g_u * np.sin(gamma)
    b = 2 * z * V_0 * V_0
    V = a / b * np.exp(-b * t) + (2 * z * V_0 * V_0 * V_0 - a) / b
    return V  # m/s


def range(t, gamma):
    z = rho * S * C_D / (2 * m)
    a = z * V_0 - g_u * np.sin(gamma)
    b = 2 * z * V_0 * V_0
    x = np.cos(gamma) * np.cos(chi) * (-a / (b * b) * np.exp(-b * t) + (2 * z * V_0 * V_0 * V_0 - a) / b) + V_x_wind * t
    return x  # m


def height(t, gamma):
    z = rho * S * C_D / (2 * m)
    a = z * V_0 - g_u * np.sin(gamma)
    b = 2 * z * V_0 * V_0
    h = (-a / (b * b) * np.exp(-b * t) + (2 * z * V_0 * V_0 * V_0 - a) / b) * np.sin(gamma)
    return h  # m


def gamma(rho, V_0):
    z = rho * S * C_D / (2 * m)
    gamma = np.arcsin(z * V_0 / g_u)
    return gamma


t = [0]
V = []
x = []
h = [0]
rho_0 = rho()
density = [rho_0]

if __name__ == "__main__":
    while h[-1] < 170 * 10 ** 3:
        if rho / rho_0 < 0.01:
            t.append(t[-1] + 0.1)
            gamma = gamma(rho_0, V_0)
            V = velocity(t[-1], gamma)
            x = range(t[-1], gamma)
            z = height(t[-1], gamma)
            density = rho(z)

            V.append(V)
            x.append(x)
            density.append(density)

            if h[0] == 0:
                h[0] = z
            else:
                h.append(z)

        else:
            V_0 = V[-1]
            rho_0 = density[-1]
            gamma = gamma(rho_0, V_0)

    plt.plot(t, V, label="Velocity")
    plt.plot(t, x, label="Range")
    plt.plot(t, h, label="Depth")
    plt.legend()
    plt.show()
