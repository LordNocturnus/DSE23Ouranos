from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt
import numpy as np
import scipy as sp

from atmospheric_model import GRAM


class Core:

    def __init__(self):
        self.mass = 0.0
        self.c_d = 0.0
        self.c_l = 0.0
        self.beta = 0.0
        self.gamma = 0.0

        self.radius = 0.0
        self.heading = 0.0

        self.alt = 0.0
        self.lat = 0.0
        self.lon = 0.0
        self.heat = 0.0
        self._vx = 0.0
        self.dlat = 0.0
        self.dlon = 0.0
        self.vy = 0.0
        self.q = 0.0
        self.ax = 0.0
        self.ay = 0.0

        self.t = 0.0

    @property
    def vx(self):
        return self._vx

    @property
    def y(self):
        return np.array([self.alt,
                         self.lat,
                         self.lon,
                         self.vy,
                         self._vx,
                         self.heat])

    @property
    def dy(self):
        return np.array([self.vy,
                         self.dlat,
                         self.dlon,
                         self.ay,
                         self.ax,
                         self.q])

    @y.setter
    def y(self, y):
        self.alt = y[0]
        self.lat = np.arcsin(np.sin(y[1]))
        self.lon = y[2] % np.pi
        if self.lon < 0.0:
            self.lon += 2 * np.pi
        self.vy = y[3]
        self.vx = y[4]
        self.heat = y[5]

    @vx.setter
    def vx(self, other):
        self._vx = other
        self.dlat = self.vx / (self.alt + self.radius) * np.sin(self.heading)
        self.dlon = self.vx / (self.alt + self.radius) * np.cos(self.heading)


def uranus_entry_sim(alt, lat, lon, vel, gamma, heading, mass, beta, c_l, c_d):
    gram =GRAM.GRAM()
    gram.altitudes = np.array([0.0])
    gram.lat = np.array([np.rad2deg(lat)])
    gram.long = np.array([np.rad2deg(lon)])
    gram.time = np.array([0.0])
    gram.run()

    core = Core()
    core.radius = gram.data.LatitudeRadius_km[0]

    core.mass = mass
    core.c_d = c_d
    core.c_l = c_l
    core.beta = beta
    core.gamma = gamma

    core._lat = lat
    core._lon = lon
    core.heading = heading

    vx = vel * np.sin(gamma)
    vy = vel * np.cos(gamma)
    y0 = np.asarray([alt - core.radius * 1000, lat, lon, vx, vy, 0.0])

    def rhs(t, y):
        core.y = y
        core.t = t

        alpha = np.arctan2(core.vy, core.vx)
        gram.altitudes = np.array([core.alt / 1000])
        gram.lat = np.array([np.rad2deg(core.lat)])
        gram.long = np.array([np.rad2deg(core.lon)])
        gram.run()
        core.radius = gram.data.LatitudeRadius_km[0]
        core.ax = -1 / 2 * gram.data.Density_kgm3[0] / (core.beta * np.cos(alpha)) * \
                  core.vx ** 2 * (core.c_l / core.c_d * np.tan(alpha) + 1)
        core.ay = 1 / 2 * gram.data.Density_kgm3[0] / (core.beta * np.sin(alpha)) * \
                  core.vy ** 2 * (core.c_l / core.c_d * 1 / np.tan(alpha) - 1) - gram.data.Gravity_ms2[0]
        core.q = (core.vx ** 2 + core.vy ** 2) ** (3/2) * gram.data.Density_kgm3[0]
        return core.dy

    def event_alt(t, y):
        print("t: ", t)
        print("alt: ", y[0])
        print("vx: ", y[4])
        return y[0]

    event_alt.terminal = True

    sol = solve_ivp(rhs, [core.t, 18000], y0,
                    events=event_alt,
                    max_step=1, method="LSODA",
                    rtol=1e-3, atol=1e-6)

    return sol


if __name__ == "__main__":
    example_5000 = [3.02877105e+07, -6.40748300e-02, -1.63500310e+00,  1.93919454e+04, -4.15342314e-01, -2.35413606e+00]

    example_2500 = [2.78349356e+07, -2.20565560e-01, -2.91545691e+00,  2.02425574e+04, -3.05261011e-01, 2.33104808e+00]

    test = uranus_entry_sim(example_2500[0]-1700000, example_2500[1], example_2500[2] + 2 * np.pi, example_2500[3],
                            np.deg2rad(-10), example_2500[5], 500, 125, 0, 1)

    gram = GRAM.GRAM()
    gram.altitudes = test.y[0] / 1000
    gram.lat = np.rad2deg(np.arcsin(np.sin(test.y[1])))
    gram.long = np.rad2deg(test.y[2] % np.pi)
    gram.long[gram.long < 0.0] += 2 * np.pi
    gram.time = test.t
    gram.run()

    plt.plot(test.y[3] / 1000, test.y[0] / 1000, label="Vx")
    plt.plot(test.y[4] / 1000, test.y[0] / 1000, label="Vy")
    plt.plot(np.sqrt(test.y[3] ** 2 + test.y[4] ** 2) / 1000, test.y[0] / 1000, label="V_total")
    plt.legend()
    plt.show()

    plt.plot(test.t, test.y[0] / 1000, label="alt")
    plt.legend()
    plt.show()

    plt.plot(test.t, (test.y[3] ** 2 + test.y[4] ** 2) ** (3/2) * np.asarray(gram.data.Density_kgm3), label="heat flux")
    plt.legend()
    plt.show()

    plt.plot(test.t, test.y[5], label="heat load")
    plt.legend()
    plt.show()

    alpha = np.arctan2(test.y[4], test.y[3])
    ax = -1 / 2 * np.asarray(gram.data.Density_kgm3) / (125 * np.cos(alpha)) * \
         test.y[3] ** 2 * (0 / 1 * np.tan(alpha) + 1)
    ay = 1 / 2 * np.asarray(gram.data.Density_kgm3) / (125 * np.sin(alpha)) * \
         test.y[4] ** 2 * (0 / 1 * 1 / np.tan(alpha) - 1) - np.asarray(gram.data.Gravity_ms2)
    plt.plot(ax, test.y[0] / 1000, label="ax")
    plt.plot(ay, test.y[0] / 1000, label="ay")
    plt.plot(np.sqrt(ax ** 2 + ay ** 2), test.y[0] / 1000, label="a_total")
    plt.legend()
    plt.show()
