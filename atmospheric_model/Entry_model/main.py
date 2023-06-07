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

        self._lat = 0.0
        self._lon = 0.0
        self.heading = 0.0

        self.alt = 0.0
        self.dis = 0.0
        self.vx = 0.0
        self.vy = 0.0
        self.ax = 0.0
        self.ay = 0.0

        self.t = 0.0

    @property
    def y(self):
        return np.array([self.alt,
                         self.dis,
                         self.vy,
                         self.vx])

    @property
    def dy(self):
        return np.array([self.vy,
                         self.vx,
                         self.ay,
                         self.ax])

    @y.setter
    def y(self, y):
        self.alt = y[0]
        self.dis = y[1]
        self.vy = y[2]
        self.vx = y[3]


def uranus_entry_sim(alt, lat, lon, vel, gamma, heading, mass, beta, c_l, c_d):
    gram =GRAM.GRAM()
    gram.altitudes = np.array(0.0)
    gram.lat = np.array(lat)
    gram.long = np.array(lon)
    gram.time = np.array(0.0)
    gram.run()
    radius = gram.data.LatitudeRadius_km[0]
    core = Core()

    core.mass = mass
    core.c_d = c_d
    core.c_l = c_l
    core.beta = beta
    core.gamma = gamma

    core._lat = lat
    core._lon = lon
    core.heading = heading

    vx = vel * np.cos(gamma)
    vy = vel * np.sin(gamma)
    y0 = np.asarray([alt - radius, 0.0, vx, vy])

    def rhs(t, y):
        core.y = y
        core.t = t

        alpha = np.arctan2(core.vy, core.vx)
        gram.altitudes = np.array(core.alt)
        gram.run()
        core.ax = -1 / 2 * gram.data.Density_kgm3[0] / (np.cos(alpha) * core.beta) * \
                  (core.vx ** 2 - core.vx * core.vy * core.c_l / core.c_d)
        core.ay = 1 / 2 * gram.data.Density_kgm3[0] / (np.sin(alpha) * core.beta) * \
                  (core.vy ** 2 + core.vx * core.vy * core.c_l / core.c_d) + gram.data.Gravity_ms2[0]
        return core.dy

    def event_dis(t, y):
        return 50000 - y[1]

    event_dis.terminal = True

    def event_alt(t, y):
        return y[0]

    event_alt.terminal = True

    sol = solve_ivp(rhs, [core.t, 18000], y0,
                    events=[event_dis, event_alt],
                    max_step=0.001, method="LSODA",
                    rtol=1e-3, atol=1e-6)

    return sol

if __name__ == "__main__":
    example_data = [3.02877105e+07, -6.40748300e-02, -1.63500310e+00,  1.93919454e+04, -4.15342314e-01, -2.35413606e+00]

    print(np.arctan2(0, 1))
    print(np.arctan2(1, 0))
    print(np.arctan2(0, -1))
    print(np.arctan2(-1, 0))
