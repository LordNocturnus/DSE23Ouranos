import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
ve_2 = 9.80665 * 454
m_shell = 1700


def func(x, dv_2, m_str2, m_prop2):
    return - m_shell - m_str2 + m_prop2/(np.exp(dv_2/ve_2)-1) - x


m_pl_data = [27200, 24600, 19000, 10800]  # payload to LEO
y = [0, 0, 0, 0]


popt, pcov = curve_fit(func, m_pl_data, y, bounds=(0, [4000, 100000, 400000]))
y_plot = []

for x in m_pl_data:
    f = func(x, popt[0], popt[1], popt[2])
    y_plot.append(f)

t = np.linspace(5000, 30000, 200)
y = []
for i in t:
    y.append(func(i, popt[0], popt[1], popt[2]))
plt.plot(t, y)
plt.show()

