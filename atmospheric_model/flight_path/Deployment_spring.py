import numpy as np
import matplotlib.pyplot as plt


def spring_constant(x, m=162.1, g=8.7, alpha=np.deg2rad(20), v=200):
    a = m * (v / x) ** 2
    b = 3 * m * g / x * np.sin(alpha)
    k = a - b
    return k


x = np.linspace(4.5, 100, 1000)
n = np.linspace(1, 500, 1000)
k = spring_constant(x)
k_n = k / n

plt.subplot(121)
plt.plot(x, k / 1000)
plt.xlabel("Ramp Length [m]")
plt.ylabel("Spring Constant [kN/m]")
plt.ylim(0, 10)
plt.subplot(122)
plt.plot(n, k_n)
plt.xlabel("Number of Springs")
plt.ylabel("Spring Constant per Spring [kN/m]")
plt.ylim(0, 10)
plt.show()
