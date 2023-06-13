import numpy as np
import matplotlib.pyplot as plt
v_leo = 7800


def earth_departure(pericentre, dv, grav_para_earth=3.986e14):
    v_leo = np.sqrt(grav_para_earth/(pericentre))
    v_inf = np.sqrt((v_leo + dv)**2 - 2*grav_para_earth/pericentre)
    theta_inf = (2 * np.arcsin(1 / (1 + pericentre * v_inf * v_inf / grav_para_earth)) + np.pi) / 2
    e = -1/np.cos(theta_inf)
    a = -grav_para_earth/(v_inf*v_inf)
    return a, e


def uranus_arrival(pericentre, v_inf, grav_para_uranus: float = 5.793939e15):
    theta_inf = (2 * np.arcsin(1 / (1 + pericentre * v_inf * v_inf / grav_para_uranus)) + np.pi) / 2


def visual_e():
    theta = np.arange(95, 180)
    e_theta = []
    for angle in theta:
        e_theta.append(-1/np.cos(np.radians(angle)))

    mu_earth = 3.986e14
    r_p = 6600000
    v_inf = np.linspace(100, 10000, 100)
    arcsin_arg = []
    theta_v = []
    e_theta_v = []


    for v in v_inf:
        arcsin_arg.append(1 / (1 + ((r_p * v * v) / mu_earth)))
        theta_v.append((2 * np.arcsin(arcsin_arg[-1]) + np.pi) / 2)
        e_theta_v.append(-1/np.cos(theta_v[-1]))

    fig, ax1 = plt.subplots()
    ax1.plot(v_inf, e_theta_v, color='purple')
    ax1.set_xlabel("V_infty")
    ax1.set_ylabel("Eccentricity e")

    ax2 = ax1.twiny()
    ax2.plot(theta, e_theta, color='pink')
    ax2.set_xlabel("Theta_infty")


    plt.title("Eccentricity vs v_infty, theta_infty")

    plt.ylim(0.9, 4)
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()
