from matplotlib import pyplot as plt
import numpy as np
import scipy as sp

_dia = 4.5
_radius1 = 1.125
_radius2 = 0.126
_bottom_angle = np.deg2rad(20)
_top_angle = np.deg2rad(59.73)


def heatshield_sizing(diameter, heat_load, interface_velocity, peak_heat_flux):
    radius1 = _radius1 / _dia * diameter
    radius2 = _radius2 / _dia * diameter

    PICA_thickness = 1.8686 * (heat_load / 10000 / (interface_velocity / 1000)) ** 0.1879 / 100
    PICA_density = 0.352 * 100**3 / 1000  # 0.352 - 0.701 g/cm^3
    CP_thickness = 1.1959 * (heat_load / 10000 / (interface_velocity / 1000)) ** 0.2102 / 100
    CP_density = 1.31 * 100**3 / 1000  # 1.31 - 1.55 g/cm^3
    HT_424_thickness = 0.0381 / 100
    HT_424_density = 0.66 / HT_424_thickness
    ACC6_thickness = 0.25 / 100
    ACC6_density = 1.6 * 100**3 / 1000  # 1.6 - 2.1 g/cm^3

    PICA = False
    CP = False
    if peak_heat_flux <= 1200 and 3.27 / 100 <= PICA_thickness <= radius2:
        PICA = True
        PICA_volume = volume(diameter, radius1, radius2, 0.0, PICA_thickness)
        PICA_weight = PICA_volume * PICA_density

    if 2.266 / 100 <= CP_thickness and CP_thickness + HT_424_thickness + ACC6_thickness <= radius2:
        CP = True
        CP_volume = volume(diameter, radius1, radius2, 0.0, CP_thickness)
        HT_424_volume = volume(diameter, radius1, radius2, CP_thickness, HT_424_thickness)
        ACC6_volume = volume(diameter, radius1, radius2, CP_thickness + HT_424_thickness, ACC6_thickness)
        CP_weight = CP_volume * CP_density
        HT_424_weight = HT_424_volume * HT_424_density
        ACC6_weight = ACC6_volume * ACC6_density
        CP_total_weight = CP_weight + HT_424_weight + ACC6_weight

    if CP and PICA:
        if CP_total_weight < PICA_weight:
            print("choose PICA")
            return PICA_weight
        else:
            print("choose CP")
        return CP_total_weight
    elif PICA:
        print("choose PICA")
        return PICA_weight
    elif CP:
        print("choose CP")
        return CP_total_weight
    else:
        raise ValueError("No heatshield feasible")


def volume(diameter, radius1, radius2, depth, thickness):
    l0 = -(radius1 - depth)
    l1 = -(radius1 - depth) * np.cos(_bottom_angle)
    l2 = (diameter / 2 - radius2 + (radius2 - depth) * np.cos(np.pi / 2 - _bottom_angle) -
          (radius1 - depth) / np.sin(_bottom_angle)) / np.tan(np.pi / 2 - _bottom_angle)
    l3 = (diameter / 2 - radius2 + (radius2 - depth) * np.cos(np.pi / 2 - _bottom_angle) -
          (radius1 - depth) / np.sin(_bottom_angle)) / np.tan(np.pi / 2 - _bottom_angle) + \
         (radius2 - depth) * (np.sin(np.pi / 2 - _bottom_angle) + np.sin(_top_angle))

    l0d = -(radius1 - depth - thickness)
    l1d = -(radius1 - depth - thickness) * np.cos(_bottom_angle)
    l2d = (diameter / 2 - radius2 + (radius2 - depth - thickness) * np.cos(np.pi / 2 - _bottom_angle) -
           (radius1 - depth - thickness) / np.sin(_bottom_angle)) / np.tan(np.pi / 2 - _bottom_angle)
    l3d = (diameter / 2 - radius2 + (radius2 - depth - thickness) * np.cos(np.pi / 2 - _bottom_angle) -
           (radius1 - depth - thickness) / np.sin(_bottom_angle)) / np.tan(np.pi / 2 - _bottom_angle) + \
          (radius2 - depth - thickness) * (np.sin(np.pi / 2 - _bottom_angle) + np.sin(_top_angle))


    vpos = sp.integrate.quad(lambda x: np.pi * s1(x, depth, radius1, radius2, _bottom_angle, _top_angle, diameter) ** 2,
                             l0, l1)[0]
    vneg = sp.integrate.quad(lambda x: np.pi * s1(x, depth + thickness, radius1, radius2, _bottom_angle, _top_angle,
                                                  diameter) ** 2, l0d, l1d)[0]

    vpos += sp.integrate.quad(lambda x: np.pi * s2(x, depth, radius1, radius2, _bottom_angle, _top_angle, diameter) ** 2,
                              l1, l2)[0]
    vneg += sp.integrate.quad(lambda x: np.pi * s2(x, depth + thickness, radius1, radius2, _bottom_angle, _top_angle,
                                                   diameter) ** 2, l1d, l2d)[0]

    vpos += sp.integrate.quad(lambda x: np.pi * s3(x, depth, radius1, radius2, _bottom_angle, _top_angle, diameter) ** 2,
                              l2, l3)[0]
    vneg += sp.integrate.quad(lambda x: np.pi * s3(x, depth + thickness, radius1, radius2, _bottom_angle, _top_angle,
                                                   diameter) ** 2, l2d, l3d)[0]

    edge = s3(l3, depth, radius1, radius2, _bottom_angle, _top_angle, diameter)
    edged = s3(l3d, depth + thickness, radius1, radius2, _bottom_angle, _top_angle, diameter)
    vneg += sp.integrate.quad(lambda x: np.pi * ((edge - edged) / (l3 - l3d) * (-l3 + x) + edge) ** 2, l3d, l3)[0]
    print(vpos, vneg)

    plt.plot(np.linspace(l0, l1, 1000), s1(np.linspace(l0, l1, 1000), depth, radius1, radius2, _bottom_angle,
                                           _top_angle, diameter))
    plt.plot(np.linspace(l0d, l1d, 1000), s1(np.linspace(l0d, l1d, 1000), depth + thickness, radius1, radius2,
                                             _bottom_angle, _top_angle, diameter))

    plt.plot(np.linspace(l1, l2, 1000), s2(np.linspace(l1, l2, 1000), depth, radius1, radius2, _bottom_angle,
                                           _top_angle, diameter))
    plt.plot(np.linspace(l1d, l2d, 1000), s2(np.linspace(l1d, l2d, 1000), depth + thickness, radius1, radius2,
                                             _bottom_angle, _top_angle, diameter))

    plt.plot(np.linspace(l2, l3, 1000), s3(np.linspace(l2, l3, 1000), depth, radius1, radius2, _bottom_angle,
                                           _top_angle, diameter))
    plt.plot(np.linspace(l2d, l3d, 1000), s3(np.linspace(l2d, l3d, 1000), depth + thickness, radius1, radius2,
                                             _bottom_angle, _top_angle, diameter))

    plt.plot(np.linspace(l3d, l3, 1000), (edge - edged) / (l3 - l3d) * (-l3 + np.linspace(l3d, l3, 1000)) + edge)
    plt.grid()
    plt.xlim(-2.5, 3)
    plt.ylim(-0.5, 5)
    plt.show()

    return vpos - vneg


def s1(x, d, r1, r2, a1, a2, dia):
    return (r1 - d) * np.sin(np.arccos(x/(r1 - d)))


def s2(x, d, r1, r2, a1, a2, dia):
    return np.tan(np.pi / 2 - a1) * x + (r1 - d) / np.sin(a1)


def s3(x, d, r1, r2, a1, a2, dia):
    x -= (dia / 2 - r2 + (r2 - d) * np.cos(np.pi / 2 - a1) -
          (r1 - d) / np.sin(a1)) / np.tan(np.pi / 2 - a1) + (r2 - d) * np.sin(np.pi / 2 - a1)
    return dia / 2 - r2 + (r2 - d) * np.cos(np.arcsin(x / (r2 - d)))


if __name__ == "__main__":
    print(heatshield_sizing(4.5, 10**8, 20000, 1100))
