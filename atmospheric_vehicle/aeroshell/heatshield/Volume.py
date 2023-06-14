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
    """print(vpos, vneg)

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
    plt.show() #"""

    return vpos - vneg