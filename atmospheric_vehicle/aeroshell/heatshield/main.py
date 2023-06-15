from matplotlib import pyplot as plt
from tudatpy.kernel.numerical_simulation import environment_setup, environment, propagation_setup
import numpy as np
import scipy as sp

import atmospheric_model.Entry_model
from atmospheric_model.GRAM import GRAM

_dia = 4.5
_radius1 = 1.125
_radius2 = 0.126
_bottom_angle = np.deg2rad(20)
_top_angle = np.deg2rad(59.73)


def heatshield_sizing(diameter, heat_load, interface_velocity, peak_heat_flux):
    radius1 = _radius1 / _dia * diameter
    radius2 = _radius2 / _dia * diameter

    heat_load = heat_load / 10000 # convert to J/cm**2
    interface_velocity = interface_velocity / 1000 # convert to km/s
    peak_heat_flux = peak_heat_flux / 10000 # convert to W/cm**2

    PICA_thickness = 1.8686 * (heat_load / (interface_velocity ** 2)) ** 0.1879 / 100
    PICA_density = 0.352 * 100**3 / 1000  # 0.352 - 0.701 g/cm^3

    if 3.27 / 100 >= PICA_thickness:
        print("PICA thin")
        PICA_thickness = 3.27 / 100
    PICA_volume = volume(diameter, radius1, radius2, 0.0, PICA_thickness)
    PICA_weight = PICA_volume * PICA_density

    CP_thickness = 1.1959 * (heat_load / (interface_velocity ** 2)) ** 0.2102 / 100
    CP_density = 1.31 * 100**3 / 1000  # 1.31 - 1.55 g/cm^3
    HT_424_thickness = 0.0381 / 100
    HT_424_density = 0.66 / HT_424_thickness
    ACC6_thickness = 0.25 / 100
    ACC6_density = 1.6 * 100**3 / 1000  # 1.6 - 2.1 g/cm^3

    if 2.266 / 100 >= CP_thickness:
        print("CP thin")
        CP_thickness = 2.266 / 100

    CP_volume = volume(diameter, radius1, radius2, 0.0, CP_thickness)
    HT_424_volume = volume(diameter, radius1, radius2, CP_thickness, HT_424_thickness)
    ACC6_volume = volume(diameter, radius1, radius2, CP_thickness + HT_424_thickness, ACC6_thickness)
    CP_weight = CP_volume * CP_density
    HT_424_weight = HT_424_volume * HT_424_density
    ACC6_weight = ACC6_volume * ACC6_density
    CP_total_weight = CP_weight + HT_424_weight + ACC6_weight

    PICA = False
    CP = False
    if PICA_thickness <= radius2:
        PICA = True

    if CP_thickness + HT_424_thickness + ACC6_thickness <= radius2:
        CP = True

    if CP and PICA:
        if CP_total_weight > PICA_weight:
            print("choose PICA")
            return PICA_weight, PICA_thickness
        else:
            print("choose CP")
        return CP_total_weight, CP_thickness + HT_424_thickness + ACC6_thickness
    elif PICA:
        print("choose PICA")
        return PICA_weight, PICA_thickness
    elif CP:
        print("choose CP")
        return CP_total_weight, CP_thickness + HT_424_thickness + ACC6_thickness
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


def s1(x, d, r1, r2, a1, a2, dia):
    return (r1 - d) * np.sin(np.arccos(x/(r1 - d)))


def s2(x, d, r1, r2, a1, a2, dia):
    return np.tan(np.pi / 2 - a1) * x + (r1 - d) / np.sin(a1)


def s3(x, d, r1, r2, a1, a2, dia):
    x -= (dia / 2 - r2 + (r2 - d) * np.cos(np.pi / 2 - a1) -
          (r1 - d) / np.sin(a1)) / np.tan(np.pi / 2 - a1) + (r2 - d) * np.sin(np.pi / 2 - a1)
    return dia / 2 - r2 + (r2 - d) * np.cos(np.arcsin(x / (r2 - d)))


def simulate_entry_heating(mass, diameter, alt, lat, lon, speed, flight_path_angle, heading_angle, acc=1):

    drag = atmospheric_model.Entry_model.CapsuleDrag(diameter, 1.125, np.deg2rad(20), lat, lon, acc)

    aero_coefficient_setting = environment_setup.aerodynamic_coefficients.custom_aerodynamic_force_coefficients(
        force_coefficient_function=drag.drag_coefficient,
        reference_area=np.pi * (drag.diameter / 2) ** 2,
        independent_variable_names=[environment.AerodynamicCoefficientsIndependentVariables.mach_number_dependent,
                                    environment.AerodynamicCoefficientsIndependentVariables.altitude_dependent])

    termination_mach_setting = propagation_setup.propagator.dependent_variable_termination(
        dependent_variable_settings=propagation_setup.dependent_variable.mach_number("Capsule", "Uranus"),
        limit_value=2.0,
        use_as_lower_limit=True)
    dependent_variables_array = atmospheric_model.Entry_model.entry_sim(mass, aero_coefficient_setting, alt, lat, lon,
                                                                        speed, flight_path_angle, heading_angle,
                                                                        [termination_mach_setting], acc=1)
    gram = GRAM()
    gram.altitudes = dependent_variables_array[:, 1] / 1000
    gram.time = dependent_variables_array[:, 0]
    gram.lat = np.rad2deg(dependent_variables_array[:, 2])
    gram.long = np.rad2deg((dependent_variables_array[:, 3] + 2 * np.pi) % (2 * np.pi))
    gram.run()

    k = 1 / (np.asarray(gram.data.H2mass_pct) / 0.0395 + np.asarray(gram.data.Hemass_pct) / 0.0797)
    q_c = k * dependent_variables_array[:, 4] ** 3 * np.sqrt(np.asarray(gram.data.Density_kgm3) /
                                                             (np.pi * (diameter/2) ** 2))
    q_r = 9.7632379 ** (-40) * diameter ** (-0.17905) * np.asarray(gram.data.Density_kgm3) ** 1.763827469 * \
          dependent_variables_array[:, 4] ** 10.993852

    q_func = sp.interpolate.interp1d(dependent_variables_array[:, 0], q_c + q_r)
    h = sp.integrate.quad(lambda x: q_func(x), dependent_variables_array[0, 0],
                          dependent_variables_array[-1, 0])[0]


    print(h, max(q_c + q_r))
    drag.gamma = np.asarray(gram.data.SpecificHeatRatio)
    drag.mach = dependent_variables_array[:, 6]
    pressure = drag.p_0_stag() * np.asarray(gram.data.Pressure_Pa)
    return h, max(q_c + q_r), max(pressure), max(dependent_variables_array[:, 5])


def itterate_heatshield(mass, diameter, alt, lat, lon, speed, flight_path_angle,
                        heading_angle, acc=1, steps=5):
    hmass = 0
    for _ in range(0, steps):
        h, q, p, a = simulate_entry_heating(mass + hmass, diameter, alt, lat, lon, speed, flight_path_angle,
                                            heading_angle, acc)
        hmass, t = heatshield_sizing(diameter, h, speed, q)

    return hmass, p, a, t


if __name__ == "__main__":
    print(itterate_heatshield(400, 4.5, 3.03327727e+07, 5.45941114e-01, -2.33346601e-02, 2.65992642e+04,
                              -5.91036848e-01, -2.96367147e+00, 1, 3))
