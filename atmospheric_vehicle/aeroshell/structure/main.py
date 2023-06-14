"""
Script to size the backshell structure. It implements calculations for pressure loads,
buckling loads and aerodynamic loads on the entire structure, plus the peak axial loads on
the top surface of the backshell due to parachute deployment.
"""
import numpy as np
import matplotlib.pyplot as plt
import atmospheric_vehicle.aeroshell.heatshield as heat_sh

def angle_cone(r_big, r_small, height):
    """
    Function that determined the angle between the base and the hypotenuse of the cone
    :param r_big: Bottom radius of the cone
    :param r_small: Top radius of the cone
    :param height: height of the cone
    :return: angle between hypotenuse and base
    """
    return np.arctan2(r_big - r_small, height)

def t_pressure(p, r_big, a, sigma_y):
    """
    Function to determine the thickness required to withstand pressure loads during entry. Formulas were derived from
    the following website: https://123sanat.com/d/catalogue/ASME-VIII-_8_-div.1-2019.pdf
    :param p: External pressure
    :param r_big: Top radius of cone
    :param a: angle between hypotenuse and base
    :param sigma_y: Yield stress of selected material
    :return: Required minimum thickness
    """
    return p * r_big * 2 / (2 * np.cos(a) * (sigma_y - 0.6 * p))

def t_hoop(p, r_big, sigma_y):
    """
    Function to calculate the required thickness to withstand pressure loads during entry. Implemented the classical
    formula to compute hoop stress areound a cylindrical tank. Use the function before, in case no reference can be
    found, then this function should be placed in the code at line 54-55
    :param p: External pressure
    :param r_big: Bottom radius of cone
    :param sigma_y: yield strength
    :return: Required minimum thickness
    """
    return p * r_big / sigma_y

def volume_truncated(r_big, r_small, h):
    """
    Function to calculate the volume of a truncated cone (shape of the backshell elements)
    :param r_big: Bottom radius of cone
    :param r_small: Top radius of cone
    :param h: Height of cone
    :return: Volume of truncated cone
    """
    return (1 / 3) * np.pi * h * (r_small ** 2 + r_small * r_big + r_big ** 2)


def buckling(E, l, I, p_load):
    """
    Function to calculate the buckling limit load
    :param E: Young's Modulus
    :param l: Length of plate
    :param I: Moment of inertia plate
    :param sigma_y: Yield stress
    :return: bool if buckling criteria is respected
    """
    return p_load < np.pi**2 * E * I / (l**2)


def aeroshell_geometry(peak_load, p_load, sigma_y, E, taper, r_thermal, h_glider, h_parachute, r_parachute, delta_t, alpha):
    """
    Function that determines all geometrical properties and mass of the backshell. Multiple things are computed,
    during integration, it should be decided what is the most useful return, if one is actually needed. As of now the
    function prints the main geometrical characteristics and mass of the backshell. The backshell is assumed to be
    formed by two truncated cones, a bottom one that has the bottom radius equal to the thermal shield radius and a
    top one that has the top radius equal to the parachute radius (in folded configuration). Taper ratio is selected
    manually to determine the transition from one truncated cone to the other.
    :param peak_load: Peak load experienced during entry (from decelerator analysis)
    :param p_load: Pressure loads experienced during entry
    :param sigma_y: Yield stress of the chosen material
    :param taper: Taper ratio of the different elements of the backshell
    :param r_thermal: Thermal shield radius
    :param h_glider: Height of the glider in the folded configuration
    :param h_parachute: Height of the parachute in the folded configuration
    :param r_parachute: Radius of the parachute in the folded configuration
    :return: Backshell mass !!!! ADD MASS OF THERMAL SHIELD SUPPORTING STRUCTURE !!!!! (Use volume function from heatshild code)
    """

    # Calculate the big and small radius of the top truncated cone. Limit case is the parachute radius
    r_top_small = np.sqrt(max(peak_load / sigma_y, np.pi * r_parachute**2) / np.pi)
    r_top_big = r_top_small / (1 - taper)

    # Calculate the big and small radius of the bottom truncated cone. Limit case is the radius of the thermal shield
    r_bottom_big = r_thermal
    r_bottom_small = r_bottom_big * (1 - taper)

    # Calculate the angle at the base of the truncated cones
    a_top = angle_cone(r_top_big, r_top_small, h_parachute)
    a_bottom = angle_cone(r_bottom_big, r_bottom_small, h_glider)

    # Calculate thickness based on pressure loads. Use personalised formula
    t_top = t_pressure(p_load, r_top_big, a_top, sigma_y)
    t_bottom = t_pressure(p_load, r_bottom_big, a_bottom, sigma_y)

    # Calculate thickness based on pressure loads. Use traditional formula for thin walled cylinders
    # t_top = t_hoop(p_load, r_top_big, sigma_y)
    # t_bottom = t_hoop(p_load, r_bottom_big, sigma_y)

    # Check buckling
    I_top = 1/12 * t_top * (h_parachute / np.cos(a_top))**3
    buck_top = buckling(E, I_top, h_parachute / np.cos(a_top), p_load)
    I_bottom = 1/12 * t_bottom * h_glider / np.cos(a_bottom)
    buck_bottom = buckling(E, I_bottom, h_glider / np.cos(a_bottom), p_load)

    if buck_bottom and buck_top:
        # Calculate the volume of the thin walled structure by subtraction
        volume_top = volume_truncated(r_top_big, r_top_small, h_parachute) - volume_truncated(r_top_big - t_top, r_top_small - t_top, h_parachute - 2 * t_top)
        volume_bottom = volume_truncated(r_bottom_big, r_bottom_small, h_glider) - volume_truncated(r_bottom_big - t_bottom, r_bottom_small - t_bottom, h_glider - 2 * t_bottom)

        # Calculate backshell mass
        mass_backshell = (volume_top + volume_bottom) * rho_backshell
    else:
        print(f'Buckling requirements is not satisfied')
    return mass_backshell


def thermal_loads(alpha, delta_T, E, sigma_y):
    return sigma_y > alpha * delta_T * E

if __name__ == "__main__":
    # Loads values
    load_peak_para = 100
    load_peak_entry = 300
    p_load = 2 * 10**6
    delta_T = 100

    # Mass Budget
    m_thermal_para = 330 + 30
    m_glider = 300


    # Size Constraints
    h_glider = 1.8
    h_parachute = 0.5
    r_thermal = 4.5
    r_parachute = 0.5
    _radius1 = 1.125
    _radius2 = 0.126

    # Capsule Properties
    # cd_capsule = 1.3
    # -- CD research -- #
    # Apollo --> cd = 0.62 (Based on ballistic coefficient, surface area and mass)
    # Pioneer Venus --> cd = 0.63 (Based on ballistic coefficient, geometry and mass)
    # Huygens Titan --> cd = 1.34 (Based on ballistic coefficient, geometry and mass)(Made an average of the ballistic coeff.)
    # Galileo Jupiter -> cd = 0.74 (Based on ballistic coefficient, surface area and mass)
    # theta_capsule = 60 * np.pi / 180
    taper_ratio = 0.2

    # Material Properties
    # -- Alluminum honeycomb with graphite eopxy -- Mars rover (https://spaceflight101.com/msl/msl-aeroshell-and-heat-shield/)
    rho_backshell = 49.7  # https://journals.sagepub.com/doi/pdf/10.1177/0021998313499949 (Q2 selected because strongest)
    sigma_y_backshell = 450 * 10 ** 6  # https://journals.sagepub.com/doi/pdf/10.1177/0021998313499949 (Q2 selected because strongest)
    E = 130 * 10**9  # https://www.azom.com/article.aspx?ArticleID=6618
    alpha = 23.6 * 10**-6

    # -- Alluminum Honeycomb --
    # rho_backshell = 54.4  # https://pdf.sciencedirectassets.com/271493/1-s2.0-S0263823100X00456/1-s2.0-S0263823199000269/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEPf%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJHMEUCIEYcpzxEXVmsIPQ7dXYvc%2F7KIwF%2Bvm5xmju2mp16E5B%2FAiEAoLFSr%2FQD8DE61bpasIYYN0CevERM48NxnWNw2HogFd8qsgUIQBAFGgwwNTkwMDM1NDY4NjUiDEg%2FA3swjHWYTtrDtiqPBVPXitudSZ5okfgxfWqURmNRI5QDeaJQIbUUnHBTXqKEUxQYB7TvbGceTSw1UWYH4rmzkgOs8YD56a5EMWyLQcvjck6L42LfuNEnLzlVmrwaY46SEpmtSi5W3D3gYZ7BYA0vKBn%2FfjiqR197rzZYKtLFuT8EifjvL%2FY6bW3rp%2FFkGNgGDryYeUSfosFpGG5LNEbk7tW4QTVWtn3RUkqSZr%2Bsp8oZKvPfc%2BjVoBeWP77d9j6jONRrIbdn%2F6UEeV0pzaNyLGPTwu9B3tVKHOkCdNxR5c5dWgFDmSrJSRNQpR69%2BpQAbVOTMEredyBZyhFLnKQ9pj%2FVjLEu3uL2AeQTW0wVr8GTqdXfjjmF440AqHckWg2NO0GbHAbrCxIXauIIAae3HOvFzGS2gY%2FT%2FwPVIiDn2WGIg4VD6urxcPGjB3jqOGuX%2FLvpb3MIJZ8xR67vAfAHQ3Yo08E7RE2eAkroglZPtbn7Uqnx5IC96d%2BEsvFPJ2hWjWQM0Yy6HQhT8wPMrFpD%2FM93cDzeIwaGbW5AwAG4D9fCbOXzqqdS5T%2BwAi9mE%2F1CfFmtGJWNiDTehe9d9DYw2wt17c8GBG%2FAPB4nsA355RWBbIpHNc%2BBgkS4OszazUfBvUITlybBWdIetY6%2BmeckwcfMl%2FEpiqemb4mY%2BIBfY8ZU0wA%2BPDxGxw3FtA6GXxLvxMca0nkjcHq85d7vB4SAR7XBe41wQWz329ss7AciFEX8ZZ78vTeXUpEtwEOAZxONCYfUB82pEq7lFfQ%2FI%2F5qOweX4cbVfn%2B2WA10rUwA7A74un3pOvquXH%2Fop88mTfUeZ5aqRnv9WjWDCNj53oTEaENu6JRXGQzgyKPvxmOVRh4CZqy96WVf5xGaXoEwyq77owY6sQGhBiWUaPlV51caDzPYCTJrqJ%2Fez6twk6Zg05oosFK7a8VlmXhnARUtOITHZ3xl2d0JZFlcX9waYOey1XoKz08tBQL9ClOv8c1lto7aXpxn6zRX5KYr1RZR%2BeqEASArmRdTRpVWXsI3GOYF0OGTGqDcP3rJSeHN6HChK%2F78d2ZHfmvwCGmCkfluewj%2Ftu5fhIqv1lwBgQmrg7UrpH3J5XXPJ0O77sqtwV1x2ii9mBnXd2s%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20230606T081127Z&X-Amz-SignedHeaders=host&X-Amz-Expires=299&X-Amz-Credential=ASIAQ3PHCVTYQXRRMFP4%2F20230606%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=bf7dc076dbd173ecbcd96ffa72c3ccf458b36f2e87be3cad24c9d2e26d0566e0&hash=381e99d39fab924be7c214f45491d50c10cf7169b2124ee73b7fc356ba479251&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0263823199000269&tid=spdf-29eaeca5-92ad-4c84-a236-8e8679c5e771&sid=38c3df1f57d6e84e14099c90f14e886d3390gxrqb&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=0b0a520b560453035557&rr=7d2f2f459d5c0bb3&cc=nl
    # sigma_y_backshell = 190 * 10**6  # https://pdf.sciencedirectassets.com/271493/1-s2.0-S0263823100X00456/1-s2.0-S0263823199000269/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEPf%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJHMEUCIEYcpzxEXVmsIPQ7dXYvc%2F7KIwF%2Bvm5xmju2mp16E5B%2FAiEAoLFSr%2FQD8DE61bpasIYYN0CevERM48NxnWNw2HogFd8qsgUIQBAFGgwwNTkwMDM1NDY4NjUiDEg%2FA3swjHWYTtrDtiqPBVPXitudSZ5okfgxfWqURmNRI5QDeaJQIbUUnHBTXqKEUxQYB7TvbGceTSw1UWYH4rmzkgOs8YD56a5EMWyLQcvjck6L42LfuNEnLzlVmrwaY46SEpmtSi5W3D3gYZ7BYA0vKBn%2FfjiqR197rzZYKtLFuT8EifjvL%2FY6bW3rp%2FFkGNgGDryYeUSfosFpGG5LNEbk7tW4QTVWtn3RUkqSZr%2Bsp8oZKvPfc%2BjVoBeWP77d9j6jONRrIbdn%2F6UEeV0pzaNyLGPTwu9B3tVKHOkCdNxR5c5dWgFDmSrJSRNQpR69%2BpQAbVOTMEredyBZyhFLnKQ9pj%2FVjLEu3uL2AeQTW0wVr8GTqdXfjjmF440AqHckWg2NO0GbHAbrCxIXauIIAae3HOvFzGS2gY%2FT%2FwPVIiDn2WGIg4VD6urxcPGjB3jqOGuX%2FLvpb3MIJZ8xR67vAfAHQ3Yo08E7RE2eAkroglZPtbn7Uqnx5IC96d%2BEsvFPJ2hWjWQM0Yy6HQhT8wPMrFpD%2FM93cDzeIwaGbW5AwAG4D9fCbOXzqqdS5T%2BwAi9mE%2F1CfFmtGJWNiDTehe9d9DYw2wt17c8GBG%2FAPB4nsA355RWBbIpHNc%2BBgkS4OszazUfBvUITlybBWdIetY6%2BmeckwcfMl%2FEpiqemb4mY%2BIBfY8ZU0wA%2BPDxGxw3FtA6GXxLvxMca0nkjcHq85d7vB4SAR7XBe41wQWz329ss7AciFEX8ZZ78vTeXUpEtwEOAZxONCYfUB82pEq7lFfQ%2FI%2F5qOweX4cbVfn%2B2WA10rUwA7A74un3pOvquXH%2Fop88mTfUeZ5aqRnv9WjWDCNj53oTEaENu6JRXGQzgyKPvxmOVRh4CZqy96WVf5xGaXoEwyq77owY6sQGhBiWUaPlV51caDzPYCTJrqJ%2Fez6twk6Zg05oosFK7a8VlmXhnARUtOITHZ3xl2d0JZFlcX9waYOey1XoKz08tBQL9ClOv8c1lto7aXpxn6zRX5KYr1RZR%2BeqEASArmRdTRpVWXsI3GOYF0OGTGqDcP3rJSeHN6HChK%2F78d2ZHfmvwCGmCkfluewj%2Ftu5fhIqv1lwBgQmrg7UrpH3J5XXPJ0O77sqtwV1x2ii9mBnXd2s%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20230606T081127Z&X-Amz-SignedHeaders=host&X-Amz-Expires=299&X-Amz-Credential=ASIAQ3PHCVTYQXRRMFP4%2F20230606%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=bf7dc076dbd173ecbcd96ffa72c3ccf458b36f2e87be3cad24c9d2e26d0566e0&hash=381e99d39fab924be7c214f45491d50c10cf7169b2124ee73b7fc356ba479251&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0263823199000269&tid=spdf-29eaeca5-92ad-4c84-a236-8e8679c5e771&sid=38c3df1f57d6e84e14099c90f14e886d3390gxrqb&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=0b0a520b560453035557&rr=7d2f2f459d5c0bb3&cc=nl

    aeroshell_geometry(load_peak_entry * (m_glider + m_thermal_para), p_load, sigma_y_backshell, E, taper_ratio, r_thermal, h_glider, h_parachute, r_parachute, delta_T, alpha)