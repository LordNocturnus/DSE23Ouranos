from atmospheric_vehicle.aeroshell.structure import aeroshell_geometry

# Loads values
load_peak = 100 * 9.81
p_load = 1 * 10 ** 6
delta_T = 100

# Mass Budget
m_thermal = 1100
m_glider = 114

# Size Constraints
h_glider = 1.8
h_parachute = 0.5
r_thermal = 4.5
r_parachute = 0.25

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
E = 130 * 10 ** 9  # https://www.azom.com/article.aspx?ArticleID=6618
alpha = 23.6 * 10 ** -6

# -- Alluminum Honeycomb --
# rho_backshell = 54.4  # https://pdf.sciencedirectassets.com/271493/1-s2.0-S0263823100X00456/1-s2.0-S0263823199000269/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEPf%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJHMEUCIEYcpzxEXVmsIPQ7dXYvc%2F7KIwF%2Bvm5xmju2mp16E5B%2FAiEAoLFSr%2FQD8DE61bpasIYYN0CevERM48NxnWNw2HogFd8qsgUIQBAFGgwwNTkwMDM1NDY4NjUiDEg%2FA3swjHWYTtrDtiqPBVPXitudSZ5okfgxfWqURmNRI5QDeaJQIbUUnHBTXqKEUxQYB7TvbGceTSw1UWYH4rmzkgOs8YD56a5EMWyLQcvjck6L42LfuNEnLzlVmrwaY46SEpmtSi5W3D3gYZ7BYA0vKBn%2FfjiqR197rzZYKtLFuT8EifjvL%2FY6bW3rp%2FFkGNgGDryYeUSfosFpGG5LNEbk7tW4QTVWtn3RUkqSZr%2Bsp8oZKvPfc%2BjVoBeWP77d9j6jONRrIbdn%2F6UEeV0pzaNyLGPTwu9B3tVKHOkCdNxR5c5dWgFDmSrJSRNQpR69%2BpQAbVOTMEredyBZyhFLnKQ9pj%2FVjLEu3uL2AeQTW0wVr8GTqdXfjjmF440AqHckWg2NO0GbHAbrCxIXauIIAae3HOvFzGS2gY%2FT%2FwPVIiDn2WGIg4VD6urxcPGjB3jqOGuX%2FLvpb3MIJZ8xR67vAfAHQ3Yo08E7RE2eAkroglZPtbn7Uqnx5IC96d%2BEsvFPJ2hWjWQM0Yy6HQhT8wPMrFpD%2FM93cDzeIwaGbW5AwAG4D9fCbOXzqqdS5T%2BwAi9mE%2F1CfFmtGJWNiDTehe9d9DYw2wt17c8GBG%2FAPB4nsA355RWBbIpHNc%2BBgkS4OszazUfBvUITlybBWdIetY6%2BmeckwcfMl%2FEpiqemb4mY%2BIBfY8ZU0wA%2BPDxGxw3FtA6GXxLvxMca0nkjcHq85d7vB4SAR7XBe41wQWz329ss7AciFEX8ZZ78vTeXUpEtwEOAZxONCYfUB82pEq7lFfQ%2FI%2F5qOweX4cbVfn%2B2WA10rUwA7A74un3pOvquXH%2Fop88mTfUeZ5aqRnv9WjWDCNj53oTEaENu6JRXGQzgyKPvxmOVRh4CZqy96WVf5xGaXoEwyq77owY6sQGhBiWUaPlV51caDzPYCTJrqJ%2Fez6twk6Zg05oosFK7a8VlmXhnARUtOITHZ3xl2d0JZFlcX9waYOey1XoKz08tBQL9ClOv8c1lto7aXpxn6zRX5KYr1RZR%2BeqEASArmRdTRpVWXsI3GOYF0OGTGqDcP3rJSeHN6HChK%2F78d2ZHfmvwCGmCkfluewj%2Ftu5fhIqv1lwBgQmrg7UrpH3J5XXPJ0O77sqtwV1x2ii9mBnXd2s%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20230606T081127Z&X-Amz-SignedHeaders=host&X-Amz-Expires=299&X-Amz-Credential=ASIAQ3PHCVTYQXRRMFP4%2F20230606%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=bf7dc076dbd173ecbcd96ffa72c3ccf458b36f2e87be3cad24c9d2e26d0566e0&hash=381e99d39fab924be7c214f45491d50c10cf7169b2124ee73b7fc356ba479251&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0263823199000269&tid=spdf-29eaeca5-92ad-4c84-a236-8e8679c5e771&sid=38c3df1f57d6e84e14099c90f14e886d3390gxrqb&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=0b0a520b560453035557&rr=7d2f2f459d5c0bb3&cc=nl
# sigma_y_backshell = 190 * 10**6  # https://pdf.sciencedirectassets.com/271493/1-s2.0-S0263823100X00456/1-s2.0-S0263823199000269/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEPf%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJHMEUCIEYcpzxEXVmsIPQ7dXYvc%2F7KIwF%2Bvm5xmju2mp16E5B%2FAiEAoLFSr%2FQD8DE61bpasIYYN0CevERM48NxnWNw2HogFd8qsgUIQBAFGgwwNTkwMDM1NDY4NjUiDEg%2FA3swjHWYTtrDtiqPBVPXitudSZ5okfgxfWqURmNRI5QDeaJQIbUUnHBTXqKEUxQYB7TvbGceTSw1UWYH4rmzkgOs8YD56a5EMWyLQcvjck6L42LfuNEnLzlVmrwaY46SEpmtSi5W3D3gYZ7BYA0vKBn%2FfjiqR197rzZYKtLFuT8EifjvL%2FY6bW3rp%2FFkGNgGDryYeUSfosFpGG5LNEbk7tW4QTVWtn3RUkqSZr%2Bsp8oZKvPfc%2BjVoBeWP77d9j6jONRrIbdn%2F6UEeV0pzaNyLGPTwu9B3tVKHOkCdNxR5c5dWgFDmSrJSRNQpR69%2BpQAbVOTMEredyBZyhFLnKQ9pj%2FVjLEu3uL2AeQTW0wVr8GTqdXfjjmF440AqHckWg2NO0GbHAbrCxIXauIIAae3HOvFzGS2gY%2FT%2FwPVIiDn2WGIg4VD6urxcPGjB3jqOGuX%2FLvpb3MIJZ8xR67vAfAHQ3Yo08E7RE2eAkroglZPtbn7Uqnx5IC96d%2BEsvFPJ2hWjWQM0Yy6HQhT8wPMrFpD%2FM93cDzeIwaGbW5AwAG4D9fCbOXzqqdS5T%2BwAi9mE%2F1CfFmtGJWNiDTehe9d9DYw2wt17c8GBG%2FAPB4nsA355RWBbIpHNc%2BBgkS4OszazUfBvUITlybBWdIetY6%2BmeckwcfMl%2FEpiqemb4mY%2BIBfY8ZU0wA%2BPDxGxw3FtA6GXxLvxMca0nkjcHq85d7vB4SAR7XBe41wQWz329ss7AciFEX8ZZ78vTeXUpEtwEOAZxONCYfUB82pEq7lFfQ%2FI%2F5qOweX4cbVfn%2B2WA10rUwA7A74un3pOvquXH%2Fop88mTfUeZ5aqRnv9WjWDCNj53oTEaENu6JRXGQzgyKPvxmOVRh4CZqy96WVf5xGaXoEwyq77owY6sQGhBiWUaPlV51caDzPYCTJrqJ%2Fez6twk6Zg05oosFK7a8VlmXhnARUtOITHZ3xl2d0JZFlcX9waYOey1XoKz08tBQL9ClOv8c1lto7aXpxn6zRX5KYr1RZR%2BeqEASArmRdTRpVWXsI3GOYF0OGTGqDcP3rJSeHN6HChK%2F78d2ZHfmvwCGmCkfluewj%2Ftu5fhIqv1lwBgQmrg7UrpH3J5XXPJ0O77sqtwV1x2ii9mBnXd2s%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20230606T081127Z&X-Amz-SignedHeaders=host&X-Amz-Expires=299&X-Amz-Credential=ASIAQ3PHCVTYQXRRMFP4%2F20230606%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=bf7dc076dbd173ecbcd96ffa72c3ccf458b36f2e87be3cad24c9d2e26d0566e0&hash=381e99d39fab924be7c214f45491d50c10cf7169b2124ee73b7fc356ba479251&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0263823199000269&tid=spdf-29eaeca5-92ad-4c84-a236-8e8679c5e771&sid=38c3df1f57d6e84e14099c90f14e886d3390gxrqb&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=0b0a520b560453035557&rr=7d2f2f459d5c0bb3&cc=nl


class Aeroshell:
    def __init__(self):
        self.r_top_small, self.r_top_big, self.r_bottom_small, self.r_bottom_big, self.mass, self.t_top, self.t_bot, \
            self.rho = aeroshell_geometry(load_peak * (m_glider + m_thermal), p_load, sigma_y_backshell, E, taper_ratio,
                                          r_thermal, h_glider, h_parachute, r_parachute, delta_T, alpha)


if __name__ == "__main__":
    print("Hello World")