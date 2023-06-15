import numpy as np

# constants:
k = 1.38 * 10 ** (-23)  # boltzmann constant
c = 300000000  # speed of light in m/s
earthRadius = 6371000.  # radius Earth in m
L_a = -0.5  # atmospheric attenuation in dB
AU = 149597870691  # one AU in m
d_EarthSC = 20.8  # max distance between SC and Earth in AU
Tnoisedown = 424  # Noise temperature in K
Tnoiseup = 763  # Noise temperature in K
TurnAroundRatio = 3599 / 3344

# antenna spacecraft:
d_antenna = 5  # antenna diameter in m
eta_antenna = 0.55
P = 60  # transmitting power in W
L_l = 0.9  # loss factor spacecraft
PointingAccuracy = 0.0572958  # pointing accuracy in deg
DR = 8000  # downlink data rate in bps
f = 32  # Downlink frequency in GHz
wavelengthdown = c / (f * 10 ** 9)

# antenna ground station:
P_gs = 800  # power ground station in W
d_gs = 70  # antenna diameter in m
L_r = 0.75  # loss factor ground station
uplinkDR = 25000  # uplink data rate in bps
f_gs = f * TurnAroundRatio  # uplink frequency in Ghz
wavelenghtup = c / (f_gs * 10 ** 9)


# calculate space loss:
def SpaceLoss(wavelength, d, AU):
    R = d * AU
    spaceloss = (wavelength / (4 * np.pi * R)) ** 2
    spaceloss = 10 * np.log10(spaceloss)
    return spaceloss


# calculate the gain:
def gain(D, wavelength, eta):
    G = np.pi ** 2 * D ** 2 / wavelength ** 2 * eta
    G = 10 * np.log10(G)
    return G


# calculate the half-power beamwidth:
def halfpowerbeamwidth(D, f):
    alpha = 21 / (D * f)
    return alpha


# calculate the antenna pointing loss:
# should be used for the transmitting and receiving antenna
def pointingLoss(alpha, pointingAccuracy):
    pointingLoss = -12. * (pointingAccuracy / alpha) ** 2
    return pointingLoss


def downlink(P, L_l, L_r, L_a, DR, Tnoise, k, d_ant=d_antenna, d_ground=d_gs, wavelength=wavelengthdown,
             eta_ant=eta_antenna, PointAccSC=PointingAccuracy, freq=f, d=d_EarthSC, AU=AU):
    P = 10 * np.log10(P)
    G_t = gain(d_ant, wavelength, eta_ant)
    G_r = gain(d_ground, wavelength, eta_ant)
    L_l = 10 * np.log10(L_l)
    L_r = 10 * np.log10(L_r)
    L_s = SpaceLoss(wavelength, d, AU)
    alphaSC = halfpowerbeamwidth(d_ant, freq)
    alphaGS = halfpowerbeamwidth(d_ground, freq)
    PointingAccuracyGS = 0.1 * alphaGS  # pointing accuracy in deg
    L_prSC = pointingLoss(alphaSC, PointAccSC)  # in [dB]
    L_prGS = pointingLoss(alphaGS, PointingAccuracyGS)  # in [dB]
    L_pr = L_prSC + L_prGS  # in [dB]
    DR = 10 * np.log10(1 / DR)  # [dB]
    Tnoise = 10 * np.log10(1 / Tnoise)  # [dB]
    k = 10 * np.log10(1 / k)
    EbN0 = P + G_t + G_r + L_l + L_r + L_s + L_pr + L_a + DR + Tnoise + k
    Eb = EbN0 - Tnoise - k
    Eb = 10 ** (Eb / 10)
    return EbN0


def uplink(f_gs, P_gs, L_l, L_r, L_a, uplinkDR, Tnoise, k, wavelength=wavelenghtup, eta_ant=eta_antenna, d_ground=d_gs,
           d_ant=d_antenna, d=d_EarthSC, AU=AU, PointAccSC=PointingAccuracy):
    P = 10 * np.log10(P_gs)  # in [dB]
    G_t = gain(d_ground, wavelength, eta_ant)  # in [dB]
    G_r = gain(d_ant, wavelength, eta_ant)  # in [dB]
    L_l = 10 * np.log10(L_l)  # in [dB]
    L_r = 10 * np.log10(L_r)  # in [dB]
    L_s = SpaceLoss(wavelength, d, AU)
    alphaSC = halfpowerbeamwidth(d_ant, f_gs)
    alphaGS = halfpowerbeamwidth(d_ground, f_gs)
    PointingAccuracyGS = 0.1 * alphaGS  # pointing accuracy in deg
    L_prSC = pointingLoss(alphaSC, 0.1 * PointAccSC)  # in [dB]
    L_prGS = pointingLoss(alphaGS, PointingAccuracyGS)  # in [dB]
    L_pr = L_prSC + L_prGS  # in [dB]
    DR = 10 * np.log10(1 / uplinkDR)  # [dB]
    Tnoise = 10 * np.log10(1 / Tnoise)  # [dB]
    k = 10 * np.log10(1 / k)
    EbN0 = P + G_t + G_r + L_l + L_r + L_s + L_pr + L_a + DR + Tnoise + k
    Eb = EbN0 - Tnoise - k
    Eb = 10 ** (Eb / 10)
    return EbN0

def total_cost(total_mass):
    return 237.93 * 1000 * total_mass * 1.43 * 0.951


if __name__ == "__main__":
    ...
    # print('**** DOWNLINK ****')
    # print('P:', downlink(P, L_l, L_r, L_a, DR, Tnoisedown, k)[0])
    # print('G_t:', downlink(P, L_l, L_r, L_a, DR, Tnoisedown, k)[1])
    # print('G_r:', downlink(P, L_l, L_r, L_a, DR, Tnoisedown, k)[2])
    # print('L_l:', downlink(P, L_l, L_r, L_a, DR, Tnoisedown, k)[3])
    # print('L_r:', downlink(P, L_l, L_r, L_a, DR, Tnoisedown, k)[4])
    # print('L_s:', downlink(P, L_l, L_r, L_a, DR, Tnoisedown, k)[5])
    # print('L_pr:', downlink(P, L_l, L_r, L_a, DR, Tnoisedown, k)[6])
    # print('DR:', downlink(P, L_l, L_r, L_a, DR, Tnoisedown, k)[7])
    # print('Tnoise:', downlink(P, L_l, L_r, L_a, DR, Tnoisedown, k)[8])
    # print('k:', downlink(P, L_l, L_r, L_a, DR, Tnoisedown, k)[9])
    # print('EbN0:', downlink(P, L_l, L_r, L_a, DR, Tnoisedown, k)[10])
    # print('Eb:', downlink(P, L_l, L_r, L_a, DR, Tnoisedown, k)[11])
    # print()
    # print('**** UPLINK ****')
    # print('P:', uplink(f_gs, P_gs, L_l, L_r, L_a,uplinkDR, Tnoiseup, k)[0])
    # print('G_t:', uplink(f_gs, P_gs, L_l, L_r, L_a, uplinkDR, Tnoiseup, k)[1])
    # print('G_r:', uplink(f_gs, P_gs, L_l, L_r, L_a, uplinkDR, Tnoiseup, k)[2])
    # print('L_l:', uplink(f_gs, P_gs, L_l, L_r, L_a, uplinkDR, Tnoiseup, k)[3])
    # print('L_r:', uplink(f_gs, P_gs, L_l, L_r, L_a, uplinkDR, Tnoiseup, k)[4])
    # print('L_s:', uplink(f_gs, P_gs, L_l, L_r, L_a, uplinkDR, Tnoiseup, k)[5])
    # print('L_pr:', uplink(f_gs, P_gs, L_l, L_r, L_a, uplinkDR, Tnoiseup, k)[6])
    # print('DR:', uplink(f_gs, P_gs, L_l, L_r, L_a, uplinkDR, Tnoiseup, k)[7])
    # print('Tnoise:', uplink(f_gs, P_gs, L_l, L_r, L_a, uplinkDR, Tnoiseup, k)[8])
    # print('k:', uplink(f_gs, P_gs, L_l, L_r, L_a, uplinkDR, Tnoiseup, k)[9])
    # print('EbN0:', uplink(f_gs, P_gs, L_l, L_r, L_a, uplinkDR, Tnoiseup, k)[10])
    # print('Eb:', uplink(f_gs, P_gs, L_l, L_r, L_a, uplinkDR, Tnoiseup, k)[11])
