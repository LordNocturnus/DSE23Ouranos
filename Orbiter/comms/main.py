import numpy as np

# constants:
c = 300000000  # speed of light in m/s
earthRadius = 6371000.  # radius Earth in m
AU = 149597870691  # one AU in m
d_EarthSC = 20.8  # max distance between SC and Earth in AU
TurnAroundRatio = 3599 / 3344

# antenna spacecraft:
d_antenna = 5  # antenna diameter in m
eta_antenna = 0.55
PointingAccuracy = 0.0572958  # pointing accuracy in deg
f = 32  # Downlink frequency in GHz
wavelengthdown = c / (f * 10 ** 9)

# antenna ground station:
d_gs = 70  # antenna diameter in m
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


def downlink(P, L_l, L_r, L_a, DR, Tnoise, k):
    P = 10 * np.log10(P)
    G_t = gain(d_antenna, wavelengthdown, eta_antenna)
    G_r = gain(d_gs, wavelengthdown, eta_antenna)
    L_l = 10 * np.log10(L_l)
    L_r = 10 * np.log10(L_r)
    L_s = SpaceLoss(wavelengthdown, d_EarthSC, AU)
    alphaSC = halfpowerbeamwidth(d_antenna, f)
    alphaGS = halfpowerbeamwidth(d_gs, f)
    PointingAccuracyGS = 0.1 * alphaGS  # pointing accuracy in deg
    L_prSC = pointingLoss(alphaSC, PointingAccuracy) # in [dB]
    L_prGS = pointingLoss(alphaGS, PointingAccuracyGS) # in [dB]
    L_pr = L_prSC + L_prGS # in [dB]
    DR = 10 * np.log10(1 / DR)  # [dB]
    Tnoise = 10 * np.log10(1 / Tnoise)  # [dB]
    k = 10 * np.log10(1 / k)
    EbN0 = P + G_t + G_r + L_l + L_r + L_s + L_pr + L_a + DR + Tnoise + k
    Eb = EbN0 - Tnoise - k
    Eb = 10 ** (Eb / 10)
    return EbN0

def uplink(f_gs, P_gs, L_l, L_r, L_a, uplinkDR, Tnoise, k):
    P = 10 * np.log10(P_gs)  # in [dB]
    G_t = gain(d_gs, wavelenghtup, eta_antenna)  # in [dB]
    G_r = gain(d_antenna, wavelenghtup, eta_antenna)  # in [dB]
    L_l = 10 * np.log10(L_l)  # in [dB]
    L_r = 10 * np.log10(L_r)  # in [dB]
    L_s = SpaceLoss(wavelenghtup, d_EarthSC, AU)
    alphaSC = halfpowerbeamwidth(d_antenna, f_gs)
    alphaGS = halfpowerbeamwidth(d_gs, f_gs)
    PointingAccuracyGS = 0.1 * alphaGS  # pointing accuracy in deg
    L_prSC = pointingLoss(alphaSC, 0.1 * PointingAccuracy)  # in [dB]
    L_prGS = pointingLoss(alphaGS, PointingAccuracyGS)  # in [dB]
    L_pr = L_prSC + L_prGS  # in [dB]
    DR = 10 * np.log10(1 / uplinkDR)  # [dB]
    Tnoise = 10 * np.log10(1 / Tnoise)  # [dB]
    k = 10 * np.log10(1 / k)
    EbN0 = P + G_t + G_r + L_l + L_r + L_s + L_pr + L_a + DR + Tnoise + k
    Eb = EbN0 - Tnoise - k
    Eb = 10 ** (Eb / 10)
    return EbN0


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
