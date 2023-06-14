import numpy as np

# https://nuke.fas.org/space/gphs.pdf
# https://openlearninglibrary.mit.edu/courses/course-v1:MITx+22.011x+3T2018/courseware/df988ecc766d475e8010e172cd81426b/e6493213bbd74be6a0f229a425e7af37/?activate_block_id=block-v1%3AMITx%2B22.011x%2B3T2018%2Btype%40sequential%2Bblock%40e6493213bbd74be6a0f229a425e7af37
# https://ebookcentral-proquest-com.tudelft.idm.oclc.org/lib/delft/reader.action?docID=693314
# Cost: https://inldigitallibrary.inl.gov/sites/sti/sti/7267852.pdf

P_start = 300  # Begin of life power one GPHS-RTG in W
tau1 = 87.7  # half life fuel in years
mass_RTG = 55.9  # mass of one GPHS-RTG in kg
costRTG1 = 145699633.36  # cost of one GPHS-RTG in FY$2022, This is the highest value. It could be around 130 million as well
missiontime = 20  # in years
#P_req = 500  # total power required for all subsystems in W


# Power at the end of life of one RTG calculations:
def powerdecay(t=missiontime, P_0=P_start, tau=tau1):
    P = P_0 * np.exp(-np.log(2) / (tau * 365 * 24 * 60 * 60) * (t * 365 * 24 * 60 * 60))
    return P


# Calculate the number of RTGs needed:
# 0.85 is path efficiency coming from SMAD
def numberRTG(P_req, t=missiontime, P_0=P_start, tau=tau1):
    M = P_req / (powerdecay(t, P_0, tau) * 0.85)
    N = np.ceil(M)
    return M, N


# Calculate mass of RTGs:
def massRTG(P_req, missiontime, m_rtg=mass_RTG, P_0=P_start, tau=tau1):
    m_RTGs = m_rtg * numberRTG(P_req, missiontime, P_0, tau)[1]
    return m_RTGs

def costRTG(P_req, missiontime, cost_rtg=costRTG1, P_0=P_start, tau=tau1):
    cost_RTG = cost_rtg * numberRTG(P_req, missiontime, P_0, tau)[1]
    return cost_RTG


if __name__ == "__main__":
    # https://nuke.fas.org/space/gphs.pdf
    # https://openlearninglibrary.mit.edu/courses/course-v1:MITx+22.011x+3T2018/courseware/df988ecc766d475e8010e172cd81426b/e6493213bbd74be6a0f229a425e7af37/?activate_block_id=block-v1%3AMITx%2B22.011x%2B3T2018%2Btype%40sequential%2Bblock%40e6493213bbd74be6a0f229a425e7af37
    # https://ebookcentral-proquest-com.tudelft.idm.oclc.org/lib/delft/reader.action?docID=693314
    # Cost: https://inldigitallibrary.inl.gov/sites/sti/sti/7267852.pdf



    print("Power at end of life of one RTG=", powerdecay(P_0, tau1, missiontime), "W")
    print("number of RTGS unrounded =", numberRTG(P_req)[0])
    print("number of RTGs =", numberRTG(P_req)[1])
    print("Total mass of RTGs =", massRTG(mass_RTG), "kg")
    print("Total cost of RTGs =", costRTG(costRTG1), "$")
