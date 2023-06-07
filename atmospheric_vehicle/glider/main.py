dod = 1
n_bat = 0.9
n_cab = 0.9
spec_energy = 1000

def tot_p_req(dh_p_req, ttc_p_req,
        adcs_p_req, pl_p_req)
    return dh_p_req+ttc_p_req+adcs_p_req+pl_p_req

def m_bat(dh_p_req, ttc_p_req, adcs_p_req,
        pl_p_req, dod, n_bat spec_energy)
    tot_p_req() = tot_p_req(dh_p_req, ttc_p_req, adcs_p_req, pl_p_req)
    return tot_p_req / dod / n_bat / n_cab / spec_energy


if __name__ == "__main__":
    print("Hello World")