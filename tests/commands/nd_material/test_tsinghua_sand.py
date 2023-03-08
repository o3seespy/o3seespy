import o3seespy as o3  # for testing only


def test_cyc_liq_cp():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.nd_material.CycLiqCP(osi, g0=1.0, kappa=1.0, h=1.0, mfc=1.0, dre1=1.0, mdc=1.0, dre2=1.0, gamma_dr=1.0,
                            alpha=1.0, d_ir=1.0, e_init=1.0, den=1.0)


def test_cyc_liq_cpsp():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.nd_material.CycLiqCPSP(osi, g0=200, kappa=0.008, h=1.8, big_m=1.25, dre1=0.6, dre2=30, gamma_dr=0.05,
                              alpha=20, d_ir=1.4, lambda_c=0.019, ksi=0.7, e_0=0.934, n_p=1.1, n_d=7.8, e_init=0.87,
                              den=1.6)
