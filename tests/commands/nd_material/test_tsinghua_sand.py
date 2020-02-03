import o3seespy as o3  # for testing only


def test_cyc_liq_cp():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.nd_material.CycLiqCP(osi, g0=1.0, kappa=1.0, h=1.0, mfc=1.0, dre1=1.0, mdc=1.0, dre2=1.0, rdr=1.0, alpha=1.0, dir=1.0, ein=1.0, rho=1.0)


def test_cyc_liq_cpsp():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.nd_material.CycLiqCPSP(osi, g0=1.0, kappa=1.0, h=1.0, big_m=1.0, dre1=1.0, dre2=1.0, rdr=1.0, alpha=1.0, dir=1.0, lambdac=1.0, ksi=1.0, e0=1.0, np=1.0, nd=1.0, ein=1.0, rho=1.0)

