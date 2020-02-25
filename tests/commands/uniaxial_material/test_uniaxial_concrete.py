import o3seespy as o3  # for testing only


def test_concrete01():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.Concrete01(osi, fpc=1.0, epsc0=1.0, fpcu=1.0, eps_u=1.0)


def test_concrete02():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.Concrete02(osi, fpc=1.0, epsc0=1.0, fpcu=1.0, eps_u=1.0, lamb=1.0, ft=1.0, ets=1.0)


def test_concrete04():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.Concrete04(osi, fc=1.0, epsc=1.0, epscu=1.0, ec=1.0, fct=1.0, et=1.0, beta=1.0)


def test_concrete06():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.Concrete06(osi, fc=1.0, e0=1.0, n=1.0, k=1.0, alpha1=1.0, fcr=1.0, ecr=1.0, b=1.0, alpha2=1.0)


def test_concrete07():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.Concrete07(osi, fc=1.0, epsc=1.0, ec=1.0, ft=1.0, et=1.0, xp=1.0, xn=1.0, r=1.0)


def test_concrete01with_sitc():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.Concrete01WithSITC(osi, fpc=1.0, epsc0=1.0, fpcu=1.0, eps_u=1.0, end_strain_sitc=0.01)


def test_confined_concrete01():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.ConfinedConcrete01(osi, sec_type='R', fpc=1.0, ec=1.0, epscu_type=1, epscu_val=1.0, nu=1, l1=1.0, l2=1.0, l3=1.0, phis=1.0, big_s=1.0, fyh=1.0, es0=1.0, ha_ratio=1.0, mu=1.0, phi_lon=1.0, gravel=1, tol=1.0, max_num_iter=1, epscu_limit=1.0, st_ratio=1)


def test_concrete_d():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.ConcreteD(osi, fc=1.0, epsc=1.0, ft=1.0, epst=1.0, ec=1.0, alphac=1.0, alphat=1.0, cesp=0.25, etap=1.15)


def test_frp_confined_concrete():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.FRPConfinedConcrete(osi, fpc1=1.0, fpc2=1.0, epsc0=1.0, big_d=1.0, c=1.0, ej=1.0, sj=1.0, tj=1.0, eju=1.0, big_s=1.0, fyl=1.0, fyh=1.0, dlong=1.0, dtrans=1.0, es=1.0, nu0=1.0, k=1.0, use_buck=1.0)


def test_concrete_cm():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.ConcreteCM(osi, fpcc=1.0, epcc=1.0, ec=1.0, rc=1.0, xcrn=1.0, ft=1.0, et=1.0, rt=1.0, xcrp=1.0, gap_close=0)

