import o3seespy as o3  # for testing only
import pytest

def test_elastic_isotropic():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.nd_material.ElasticIsotropic(osi, e_mod=1.0, nu=1.0, rho=0.0)


def test_elastic_isotropic_prebuild():
    mat = o3.nd_material.ElasticIsotropic(None, e_mod=1.0, nu=1.0, rho=0.0)
    osi = o3.OpenSeesInstance(ndm=2)
    mat.build(osi)


def test_elastic_orthotropic():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.nd_material.ElasticOrthotropic(osi, ex=1.0, ey=1.0, ez=1.0, nu_xy=1.0, nu_yz=1.0, nu_zx=1.0, gxy=1.0, gyz=1.0, gzx=1.0, rho=0.0)


def test_j2plasticity():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.nd_material.J2Plasticity(osi, k_mod=1.0, g_mod=1.0, sig0=1.0, sig_inf=1.0, delta=1.0, big_h=1.0)


def test_drucker_prager():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.nd_material.DruckerPrager(osi, k_mod=1.0, g_mod=1.0, sigma_y=1.0, rho=1.0, rho_bar=1.0, kinf=1.0, ko=1.0, delta1=1.0, delta2=1.0, big_h=1.0, theta=1.0, density=1.0, atm_pressure=101e3)


@pytest.mark.skip()
def test_damage2p():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.nd_material.Damage2p(osi, fcc=1.0, fct=1.0, e_mod=1.0, ni=1.0, gt=1.0, gc=1.0, rho_bar=1.0, big_h=1.0, theta=1.0, tangent=1.0)


def test_plane_stress():
    osi = o3.OpenSeesInstance(ndm=3)
    mat_3d = o3.nd_material.ElasticIsotropic(osi, e_mod=1.0, nu=1.0, rho=0.0)
    o3.nd_material.PlaneStress(osi, mat3d=mat_3d)


def test_plane_strain():
    osi = o3.OpenSeesInstance(ndm=2)
    mat_3d = o3.nd_material.ElasticIsotropic(osi, e_mod=1.0, nu=1.0, rho=0.0)
    o3.nd_material.PlaneStrain(osi, mat3d=mat_3d)


def test_multiaxial_cyclic_plasticity():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.nd_material.MultiaxialCyclicPlasticity(osi, rho=1.0, k_mod=1.0, g_mod=1.0, su=1.0, ho=1.0, h=1.0, m=1.0, beta=1.0, k_coeff=1.0)


def test_bounding_cam_clay():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.nd_material.BoundingCamClay(osi, mass_density=1.0, big_c=1.0, bulk_mod=1.0, ocr=1.0, mu_o=1.0, alpha=1.0, lamb=1.0, h=1.0, m=1.0)


def test_plate_fiber():
    osi = o3.OpenSeesInstance(ndm=3)
    mat = o3.nd_material.ElasticIsotropic(osi, e_mod=1.0, nu=1.0, rho=0.0)
    o3.nd_material.PlateFiber(osi, three_d=mat)


def test_fsam():
    osi = o3.OpenSeesInstance(ndm=3)
    s_x = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    s_y = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    conc = o3.uniaxial_material.ConcreteCM(osi, fpcc=1.0, epcc=1.0, ec=1.0, rc=1.0, xcrn=1.0, ft=1.0, et=1.0, rt=1.0, xcrp=1.0,
                                    gap_close=0)
    o3.nd_material.FSAM(osi, rho=1.0, s_x=s_x, s_y=s_y, conc=conc, rou_x=1.0, rou_y=1.0, nu=1.0, alfadow=1.0)


def test_manzari_dafalias():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.nd_material.ManzariDafalias(osi, g0=1.0, nu=1.0, e_init=1.0, mc=1.0, c=1.0, lambda_c=1.0, e0=1.0, ksi=1.0, p_atm=1.0, m=1.0, h0=1.0, ch=1.0, nb=1.0, a0=1.0, nd=1.0, z_max=1.0, cz=1.0, den=1.0)


def test_pm4sand():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.nd_material.PM4Sand(osi, d_r=1.0, g_o=1.0, h_po=1.0, den=1.0, p_atm=1.0, h_o=1.0, e_max=1.0, e_min=1.0, n_b=1.0, n_d=1.0, a_do=1.0, z_max=1.0, c_z=1.0, c_e=1.0, phi_cv=1.0, nu=1.0, g_degr=1.0, c_dr=1.0, c_kaf=1.0, q_bolt=1.0, r_bolt=1.0, m_par=1.0, f_sed=1.0, p_sed=1.0)


@pytest.mark.skip()
def test_stress_density_model():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.nd_material.StressDensityModel(osi, den=1.0, e_init=1.0, big_a=1.0, n=1.0, nu=1.0, a1=1.0, b1=1.0, a2=1.0, b2=1.0, a3=1.0, b3=1.0, fd=1.0, mu_0=1.0, mu_cyc=1.0, sc=1.0, big_m=1.0, p_atm=1.0)


def test_acoustic_medium():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.nd_material.AcousticMedium(osi, k_mod=1.0, rho=1.0)

