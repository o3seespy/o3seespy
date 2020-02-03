import o3seespy as o3  # for testing only
import pytest


def test_hardening():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.Hardening(osi, e_mod=1.0, sigma_y=1.0, h_iso=1.0, h_kin=1.0, eta=0.0)


@pytest.mark.skip()
def test_cast():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.Cast(osi, n=1, bo=1.0, h=1.0, fy=1.0, e_mod=1.0, big_l=1.0, b=1.0, ro=1.0, c_r1=1.0, c_r2=1.0, a1=None, a2=1.0, a3=None, a4=1.0)


def test_viscous_damper():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.ViscousDamper(osi, k_el=1.0, cd=1.0, alpha=1.0, l_gap=0.0, nm=1, rel_tol=1e-6, abs_tol=1e-10, max_half=15)


def test_bilinear_oil_damper():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.BilinearOilDamper(osi, k_el=1.0, cd=1.0, fr=1.0, p=1.0, l_gap=0.0, nm=1, rel_tol=1e-6, abs_tol=1e-10, max_half=15)


def test_bilin():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.Bilin(osi, k0=1.0, as__plus=1.0, as__neg=1.0, my__plus=1.0, my__neg=1.0, lamda_s=1.0, lamda_c=1.0, lamda_a=1.0, lamda_k=1.0, c_s=1.0, c_c=1.0, c_a=1.0, c_k=1.0, theta_p__plus=1.0, theta_p__neg=1.0, theta_pc__plus=1.0, theta_pc__neg=1.0, res__pos=1.0, res__neg=1.0, theta_u__plus=1.0, theta_u__neg=1.0, d__plus=1.0, d__neg=1.0, n_factor=0.0)


def test_mod_imk_peak_oriented():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.ModIMKPeakOriented(osi, k0=1.0, as__plus=1.0, as__neg=1.0, my__plus=1.0, my__neg=1.0, lamda_s=1.0, lamda_c=1.0, lamda_a=1.0, lamda_k=1.0, c_s=1.0, c_c=1.0, c_a=1.0, c_k=1.0, theta_p__plus=1.0, theta_p__neg=1.0, theta_pc__plus=1.0, theta_pc__neg=1.0, res__pos=1.0, res__neg=1.0, theta_u__plus=1.0, theta_u__neg=1.0, d__plus=1.0, d__neg=1.0)


def test_mod_imk_pinching():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.ModIMKPinching(osi, k0=1.0, as__plus=1.0, as__neg=1.0, my__plus=1.0, my__neg=1.0, fpr_pos=1.0, fpr_neg=1.0, a_pinch=1.0, lamda_s=1.0, lamda_c=1.0, lamda_a=1.0, lamda_k=1.0, c_s=1.0, c_c=1.0, c_a=1.0, c_k=1.0, theta_p__plus=1.0, theta_p__neg=1.0, theta_pc__plus=1.0, theta_pc__neg=1.0, res__pos=1.0, res__neg=1.0, theta_u__plus=1.0, theta_u__neg=1.0, d__plus=1.0, d__neg=1.0)


def test_saws():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.SAWS(osi, f0=1.0, fi=1.0, du=1.0, s0=1.0, r1=1.0, r2=1.0, r3=1.0, r4=1.0, alpha=1.0, beta=1.0)


@pytest.mark.skip()
def test_bar_slip():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.BarSlip(osi, fc=1.0, fy=1.0, es=1.0, fu=1.0, eh=1.0, db=1.0, ld=1.0, nb=1, depth=1.0, height=1.0, anc_lratio=1.0, bs_flag="Strong", otype="beamtop", damage='Damage', unit='psi')


def test_bond_sp01():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.BondSP01(osi, fy=1.0, sy=1.0, fu=1.0, su=1.0, b=1.0, big_r=1.0)


def test_fatigue():
    osi = o3.OpenSeesInstance(ndm=2)
    other = o3.uniaxial_material.Hardening(osi, e_mod=1.0, sigma_y=1.0, h_iso=1.0, h_kin=1.0, eta=0.0)
    o3.uniaxial_material.Fatigue(osi, other=other, e0=0.191, m=-0.458, min=-1e16, max=1e16)


@pytest.mark.skip()
def test_impact_material():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.ImpactMaterial(osi, k1=1.0, k2=1.0, sigy=1.0, gap=1.0)


@pytest.mark.skip()
def test_hyperbolic_gap_material():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.HyperbolicGapMaterial(osi, kmax=1.0, kur=1.0, rf=1.0, fult=1.0, gap=1.0)


@pytest.mark.skip()
def test_limit_state():
    osi = o3.OpenSeesInstance(ndm=2)
    curve = 1
    o3.uniaxial_material.LimitState(osi, s1p=1.0, e1p=1.0, s2p=1.0, e2p=1.0, s3p=1.0, e3p=1.0, s1n=1.0, e1n=1.0, s2n=1.0, e2n=1.0, s3n=1.0, e3n=1.0, pinch_x=1.0, pinch_y=1.0, damage1=1.0, damage2=1.0, beta=1.0, curve=curve, curve_type=1)


def test_min_max():
    osi = o3.OpenSeesInstance(ndm=2)
    other = o3.uniaxial_material.Hardening(osi, e_mod=1.0, sigma_y=1.0, h_iso=1.0, h_kin=1.0, eta=0.0)
    o3.uniaxial_material.MinMax(osi, other=other, min_strain=1e-16, max_strain=1e16)


def test_elastic_bilin():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.ElasticBilin(osi, ep1=1.0, ep2=1.0, eps_p2=1.0, en1=None, en2=None, eps_n2=None)


def test_elastic_multi_linear():
    osi = o3.OpenSeesInstance(ndm=2)
    strain = [1.0, 1.0]
    stress = [1.0, 1.0]
    o3.uniaxial_material.ElasticMultiLinear(osi, eta=0.0, strain=strain, stress=stress)


def test_multi_linear():
    osi = o3.OpenSeesInstance(ndm=2)
    pts = [0.01, 1.0, 0.02, 2.0]
    o3.uniaxial_material.MultiLinear(osi, pts=pts)


def test_init_strain_material():
    osi = o3.OpenSeesInstance(ndm=2)
    other = o3.uniaxial_material.ModIMKPinching(osi, k0=1.0, as__plus=1.0, as__neg=1.0, my__plus=1.0, my__neg=1.0, fpr_pos=1.0, fpr_neg=1.0, a_pinch=1.0, lamda_s=1.0, lamda_c=1.0, lamda_a=1.0, lamda_k=1.0, c_s=1.0, c_c=1.0, c_a=1.0, c_k=1.0, theta_p__plus=1.0, theta_p__neg=1.0, theta_pc__plus=1.0, theta_pc__neg=1.0, res__pos=1.0, res__neg=1.0, theta_u__plus=1.0, theta_u__neg=1.0, d__plus=1.0, d__neg=1.0)
    o3.uniaxial_material.InitStrainMaterial(osi, other=other, init_strain=1.0)


def test_init_stress_material():
    osi = o3.OpenSeesInstance(ndm=2)
    other = o3.uniaxial_material.ModIMKPinching(osi, k0=1.0, as__plus=1.0, as__neg=1.0, my__plus=1.0, my__neg=1.0,
                                                fpr_pos=1.0, fpr_neg=1.0, a_pinch=1.0, lamda_s=1.0, lamda_c=1.0,
                                                lamda_a=1.0, lamda_k=1.0, c_s=1.0, c_c=1.0, c_a=1.0, c_k=1.0,
                                                theta_p__plus=1.0, theta_p__neg=1.0, theta_pc__plus=1.0,
                                                theta_pc__neg=1.0, res__pos=1.0, res__neg=1.0, theta_u__plus=1.0,
                                                theta_u__neg=1.0, d__plus=1.0, d__neg=1.0)
    o3.uniaxial_material.InitStressMaterial(osi, other=other, init_stress=1.0)


def test_path_independent():
    osi = o3.OpenSeesInstance(ndm=2)
    other = o3.uniaxial_material.Hardening(osi, e_mod=1.0, sigma_y=1.0, h_iso=1.0, h_kin=1.0, eta=0.0)
    o3.uniaxial_material.PathIndependent(osi, other=other)


def test_ecc01():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.ECC01(osi, sigt0=1.0, epst0=1.0, sigt1=1.0, epst1=1.0, epst2=1.0, sigc0=1.0, epsc0=1.0, epsc1=1.0, alpha_t1=1.0, alpha_t2=1.0, alpha_c=1.0, alpha_cu=1.0, beta_t=1.0, beta_c=1.0)


def test_self_centering():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.SelfCentering(osi, k1=1.0, k2=1.0, sig_act=1.0, beta=1.0, eps_slip=0, eps_bear=0, r_bear=None)


def test_viscous():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.Viscous(osi, big_c=1.0, alpha=1.0)


def test_bouc_wen():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.BoucWen(osi, alpha=1.0, ko=1.0, n=1.0, gamma=1.0, beta=1.0, ao=1.0, delta_a=1.0, delta_nu=1.0, delta_eta=1.0)


def test_bwbn():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.BWBN(osi, alpha=1.0, ko=1.0, n=1.0, gamma=1.0, beta=1.0, ao=1.0, q=1.0, zetas=1.0, p=1.0, shi=1.0, delta_shi=1.0, lamb=1.0, tol=1.0, max_iter=1.0)


def test_axial_sp():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.AxialSp(osi, sce=1.0, fty=1.0, fcy=1.0, bte=1.0, bty=1.0, bcy=1.0, fcr=1.0)


def test_axial_sp_hd():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.AxialSpHD(osi, sce=1.0, fty=1.0, fcy=1.0, bte=1.0, bty=1.0, bth=1.0, bcy=1.0, fcr=1.0, ath=1.0)


def test_cfswswp():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.CFSWSWP(osi, height=1.0, width=1.0, fut=1.0, tf=1.0, ife=1.0, ifi=1.0, ts=1.0, np=1.0, ds=1.0, vs=1.0, sc=1.0, nc=1.0, otype=1, opening_area=1.0, opening_length=1.0)


def test_cfssswp():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.CFSSSWP(osi, height=1.0, width=1.0, fuf=1.0, fyf=1.0, tf=1.0, af=1.0, fus=1.0, fys=1.0, ts=1.0, np=1.0, ds=1.0, vs=1.0, sc=1.0, dt=1.0, opening_area=1.0, opening_length=1.0)

