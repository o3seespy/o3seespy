import o3seespy as o3  # for testing only
import pytest


def test_elastic2d():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.section.Elastic2D(osi, e_mod=1.0, area=1.0, iz=1.0, g_mod=0.0, alpha_y=0.0)


def test_elastic3d():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.section.Elastic3D(osi, e_mod=1.0, area=1.0, iz=1.0, iy=1.0, g_mod=1.0, jxx=1.0, alpha_y=0.0, alpha_z=0.0)


def test_fiber():
    osi = o3.OpenSeesInstance(ndm=2)
    mat = o3.uniaxial_material.ElasticPP(osi, e_mod=1.0, epsy_p=1.0, epsy_n=None, eps0=0.0)
    o3.section.Fiber(osi, torsion_mat=mat)


def test_fiber_thermal():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.section.FiberThermal(osi, gj=0.0)


def test_nd_fiber():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.section.NDFiber(osi)


def test_wf_section2d():
    osi = o3.OpenSeesInstance(ndm=2)
    mat = o3.uniaxial_material.Elastic(osi, 1.0)
    o3.section.WFSection2D(osi, mat=mat, d=1.0, tw=1.0, bf=1.0, tf=1.0, nfw=1.0, nff=1.0)


def test_rc_section2d():
    osi = o3.OpenSeesInstance(ndm=2)
    core_mat = o3.uniaxial_material.Concrete01(osi, -6.0, -0.004, -5.0, -0.014)
    cover_mat = o3.uniaxial_material.Concrete01(osi, -5.0, -0.002, 0.0, -0.006)
    steel_mat = o3.uniaxial_material.Steel01(osi, 60.0, 30000.0, 0.01)
    o3.section.RCSection2D(osi, core_mat=core_mat, cover_mat=cover_mat, steel_mat=steel_mat, d=1.0, b=1.0, cover_depth=1.0, atop=1.0, abot=1.0, aside=1.0, nfcore=1.0, nfcover=1.0, nfs=1.0)


def test_rc_circular_section():
    osi = o3.OpenSeesInstance(ndm=2)
    core_mat = o3.uniaxial_material.Concrete01(osi, -6.0, -0.004, -5.0, -0.014)
    cover_mat = o3.uniaxial_material.Concrete01(osi, -5.0, -0.002, 0.0, -0.006)
    steel_mat = o3.uniaxial_material.Steel01(osi, 60.0, 30000.0, 0.01)
    o3.section.RCCircularSection(osi, core_mat=core_mat, cover_mat=cover_mat, steel_mat=steel_mat, d=1.0,
                                 cover_depth=0.10, a_s=0.1, nrings_core=2, nrings_cover=2, newedges=2, nsteel=4, gj=0.0)


def test_parallel():
    osi = o3.OpenSeesInstance(ndm=2)
    secs = [o3.section.Elastic2D(osi, e_mod=1.0, area=1.0, iz=1.0, g_mod=0.0, alpha_y=0.0),
             o3.section.Elastic2D(osi, e_mod=1.0, area=1.0, iz=1.0, g_mod=0.0, alpha_y=0.0)]
    o3.section.Parallel(osi, secs)


def test_aggregator():
    osi = o3.OpenSeesInstance(ndm=2)
    section = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
    mats = [[o3.uniaxial_material.Elastic(osi, 1.0, 1.0), o3.cc.P], [o3.uniaxial_material.Elastic(osi, 1.0, 1.0), o3.cc.M_Z]]
    o3.section.Aggregator(osi, mats=mats, section=section)


def test_uniaxial():
    osi = o3.OpenSeesInstance(ndm=2)
    mat = o3.uniaxial_material.Elastic(osi, 1.0)
    o3.section.Uniaxial(osi, mat=mat, quantity='P')


def test_elastic_membrane_plate_section():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.section.ElasticMembranePlateSection(osi, e_mod=1.0, nu=1.0, h=1.0, rho=1.0)


@pytest.mark.skip()  # needs update to opensees
def test_plate_fiber():
    osi = o3.OpenSeesInstance(ndm=2)
    mat = o3.uniaxial_material.Elastic(osi, 1.0)
    o3.section.PlateFiber(osi, mat=mat, h=1.0)


def test_bidirectional():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.section.Bidirectional(osi, e_mod=1.0, fy=1.0, hiso=1.0, hkin=1.0, code1='Vy', code2='P')


def test_iso2spring():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.section.Isolator2spring(osi, tol=1.0, k1=1.0, fyo=1.0, k2o=1.0, kvo=1.0, hb=1.0, pe=1.0, po=0.0)


@pytest.mark.skip()  # unexpected issue with opensees implementation
def test_layered_shell():
    osi = o3.OpenSeesInstance(ndm=2)
    mats = [[o3.nd_material.ElasticIsotropic(osi, e_mod=1.0, v=1.0, rho=0.0), 1.0],
            [o3.nd_material.ElasticIsotropic(osi, e_mod=1.0, v=1.0, rho=0.0), 1.0],
            [o3.nd_material.ElasticIsotropic(osi, e_mod=1.0, v=1.0, rho=0.0), 1.0]]
    o3.section.LayeredShell(osi, mats=mats)

