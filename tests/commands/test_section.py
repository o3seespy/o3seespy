import o3seespy as o3  # for testing only
import pytest


def test_elastic2d():
    osi = o3.OpenseesInstance(dimensions=2)
    o3.section.Elastic2D(osi, big_e=1.0, big_a=1.0, iz=1.0, big_g=0.0, alpha_y=0.0)

def test_elastic3d():
    osi = o3.OpenseesInstance(dimensions=2)
    o3.section.Elastic3D(osi, big_e=1.0, big_a=1.0, iz=1.0, iy=1.0, big_g=1.0, big_j=1.0, alpha_y=0.0, alpha_z=0.0)


def test_fiber():
    osi = o3.OpenseesInstance(dimensions=2)
    o3.section.Fiber(osi, gj=0.0)

def test_fiber_thermal():
    osi = o3.OpenseesInstance(dimensions=2)
    o3.section.FiberThermal(osi, gj=0.0)


def test_nd_fiber():
    osi = o3.OpenseesInstance(dimensions=2)
    o3.section.NDFiber(osi)


def test_wf_section2d():
    osi = o3.OpenseesInstance(dimensions=2)
    mat = o3.uniaxial_material.Elastic(osi, 1.0)
    o3.section.WFSection2D(osi, mat=mat, d=1.0, tw=1.0, bf=1.0, tf=1.0, nfw=1.0, nff=1.0)


@pytest.mark.skip()
def test_rc_section2d():
    osi = o3.OpenseesInstance(dimensions=2)
    o3.section.RCSection2D(osi, core=1, cover=1, steel=1, d=1.0, b=1.0, cover_depth=1.0, atop=1.0, abot=1.0, aside=1.0, nfcore=1.0, nfcover=1.0, nfs=1.0)


@pytest.mark.skip()
def test_rc_circular_section():
    osi = o3.OpenseesInstance(dimensions=2)
    o3.section.RCCircularSection(osi, core=1, cover=1, steel=1, d=1.0, cover_depth=1.0, a_s=1.0, nrings_core=1, nrings_cover=1, newedges=1, nsteel=1)


def test_parallel():
    osi = o3.OpenseesInstance(dimensions=2)
    sects = [o3.section.Elastic2D(osi, big_e=1.0, big_a=1.0, iz=1.0, big_g=0.0, alpha_y=0.0),
             o3.section.Elastic2D(osi, big_e=1.0, big_a=1.0, iz=1.0, big_g=0.0, alpha_y=0.0)]
    o3.section.Parallel(osi, sects)


@pytest.mark.skip()
def test_aggregator():
    osi = o3.OpenseesInstance(dimensions=2)
    section = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
    mats = [o3.uniaxial_material.Elastic(osi, 1.0, 1.0), o3.uniaxial_material.Elastic(osi, 1.0, 1.0)]
    o3.section.Aggregator(osi, mats=mats, section=section)


def test_uniaxial():
    osi = o3.OpenseesInstance(dimensions=2)
    mat = o3.uniaxial_material.Elastic(osi, 1.0)
    o3.section.Uniaxial(osi, mat=mat, quantity='P')


def test_elastic_membrane_plate_section():
    osi = o3.OpenseesInstance(dimensions=2)
    o3.section.ElasticMembranePlateSection(osi, big_e=1.0, nu=1.0, h=1.0, rho=1.0)


@pytest.mark.skip()
def test_plate_fiber():
    osi = o3.OpenseesInstance(dimensions=2)
    mat = o3.uniaxial_material.Elastic(osi, 1.0)
    o3.section.PlateFiber(osi, mat=mat, h=1.0)


def test_bidirectional():
    osi = o3.OpenseesInstance(dimensions=2)
    o3.section.Bidirectional(osi, big_e=1.0, fy=1.0, hiso=1.0, hkin=1.0, code1='Vy', code2='P')


@pytest.mark.skip()
def test_iso2spring():
    osi = o3.OpenseesInstance(dimensions=2)
    o3.section.Iso2spring(osi, tol=1.0, k1=1.0, fyo=1.0, k2o=1.0, kvo=1.0, hb=1.0, pe=1.0, po=0.0)


@pytest.mark.skip()
def test_layered_shell():
    osi = o3.OpenseesInstance(dimensions=2)
    o3.section.LayeredShell(osi, n_layers=1, mats=1)

