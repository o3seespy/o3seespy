import o3seespy as o3  # for testing only
import pytest


def test_zero_length():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [0, 0]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(2)]
    bilinear_mat = o3.uniaxial_material.Steel01(osi, fy=300., e0=200.0e3, b=0.01)
    o3.element.ZeroLength(osi, ele_nodes, mats=[bilinear_mat], dirs=[o3.cc.DOF2D_X], r_flag=1)


@pytest.mark.skip()
def test_zero_length_nd():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [0, 0]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(2)]
    mat = o3.nd_material.ElasticIsotropic(osi, 1.0, 0.3)
    uni = o3.uniaxial_material.Elastic(osi, 1.0)
    o3.element.ZeroLengthND(osi, ele_nodes=ele_nodes, mat=mat, uni=uni, orient=[1, 2, 3, 4, 5, 6])


def test_zero_length_section():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [0, 0]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(2)]
    sec = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
    o3.element.ZeroLengthSection(osi, ele_nodes=ele_nodes, sec=sec, r_flag=1.0, orient=[1])


def test_coupled_zero_length():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [0, 0]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(2)]
    mat = o3.uniaxial_material.Elastic(osi, 1.0)
    o3.element.CoupledZeroLength(osi, ele_nodes=ele_nodes, dirn1=1, dirn2=1, mat=mat, r_flag=1)


@pytest.mark.skip()
def test_zero_length_contact2dnormal():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [0, 0]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(2)]
    o3.element.ZeroLengthContact2Dnormal(osi, ele_nodes=ele_nodes, kn=1.0, kt=1.0, mu=1.0, nx=1, ny=0)


def test_zero_length_contact3d():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(2)]
    o3.element.ZeroLengthContact3D(osi, ele_nodes=ele_nodes, kn=1.0, kt=1.0, mu=1.0, c=1.0, dir=1)


@pytest.mark.skip()
def test_zero_length_contact_nts2d():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [0, 0]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(2)]
    o3.element.ZeroLengthContactNTS2D(osi, s_nd_num=1, m_nd_num=1, nodes=ele_nodes, kn=1.0, kt=1.0, phi=1.0)


@pytest.mark.skip()
def test_zero_length_interface2ddof():
    osi = o3.OpenSeesInstance(ndm=2)
    nodes = [1, 1]
    o3.element.ZeroLengthInterface2Ddof(osi, s_nd_num=1, m_nd_num=1, sdof=1, mdof=1, nodes=nodes, kn=1.0, kt=1.0, phi=1.0)


def test_zero_length_impact3d():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
    o3.element.ZeroLengthImpact3D(osi, ele_nodes=ele_nodes, direction=1, init_gap=1.0, friction_ratio=1.0, kt=1.0, kn=1.0, kn2=1.0, delta_y=1.0, cohesion=1.0)

