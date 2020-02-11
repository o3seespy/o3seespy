import o3seespy as o3  # for testing only
import pytest

def test_quad_up():
    osi = o3.OpenSeesInstance(ndm=2)
    mat = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
    coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    o3.element.QuadUP(osi, ele_nodes=ele_nodes, thick=1.0, mat=mat, bulk=1.0, fmass=1.0, h_perm=1.0, v_perm=1.0, b1=0, b2=0, t=0)


@pytest.mark.skip()
def test_brick_up():
    osi = o3.OpenSeesInstance(ndm=3, ndf=6)
    coords = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    mat = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
    o3.element.BrickUP(osi, ele_nodes=ele_nodes, mat=mat, bulk=1.0, fmass=1.0, perm_x=1.0, perm_y=1.0, perm_z=1.0, b_x=0, b_y=0, b_z=0)


def test_bbar_quad_up():
    osi = o3.OpenSeesInstance(ndm=2)
    mat = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
    coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    o3.element.BbarQuadUP(osi, ele_nodes=ele_nodes, thick=1.0, mat=mat, bulk=1.0, fmass=1.0, h_perm=1.0, v_perm=1.0, b1=0, b2=0, t=0)


@pytest.mark.skip()
def test_bbar_brick_up():
    osi = o3.OpenSeesInstance(ndm=3)
    coords = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    mat = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
    o3.element.BbarBrickUP(osi, ele_nodes=ele_nodes, mat=mat, bulk=1.0, fmass=1.0, perm_x=1.0, perm_y=1.0, perm_z=1.0, b_x=0, b_y=0, b_z=0)


@pytest.mark.skip()
def test_n94quad_up():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1], [0.5, 0.5, 0.5]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    mat = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
    o3.element.N94QuadUP(osi, ele_nodes=ele_nodes, thick=1.0, mat=mat, bulk=1.0, fmass=1.0, h_perm=1.0, v_perm=1.0, b1=0, b2=0)


@pytest.mark.skip()
def test_n208brick_up():
    osi = o3.OpenSeesInstance(ndm=3)
    coords = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    mat = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
    o3.element.N208BrickUP(osi, ele_nodes=ele_nodes, mat=mat, bulk=1.0, fmass=1.0, perm_x=1.0, perm_y=1.0, perm_z=1.0, b_x=0, b_y=0, b_z=0)

