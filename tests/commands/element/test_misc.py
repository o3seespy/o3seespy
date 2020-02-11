import o3seespy as o3  # for testing only
import pytest


def test_surface_load():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
    o3.element.SurfaceLoad(osi, ele_nodes=ele_nodes, p=1.0)


def test_vs3d4():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
    o3.element.VS3D4(osi, ele_nodes=ele_nodes, big_e=1.0, big_g=1.0, rho=1.0, big_r=1.0, alpha_n=1.0, alpha_t=1.0)


@pytest.mark.skip()
def test_ac3d8():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
    mat = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
    o3.element.AC3D8(osi, ele_nodes=ele_nodes, mat=mat)


@pytest.mark.skip()
def test_asi3d8():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.element.ASI3D8(osi, ele_nodes1=1, ele_nodes2=1)


@pytest.mark.skip()
def test_av3d4():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
    mat = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
    o3.element.AV3D4(osi, ele_nodes=ele_nodes, mat=mat)

