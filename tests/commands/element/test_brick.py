import o3seespy as o3  # for testing only
import pytest

def test_std_brick():
    osi = o3.OpenSeesInstance(ndm=3)
    coords = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    mat = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
    o3.element.StdBrick(osi, ele_nodes=ele_nodes, mat=mat, b1=1.0, b2=1.0, b3=1.0)


def test_bbar_brick():
    osi = o3.OpenSeesInstance(ndm=3)
    coords = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    mat = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
    o3.element.BbarBrick(osi, ele_nodes=ele_nodes, mat=mat, b1=1.0, b2=1.0, b3=1.0)


@pytest.mark.skip()
def test_n20node_brick():
    osi = o3.OpenSeesInstance(ndm=3, ndf=6)
    coords = [[0, 0, 0], [-1, 0, 0], [-1, -1, 0], [0, -1, 0],
              [0, 0, -1], [-1, 0, -1], [-1, -1, -1], [0, -1, -1],
              [-0.5, 0, 0], [-1, -0.5, 0], [-0.5, -1, 0], [0, -0.5, 0],
              [-0.5, 0, -1], [-1, -0.5, -1], [-0.5, -1, -1], [0, -0.5, -1],
              [0, 0, -0.5], [-1, 0, -0.5], [-1, -1, -0.5], [0, -1, -0.5],
              ]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    mat = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
    o3.element.N20NodeBrick(osi, ele_nodes=ele_nodes, mat=mat, bf1=1.0, bf2=1.0, bf3=1.0, mass_den=1.0)


def test_ss_pbrick():
    osi = o3.OpenSeesInstance(ndm=3)
    coords = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    mat = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
    o3.element.SSPbrick(osi, ele_nodes=ele_nodes, mat=mat, b1=1.0, b2=1.0, b3=1.0)

