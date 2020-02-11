import o3seespy as o3  # for testing only
import pytest


def test_beam_column_joint():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
    mats = [o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None) for x in range(13)]
    o3.element.BeamColumnJoint(osi, ele_nodes, *mats, ele_height_fac=1.0, ele_width_fac=1.0)


def test_elastic_tubular_joint():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
    o3.element.ElasticTubularJoint(osi, ele_nodes=ele_nodes, brace__diameter=1.0, brace__angle=1.0, big_e=1.0, chord__diameter=1.0, chord__thickness=1.0, chord__angle=1.0)


@pytest.mark.skip()
def test_joint2d():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0], [1, 1], [0, 1], [0.5, 0.5]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    o3.element.Joint2D(osi, ele_nodes=ele_nodes, mat1=1, mat2=1, mat3=1, mat4=1, mat_c=1, lrg_dsp='', dmg='', dmg1dmg2dmg3dmg4dmg_c=1)

