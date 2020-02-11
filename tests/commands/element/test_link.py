import o3seespy as o3  # for testing only
import pytest


@pytest.mark.skip()
def test_two_node_link():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [0, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(2)]
    mats = [o3.uniaxial_material.Elastic(osi, 1.0),
            o3.uniaxial_material.Elastic(osi, 1.0)]
    p_delta_vals = [1.0, 1.0]
    o3.element.TwoNodeLink(osi, ele_nodes=ele_nodes, mats=mats, dir=[1, 1], p_delta_vals=p_delta_vals)

