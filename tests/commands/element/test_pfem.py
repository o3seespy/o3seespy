import o3seespy as o3  # for testing only
import pytest


@pytest.mark.skip()
def test_pfem_element_bubble():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
    o3.element.PFEMElementBubble(osi, ele_nodes, rho=1.0, mu=1.0, b1=1.0, b2=1.0, b3=1.0, thickness=1.0, kappa=1.0)


@pytest.mark.skip()
def test_pfem_element_compressible():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
    o3.element.PFEMElementCompressible(osi, ele_nodes, rho=1.0, mu=1.0, b1=1.0, b2=1.0, thickness=1.0, kappa=1.0)

