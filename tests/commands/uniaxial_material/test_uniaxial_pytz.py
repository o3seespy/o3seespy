import o3seespy as o3  # for testing only
import pytest

def test_py_simple1():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.PySimple1(osi, soil_type=1, pult=1.0, y50=1.0, cd=1.0, c=0.0)


def test_tz_simple1():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.TzSimple1(osi, soil_type=1, tult=1.0, z50=1.0, c=0.0)


def test_qz_simple1():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.QzSimple1(osi, qz_type=1, qult=1.0, z50=1.0, suction=0.0, c=0.0)


@pytest.mark.skip('need to define bulk elements')
def test_py_liq1():
    osi = o3.OpenSeesInstance(ndm=2)
    # Define nodes
    node1 = o3.node.Node(osi, 0.0, 0.0)
    node2 = o3.node.Node(osi, 1.0, 0.0)
    node3 = o3.node.Node(osi, 1.0, 1.0)
    node4 = o3.node.Node(osi, 0.0, 1.0)
    # Define a simple elastic material
    mat = o3.uniaxial_material.Elastic(osi, 1.0)
    # Define an SSPquad element using the nodes
    ssp_ele = o3.element.SSPquad(osi, [node1, node2, node3, node4], mat, 1.0, 1.0)
    # Define the PyLiq1 material
    o3.uniaxial_material.PyLiq1(osi, soil_type=1, pult=1.0, y50=1.0, cd=1.0, c=1.0, p_res=1.0, ele1=ssp_ele, ele2=ssp_ele)


@pytest.mark.skip('need to define bulk elements')
def test_tz_liq1():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.TzLiq1(osi, tz_type=1, tult=1.0, z50=1.0, c=1.0, ele1=1.0, ele2=1.0)

