import o3seespy as o3  # for testing only
import pytest


@pytest.mark.skip()
def test_simple_contact2d():
    osi = o3.OpenSeesInstance(ndm=2)
    i_node = o3.node.Node(osi, 0.0, 0.0)
    j_node = o3.node.Node(osi, 0.0, 1.0)
    mat = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
    o3.element.SimpleContact2D(osi, i_node=i_node, j_node=j_node, s_node=1, l_node=1, mat=mat, g_tol=1.0, f_tol=1.0)


@pytest.mark.skip()
def test_simple_contact3d():
    osi = o3.OpenSeesInstance(ndm=2)
    i_node = o3.node.Node(osi, 0.0, 0.0)
    j_node = o3.node.Node(osi, 0.0, 1.0)
    mat = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
    o3.element.SimpleContact3D(osi, i_node=i_node, j_node=j_node, k_node=1, l_node=1, s_node=1, lagr_node=1, mat=mat, g_tol=1.0, f_tol=1.0)


@pytest.mark.skip()
def test_beam_contact2d():
    osi = o3.OpenSeesInstance(ndm=2)
    i_node = o3.node.Node(osi, 0.0, 0.0)
    j_node = o3.node.Node(osi, 0.0, 1.0)
    mat = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
    o3.element.BeamContact2D(osi, i_node=i_node, j_node=j_node, s_node=1, l_node=1, mat=mat, width=1.0, g_tol=1.0, f_tol=1.0, c_flag=1)


@pytest.mark.skip()
def test_beam_contact3d():
    osi = o3.OpenSeesInstance(ndm=2)
    i_node = o3.node.Node(osi, 0.0, 0.0)
    j_node = o3.node.Node(osi, 0.0, 1.0)
    mat = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
    o3.element.BeamContact3D(osi, i_node=i_node, j_node=j_node, s_node=1, l_node=1, radius=1.0, crd_transf=1, mat=mat, g_tol=1.0, f_tol=1.0, c_flag=1)


@pytest.mark.skip()
def test_beam_end_contact3d():
    osi = o3.OpenSeesInstance(ndm=2)
    i_node = o3.node.Node(osi, 0.0, 0.0)
    j_node = o3.node.Node(osi, 0.0, 1.0)
    o3.element.BeamEndContact3D(osi, i_node=i_node, j_node=j_node, s_node=1, l_node=1, radius=1.0, g_tol=1.0, f_tol=1.0, c_flag=1.0)

