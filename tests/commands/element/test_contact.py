import o3seespy as o3  # for testing only
import pytest


def test_simple_contact2d():
    osi = o3.OpenSeesInstance(ndm=2, ndf=2)
    i_node = o3.node.Node(osi, 0.0, 0.0)
    j_node = o3.node.Node(osi, 0.0, 1.0)
    c_node = o3.node.Node(osi, 0.0, 0.0)
    l_node = o3.node.Node(osi, 0.0, 1.0)
    mat = o3.nd_material.ContactMaterial2D(osi, mu=1.0, g_mod=1.0, c=1.0, t=1.0)
    o3.element.SimpleContact2D(osi, i_node=i_node, j_node=j_node, c_node=c_node, l_node=l_node, mat=mat,
                               g_tol=1.0, f_tol=1.0)


@pytest.mark.skip()
def test_simple_contact3d():
    osi = o3.OpenSeesInstance(ndm=3, ndf=3)
    i_node = o3.node.Node(osi, 0.0, 0.0, 0.0)
    j_node = o3.node.Node(osi, 0.0, 1.0, 0.0)
    k_node = o3.node.Node(osi, 0.0, 0.0, 0.0)
    c_node = o3.node.Node(osi, 0.0, 0.0, 0.0)
    l_node = o3.node.Node(osi, 0.0, 1.0, 0.0)
    lagr_node = o3.node.Node(osi, 0.0, 1.0, 0.0)
    mat = o3.nd_material.ContactMaterial3D(osi, mu=1.0, g_mod=1.0, c=1.0, t=1.0)
    o3.element.SimpleContact3D(osi, i_node=i_node, j_node=j_node, k_node=k_node, l_node=l_node, c_node=c_node,
                               lagr_node=lagr_node, mat=mat, g_tol=1.0, f_tol=1.0)


def test_beam_contact2d():
    osi = o3.OpenSeesInstance(ndm=2)
    i_node = o3.node.Node(osi, 0.0, 0.0)
    j_node = o3.node.Node(osi, 0.0, 1.0)
    s_node = o3.node.Node(osi, 0.0, 1.0)
    l_node = o3.node.Node(osi, 0.0, 1.0)
    mat = o3.nd_material.ContactMaterial2D(osi, mu=1.0, g_mod=1.0, c=1.0, t=1.0)
    o3.element.BeamContact2D(osi, i_node=i_node, j_node=j_node, s_node=s_node, l_node=l_node, mat=mat, width=1.0,
                             g_tol=1.0, f_tol=1.0, c_flag=1)


@pytest.mark.skip()
def test_beam_contact3d():
    osi = o3.OpenSeesInstance(ndm=3, ndf=6)
    i_node = o3.node.Node(osi, 0.0, 0.0, 0.0)
    j_node = o3.node.Node(osi, 0.0, 1.0, 0.0)
    osi.reset_model_params(ndm=3, ndf=3)
    c_node = o3.node.Node(osi, 0.0, 1.0, 0.0)
    l_node = o3.node.Node(osi, 0.0, 1.0, 0.0)
    osi.reset_model_params(ndm=3, ndf=6)
    mat = o3.nd_material.ContactMaterial3D(osi, mu=1.0, g_mod=1.0, c=1.0, t=1.0)
    crd_transf = o3.geom_transf.Linear3D(osi, vecxz=[1.0, 1.0], d_i=[1.0, 1.0], d_j=[1.0, 1.0])
    o3.element.BeamContact3D(osi, i_node=i_node, j_node=j_node, c_node=c_node, l_node=l_node, radius=1.0,
                             crd_transf=crd_transf, mat=mat, g_tol=1.0, f_tol=1.0, c_flag=1)


def test_beam_end_contact3d():
    osi = o3.OpenSeesInstance(ndm=3, ndf=6)
    i_node = o3.node.Node(osi, 0.0, 0.0, 0.0)
    j_node = o3.node.Node(osi, 0.0, 1.0, 0.0)
    osi.reset_model_params(ndm=3, ndf=3)
    c_node = o3.node.Node(osi, 0.0, 1.0, 0.0)
    l_node = o3.node.Node(osi, 0.0, 1.0, 0.0)
    o3.element.BeamEndContact3D(osi, i_node=i_node, j_node=j_node, c_node=c_node, l_node=l_node, radius=1.0,
                                g_tol=1.0, f_tol=1.0, c_flag=1.0)

