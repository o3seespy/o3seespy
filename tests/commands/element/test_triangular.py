import o3seespy as o3  # for testing only


def test_tri31():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0], [1, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(3)]
    mat = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
    o3.element.Tri31(osi, ele_nodes=ele_nodes, thick=1.0, otype=o3.cc.PLANE_STRAIN, mat=mat, pressure=1.0, rho=1.0, b1=1.0, b2=1.0)

