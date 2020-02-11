import o3seespy as o3  # for testing only


def test_four_node_tetrahedron():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
    mat = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
    o3.element.FourNodeTetrahedron(osi, ele_nodes=ele_nodes, mat=mat, b1=1.0, b2=1.0, b3=1.0)

