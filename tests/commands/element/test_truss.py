import o3seespy as o3  # for testing only


def test_truss():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    mat = o3.uniaxial_material.Elastic(osi, 1.0)
    o3.element.Truss(osi, ele_nodes=ele_nodes, big_a=1.0, mat=mat, rho=1.0, c_flag=1.0, r_flag=1.0)


def test_corot_truss():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    mat = o3.uniaxial_material.Elastic(osi, 1.0)
    o3.element.CorotTruss(osi, ele_nodes=ele_nodes, big_a=1.0, mat=mat, rho=1.0, c_flag=1.0, r_flag=1.0)

