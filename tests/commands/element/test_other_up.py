import o3seespy as o3  # for testing only


def test_ss_pquad_up():
    osi = o3.OpenSeesInstance(ndm=2)
    obj = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
    ele_nodes = [o3.node.Node(osi, 0.0, x) for x in range(4)]
    o3.element.SSPquadUP(osi, ele_nodes=ele_nodes, mat=obj, thick=1.0, f_bulk=1.0, f_den=1.0, k1=1.0, k2=1.0, void=1.0, alpha=1.0, b1=0.0, b2=0.0)


def test_ss_pbrick_up():
    osi = o3.OpenSeesInstance(ndm=2)
    obj = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
    ele_nodes = [o3.node.Node(osi, 0.0, x) for x in range(8)]
    o3.element.SSPbrickUP(osi, ele_nodes=ele_nodes, mat=obj, f_bulk=1.0, f_den=1.0, k1=1.0, k2=1.0, k3=1.0, void=0.5, alpha=1.0e-5, b1=1.0, b2=1.0, b3=1.0)


if __name__ == '__main__':
    test_ss_pquad_up()
