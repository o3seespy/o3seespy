import o3seespy as o3  # for testing only


def test_catenary_cable():
    osi = o3.OpenSeesInstance(ndm=2)
    i_node = o3.node.Node(osi, 0.0, 0.0)
    j_node = o3.node.Node(osi, 0.0, 1.0)
    o3.element.CatenaryCable(osi, i_node=i_node, j_node=j_node, weight=1.0, big_e=1.0, big_a=1.0, l0=1.0, alpha=1.0, temperature_change=1.0, rho=1.0, error_tol=1.0, nsubsteps=1, mass_type=1)

