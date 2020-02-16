import numpy as np
import o3seespy.extensions
import o3seespy as o3


def test_2d_site_period():

    osi = o3.OpenSeesInstance(ndm=2, ndf=2, state=3)
    # Establish nodes
    node_depths = np.arange(0, 10, 1)
    n_node_rows = len(node_depths)
    x_nodes = np.arange(0, 10, 1)
    nx = len(x_nodes)
    nd = {}

    for yy in range(0, len(node_depths)):
        for xx in range(nx):
            # Establish left and right nodes
            nd[f"X{xx}Y{yy}"] = o3.node.Node(osi, x_nodes[xx], -node_depths[yy])
            # set x and y dofs equal for left and right nodes
        o3.EqualDOF(osi, nd[f"X0Y{yy}"], nd[f"X{nx - 1}Y{yy}"], [o3.cc.X])

    # Fix base nodes
    for xx in range(nx):
        o3.Fix2DOF(osi, nd[f"X{xx}Y{n_node_rows -1}"], o3.cc.FIXED, o3.cc.FIXED)
    for yy in range(0, len(node_depths) - 1):
        for xx in range(nx):
            o3.Fix2DOF(osi, nd[f"X{xx}Y{yy}"], o3.cc.FREE, o3.cc.FIXED)

    vs = 150.0
    rho = 1.8
    g_mod = vs ** 2 * rho
    poissons_ratio = 0.3
    e_mod = 2 * g_mod * (1 + poissons_ratio)
    ele_thick = 1.0

    soil_mat = o3.nd_material.ElasticIsotropic(osi, e_mod=e_mod, v=poissons_ratio, rho=rho)
    eles = []
    for yy in range(0, len(node_depths) - 1):
        for xx in range(nx - 1):
            # def element
            nodes = [nd[f"X{xx}Y{yy + 1}"], nd[f"X{xx + 1}Y{yy + 1}"], nd[f"X{xx + 1}Y{yy}"], nd[f"X{xx}Y{yy}"]]
            eles.append(o3.element.SSPquad(osi, nodes, soil_mat, o3.cc.PLANE_STRAIN, ele_thick, 0.0, 0.0))

    # set damping based on first eigen mode
    angular_freq_sqrd = o3.get_eigen(osi, solver='fullGenLapack', n=1)
    if hasattr(angular_freq_sqrd, '__len__'):
        angular_freq = angular_freq_sqrd[0] ** 0.5
    else:
        angular_freq = angular_freq_sqrd ** 0.5
    # o3.extensions.to_py_file(osi, 'many_eles_2d.py')
    response_period = 2 * np.pi / angular_freq
    expected_period = 4 * max(node_depths) / vs
    print('response_period: ', response_period, expected_period, response_period / expected_period)
    assert np.isclose(response_period, expected_period, rtol=0.01)

