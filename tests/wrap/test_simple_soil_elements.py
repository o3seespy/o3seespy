import numpy as np
import o3seespy.extensions
import o3seespy as o3
import pytest

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

    soil_mat = o3.nd_material.ElasticIsotropic(osi, e_mod=e_mod, nu=poissons_ratio, rho=rho)
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


def _run_dyn_1d_site_response(region_based=None):

    osi = o3.OpenSeesInstance(ndm=2, ndf=2, state=3)
    # Establish nodes
    node_depths = np.arange(0, 5, 1)
    n_node_rows = len(node_depths)
    x_nodes = np.arange(0, 2, 1)
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

    soil_mat = o3.nd_material.ElasticIsotropic(osi, e_mod=e_mod, nu=poissons_ratio, rho=rho)
    eles = []
    for yy in range(0, len(node_depths) - 1):
        for xx in range(nx - 1):
            # def element
            nodes = [nd[f"X{xx}Y{yy + 1}"], nd[f"X{xx + 1}Y{yy + 1}"], nd[f"X{xx + 1}Y{yy}"], nd[f"X{xx}Y{yy}"]]
            eles.append(o3.element.SSPquad(osi, nodes, soil_mat, o3.cc.PLANE_STRAIN, ele_thick, 0.0, 0.0))

    freqs = np.array([0.5, 15.0])
    xi = 0.1  # high damping to see effects
    ang_f = np.pi / freqs
    alpha_m = xi * 2.0 * ang_f[0] * ang_f[1] / (ang_f[0] + ang_f[1])  # mass proportional
    beta_k = xi * 2.0 / (ang_f[0] + ang_f[1])  # stiffness proportional
    if region_based == 'node':
        o3.region.NodeRegion(osi, nodes='all', rayleigh={'alpha_m': alpha_m, 'beta_k_init': beta_k})
    if region_based == 'ele':
        o3.region.ElementRegion(osi, eles='all', rayleigh={'alpha_m': alpha_m, 'beta_k_init': beta_k})
    else:
        o3.rayleigh.Rayleigh(osi, alpha_m=alpha_m, beta_k=0.0, beta_k_init=beta_k, beta_k_comm=0.0)

    # Define the dynamic analysis
    dt = 0.01
    vals = np.sin(np.linspace(0, np.pi, 50))
    acc_series = o3.time_series.Path(osi, dt=dt, values=-vals)  # should be negative
    o3.pattern.UniformExcitation(osi, dir=o3.cc.X, accel_series=acc_series)
    # Run the dynamic analysis
    o3.wipe_analysis(osi)

    o3.algorithm.Newton(osi)
    o3.system.SparseGeneral(osi)
    o3.numberer.RCM(osi)
    o3.constraints.Transformation(osi)
    o3.integrator.Newmark(osi, 0.5, 0.25)
    o3.analysis.Transient(osi)

    o3.test_check.EnergyIncr(osi, tol=1.0e-3, max_iter=10)
    analysis_time = 1.0
    analysis_dt = 0.001
    rec_dt = 0.01
    tn = o3.recorder.NodeToArrayCache(osi, nd["X0Y0"], [o3.cc.DOF2D_X], 'accel', dt=rec_dt)

    if region_based:
        o3.extensions.to_py_file(osi, 'wr.py')
    else:
        o3.extensions.to_py_file(osi, 'wor.py')
    curr_time = o3.get_time(osi)
    while curr_time < analysis_time:
        status = o3.analyze(osi, 1, analysis_dt)
        curr_time = o3.get_time(osi)
        if status != 0:
            print('Not ok')
            print(status)

    o3.wipe(osi)

    x_acc = tn.collect()
    return x_acc


def test_region_based_damping():
    wnr = _run_dyn_1d_site_response(region_based='node')
    wer = _run_dyn_1d_site_response(region_based='ele')
    wor = _run_dyn_1d_site_response(region_based=None)
    assert np.isclose(wor, wnr).all()
    assert np.isclose(wor, wer).all()


if __name__ == '__main__':
    test_2d_site_period()
    test_region_based_damping()