import eqsig
from eqsig import sdof
import numpy as np
import o3seespy.extensions
import o3seespy as o3


def test_1d_period(mass, k_spring, f_yield, motion, dt, xi=0.05, r_post=0.0):

    osi = o3.OpenSeesInstance(ndm=2, ndf=2, state=3)

    ele_width = 1.0
    # Establish nodes
    node_depths = np.arange(0, 2, 1)
    n_node_rows = len(node_depths)
    nd = {}
    for i in range(0, len(node_depths)):
        # Establish left and right nodes
        nd[f"R{i}L"] = o3.node.Node(osi, 0, -node_depths[i])
        nd[f"R{i}R"] = o3.node.Node(osi, ele_width, -node_depths[i])
        # set x and y dofs equal for left and right nodes
        o3.EqualDOF(osi, nd[f"R{i}L"], nd[f"R{i}R"], [o3.cc.X])

    # Fix base nodes
    o3.Fix2DOF(osi, nd[f"R{n_node_rows - 1}L"], o3.cc.FIXED, o3.cc.FIXED)
    o3.Fix2DOF(osi, nd[f"R{n_node_rows - 1}R"], o3.cc.FIXED, o3.cc.FIXED)
    o3.Fix2DOF(osi, nd[f"R0L"], o3.cc.FREE, o3.cc.FIXED)
    o3.Fix2DOF(osi, nd[f"R0R"], o3.cc.FREE, o3.cc.FIXED)

    soil_mat = o3.nd_material.ElasticIsotropic(osi, e_mod=20.0, v=0.3)
    eles = []
    for i in range(n_node_rows - 1):
        # def element
        nodes = [nd[f"R{i + 1}L"], nd[f"R{i + 1}R"], nd[f"R{i}R"], nd[f"R{i}L"]]
        eles.append(o3.element.SSPquad(osi, nodes, soil_mat, 'PlaneStrain', 1.0, 0.0, 0.0))

    # nodal mass (weight / g):
    o3.Mass(osi, nd["R0L"], mass, mass, 0.)
    o3.Mass(osi, nd["R0R"], mass, mass, 0.)

    # Define the dynamic analysis
    acc_series = o3.time_series.Path(osi, dt=dt, values=-motion)  # should be negative
    o3.pattern.UniformExcitation(osi, dir=o3.cc.X, accel_series=acc_series)

    # set damping based on first eigen mode
    angular_freq = o3.get_eigen(osi, solver='fullGenLapack', n=1) ** 0.5
    response_period = 2 * np.pi / angular_freq
    print('response_period: ', response_period)
    beta_k = 2 * xi / angular_freq
    o3.rayleigh.Rayleigh(osi, alpha_m=0.0, beta_k=beta_k, beta_k_init=0.0, beta_k_comm=0.0)

    # Run the dynamic analysis

    o3.wipe_analysis(osi)

    o3.algorithm.Newton(osi)
    o3.system.SparseGeneral(osi)
    o3.numberer.RCM(osi)
    o3.constraints.Transformation(osi)
    o3.integrator.Newmark(osi, 0.5, 0.25)
    o3.analysis.Transient(osi)

    o3.test_check.EnergyIncr(osi, tol=1.0e-10, max_iter=10)
    analysis_time = (len(motion) - 1) * dt
    analysis_dt = 0.001
    outputs = {
        "time": [],
        "rel_disp": [],
        "rel_accel": [],
        "rel_vel": [],
        "force": []
    }
    o3.extensions.to_py_file(osi, 'many_eles.py')


    while o3.get_time(osi) < analysis_time:
        o3.analyze(osi, 1, analysis_dt)
        curr_time = o3.get_time(osi)
        outputs["time"].append(curr_time)
        outputs["rel_disp"].append(o3.get_node_disp(osi, nd["R0L"], o3.cc.X))
        outputs["rel_vel"].append(o3.get_node_vel(osi, nd["R0L"], o3.cc.X))
        outputs["rel_accel"].append(o3.get_node_accel(osi, nd["R0L"], o3.cc.X))
        o3.gen_reactions(osi)
        outputs["force"].append(o3.get_ele_response(osi, eles[0], 'stress'))  # Negative since diff node
    o3.wipe(osi)
    for item in outputs:
        outputs[item] = np.array(outputs[item])

    return outputs


def test_sdof():
    """
    Create a plot of an elastic analysis, nonlinear analysis and closed form elastic

    :return:
    """

    from tests.conftest import TEST_DATA_DIR

    record_path = ''
    record_filename = 'test_motion_dt0p01.txt'
    dt = 0.01
    rec = np.loadtxt(record_path + record_filename)
    acc_signal = eqsig.AccSignal(rec, dt)
    period = 1.0
    xi = 0.05
    mass = 1.0
    f_yield = 1.5  # Reduce this to make it nonlinear
    r_post = 0.0

    periods = np.array([period])
    resp_u, resp_v, resp_a = sdof.response_series(motion=rec, dt=dt, periods=periods, xi=xi)

    k_spring = 4 * np.pi ** 2 * mass / period ** 2
    outputs = get_inelastic_response(mass, k_spring, f_yield, rec, dt, xi=xi, r_post=r_post)
    ux_opensees = outputs["rel_disp"]


    time = acc_signal.time

    run = 1
    if run:
        import matplotlib.pyplot as plt
        bf, sps = plt.subplots(nrows=3)

        sps[0].plot(time, resp_u[0], lw=0.7, c='r')
        sps[0].plot(outputs["time"], outputs["rel_disp"], ls='--')
        sps[1].plot(outputs['rel_disp'], outputs['force'])
        sps[2].plot(time, resp_a[0], lw=0.7, c='r')
        plt.show()


if __name__ == '__main__':
    test_sdof()
