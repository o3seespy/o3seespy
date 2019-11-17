import eqsig
from eqsig import sdof
import numpy as np

import openseespy.opensees as opy
from o3seespy import cc as opc


def get_inelastic_response(mass, k_spring, f_yield, motion, dt, xi=0.05, r_post=0.0):
    """
    Run seismic analysis of a nonlinear SDOF

    :param mass: SDOF mass
    :param k_spring: spring stiffness
    :param f_yield: yield strength
    :param motion: list, acceleration values
    :param dt: float, time step of acceleration values
    :param xi: damping ratio
    :param r_post: post-yield stiffness
    :return:
    """

    opy.wipe()
    opy.model('basic', '-ndm', 2, '-ndf', 3)  # 2 dimensions, 3 dof per node

    # Establish nodes
    bot_node = 1
    top_node = 2
    opy.node(bot_node, 0., 0.)
    opy.node(top_node, 0., 0.)

    # Fix bottom node
    opy.fix(top_node, opc.FREE, opc.FIXED, opc.FIXED)
    opy.fix(bot_node, opc.FIXED, opc.FIXED, opc.FIXED)
    # Set out-of-plane DOFs to be slaved
    opy.equalDOF(1, 2, *[2, 3])

    # nodal mass (weight / g):
    opy.mass(top_node, mass, 0., 0.)

    # Define material
    bilinear_mat_tag = 1
    mat_type = "Steel01"
    mat_props = [f_yield, k_spring, r_post]
    opy.uniaxialMaterial(mat_type, bilinear_mat_tag, *mat_props)

    # Assign zero length element
    beam_tag = 1
    opy.element('zeroLength', beam_tag, bot_node, top_node, "-mat", bilinear_mat_tag, "-dir", 1, '-doRayleigh', 1)

    # Define the dynamic analysis
    load_tag_dynamic = 1
    pattern_tag_dynamic = 1

    values = list(-1 * motion)  # should be negative
    opy.timeSeries('Path', load_tag_dynamic, '-dt', dt, '-values', *values)
    opy.pattern('UniformExcitation', pattern_tag_dynamic, opc.X, '-accel', load_tag_dynamic)

    # set damping based on first eigen mode
    angular_freq2 = opy.eigen('-fullGenLapack', 1)
    if hasattr(angular_freq2, '__len__'):
        angular_freq2 = angular_freq2[0]
    angular_freq = angular_freq2 ** 0.5
    alpha_m = 0.0
    beta_k = 2 * xi / angular_freq
    beta_k_comm = 0.0
    beta_k_init = 0.0

    opy.rayleigh(alpha_m, beta_k, beta_k_init, beta_k_comm)

    # Run the dynamic analysis

    opy.wipeAnalysis()

    opy.algorithm('Newton')
    opy.system('SparseGeneral')
    opy.numberer('RCM')
    opy.constraints('Transformation')
    opy.integrator('Newmark', 0.5, 0.25)
    opy.analysis('Transient')

    tol = 1.0e-10
    iterations = 10
    opy.test('EnergyIncr', tol, iterations, 0, 2)
    analysis_time = (len(values) - 1) * dt
    analysis_dt = 0.001
    outputs = {
        "time": [],
        "rel_disp": [],
        "rel_accel": [],
        "rel_vel": [],
        "force": []
    }

    while opy.getTime() < analysis_time:
        curr_time = opy.getTime()
        opy.analyze(1, analysis_dt)
        outputs["time"].append(curr_time)
        outputs["rel_disp"].append(opy.nodeDisp(top_node, 1))
        outputs["rel_vel"].append(opy.nodeVel(top_node, 1))
        outputs["rel_accel"].append(opy.nodeAccel(top_node, 1))
        opy.reactions()
        outputs["force"].append(-opy.nodeReaction(bot_node, 1))  # Negative since diff node
    opy.wipe()
    for item in outputs:
        outputs[item] = np.array(outputs[item])

    return outputs


def test_sdof():
    """
    Create a plot of an elastic analysis, nonlinear analysis and closed form elastic

    :return:
    """

    from tests.conftest import TEST_DATA_DIR

    record_path = TEST_DATA_DIR
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
    outputs_elastic = get_inelastic_response(mass, k_spring, f_yield * 100, rec, dt, xi=xi, r_post=r_post)
    ux_opensees = outputs["rel_disp"]
    disp_inelastic_final = ux_opensees[-1]

    time = acc_signal.time
    acc_opensees_elastic = np.interp(time, outputs_elastic["time"], outputs_elastic["rel_accel"]) - rec
    ux_opensees_elastic = np.interp(time, outputs_elastic["time"], outputs_elastic["rel_disp"])
    diff_disp = abs(np.sum(ux_opensees_elastic - resp_u[0]))
    diff_acc = abs(np.sum(acc_opensees_elastic - resp_a[0]))
    assert diff_disp < 1.0e-4
    assert diff_acc < 5.0e-4
    assert np.isclose(disp_inelastic_final, 0.0186556)


if __name__ == '__main__':
    test_sdof()
