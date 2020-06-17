from eqsig import sdof
import numpy as np
import os

import o3seespy as o3


def get_inelastic_response(mass, k_spring, f_yield, motion, dt, xi=0.05, r_post=0.0):
    """
    Run seismic analysis of a nonlinear SDOF

    Parameters
    ----------
    mass: float
        SDOF mass
    k_spring: float
        Spring stiffness
    f_yield: float
        Yield strength
    motion: array_like,
        Acceleration values
    dt: float, time step of acceleration values
    xi: damping ratio
    r_post: post-yield stiffness

    Returns
    -------
    outputs: dict
        Dictionary containing time series from analysis
    """
    osi = o3.OpenSeesInstance(ndm=2, state=0)

    # Establish nodes
    bot_node = o3.node.Node(osi, 0, 0)
    top_node = o3.node.Node(osi, 0, 0)

    # Fix bottom node
    o3.Fix3DOF(osi, top_node, o3.cc.FREE, o3.cc.FIXED, o3.cc.FIXED)
    o3.Fix3DOF(osi, bot_node, o3.cc.FIXED, o3.cc.FIXED, o3.cc.FIXED)
    # Set out-of-plane DOFs to be constrained
    o3.EqualDOF(osi, top_node, bot_node, [o3.cc.Y, o3.cc.ROTZ])

    # nodal mass (weight / g):
    o3.Mass(osi, top_node, mass, 0., 0.)

    # Define material
    bilinear_mat = o3.uniaxial_material.Steel01(osi, fy=f_yield, e0=k_spring, b=r_post)

    # Assign zero length element, # Note: pass actual node and material objects into element
    o3.element.ZeroLength(osi, [bot_node, top_node], mats=[bilinear_mat], dirs=[o3.cc.DOF2D_X], r_flag=1)

    # Define the dynamic analysis
    acc_series = o3.time_series.Path(osi, dt=dt, values=-motion)  # should be negative
    o3.pattern.UniformExcitation(osi, dir=o3.cc.X, accel_series=acc_series)

    # set damping based on first eigen mode
    angular_freq_sqrd = o3.get_eigen(osi, solver='fullGenLapack', n=1)
    if hasattr(angular_freq_sqrd, '__len__'):
        angular_freq = angular_freq_sqrd[0] ** 0.5
    else:
        angular_freq = angular_freq_sqrd ** 0.5
    response_period = 2 * np.pi / angular_freq
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

    while o3.get_time(osi) < analysis_time:

        o3.analyze(osi, 1, analysis_dt)
        curr_time = o3.get_time(osi)
        outputs["time"].append(curr_time)
        outputs["rel_disp"].append(o3.get_node_disp(osi, top_node, o3.cc.X))
        outputs["rel_vel"].append(o3.get_node_vel(osi, top_node, o3.cc.X))
        outputs["rel_accel"].append(o3.get_node_accel(osi, top_node, o3.cc.X))
        o3.gen_reactions(osi)
        outputs["force"].append(-o3.get_node_reaction(osi, bot_node, o3.cc.X))  # Negative since diff node
    o3.wipe(osi)
    for item in outputs:
        outputs[item] = np.array(outputs[item])

    return outputs


def test_sdof(show=0):
    """
    Create a plot of an elastic analysis, nonlinear analysis and closed form elastic
    """

    folder_path = os.path.dirname(os.path.abspath(__file__))
    record_filename = folder_path + '/test_motion_dt0p01.txt'
    dt = 0.01
    rec = np.loadtxt(record_filename)
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

    time = np.arange(len(rec)) * dt

    acc_opensees_elastic = np.interp(time, outputs_elastic["time"], outputs_elastic["rel_accel"]) - rec
    acc_opensees_inelastic = np.interp(time, outputs["time"], outputs["rel_accel"]) - rec
    ux_opensees_elastic = np.interp(time, outputs_elastic["time"], outputs_elastic["rel_disp"])
    diff_disp = np.sum(abs(ux_opensees_elastic - resp_u[0]))
    diff_acc = np.sum(abs(acc_opensees_elastic - resp_a[0]))
    assert diff_disp < 1.0e-2, diff_disp
    assert diff_acc < 5.0e-1, diff_acc
    assert np.isclose(disp_inelastic_final, 0.0186556)
    if show:
        import matplotlib.pyplot as plt
        bf, sps = plt.subplots(nrows=3, sharex='col')
        sps[0].plot(outputs_elastic["time"], outputs_elastic["rel_disp"], label='O3 Elastic', lw=0.7, c='b')
        sps[0].plot(outputs["time"], outputs["rel_disp"], label='O3 Inelastic', lw=0.7, c='r')
        sps[0].plot(time, resp_u[0], label='Closed form', lw=1.4, c='g', ls='--')
        sps[1].plot(outputs_elastic["time"], outputs_elastic["rel_vel"])
        sps[1].plot(time, resp_v[0], lw=1.4, c='g', ls='--')
        sps[2].plot(time, acc_opensees_elastic, c='b')
        sps[2].plot(time, resp_a[0], lw=1.4, c='g', ls='--')
        sps[2].plot(time, acc_opensees_inelastic, c='r')
        sps[2].plot(time, rec, lw=0.7, c='k', label='Input', zorder=0)
        sps[2].axhline(f_yield, c=(0.4, 0.4, 0.4), label='$f_{yield}$', ls='--', lw=0.5, zorder=-1)
        sps[2].axhline(-f_yield, c=(0.4, 0.4, 0.4), ls='--', lw=0.5, zorder=-1)
        sps[0].set_ylabel('Disp. [m]')
        sps[1].set_ylabel('Vel. [m/s]')
        sps[2].set_ylabel('Accel. [m/s2]')
        sps[-1].set_xlabel('Time [s]')
        sps[-1].set_xlim([0, time[-1]])
        sps[0].legend(loc='upper right')
        sps[2].legend(loc='upper right')
        plt.show()


if __name__ == '__main__':
    test_sdof(show=1)
