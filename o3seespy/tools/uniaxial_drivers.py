import numpy as np
import o3seespy as o3


def run_uniaxial_disp_driver(osi, mat_obj, disps, target_d_inc=1.0e-5):
    """
    A Uniaxial material displacement controlled driver

    Parameters
    ----------
    osi: o3.OpenSeesInstance()
        An Opensees instance
    mat_obj: o3.uniaxial_material.UniaxialMaterialBase()
        An instance of uniaxial material
    disps: array_like
        Target displacements
    target_d_inc: float
        Target displacement increment

    Returns
    -------
    disp: array_like
        Actual displacements
    react: array_like
        Reactions at each displacement
    """
    left_node = o3.node.Node(osi, 0, 0)
    right_node = o3.node.Node(osi, 0, 0)
    o3.Fix1DOF(osi, left_node, o3.cc.FIXED)
    o3.Fix1DOF(osi, right_node, o3.cc.FREE)
    ele = o3.element.ZeroLength(osi, [left_node, right_node], mats=[mat_obj], dirs=[o3.cc.DOF2D_X], r_flag=1)

    disp = []
    react = []

    o3.constraints.Plain(osi)
    o3.numberer.RCM(osi)
    o3.system.BandGeneral(osi)
    o3.test_check.NormDispIncr(osi, 0.002, 10, p_flag=0)
    o3.algorithm.Newton(osi)
    o3.integrator.DisplacementControl(osi, right_node, o3.cc.X, -target_d_inc)
    o3.analysis.Static(osi)
    ts_po = o3.time_series.Linear(osi, factor=1)
    o3.pattern.Plain(osi, ts_po)
    o3.Load(osi, right_node, [1.0])

    disp.append(0)
    react.append(0)
    d_incs = np.diff(disps, prepend=0)
    for i in range(len(disps)):
        d_inc_i = d_incs[i]
        if target_d_inc < abs(d_inc_i):
            n = int(abs(d_inc_i / target_d_inc))
            d_step = d_inc_i / n
        else:
            n = 1
            d_step = d_inc_i
        o3.integrator.DisplacementControl(osi, right_node, o3.cc.X, d_step)
        o3.analyze(osi, n)
        o3.gen_reactions(osi)
        react.append(o3.get_ele_response(osi, ele, 'force')[0])
        end_disp = -o3.get_node_disp(osi, right_node, dof=o3.cc.X)
        disp.append(end_disp)
    return np.array(disp), np.array(react)


def run_uniaxial_force_driver(osi, mat_obj, forces, d_step=0.001, max_steps=10000, handle='silent'):
    """
    A Uniaxial material force-defined driver

    Parameters
    ----------
    osi: o3.OpenSeesInstance()
        An Opensees instance
    mat_obj: o3.uniaxial_material.UniaxialMaterialBase()
        An instance of uniaxial material
    forces: array_like
        Target forces
    d_step: float
        Displacement increment
    max_steps: int
        Maximum number of steps to take to achieve target force
    handle: str
        Behaviour if target force not reached, If 'silent' then  change to next target force,
        if 'warn' then print warning and go to next force,
        else raise error.

    Returns
    -------
    disp: array_like
        Actual displacements
    react: array_like
        Reactions at each displacement
    """
    left_node = o3.node.Node(osi, 0, 0)
    right_node = o3.node.Node(osi, 0, 0)
    o3.Fix1DOF(osi, left_node, o3.cc.FIXED)
    o3.Fix1DOF(osi, right_node, o3.cc.FREE)
    ele = o3.element.ZeroLength(osi, [left_node, right_node], mats=[mat_obj], dirs=[o3.cc.DOF2D_X], r_flag=1)

    o3.constraints.Plain(osi)
    o3.numberer.RCM(osi)
    o3.system.BandGeneral(osi)
    o3.test_check.NormDispIncr(osi, 0.002, 10, p_flag=0)
    o3.algorithm.Newton(osi)
    o3.integrator.DisplacementControl(osi, right_node, o3.cc.X, 0.0001)
    o3.analysis.Static(osi)
    ts_po = o3.time_series.Linear(osi, factor=1)
    o3.pattern.Plain(osi, ts_po)
    o3.Load(osi, right_node, [1.0])

    react = 0
    disp = [0]
    reacts = [react]

    diffs = np.diff(forces, prepend=0)
    orys = np.where(diffs >= 0, 1, -1)
    for i in range(len(forces)):
        ory = orys[i]
        o3.integrator.DisplacementControl(osi, right_node, o3.cc.X, -d_step * ory)
        for j in range(max_steps):
            if react * ory < forces[i] * ory:
                o3.analyze(osi, 1)
            else:
                break
            o3.gen_reactions(osi)
            react = o3.get_ele_response(osi, ele, 'force')[0]
            reacts.append(react)
            end_disp = -o3.get_node_disp(osi, right_node, dof=o3.cc.X)
            disp.append(end_disp)
        if j == max_steps - 1:
            if handle == 'silent':
                break
            if handle == 'warn':
                print(f'Target force not reached: force={react:.4g}, target: {forces[i]:.4g}')
            else:
                raise ValueError()
    return np.array(disp), np.array(reacts)


def example_run_disp_gen():
    osi = o3.OpenSeesInstance(ndm=1, ndf=1, state=0)
    mat_obj1 = o3.uniaxial_material.PySimple1(osi, 1, 1e3, 0.05, 1.0, 0.0)
    t = np.arange(0, 20, 0.01)
    umax = 0.6
    disps = np.sin(t) * np.arange(len(t)) / len(t) * umax
    disp, react = o3.tools.run_uniaxial_disp_driver(osi, mat_obj1, disps)

    import matplotlib.pyplot as plt
    plt.plot(disp, react)
    plt.show()


def example_run_force_gen():
    osi = o3.OpenSeesInstance(ndm=1, ndf=1, state=0)
    mat_obj1 = o3.uniaxial_material.PySimple1(osi, 1, 1e3, 0.05, 1.0, 0.0)
    forces = [400, 50, 550, -10, 800, -800, 800]
    disp, react = o3.tools.run_uniaxial_force_driver(osi, mat_obj1, forces)

    import matplotlib.pyplot as plt
    plt.plot(disp, react)
    plt.show()



if __name__ == '__main__':
    # example_run_force_gen()
    example_run_disp_gen()