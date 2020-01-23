import numpy as np
import o3seespy as o3


def uniaxial_disp_cont_driver(mat_obj, mat_inputs, disps, dmax_lim=0.01):
    """
    A uniaxial displacement controlled driver
    """

    osi = o3.OpenseesInstance(dimensions=2, state=3)

    # Establish nodes
    left_node = o3.node.Node(osi, 0, 0)
    right_node = o3.node.Node(osi, 0, 0)

    # Fix bottom node
    o3.Fix(osi, left_node, o3.cc.FIXED, o3.cc.FIXED, o3.cc.FIXED)
    o3.Fix(osi, right_node, o3.cc.FREE, o3.cc.FIXED, o3.cc.FIXED)
    mat = mat_obj(osi, *mat_inputs)
    ele = o3.element.ZeroLength(osi, left_node, right_node, mat_x=mat, r_flag=1)

    disp = []
    react = []

    d_inc = 0.001
    d_max = 0.2
    n = int(d_max / d_inc)
    o3.constraints.Plain(osi)
    o3.numberer.RCM(osi)
    o3.system.BandGeneral(osi)
    o3.test_check.NormDispIncr(osi, 0.002, 10, p_flag=0)
    o3.algorithm.Newton(osi)
    o3.integrator.DisplacementControl(osi, right_node, o3.cc.X, -d_inc)
    o3.analysis.Static(osi)
    ts_po = o3.time_series.Linear(osi, factor=1)
    o3.pattern.Plain(osi, ts_po)
    o3.Load(osi, right_node, [1.0, 0.0, 0])

    disp.append(0)
    react.append(0)
    d_incs = np.diff(disps, prepend=0)
    for i in range(len(disps)):
        d_inc = d_incs[i]
        if dmax_lim < d_inc:
            n = int(d_inc / dmax_lim)
            d_step = d_inc / n
        else:
            d_step = d_inc
        o3.integrator.DisplacementControl(osi, right_node, o3.cc.X, d_step)
        o3.analyze(osi, n)
        o3.gen_reactions(osi)
        react.append(o3.get_ele_response(osi, ele, 'force')[0])
        end_disp = -o3.get_node_disp(osi, right_node, dof=o3.cc.X)
        disp.append(end_disp)
    return np.array(disp), np.array(react)
