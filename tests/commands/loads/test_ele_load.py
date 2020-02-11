import o3seespy as o3  # for testing only
import pytest
import openseespy.opensees as opy
import numpy as np


def test_ele_load_uniform():
    osi = o3.OpenSeesInstance(ndm=2, state=3)
    ele_len = 2.0
    coords = [[0, 0], [ele_len, 0]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    transf = o3.geom_transf.Linear2D(osi, [])
    ele = o3.element.ElasticBeamColumn2D(osi, ele_nodes=ele_nodes,
                                   area=1.0, e_mod=1.0, iz=1.0, transf=transf, mass=1.0, c_mass="string")

    for i, node in enumerate(ele_nodes):
        if i == 0:
            o3.Fix3DOF(osi, node, o3.cc.FIXED, o3.cc.FIXED, o3.cc.FIXED)
        else:
            o3.Fix3DOF(osi, node, o3.cc.FREE, o3.cc.FIXED, o3.cc.FIXED)

    ts_po = o3.time_series.Linear(osi, factor=1)
    o3.pattern.Plain(osi, ts_po)
    udl = 10.
    o3.EleLoad2DUniform(osi, ele, w_y=-udl)

    tol = 1.0e-4
    o3.constraints.Plain(osi)
    o3.numberer.RCM(osi)
    o3.system.BandGeneral(osi)
    o3.test_check.NormDispIncr(osi, tol, 6)
    o3.algorithm.Newton(osi)
    n_steps_gravity = 1
    d_gravity = 1. / n_steps_gravity
    o3.integrator.LoadControl(osi, d_gravity, num_iter=10)
    o3.analysis.Static(osi)
    o3.analyze(osi, n_steps_gravity)
    opy.reactions()
    ele_loads = o3.get_ele_response(osi, ele, 'force')
    assert np.isclose(ele_loads[0], 0.0)
    assert np.isclose(ele_loads[1], udl * ele_len / 2)
    assert np.isclose(ele_loads[2], udl * ele_len ** 2 / 12)
    assert np.isclose(ele_loads[3], 0.0)
    assert np.isclose(ele_loads[4], udl * ele_len / 2)
    assert np.isclose(ele_loads[5], -udl * ele_len ** 2 / 12)
    assert np.isclose(o3.get_node_reaction(osi, ele_nodes[0], o3.cc.Y), udl * ele_len / 2)


if __name__ == '__main__':
    test_ele_load_uniform()
