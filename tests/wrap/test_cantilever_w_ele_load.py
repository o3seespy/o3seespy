import o3seespy as o3
import numpy as np


def test_cantilever_w_force_beam_column():
    def run_w_settings(ele_len, e_mod, i_sect, pload, udl, beam_in_parts, use_pload):
        osi = o3.OpenSeesInstance(ndm=2, state=3)

        # Establish nodes
        left_node = o3.node.Node(osi, 0, 0)
        right_node = o3.node.Node(osi, ele_len, 0)

        # Fix bottom node
        o3.Fix3DOF(osi, left_node, o3.cc.FIXED, o3.cc.FIXED, o3.cc.FIXED)
        o3.Fix3DOF(osi, right_node, o3.cc.FREE, o3.cc.FREE, o3.cc.FREE)

        area = 0.5
        lp_i = 0.2
        lp_j = 0.2

        if beam_in_parts:
            elastic = 1

            if elastic:
                left_sect = o3.section.Elastic2D(osi, e_mod, area, i_sect)
            else:
                m_cap = 600.0
                b = 0.05
                phi = m_cap / (e_mod * i_sect)
                # mat_flex = o3.uniaxial_material.Steel01(osi, m_cap, e0=e_mod * i_sect, b=b)
                mat_flex = o3.uniaxial_material.ElasticBilin(osi, e_mod * i_sect, e_mod * i_sect * b, phi)
                mat_axial = o3.uniaxial_material.Elastic(osi, e_mod * area)
                left_sect = o3.section.Aggregator(osi, mats=[mat_axial.tag, o3.cc.P, mat_flex.tag, o3.cc.M_Z, mat_flex.tag,
                                                             o3.cc.M_Y])
            right_sect = o3.section.Elastic2D(osi, e_mod, area, i_sect)
            centre_sect = o3.section.Elastic2D(osi, e_mod, area, i_sect)
            integ = o3.beam_integration.HingeMidpoint(osi, left_sect, lp_i, right_sect, lp_j, centre_sect)
            beam_transf = o3.geom_transf.Linear2D(osi, )
            ele = o3.element.ForceBeamColumn(osi, [left_node, right_node], beam_transf, integ)
        else:
            beam_transf = o3.geom_transf.Linear2D(osi)
            ele = o3.element.ElasticBeamColumn2D(osi, [left_node, right_node], area, e_mod, i_sect, beam_transf)
        # Apply gravity loads

        # If true then load applied along beam

        ts_po = o3.time_series.Linear(osi, factor=1)
        o3.pattern.Plain(osi, ts_po)
        if use_pload:
            o3.Load(osi, right_node, [0, -pload, 0])
        else:
            o3.EleLoad2DUniform(osi, ele, -udl)

        tol = 1.0e-3
        o3.constraints.Plain(osi)
        o3.numberer.RCM(osi)
        o3.system.BandGeneral(osi)
        n_steps_gravity = 10
        o3.integrator.LoadControl(osi, 1. / n_steps_gravity, num_iter=10)
        o3.test_check.NormDispIncr(osi, tol, 10)
        o3.algorithm.Linear(osi)
        o3.analysis.Static(osi)
        o3.analyze(osi, n_steps_gravity)
        o3.gen_reactions(osi)
        end_disp = o3.get_node_disp(osi, right_node, dof=o3.cc.Y)
        return o3.get_ele_response(osi, ele, 'force')[:3], end_disp

    ele_len_o = 4.0
    pload_o = 10.0 * ele_len_o
    udl_o = 10.0
    e_mod_o = 200.0
    i_sect_o = 0.1

    # point load
    w_pload = 1
    disp_expected = -pload_o * ele_len_o ** 3 / (3 * e_mod_o * i_sect_o)
    r, disp = run_w_settings(ele_len_o, e_mod_o, i_sect_o, pload_o, udl_o, beam_in_parts=0, use_pload=w_pload)
    assert np.isclose(disp, disp_expected, rtol=0.01)
    assert np.isclose(r[1], pload_o, rtol=0.01)
    assert np.isclose(r[2], pload_o * ele_len_o, rtol=0.01)
    r, disp = run_w_settings(ele_len_o, e_mod_o, i_sect_o, pload_o, udl_o, beam_in_parts=1, use_pload=w_pload)
    assert np.isclose(disp, disp_expected, rtol=0.01), (disp, disp_expected)
    assert np.isclose(r[1], pload_o, rtol=0.01)
    assert np.isclose(r[2], pload_o * ele_len_o, rtol=0.01)

    # UDL
    w_pload = 0
    disp_expected = -udl_o * ele_len_o ** 4 / (8 * e_mod_o * i_sect_o)
    v_expected = udl_o * ele_len_o
    m_expected = udl_o * ele_len_o ** 2 / 2
    r, disp = run_w_settings(ele_len_o, e_mod_o, i_sect_o, pload_o, udl_o, beam_in_parts=0, use_pload=w_pload)
    assert np.isclose(disp, disp_expected, rtol=0.01), (disp, disp_expected)
    assert np.isclose(r[1], v_expected, rtol=0.01)
    assert np.isclose(r[2], m_expected, rtol=0.01)
    r, disp = run_w_settings(ele_len_o, e_mod_o, i_sect_o, pload_o, udl_o, beam_in_parts=1, use_pload=w_pload)
    assert np.isclose(disp, disp_expected, rtol=0.05), (disp, disp_expected)
    assert np.isclose(r[1], v_expected, rtol=0.01)
    assert np.isclose(r[2], m_expected, rtol=0.02), (r[2], m_expected)


def test_disp_control_cantilever_nonlinear():
    osi = o3.OpenSeesInstance(ndm=2, state=3)

    ele_len = 4.0

    # Establish nodes
    left_node = o3.node.Node(osi, 0, 0)
    right_node = o3.node.Node(osi, ele_len, 0)

    # Fix bottom node
    o3.Fix3DOF(osi, left_node, o3.cc.FIXED, o3.cc.FIXED, o3.cc.FIXED)
    o3.Fix3DOF(osi, right_node, o3.cc.FREE, o3.cc.FREE, o3.cc.FREE)
    e_mod = 200.0
    i_sect = 0.1
    area = 0.5
    lp_i = 0.1
    lp_j = 0.1

    m_cap = 0.30
    b = 0.05
    phi = m_cap / (e_mod * i_sect)
    # mat_flex = o3.uniaxial_material.Steel01(osi, m_cap, e0=e_mod * i_sect, b=b)
    mat_flex = o3.uniaxial_material.ElasticBilin(osi, e_mod * i_sect, e_mod * i_sect * b, phi)
    mat_axial = o3.uniaxial_material.Elastic(osi, e_mod * area)
    left_sect = o3.section.Aggregator(osi, mats=[[mat_axial, o3.cc.P], [mat_flex, o3.cc.M_Z], [mat_flex, o3.cc.M_Y]])
    right_sect = o3.section.Elastic2D(osi, e_mod, area, i_sect)
    centre_sect = o3.section.Elastic2D(osi, e_mod, area, i_sect)
    integ = o3.beam_integration.HingeMidpoint(osi, left_sect, lp_i, right_sect, lp_j, centre_sect)
    beam_transf = o3.geom_transf.Linear2D(osi, )
    ele = o3.element.ForceBeamColumn(osi, [left_node, right_node], beam_transf, integ)

    # start displacement controlled
    d_inc = 0.01

    # opy.wipeAnalysis()
    o3.constraints.Plain(osi)
    o3.numberer.RCM(osi)
    o3.system.BandGeneral(osi)
    o3.test_check.NormUnbalance(osi, 2, max_iter=10, p_flag=0)
    # o3.test_check.FixedNumIter(osi, max_iter=10)
    # o3.test_check.NormDispIncr(osi, 0.002, 10, p_flag=0)
    o3.algorithm.Newton(osi)
    o3.integrator.DisplacementControl(osi, right_node, o3.cc.Y, d_inc)
    o3.analysis.Static(osi)
    ts_po = o3.time_series.Linear(osi, factor=1)
    o3.pattern.Plain(osi, ts_po)
    o3.Load(osi, right_node, [0.0, 1.0, 0])
    o3.analyze(osi, 4)
    end_disp = o3.get_node_disp(osi, right_node, dof=o3.cc.Y)
    r = o3.get_ele_response(osi, ele, 'force')[:3]
    k = r[1] / -end_disp
    k_elastic_expected = 1. / (ele_len ** 3 / (3 * e_mod * i_sect))
    assert np.isclose(k, k_elastic_expected)
    o3.analyze(osi, 6)
    end_disp = o3.get_node_disp(osi, right_node, dof=o3.cc.Y)
    r = o3.get_ele_response(osi, ele, 'force')[:3]
    k = r[1] / -end_disp
    assert k < 0.95 * k_elastic_expected


def run():
    osi = o3.OpenSeesInstance(ndm=2, state=3)

    ele_len = 4.0

    # Establish nodes
    left_node = o3.node.Node(osi, 0, 0)
    right_node = o3.node.Node(osi, ele_len, 0)

    # Fix bottom node
    o3.Fix(osi, left_node, o3.cc.FIXED, o3.cc.FIXED, o3.cc.FIXED)
    o3.Fix(osi, right_node, o3.cc.FREE, o3.cc.FREE, o3.cc.FREE)
    e_mod = 200.0e2
    i_sect = 0.1
    area = 0.5
    lp_i = 0.1
    lp_j = 0.1
    beam_in_parts = 1
    if beam_in_parts:
        elastic = 0

        if elastic:
            left_sect = o3.section.Elastic2D(osi, e_mod, area, i_sect)
        else:
            m_cap = 18.80
            b = 0.05
            phi = m_cap / (e_mod * i_sect)
            # mat_flex = o3.uniaxial_material.Steel01(osi, m_cap, e0=e_mod * i_sect, b=b)
            mat_flex = o3.uniaxial_material.ElasticBilin(osi, e_mod * i_sect, e_mod * i_sect * b, phi)
            mat_axial = o3.uniaxial_material.Elastic(osi, e_mod * area)
            left_sect = o3.section.Aggregator(osi, mats=[mat_axial.tag, o3.cc.P, mat_flex.tag, o3.cc.M_Z, mat_flex.tag, o3.cc.M_Y])
        right_sect = o3.section.Elastic2D(osi, e_mod, area, i_sect)
        centre_sect = o3.section.Elastic2D(osi, e_mod, area, i_sect)
        integ = o3.beam_integration.HingeMidpoint(osi, left_sect, lp_i, right_sect, lp_j, centre_sect)
        beam_transf = o3.geom_transf.Linear2D(osi, )
        ele = o3.element.ForceBeamColumn(osi, [left_node, right_node], beam_transf, integ)
    else:
        beam_transf = o3.geom_transf.Linear2D(osi)
        ele = o3.element.ElasticBeamColumn2D(osi, [left_node, right_node], area, e_mod, i_sect, beam_transf)

    w_gloads = 1
    if w_gloads:
        # Apply gravity loads
        # If true then load applied along beam
        use_pload = 0
        pload = 1.0 * ele_len
        udl = 2.0
        ts_po = o3.time_series.Linear(osi, factor=1)
        o3.pattern.Plain(osi, ts_po)
        if use_pload:
            o3.Load(osi, right_node, [0, -pload, 0])
        else:
            o3.EleLoad2DUniform(osi, ele, -udl)

        tol = 1.0e-3
        o3.constraints.Plain(osi)
        o3.numberer.RCM(osi)
        o3.system.BandGeneral(osi)
        n_steps_gravity = 10
        o3.integrator.LoadControl(osi, 1. / n_steps_gravity, num_iter=10)
        o3.test_check.NormDispIncr(osi, tol, 10)
        o3.algorithm.Linear(osi)
        o3.analysis.Static(osi)
        o3.analyze(osi, n_steps_gravity)
        o3.gen_reactions(osi)
        print('reactions: ', o3.get_ele_response(osi, ele, 'force')[:3])
        end_disp = o3.get_node_disp(osi, right_node, dof=o3.cc.Y)
        print(f'end_disp: {end_disp}')
        if use_pload:
            disp_expected = pload * ele_len ** 3 / (3 * e_mod * i_sect)
            print(f'v_expected: {pload}')
            print(f'm_expected: {pload * ele_len}')
            print(f'disp_expected: {disp_expected}')
        else:
            v_expected = udl * ele_len
            m_expected = udl * ele_len ** 2 / 2
            disp_expected = udl * ele_len ** 4 / (8 * e_mod * i_sect)
            print(f'v_expected: {v_expected}')
            print(f'm_expected: {m_expected}')
            print(f'disp_expected: {disp_expected}')

        # o3.extensions.to_py_file(osi, 'temp4.py')
    disp_load = 1
    if disp_load:
        o3.load_constant(osi, time=0.0)
        end_disp_init = o3.get_node_disp(osi, right_node, dof=o3.cc.Y)
        # start displacement controlled
        d_inc = -0.01

        # opy.wipeAnalysis()
        o3.numberer.RCM(osi)
        o3.system.BandGeneral(osi)
        o3.test_check.NormUnbalance(osi, 2, max_iter=10, p_flag=0)
        # o3.test_check.FixedNumIter(osi, max_iter=10)
        # o3.test_check.NormDispIncr(osi, 0.002, 10, p_flag=0)
        o3.algorithm.Newton(osi)
        o3.integrator.DisplacementControl(osi, right_node, o3.cc.Y, d_inc)
        o3.analysis.Static(osi)
        ts_po = o3.time_series.Linear(osi, factor=1)
        o3.pattern.Plain(osi, ts_po)
        o3.Load(osi, right_node, [0.0, 1.0, 0])
        ok = o3.analyze(osi, 10)
        end_disp = o3.get_node_disp(osi, right_node, dof=o3.cc.Y)
        print(f'end_disp: {end_disp}')
        r = o3.get_ele_response(osi, ele, 'force')[:3]
        print('reactions: ', r)
        k = r[1] / -(end_disp - end_disp_init)
        print('k: ', k)
        k_elastic_expected = 1. / (ele_len ** 3 / (3 * e_mod * i_sect))
        print('k_elastic_expected: ', k_elastic_expected)


if __name__ == '__main__':
    # test_disp_control_cantilever_nonlinear()
    run()



