import o3seespy as o3


def run(use_pload):
    # If pload=true then apply point load at end, else apply distributed load along beam

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

    elastic = 0

    if elastic:
        left_sect = o3.section.Elastic2D(osi, e_mod, area, i_sect)
    else:
        m_cap = 14.80
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

    w_gloads = 1
    if w_gloads:
        # Apply gravity loads
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
        k = r[1] / -end_disp
        print('k: ', k)
        k_elastic_expected = 1. / (ele_len ** 3 / (3 * e_mod * i_sect))
        print('k_elastic_expected: ', k_elastic_expected)


if __name__ == '__main__':
    run(use_pload=0)