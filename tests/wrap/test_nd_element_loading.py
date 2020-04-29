import o3seespy as o3


def load_element(nu_init, nu_setp):
    esig_v0 = 100.0

    osi = o3.OpenSeesInstance(ndm=2, ndf=2)

    mat = o3.nd_material.ElasticIsotropic(osi, 1.0e6, nu=nu_init)
    h_ele = 1.
    nodes = [
        o3.node.Node(osi, 0.0, 0.0),
        o3.node.Node(osi, h_ele, 0.0),
        o3.node.Node(osi, h_ele, h_ele),
        o3.node.Node(osi, 0.0, h_ele)
    ]

    # Fix bottom node
    o3.Fix2DOF(osi, nodes[0], o3.cc.FIXED, o3.cc.FIXED)
    o3.Fix2DOF(osi, nodes[1], o3.cc.FIXED, o3.cc.FIXED)
    o3.Fix2DOF(osi, nodes[2], o3.cc.FREE, o3.cc.FREE)
    o3.Fix2DOF(osi, nodes[3], o3.cc.FREE, o3.cc.FREE)
    # Set out-of-plane DOFs to be slaved
    o3.EqualDOF(osi, nodes[2], nodes[3], [o3.cc.X, o3.cc.Y])

    ele = o3.element.SSPquad(osi, nodes, mat, 'PlaneStrain', 1, 0.0, 0.0)

    o3.constraints.Transformation(osi)
    o3.test_check.NormDispIncr(osi, tol=1.0e-3, max_iter=35, p_flag=0)
    o3.algorithm.Newton(osi)
    o3.numberer.RCM(osi)
    o3.system.FullGeneral(osi)
    o3.integrator.DisplacementControl(osi, nodes[2], o3.cc.DOF2D_Y, 0.005)
    # o3.rayleigh.Rayleigh(osi, a0, a1, 0.0, 0.0)
    o3.analysis.Static(osi)
    o3.set_parameter(osi, value=nu_setp, eles=[ele], args=['nu', 1])

    # Add static vertical pressure and stress bias
    # time_series = o3.time_series.Path(osi, time=[0, 100, 1e10], values=[0, 1, 1])
    # o3.pattern.Plain(osi, time_series)
    ts0 = o3.time_series.Linear(osi, factor=1)
    o3.pattern.Plain(osi, ts0)
    o3.Load(osi, nodes[2], [0, -esig_v0 / 2])
    o3.Load(osi, nodes[3], [0, -esig_v0 / 2])

    o3.analyze(osi, num_inc=100)
    stresses = o3.get_ele_response(osi, ele, 'stress')
    print('init_stress0: ', stresses)


