import o3seespy as o3  # for testing only


def test_load_control():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.integrator.LoadControl(osi, incr=1.0, num_iter=1, min_incr=None, max_incr=None)


def test_displacement_control():
    osi = o3.OpenSeesInstance(ndm=2)
    node = o3.node.Node(osi, 0.0, 0.0)
    o3.integrator.DisplacementControl(osi, node, dof=1, incr=1.0, num_iter=1, d_umin=None, d_umax=None)


def test_parallel_displacement_control():
    osi = o3.OpenSeesInstance(ndm=2)
    node = o3.node.Node(osi, 0.0, 0.0)
    o3.integrator.ParallelDisplacementControl(osi, node, dof=1, incr=1.0, num_iter=1, d_umin=None, d_umax=None)


def test_min_unbal_disp_norm():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.integrator.MinUnbalDispNorm(osi, dlambda1=1.0, jd=1, min_lambda=None, max_lambda=None, det=None)


def test_arc_length():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.integrator.ArcLength(osi, s=1.0, alpha=1.0)


def test_central_difference():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.integrator.CentralDifference(osi)


def test_newmark():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.integrator.Newmark(osi, gamma=1.0, beta=1.0, form=1)


def test_hht():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.integrator.HHT(osi, alpha=1.0, gamma=1.0, beta=1.0)


def test_generalized_alpha():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.integrator.GeneralizedAlpha(osi, alpha_m=1.0, alpha_f=1.0, gamma=1.0, beta=1.0)


def test_trbdf2():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.integrator.TRBDF2(osi)


def test_explicit_difference():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.integrator.ExplicitDifference(osi)