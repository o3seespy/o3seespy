import o3seespy as o3


def test_fix1dof():
    osi = o3.OpenSeesInstance(ndm=2, ndf=1)
    node = o3.node.Node(osi, 3, 4)
    o3.Fix1DOF(osi, node, o3.cc.FIXED)
    sn = o3.node.build_regular_node_mesh(osi, [3, 4], [1, 6])
    o3.Fix1DOFMulti(osi, sn[0], o3.cc.FREE)


def test_fix2dof():
    osi = o3.OpenSeesInstance(ndm=2, ndf=2)
    node = o3.node.Node(osi, 3, 4)
    o3.Fix2DOF(osi, node, o3.cc.FIXED, o3.cc.FIXED)
    sn = o3.node.build_regular_node_mesh(osi, [3, 4], [1, 6])
    o3.Fix2DOFMulti(osi, sn[0], o3.cc.FREE, o3.cc.FIXED)


def test_fix3dof():
    osi = o3.OpenSeesInstance(ndm=2, ndf=3)
    node = o3.node.Node(osi, 3, 4)
    o3.Fix3DOF(osi, node, o3.cc.FIXED, o3.cc.FIXED, o3.cc.FIXED)
    sn = o3.node.build_regular_node_mesh(osi, [3, 4], [1, 6])
    o3.Fix3DOFMulti(osi, sn[0], o3.cc.FREE, o3.cc.FIXED, o3.cc.FIXED)


def test_fix6dof():
    osi = o3.OpenSeesInstance(ndm=2, ndf=6)
    node = o3.node.Node(osi, 3, 4)
    o3.Fix6DOF(osi, node, o3.cc.FIXED, o3.cc.FIXED, o3.cc.FIXED, o3.cc.FIXED, o3.cc.FIXED, o3.cc.FIXED)
    sn = o3.node.build_regular_node_mesh(osi, [3, 4], [1, 6])
    o3.Fix6DOFMulti(osi, sn[0], o3.cc.FREE, o3.cc.FIXED, o3.cc.FIXED, o3.cc.FIXED, o3.cc.FIXED, o3.cc.FIXED)


def test_equal_dof():
    osi = o3.OpenSeesInstance(ndm=2, ndf=2)
    node1 = o3.node.Node(osi, 3, 4)
    node2 = o3.node.Node(osi, 3, 4)
    node3 = o3.node.Node(osi, 3, 4)
    o3.EqualDOF(osi, node1, node2, [o3.cc.DOF2D_X, o3.cc.DOF2D_Y])  # 1-to-1
    sn = o3.node.build_regular_node_mesh(osi, [3, 4, 5], [1, 6])
    o3.EqualDOFMulti(osi, node3, sn[0], [o3.cc.DOF2D_X, o3.cc.DOF2D_Y])  # 1-to-many
    o3.EqualDOFMulti(osi, sn[1], sn[2], [o3.cc.DOF2D_Y])  # many 1-to-1s


def test_add_fixity_to_dof():
    osi = o3.OpenSeesInstance(ndm=3, ndf=3)

    n1 = o3.node.Node(osi, 1., 1., 1.)
    o3.Fix3DOF(osi, n1, o3.cc.FIXED, o3.cc.FREE, o3.cc.FIXED)
    o3.Fix3DOF(osi, n1, o3.cc.FREE, o3.cc.FREE, o3.cc.FREE)
    # o3.Fix3DOF(osi, n1, o3.cc.FREE, o3.cc.FREE, o3.cc.FIXED)  # fails with this
    o3.add_fixity_to_dof(osi, o3.cc.DOF2D_ROTZ, [n1])


