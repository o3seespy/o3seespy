import o3seespy as o3
import numpy as np


def test_build_regular_node_mesh_2dxy():
    xs = [1, 2, 3]
    ys = [3, 4]
    zs = 4.5
    osi = o3.OpenSeesInstance(ndm=3, ndf=3)
    sn = o3.node.build_regular_node_mesh(osi, xs, ys, zs)
    sn = np.array(sn)
    print(sn.shape)
    assert len(sn.shape) == 2  # two axes
    assert sn.shape[0] == 3  # y-axis
    assert sn.shape[1] == 2  # x-axis
    assert sn[2][1].x == 3.0
    assert sn[2][1].y == 4.0
    assert sn[2][1].z == 4.5, sn[2][1].z
    # No z-axis
    zs = None
    osi = o3.OpenSeesInstance(ndm=2, ndf=3)
    sn = o3.node.build_regular_node_mesh(osi, xs, ys, zs)
    sn = np.array(sn)
    assert len(sn.shape) == 2  # two axes

    assert sn.shape[0] == 3  # x-axis
    assert sn.shape[1] == 2  # y-axis
    assert sn[2][1].x == 3.0
    assert sn[2][1].y == 4.0
    assert not hasattr(sn[2][1], 'z'), sn[1][2].z


def test_build_regular_node_mesh_2dxz():
    xs = [1, 2, 3]
    ys = [3]
    zs = [4.5, 6.0]

    osi = o3.OpenSeesInstance(ndm=3, ndf=3)
    sn = o3.node.build_regular_node_mesh(osi, xs, ys, zs)

    sn = np.array(sn)
    assert len(sn.shape) == 3  # three axes
    assert sn.shape[0] == 3  # y-axis
    assert sn.shape[1] == 1  # x-axis
    assert sn.shape[2] == 2  # x-axis
    assert sn[2][0][1].x == 3.0

