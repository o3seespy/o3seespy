

def build_graph_links(mesh_eles, inactive_value=-1):
    import numpy as np
    mesh_eles = np.where(mesh_eles == inactive_value, -1, mesh_eles)
    mflat = mesh_eles.flatten()
    mflat_cleaned = mflat[np.where(mflat > -1)]
    graph_inds = np.arange(len(mflat_cleaned), dtype=int)
    mlist = np.ones(np.max(mesh_eles) + 1, dtype=int) * (len(mflat_cleaned) + 2)  # sets Null values to exceed length
    mlist[mflat_cleaned] = graph_inds
    rolled_down = np.roll(mesh_eles, 1, axis=0)
    rolled_down[0, :] = -1
    rolled_up = np.roll(mesh_eles, -1, axis=0)
    rolled_up[-1, :] = -1
    rolled_left = np.roll(mesh_eles, -1, axis=1)
    rolled_left[:, -1] = -1
    rolled_right = np.roll(mesh_eles, 1, axis=1)
    rolled_right[:, 0] = -1
    rolled_down_right = np.roll(rolled_down, 1, axis=1)
    rolled_down_right[:, 0] = -1
    rolled_down_left = np.roll(rolled_down, -1, axis=1)
    rolled_down_left[:, -1] = -1
    rolled_up_right = np.roll(rolled_up, 1, axis=1)
    rolled_up_right[:, 0] = -1
    rolled_up_left = np.roll(rolled_up, -1, axis=1)
    rolled_up_left[:, -1] = -1
    all_rolled = np.array([rolled_down, rolled_up, rolled_left, rolled_right,
                           rolled_down_left, rolled_down_right, rolled_up_left, rolled_up_right])

    g_adjs = []
    for i in range(len(mesh_eles)):
        for j in range(len(mesh_eles[i])):
            all_inds = all_rolled[:, i, j]
            inds = all_inds[np.where(all_inds > -1)]
            g_adjs.append((mlist[inds]))

    return g_adjs, mlist


def build_nodes_in_partition_2dmatrix(ele_partitions, pid, wrap=None):
    import numpy as np
    if wrap == 'x':
        
        partitions = np.pad(ele_partitions, ((0, 0), (0, 1)), 'constant', constant_values=-1)
        x_vals = partitions[0]
        ps = list(partitions)
        ps.append(x_vals)
        partitions = np.array(ps)
        # partitions = np.insert(partitions, len(partitions), x_vals)
    else:
        partitions = np.pad(ele_partitions, ((0, 1), (0, 1)), 'constant', constant_values=-1)
    ni = len(partitions)
    nj = len(partitions[0])
    rolled_down = np.roll(partitions, 1, axis=0)
    rolled_right = np.roll(partitions, 1, axis=1)
    rolled_down_right = np.roll(rolled_down, 1, axis=1)
    all_rolled = np.array([partitions, rolled_down, rolled_right, rolled_down_right]).transpose([1, 2, 0])
    all_rolled_f = all_rolled.reshape(ni * nj, 4)
    act = np.any(np.isin(all_rolled_f, [pid]), axis=1)
    act = act.astype(dtype=np.int32)
    act = act.reshape([ni, nj])
    return act


def build_nodes_in_partition_3dmatrix(ele_partitions, pid):
    import numpy as np
    partitions = np.pad(ele_partitions, ((0, 1), (0, 1), (0, 1)), 'constant', constant_values=-1)
    ni = len(partitions)
    nj = len(partitions[0])
    nk = len(partitions[0][0])
    rolled_down = np.roll(partitions, 1, axis=0)
    rolled_right = np.roll(partitions, 1, axis=1)
    rolled_down_right = np.roll(rolled_down, 1, axis=1)
    rolled_forward = np.roll(partitions, 1, axis=2)
    rolled_right_forward = np.roll(rolled_right, 1, axis=2)
    rolled_down_right_forward = np.roll(rolled_down_right, 1, axis=2)

    all_rolled = np.array([partitions,
                           rolled_down,
                           rolled_right,
                           rolled_down_right,
                           rolled_forward,
                           rolled_right_forward,
                           rolled_down_right_forward
                           ]).transpose([1, 2, 3, 0])
    all_rolled_f = all_rolled.reshape(ni * nj * nk, 7)
    act = np.any(np.isin(all_rolled_f, [pid]), axis=1)
    act = act.astype(dtype=np.int32)
    act = act.reshape([ni, nj, nk])
    return act


def is_node_in_partion(ele_partitions, node_i, node_j, pid):
    return pid in ele_partitions[node_i:min(len(ele_partitions), node_i + 1)][node_j: min(len(ele_partitions[0]), node_j)]


def build_graph_links_3d(mesh_eles, inactive_value=-1):
    import numpy as np
    mesh_eles = np.where(mesh_eles == inactive_value, -1, mesh_eles)
    mflat = mesh_eles.flatten()
    mflat_cleaned = mflat[np.where(mflat > -1)]
    graph_inds = np.arange(len(mflat_cleaned), dtype=int)
    mlist = np.ones(np.max(mesh_eles) + 1, dtype=int) * (len(mflat_cleaned) + 2)  # sets Null values to exceed length
    mlist[mflat_cleaned] = graph_inds
    # There should be 26 connections
    # 6 centres x = down/up, y = left/right, z = forward/back
    rolled_down = np.roll(mesh_eles, 1, axis=0)
    rolled_down[0, :, :] = -1
    rolled_up = np.roll(mesh_eles, -1, axis=0)
    rolled_up[-1, :, :] = -1
    rolled_left = np.roll(mesh_eles, -1, axis=1)
    rolled_left[:, -1, :] = -1
    rolled_right = np.roll(mesh_eles, 1, axis=1)
    rolled_right[:, 0, :] = -1
    rolled_forward = np.roll(mesh_eles, 1, axis=2)
    rolled_forward[:, :, 0] = -1
    rolled_back = np.roll(mesh_eles, -1, axis=2)
    rolled_back[:, :, -1] = -1

    #     __________
    #    /   /  /  /|
    #    ---------- |
    #   /  /  /  / |/
    #   ---------

    # 12 edge bits
    rolled_down_right = np.roll(rolled_down, 1, axis=1)
    rolled_down_right[:, 0, :] = -1
    rolled_down_left = np.roll(rolled_down, -1, axis=1)
    rolled_down_left[:, -1, :] = -1
    rolled_up_right = np.roll(rolled_up, 1, axis=1)
    rolled_up_right[:, 0, :] = -1
    rolled_up_left = np.roll(rolled_up, -1, axis=1)
    rolled_up_left[:, -1, :] = -1

    rolled_down_forward = np.roll(rolled_down, 1, axis=2)
    rolled_down_forward[:, :, 0] = -1
    rolled_down_back = np.roll(rolled_down, -1, axis=2)
    rolled_down_back[:, :, -1] = -1
    rolled_up_forward = np.roll(rolled_up, 1, axis=2)
    rolled_up_forward[:, :, 0] = -1
    rolled_up_back = np.roll(rolled_up, -1, axis=2)
    rolled_up_back[:, :, -1] = -1

    rolled_left_forward = np.roll(rolled_left, 1, axis=2)
    rolled_left_forward[:, :, 0] = -1
    rolled_left_back = np.roll(rolled_left, -1, axis=2)
    rolled_left_back[:, :, -1] = -1
    rolled_right_forward = np.roll(rolled_right, 1, axis=2)
    rolled_right_forward[:, :, 0] = -1
    rolled_right_back = np.roll(rolled_right, -1, axis=2)
    rolled_right_back[:, :, -1] = -1

    # 8 corners
    rolled_down_forward_right = np.roll(rolled_down_forward, 1, axis=1)
    rolled_down_forward_right[:, 0] = -1
    rolled_down_forward_left = np.roll(rolled_down_forward, -1, axis=1)
    rolled_down_forward_left[:, -1] = -1
    rolled_up_forward_right = np.roll(rolled_up_forward, 1, axis=1)
    rolled_up_forward_right[:, 0] = -1
    rolled_up_forward_left = np.roll(rolled_up_forward, -1, axis=1)
    rolled_up_forward_left[:, -1] = -1

    rolled_down_back_right = np.roll(rolled_down_back, 1, axis=1)
    rolled_down_back_right[:, 0] = -1
    rolled_down_back_left = np.roll(rolled_down_back, -1, axis=1)
    rolled_down_back_left[:, -1] = -1
    rolled_up_back_right = np.roll(rolled_up_back, 1, axis=1)
    rolled_up_back_right[:, 0] = -1
    rolled_up_back_left = np.roll(rolled_up_back, -1, axis=1)
    rolled_up_back_left[:, -1] = -1

    all_rolled = np.array([
        rolled_down, rolled_up, rolled_left, rolled_right, rolled_forward, rolled_back,  # centres
        rolled_down_left, rolled_down_right, rolled_up_left, rolled_up_right,  # 4-edges
        rolled_down_forward, rolled_down_back, rolled_up_forward, rolled_up_back,  # 4-edges
        rolled_left_forward, rolled_left_back, rolled_right_forward, rolled_right_back,  # 4-edges
        rolled_down_forward_right, rolled_down_forward_left, rolled_down_back_right, rolled_down_back_left,  # 4-corners
        rolled_up_forward_right, rolled_up_forward_left, rolled_up_back_right, rolled_up_back_left,  # 4-corners
    ])

    g_adjs = []
    for i in range(len(mesh_eles)):
        for j in range(len(mesh_eles[i])):
            for k in range(len(mesh_eles[i][j])):
                all_inds = all_rolled[:, i, j, k]
                inds = all_inds[np.where(all_inds > -1)]
                g_adjs.append((mlist[inds]))

    # consider a ele_nodes to ele-joins fn
    return g_adjs, mlist


def array_contains_the_same(a1, a2):
    assert len(a1) == len(a2), (len(a1), len(a2))
    for item in a1:
        assert item in a2, (item, a2)


def run_build_graph_links_3d():
    erows = 3
    ecols = 3
    edepth = 3
    # 0 1 2    9 10 11  18 19 20
    # 3 4 5   12(13)14  21 22 23
    # 6 7 8   15 16 17  24 25 26
    es = np.arange(0, erows * ecols * edepth, 1)
    eles = np.reshape(es, [erows, ecols, edepth])
    # eles = np.where(fem.soil_grid != fem.inactive_value, eles, -1)
    g_adjs, mlist = build_graph_links_3d(eles)
    eles0 = [1, 3, 4, 9, 10, 12, 13]
    array_contains_the_same(eles0, g_adjs[0])
    eles1 = [0, 2, 3, 4, 5, 9, 10, 11, 12, 13, 14]
    array_contains_the_same(eles1, g_adjs[1])
    eles2 = [1, 4, 5, 10, 11, 13, 14]
    array_contains_the_same(eles2, g_adjs[2])
    eles3 = [0, 1, 4, 6, 7, 9, 10, 12, 13, 15, 16]
    array_contains_the_same(eles3, g_adjs[3])
    eles12 = [0, 1, 3, 4, 6, 7, 9, 10, 13, 15, 16, 18, 19, 21, 22, 24, 25]
    array_contains_the_same(eles12, g_adjs[12])
    eles13 = list(np.arange(27))
    eles13.remove(13)
    array_contains_the_same(eles13, g_adjs[13])
    eles21 = [9, 10, 12, 13, 15, 16, 18, 19, 22, 24, 25]
    array_contains_the_same(eles21, g_adjs[21])
    eles15 = [3, 4, 6, 7, 12, 13, 16, 21, 22, 24, 25]
    array_contains_the_same(eles15, g_adjs[15])


if __name__ == '__main__':
    run_build_graph_links_3d()