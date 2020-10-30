

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


def build_nodes_in_partition_2dmatrix(ele_partitions, pid):
    import numpy as np
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


def is_node_in_partion(ele_partitions, node_i, node_j, pid):
    return pid in ele_partitions[node_i:min(len(ele_partitions), node_i + 1)][node_j: min(len(ele_partitions[0]), node_j)]
