import os
import o3seespy as o3


class Results2D(object):
    cache_path = ''
    coords = None
    ele2node_tags = None
    x_disp = None
    y_disp = None
    node_c = None
    dynamic = True
    used_r_starter = 0
    mat2ele_tags = None  # Assume 1-to-1 so if it uses a section then should be null
    sect2ele_tags = None  # Store position and tag - UNUSED
    mat2sect_tags = None  # UNUSED  # TODO: implement
    n_nodes_per_ele = [2, 4, 8]  # for 2D

    def __init__(self):
        from numpy import savetxt, loadtxt
        self.savetxt = savetxt
        self.loadtxt = loadtxt
        self.meta_files = ['node_c', 'mat2ele_tags', 'sect2ele_tags', 'mat2sect_tags']
        self.meta_fmt = [None, '%i', '%i', '%i']

    def start_recorders(self, osi, dt=None):
        self.used_r_starter = 1
        self.coords = o3.get_all_node_coords(osi)
        self.ele2node_tags = o3.get_all_ele_node_tags_as_dict(osi)
        self.dt = dt
        if self.dynamic:
            o3.recorder.NodesToFile(osi, self.cache_path + 'x_disp.txt', 'all', [o3.cc.DOF2D_X], 'disp', nsd=4)
            o3.recorder.NodesToFile(osi, self.cache_path + 'y_disp.txt', 'all', [o3.cc.DOF2D_Y], 'disp', nsd=4)

    def wipe_old_files(self):
        for node_len in self.n_nodes_per_ele:
            ffp = self.cache_path + f'ele_node_tags_{node_len}.txt'
            if os.path.exists(ffp):
                os.remove(ffp)
        if not self.used_r_starter:
            try:
                os.remove(self.cache_path + 'x_disp.txt')
            except FileNotFoundError:
                pass
            try:
                os.remove(self.cache_path + 'y_disp.txt')
            except FileNotFoundError:
                pass

        for fname in self.meta_files:
            try:
                os.remove(self.cache_path + f'{fname}.txt')
            except FileNotFoundError:
                pass

    def save_to_cache(self):
        self.wipe_old_files()
        self.savetxt(self.cache_path + 'coords.txt', self.coords)
        for node_len in self.ele2node_tags:
            oo = []
            for ele_tag in self.ele2node_tags[node_len]:
                oo.append([ele_tag] + self.ele2node_tags[node_len][ele_tag])
            self.savetxt(self.cache_path + f'ele2node_tags_{node_len}.txt', oo, fmt='%i')

        for i, fname in enumerate(self.meta_files):
            vals = getattr(self, fname)
            if vals is not None:
                self.savetxt(self.cache_path + f'{fname}.txt', vals, fmt=self.meta_fmt[i])
        if self.dynamic:
            if not self.used_r_starter:
                self.savetxt(self.cache_path + 'x_disp.txt', self.x_disp)
                self.savetxt(self.cache_path + 'y_disp.txt', self.y_disp)

    def load_from_cache(self):
        self.coords = self.loadtxt(self.cache_path + 'coords.txt')
        self.ele2node_tags = {}
        for node_len in self.n_nodes_per_ele:
            try:
                oo = self.loadtxt(self.cache_path + f'ele2node_tags_{node_len}.txt', ndmin=2)
                self.ele2node_tags[node_len] = {}
                for i in range(len(oo)):
                    self.ele2node_tags[node_len][oo[i, 0]] = oo[i, 1:]
            except OSError:
                continue
        for fname in self.meta_files:
            try:
                data = self.loadtxt(self.cache_path + f'{fname}.txt')
                if len(data) == 0:
                    data = None
                setattr(self, fname, data)
            except OSError:
                pass

        if self.dynamic:
            self.x_disp = self.loadtxt(self.cache_path + 'x_disp.txt')
            self.y_disp = self.loadtxt(self.cache_path + 'y_disp.txt')

