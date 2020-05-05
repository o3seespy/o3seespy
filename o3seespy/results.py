import os
import o3seespy as o3


class Results2D(object):
    cache_path = ''
    coords = None
    ele_node_tags = None
    x_disp = None
    y_disp = None
    node_c = None
    dynamic = True
    used_r_starter = 0

    def __init__(self):
        from numpy import savetxt, loadtxt
        self.savetxt = savetxt
        self.loadtxt = loadtxt

    def start_recorders(self, osi):
        self.used_r_starter = 1
        self.coords = o3.get_all_node_coords(osi)
        self.ele_node_tags = o3.get_all_ele_node_tags_as_dict(osi)
        if self.dynamic:
            o3.recorder.NodesToFile(osi, self.cache_path + 'x_disp.txt', 'all', [o3.cc.DOF2D_X], 'disp', nsd=4)
            o3.recorder.NodesToFile(osi, self.cache_path + 'y_disp.txt', 'all', [o3.cc.DOF2D_Y], 'disp', nsd=4)

    def wipe_old_files(self):
        for node_len in [2, 4, 8]:
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
        try:
            os.remove(self.cache_path + 'node_c.txt')
        except FileNotFoundError:
            pass

    def save_to_cache(self):
        self.wipe_old_files()
        self.savetxt(self.cache_path + 'coords.txt', self.coords)
        for node_len in self.ele_node_tags:
            self.savetxt(self.cache_path + f'ele_node_tags_{node_len}.txt', self.ele_node_tags[node_len], fmt='%i')
        if self.dynamic:
            if not self.used_r_starter:
                self.savetxt(self.cache_path + 'x_disp.txt', self.x_disp)
                self.savetxt(self.cache_path + 'y_disp.txt', self.y_disp)
            if self.node_c is not None:
                self.savetxt(self.cache_path + 'node_c.txt', self.node_c)

    def load_from_cache(self):
        self.coords = self.loadtxt(self.cache_path + 'coords.txt')
        self.ele_node_tags = {}
        for node_len in [2, 4, 8]:
            try:
                self.ele_node_tags[node_len] = self.loadtxt(self.cache_path + f'ele_node_tags_{node_len}.txt', ndmin=2)
            except OSError:
                continue

        if self.dynamic:
            self.x_disp = self.loadtxt(self.cache_path + 'x_disp.txt')
            self.y_disp = self.loadtxt(self.cache_path + 'y_disp.txt')
            try:
                self.node_c = self.loadtxt(self.cache_path + 'node_c.txt')
                if len(self.node_c) == 0:
                    self.node_c = None
            except OSError:
                pass
