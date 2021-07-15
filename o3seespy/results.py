import os
import o3seespy as o3


class Results2D(object):
    coords = None
    time = None
    x_disp = None
    y_disp = None
    node_c = None
    ele_c = None
    selected_nodes = None
    used_r_starter = 0
    mat2ele_tags = None  # Assume 1-to-1 so if it uses a section then should be null
    sect2ele_tags = None  # Store position and tag - UNUSED
    mat2sect_tags = None  # UNUSED  # TODO: implement
    n_nodes_per_ele = [2, 4, 8]  # for 2D
    ele_num_base = 1
    selected_nodes = None
    _selected_node_tags = None

    def __init__(self, cache_path='', dt=None, dynamic=False, man_nodes=False):
        self.cache_path = cache_path
        self._dt = dt
        self.dynamic = dynamic
        from numpy import savetxt, loadtxt
        self.savetxt = savetxt
        self.loadtxt = loadtxt
        self.ele2node_tags = {}
        self.meta_files = ['node_c', 'ele_c', 'mat2ele_tags', 'sect2ele_tags', 'mat2sect_tags']
        self.meta_fmt = [None, '%i', '%i', '%i']
        self.pseudo_dt = None  # use if recording steps of a static analysis
        self.man_nodes = man_nodes

    def start_recorders(self, osi, dt=None):  # TODO: handle recorder time step
        self.used_r_starter = 1
        if self.man_nodes and self.selected_node_tags is None:
            self._selected_node_tags = o3.get_node_tags(osi)
        if self.coords is None:
            if self.selected_node_tags is not None:
                self.coords = []
                for node_tag in self.selected_node_tags:
                    self.coords.append(o3.get_node_coords(osi, node_tag, node_as_tag=True))
            else:
                self.coords = o3.get_all_node_coords(osi)
        if not self.ele2node_tags:
            self.ele2node_tags = o3.get_all_ele_node_tags_as_dict(osi)
        if dt is not None:
            self._dt = dt
        if self.dynamic:
            node_tags = 'all'
            if self.selected_node_tags is not None:
                node_tags = self.selected_node_tags
            o3.recorder.NodesToFile(osi, f'{self.cache_path}x_disp.txt', node_tags, [o3.cc.DOF2D_X], 'disp', nsd=4, dt=dt, nodes_as_tags=True)
            o3.recorder.NodesToFile(osi, f'{self.cache_path}y_disp.txt', node_tags, [o3.cc.DOF2D_Y], 'disp', nsd=4, dt=dt, nodes_as_tags=True)
            if not self.pseudo_dt:
                o3.recorder.TimeToFile(osi, f'{self.cache_path}timer.txt', nsd=4, dt=dt)

    def wipe_old_files(self):
        sfiles = ['ele2node_tags', 'coords', 'selected_node_tags']
        if not self.used_r_starter:
            sfiles += ['x_disp', 'y_disp', 'timer']
        for fname in sfiles:
            try:
                os.remove(f'{self.cache_path}{fname}.txt')
            except FileNotFoundError:
                pass

        for fname in self.meta_files:
            try:
                os.remove(f'{self.cache_path}{fname}.txt')
            except FileNotFoundError:
                pass

    @property
    def selected_node_tags(self):
        if self._selected_node_tags is None:
            if self.selected_nodes is None:
                return None
            self._selected_node_tags = [x.tag for x in self.selected_nodes]
        return self._selected_node_tags

    @selected_node_tags.setter
    def selected_node_tags(self, tags):
        self._selected_node_tags = tags

    def save_to_cache(self):
        self.wipe_old_files()
        if self.coords is not None:
            self.savetxt(self.cache_path + 'coords.txt', self.coords)
        ostr = [f'{ele_tag} ' + ' '.join([str(x) for x in self.ele2node_tags[ele_tag]]) + '\n' for ele_tag in self.ele2node_tags]
        open(self.cache_path + 'ele2node_tags.txt', 'w').writelines(ostr)
        if self.selected_node_tags is not None:
            self.savetxt(self.cache_path + 'selected_node_tags.txt', self.selected_node_tags, fmt='%i')

        for i, fname in enumerate(self.meta_files):
            vals = getattr(self, fname)
            if vals is not None:
                self.savetxt(self.cache_path + f'{fname}.txt', vals, fmt=self.meta_fmt[i])
        if self.dynamic:
            if not self.used_r_starter:
                self.savetxt(self.cache_path + 'x_disp.txt', self.x_disp)
                self.savetxt(self.cache_path + 'y_disp.txt', self.y_disp)
                self.savetxt(self.cache_path + 'timer.txt', self.time)
            elif self.pseudo_dt:
                from numpy import arange
                x_disp = self.loadtxt(f'{self.cache_path}x_disp.txt', ndmin=2)
                self.time = arange(len(x_disp[:, 0])) * self.pseudo_dt
                self.savetxt(self.cache_path + 'timer.txt', self.time)

    def load_from_cache(self):
        self.coords = self.loadtxt(self.cache_path + 'coords.txt')

        try:
            self.selected_node_tags = self.loadtxt(self.cache_path + 'selected_node_tags.txt', dtype=int)
        except OSError:
            pass

        self.ele2node_tags = {}
        lines = open(self.cache_path + 'ele2node_tags.txt').read().splitlines()
        for line in lines:
            parts = [int(x) for x in line.split()]
            self.ele2node_tags[parts[0]] = parts[1:]
        for fname in self.meta_files:
            try:
                data = self.loadtxt(self.cache_path + f'{fname}.txt')
                if len(data) == 0:
                    data = None
                setattr(self, fname, data)
            except OSError:
                pass

        if self.dynamic:
            self.x_disp = self.loadtxt(f'{self.cache_path}x_disp.txt')
            self.y_disp = self.loadtxt(f'{self.cache_path}y_disp.txt')
            self.time = self.loadtxt(f'{self.cache_path}timer.txt', ndmin=2)[:, 0]

    @property
    def dt(self):
        if self._dt is None:
            if self.time is not None:
                self._dt = (self.time[-1] - self.time[0]) / (len(self.time) - 1)
        return self._dt

    @dt.setter
    def dt(self, dt):
        self._dt = dt

    def rezero_node_tags(self, osi=None):
        from numpy import array, arange, searchsorted, where
        if self.selected_node_tags is None:
            node_tags = array(o3.get_node_tags(osi))
        else:
            node_tags = self.selected_node_tags
        new_node_tags = arange(1, len(node_tags) + 1)
        sidx = node_tags.argsort()
        k = node_tags[sidx]
        v = new_node_tags[sidx]
        for ele_tag in self.ele2node_tags:
            curr_tags = self.ele2node_tags[ele_tag]
            idx = searchsorted(k, curr_tags)
            assert max(idx) < len(k)
            mask = k[idx] == curr_tags
            self.ele2node_tags[ele_tag] = where(mask, v[idx], len(k))
