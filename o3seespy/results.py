import os
import o3seespy as o3
import glob


def calc_quad_centroids(xs, ys, axis=-1):  # use axis=-2 for time domain data, since -1 is time
    import numpy as np
    x0 = np.array(xs)
    y0 = np.array(ys)
    x1 = np.roll(xs, 1, axis=axis)
    y1 = np.roll(ys, 1, axis=axis)
    a = x0 * y1 - x1 * y0
    xc = np.sum((x0 + x1) * a, axis=axis)
    yc = np.sum((y0 + y1) * a, axis=axis)

    area = 0.5 * np.sum(a, axis=axis)
    xc /= (6.0 * area)
    yc /= (6.0 * area)

    return xc, yc


class Results(object):
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
    dyn_files = None
    ndm = None

    def __init__(self, cache_path='', dt=None, dynamic=False, man_nodes=False, prefix='', mp=0):
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
        self.prefix = prefix
        self.mp = mp
        self.dyn_files = ['timer']
        if mp:
            self.pstr = f'-pid{o3.mp.get_pid()}'
        else:
            self.pstr = ''

    def start_recorders(self, osi, dt=None):  # TODO: handle recorder time step
        self.used_r_starter = 1
        if self.man_nodes and self.selected_node_tags is None:
            raise ValueError('if man_nodes then selected_node_tags must be manually set')
            # self._selected_node_tags = o3.get_node_tags(osi)
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
            o3.recorder.NodesToFile(osi, f'{self.cache_path}{self.prefix}x_disp{self.pstr}.txt', node_tags,
                                    [o3.cc.X], 'disp', nsd=5, dt=dt, nodes_as_tags=True)
            o3.recorder.NodesToFile(osi, f'{self.cache_path}{self.prefix}y_disp{self.pstr}.txt', node_tags,
                                    [o3.cc.Y], 'disp', nsd=5, dt=dt, nodes_as_tags=True)
            if self.ndm == 3:
                o3.recorder.NodesToFile(osi, f'{self.cache_path}{self.prefix}z_disp{self.pstr}.txt', node_tags,
                                        [o3.cc.DOF3D_Z], 'disp', nsd=5, dt=dt, nodes_as_tags=True)
            if not self.pseudo_dt:
                if self.mp:
                    if self.pstr == '-pid0':
                        all_node_tags = o3.get_node_tags(osi)
                        o3.recorder.TimeToFile(osi, f'{self.cache_path}{self.prefix}timer.txt', nsd=5, dt=dt,
                                               dummy_node_tag=all_node_tags[0])
                else:
                    o3.recorder.TimeToFile(osi, f'{self.cache_path}{self.prefix}timer.txt', nsd=5, dt=dt)

    def wipe_old_files(self):
        sfiles = ['ele2node_tags', 'coords', 'selected_node_tags']
        if not self.used_r_starter:
            sfiles += self.dyn_files
        for fname in sfiles:
            ffps = glob.glob(f'{self.cache_path}{self.prefix}{fname}*.txt')
            for ffp in ffps:
                try:
                    os.remove(ffp)
                except FileNotFoundError:
                    pass

        for fname in self.meta_files:
            ffps = glob.glob(f'{self.cache_path}{self.prefix}{fname}*.txt')
            for ffp in ffps:
                try:
                    os.remove(ffp)
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

    def save_to_cache(self, wipe_old_files=1):
        pf = self.prefix
        if wipe_old_files:
            self.wipe_old_files()
        if self.coords is not None:
            self.savetxt(self.cache_path + f'{pf}coords{self.pstr}.txt', self.coords, fmt='%.6g')
        ostr = [f'{ele_tag} ' + ' '.join([str(x) for x in self.ele2node_tags[ele_tag]]) + '\n' for ele_tag in self.ele2node_tags]
        open(self.cache_path + f'{pf}ele2node_tags{self.pstr}.txt', 'w').writelines(ostr)

        if self.selected_node_tags is not None and len(self.selected_node_tags):
            self.savetxt(self.cache_path + f'{pf}selected_node_tags{self.pstr}.txt', self.selected_node_tags, fmt='%i')


        for i, fname in enumerate(self.meta_files):
            vals = getattr(self, fname)
            if vals is not None:
                self.savetxt(self.cache_path + f'{pf}{fname}.txt', vals, fmt=self.meta_fmt[i])
        if self.dynamic:
            if not self.used_r_starter:
                for fname in self.dyn_files:
                    self.savetxt(self.cache_path + f'{pf}{fname}.txt', getattr(self, fname), fmt='%.5g')
                # self.savetxt(self.cache_path + f'{pf}y_disp.txt', self.y_disp, fmt='%.5g')
                # self.savetxt(self.cache_path + f'{pf}timer.txt', self.time, fmt='%.5g')
            elif self.pseudo_dt:
                from numpy import arange  # TODO: support mp, where x_disp-pid*
                x_disp = self.loadtxt(f'{self.cache_path}{pf}x_disp.txt', ndmin=2)
                self.time = arange(len(x_disp[:, 0])) * self.pseudo_dt
                self.savetxt(self.cache_path + f'{pf}timer.txt', self.time, fmt='%.5g')

    def load_from_cache(self):
        prefix = self.prefix
        from numpy import concatenate, unique
        if self.mp:

            ffps = glob.glob(self.cache_path + f'{prefix}coords*pid*.txt')  # old files used -pid
            ffps.sort()
            temp_c0 = []
            for ffp in ffps:
                temp_c0.append(self.loadtxt(ffp))
            self.coords = concatenate(temp_c0)
        else:
            self.coords = self.loadtxt(self.cache_path + f'{prefix}coords.txt')
        ffps = glob.glob(self.cache_path + f'{prefix}selected_node_tags*.txt')
        ffps.sort()
        # try:
        uindices = None
        temp_c = []
        for ffp in ffps:
            temp_c.append(self.loadtxt(ffp, dtype=int))
        if len(temp_c):
            snts = concatenate(temp_c)
            uniq, uindices = unique(snts, return_index=True)
            self.selected_node_tags = uniq
            self.coords = self.coords[uindices]
        # except OSError as e:
        #     pass
        # except ValueError as e:  # if empty list
        #     pass


        self.ele2node_tags = {}
        self.ele2pid = {}
        ffps = glob.glob(self.cache_path + f'{prefix}ele2node_tags*.txt')
        ffps.sort()
        for ffp in ffps:
            if 'pid' in ffp:
                pid = ffp.split('pid')[-1]
                pid = int(pid.split('.')[0])
            else:
                pid = 0
            lines = open(ffp).read().splitlines()
            for line in lines:
                parts = [int(x) for x in line.split()]
                self.ele2node_tags[parts[0]] = parts[1:]
                self.ele2pid[parts[0]] = pid
        for fname in self.meta_files:
            try:
                data = self.loadtxt(self.cache_path + f'{prefix}{fname}.txt')
                if len(data) == 0:
                    data = None
                setattr(self, fname, data)
            except OSError:
                pass

        if self.dynamic:
            if self.mp:
                for pp in self.dyn_files:
                    if pp == 'timer':
                        continue
                    ffps = glob.glob(self.cache_path + f'{prefix}{pp}*.txt')
                    ffps.sort()
                    temp = []
                    for ffp in ffps:
                        temp.append(self.loadtxt(ffp, ndmin=2).T)
                    if uindices is not None:
                        if len(temp) == 0:
                            print('not found: ', prefix, pp)
                        all_vals = concatenate(temp)
                        vals = all_vals[uindices]
                    else:
                        vals = concatenate(temp)
                    setattr(self, pp, vals.T)
            else:
                for fname in self.dyn_files:
                    setattr(self, fname, self.loadtxt(f'{self.cache_path}{prefix}{fname}.txt'))

                # self.y_disp = self.loadtxt(f'{self.cache_path}{prefix}y_disp.txt')
            self.time = self.loadtxt(f'{self.cache_path}{prefix}timer.txt', ndmin=2)[:, 0]

    @property
    def dt(self):
        if self._dt is None:
            if self.time is not None:
                self._dt = (self.time[-1] - self.time[0]) / (len(self.time) - 1)
        return self._dt

    @dt.setter
    def dt(self, dt):
        self._dt = dt

    def rezero_node_tags(self, osi=None, first=1):
        from numpy import array, arange, searchsorted, where, clip
        if self.selected_node_tags is None:
            node_tags = array(o3.get_node_tags(osi))
        else:
            node_tags = self.selected_node_tags
        new_node_tags = arange(first, len(node_tags) + first)
        sidx = node_tags.argsort()
        k = node_tags[sidx]
        v = new_node_tags[sidx]
        for ele_tag in self.ele2node_tags:
            curr_tags = self.ele2node_tags[ele_tag]
            idx = searchsorted(k, curr_tags)
            assert max(idx) < len(k)
            # idx = clip(idx, None, len(k) - 1)
            mask = k[idx] == curr_tags
            self.ele2node_tags[ele_tag] = where(mask, v[idx], len(k))


    def get_eles_by_n_nodes(self, n_nodes):
        eles_by_n_nodes = {2: [], 4: [], 8: []}  # TODO: extend to all 5 node and 9
        for ele in self.ele2node_tags:
            nn = len(self.ele2node_tags[ele])
            eles_by_n_nodes[nn].append(ele)
        return eles_by_n_nodes[n_nodes]


class Results2D(Results):

    def __init__(self, cache_path='', dt=None, dynamic=False, man_nodes=False, prefix='', mp=0):
        super(Results2D, self).__init__(cache_path=cache_path, dt=dt, dynamic=dynamic, man_nodes=man_nodes,
                                        prefix=prefix, mp=mp)
        self.dyn_files += ['x_disp', 'y_disp']
        self.ndm = 2

    def compute_ele_strains_and_disps(self):  # currently only available for quad elements
        import numpy as np
        rd = {}
        eles = self.get_eles_by_n_nodes(4)
        nodes = np.array([self.ele2node_tags[ele] for ele in eles]) - 1
        xd = self.x_disp[:, nodes].transpose(1, 2, 0)
        yd = self.y_disp[:, nodes].transpose(1, 2, 0)
        xc = self.coords[nodes, 0]
        yc = self.coords[nodes, 1]
        x_disp_ele, y_disp_ele = calc_quad_centroids(xc[:, :, np.newaxis] + xd, yc[:, :, np.newaxis] + yd, axis=-2)

        rd['XDISP'] = x_disp_ele - x_disp_ele[:, 0][:, np.newaxis]
        rd['YDISP'] = y_disp_ele - y_disp_ele[:, 0][:, np.newaxis]

        # return
        # nodes must be anti-clockwise
        xc0 = xc[0]
        yc0 = yc[0]
        i = 0
        for i in range(4):
            if xc0[i%4] < xc0[(i+1)%4] and xc0[(i+2)%4] > xc0[(i+3)%4] and yc0[(i+1)%4] < yc0[(i+2)%4]:  # i=bottom-left
                break
            if i == 3:
                print('WARNING could not find bottom left')
        inds = np.roll(np.arange(4), i)
        xc = xc[:, inds]
        yc = yc[:, inds]
        xd = xd[:, inds]
        yd = yd[:, inds]

        xlen = (xc[:, 1] - xc[:, 0] + xc[:, 2] - xc[:, 3]) / 2
        xdelta = (xd[:, 1] - xd[:, 0] + xd[:, 2] - xd[:, 3]) / 2
        rd['EPS_XX'] = xdelta / xlen[:, np.newaxis]
        ylen = (yc[:, 2] - yc[:, 1] + yc[:, 3] - yc[:, 0]) / 2
        ydelta = (yd[:, 2] - yd[:, 1] + yd[:, 3] - yd[:, 0]) / 2
        rd['EPS_YY'] = ydelta / ylen[:, np.newaxis]

        xd_bot = (xd[:, 1] + xd[:, 0]) / 2
        xd_top = (xd[:, 2] + xd[:, 3]) / 2
        yd_lhs = (yd[:, 3] + yd[:, 0]) / 2
        yd_rhs = (yd[:, 2] + yd[:, 1]) / 2
        rd['EPS_XY'] = ((xd_bot - xd_top) / ylen[:, np.newaxis] + (yd_rhs - yd_lhs) / xlen[:, np.newaxis]) / 2  # +ve in anti-clockwise
        rd['SSTR_MAX'] = np.sqrt(((rd['EPS_XX'] - rd['EPS_YY']) / 2) ** 2 + rd['EPS_XY'] ** 2)
        rd['EPS_VOL'] = rd['EPS_XX'] + rd['EPS_YY']
        return rd


class Results3D(Results):
    z_disp = None

    def __init__(self, cache_path='', dt=None, dynamic=False, man_nodes=False, prefix='', mp=0):
        super(Results3D, self).__init__(cache_path=cache_path, dt=dt, dynamic=dynamic, man_nodes=man_nodes,
                                        prefix=prefix, mp=mp)
        self.dyn_files += ['x_disp', 'y_disp', 'z_disp']
        self.ndm = 3


def rezero_node_tags(node_tags, ele2node_tags, first=0):
    from numpy import array, arange, searchsorted, where, clip

    new_node_tags = arange(first, len(node_tags) + first)
    sidx = node_tags.argsort()
    k = node_tags[sidx]
    v = new_node_tags[sidx]
    for ele_tag in ele2node_tags:
        curr_tags = ele2node_tags[ele_tag]
        idx = searchsorted(k, curr_tags)
        assert max(idx) < len(k)
        mask = k[idx] == curr_tags
        ele2node_tags[ele_tag] = where(mask, v[idx], len(k))
    return ele2node_tags