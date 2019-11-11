from o3seespy.command.element.base_element import ElementBase


class TwoNodeLink(ElementBase):
    op_type = 'twoNodeLink'

    def __init__(self, osi, ele_nodes, mat_tags=None, dir=None, p_delta_vals=None, shear_dist=None, do_rayleigh=False, orient=None, mass=None):
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.mat_tags = mat_tags
        self.dir = dir
        self.p_delta_vals = p_delta_vals
        self.shear_dist = shear_dist
        self.do_rayleigh = do_rayleigh
        self.orient = orient
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes]
        if getattr(self, 'mat_tags') is not None:
            self._parameters += ['-mat', *self.mat_tags]
        if getattr(self, 'dir') is not None:
            self._parameters += ['-dir', *self.dir]
        if getattr(self, 'p_delta_vals') is not None:
            self._parameters += ['-pDelta', *self.p_delta_vals]
        if getattr(self, 'shear_dist') is not None:
            self._parameters += ['-shearDist', *self.shear_dist]
        if getattr(self, 'do_rayleigh'):
            self._parameters += ['-doRayleigh']
        if getattr(self, 'orient') is not None:
            self._parameters += ['--orient', *self.orient]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        self.to_process(osi)
