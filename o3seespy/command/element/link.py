from o3seespy.command.element.base_element import ElementBase


class TwoNodeLinkdir(ElementBase):
    op_type = 'twoNodeLink'

    def __init__(self, osi, ele_nodes, mat_tags=None, dirs, mass=None):
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.mat_tags = mat_tags
        self.dirs = dirs
        self.mass = float(mass)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, '-dir', *self.dirs]
        if getattr(self, 'mat_tags') is not None:
            self._parameters += ['-mat', *self.mat_tags]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        self.to_process(osi)

class TwoNodeLinkorient(ElementBase):
    op_type = 'twoNodeLink'

    def __init__(self, osi, ele_nodes, mat_tags=None, vecx, vecy, p_delta_vals=None, mass=None):
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.mat_tags = mat_tags
        self.vecx = vecx
        self.vecy = vecy
        self.p_delta_vals = p_delta_vals
        self.mass = float(mass)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, '-orient', *self.vecx, *self.vecy]
        if getattr(self, 'mat_tags') is not None:
            self._parameters += ['-mat', *self.mat_tags]
        if getattr(self, 'p_delta_vals') is not None:
            self._parameters += ['-pDelta', *self.p_delta_vals]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        self.to_process(osi)

class TwoNodeLinkshearDist(ElementBase):
    op_type = 'twoNodeLink'

    def __init__(self, osi, ele_nodes, mat_tags=None, s_dratios, do_rayleigh=False, mass=None):
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.mat_tags = mat_tags
        self.s_dratios = s_dratios
        self.do_rayleigh = do_rayleigh
        self.mass = float(mass)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, '-shearDist', *self.s_dratios]
        if getattr(self, 'mat_tags') is not None:
            self._parameters += ['-mat', *self.mat_tags]
        if getattr(self, 'do_rayleigh'):
            self._parameters += ['-doRayleigh']
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        self.to_process(osi)
