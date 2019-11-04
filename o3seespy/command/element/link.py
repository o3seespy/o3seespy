from o3seespy.command.element.base_element import ElementBase


class TwoNodeLinkdir(ElementBase):

    def __init__(self, osi, ele_nodes, mat_tags=None, dirs):
        self.ele_nodes = ele_nodes
        self.mat_tags = mat_tags
        self.dirs = dirs
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, '-dir', *self.dirs]
        if getattr(self, 'mat_tags') is not None:
            self._parameters += ['-mat', *self.mat_tags]
        self.to_process(osi)

class TwoNodeLinkorient(ElementBase):

    def __init__(self, osi, ele_nodes, mat_tags=None, vecx, vecy, p_delta_vals=None):
        self.ele_nodes = ele_nodes
        self.mat_tags = mat_tags
        self.vecx = vecx
        self.vecy = vecy
        self.p_delta_vals = p_delta_vals
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, '-orient', *self.vecx, *self.vecy]
        if getattr(self, 'mat_tags') is not None:
            self._parameters += ['-mat', *self.mat_tags]
        if getattr(self, 'p_delta_vals') is not None:
            self._parameters += ['-pDelta', *self.p_delta_vals]
        self.to_process(osi)

class TwoNodeLinkshearDist(ElementBase):

    def __init__(self, osi, ele_nodes, mat_tags=None, s_dratios, doRayleigh=False):
        self.ele_nodes = ele_nodes
        self.mat_tags = mat_tags
        self.s_dratios = s_dratios
        self.do_rayleigh = do_rayleigh
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, '-shearDist', *self.s_dratios]
        if getattr(self, 'mat_tags') is not None:
            self._parameters += ['-mat', *self.mat_tags]
        if getattr(self, 'do_rayleigh') is not None:
            self._parameters += ['-do_rayleigh']
        self.to_process(osi)

class TwoNodeLinkmass(ElementBase):

    def __init__(self, osi, ele_nodes, mat_tags=None, m):
        self.ele_nodes = ele_nodes
        self.mat_tags = mat_tags
        self.m = float(m)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, '-mass', self.m]
        if getattr(self, 'mat_tags') is not None:
            self._parameters += ['-mat', *self.mat_tags]
        self.to_process(osi)
