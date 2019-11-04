from o3seespy.command.element.base_element import ElementBase


class SurfaceLoad(ElementBase):

    def __init__(self, osi, ele_nodes, p):
        self.ele_nodes = ele_nodes
        self.p = float(p)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.p]
        self.to_process(osi)


class VS3D4(ElementBase):

    def __init__(self, osi, ele_nodes, big_e, big_g, rho, big_r, alpha_n, alpha_t):
        self.ele_nodes = ele_nodes
        self.big_e = float(big_e)
        self.big_g = float(big_g)
        self.rho = float(rho)
        self.big_r = float(big_r)
        self.alpha_n = float(alpha_n)
        self.alpha_t = float(alpha_t)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.big_e, self.big_g, self.rho, self.big_r, self.alpha_n, self.alpha_t]
        self.to_process(osi)


class AC3D8(ElementBase):

    def __init__(self, osi, ele_nodes, mat):
        self.ele_nodes = ele_nodes
        self.mat = mat.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mat.tag]
        self.to_process(osi)


class ASI3D8(ElementBase):

    def __init__(self, osi, ele_nodes1, ele_nodes2):
        self.ele_nodes1 = ele_nodes1
        self.ele_nodes2 = ele_nodes2
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes1, *self.ele_nodes2]
        self.to_process(osi)


class AV3D4(ElementBase):

    def __init__(self, osi, ele_nodes, mat):
        self.ele_nodes = ele_nodes
        self.mat = mat.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mat.tag]
        self.to_process(osi)
