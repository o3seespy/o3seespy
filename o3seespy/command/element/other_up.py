from o3seespy.command.element.base_element import ElementBase


class SSPquadUP(ElementBase):
    op_type = 'SSPquadUP'

    def __init__(self, osi, ele_nodes, mat, thick, f_bulk, f_den, k1, k2, void, alpha, b1, b2):
        self.ele_nodes = ele_nodes
        self.mat = mat.tag
        self.thick = float(thick)
        self.f_bulk = float(f_bulk)
        self.f_den = float(f_den)
        self.k1 = float(k1)
        self.k2 = float(k2)
        self.void = float(void)
        self.alpha = float(alpha)
        self.b1 = float(b1)
        self.b2 = float(b2)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mat.tag, self.thick, self.f_bulk, self.f_den, self.k1, self.k2, self.void, self.alpha, self.b1, self.b2]
        self.to_process(osi)


class SSPbrickUP(ElementBase):
    op_type = 'SSPbrickUP'

    def __init__(self, osi, ele_nodes, mat, f_bulk, f_den, k1, k2, k3, void, alpha, b1, b2, b3):
        self.ele_nodes = ele_nodes
        self.mat = mat.tag
        self.f_bulk = float(f_bulk)
        self.f_den = float(f_den)
        self.k1 = float(k1)
        self.k2 = float(k2)
        self.k3 = float(k3)
        self.void = float(void)
        self.alpha = float(alpha)
        self.b1 = float(b1)
        self.b2 = float(b2)
        self.b3 = float(b3)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mat.tag, self.f_bulk, self.f_den, self.k1, self.k2, self.k3, self.void, self.alpha, self.b1, self.b2, self.b3]
        self.to_process(osi)
