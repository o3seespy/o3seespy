from o3seespy.command.element.base_element import ElementBase


class StdBrick(ElementBase):
    op_type = 'stdBrick'

    def __init__(self, osi, ele_nodes, mat, b1, b2, b3):
        self.ele_nodes = ele_nodes
        self.mat = mat.tag
        self.b1 = float(b1)
        self.b2 = float(b2)
        self.b3 = float(b3)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mat.tag, self.b1, self.b2, self.b3]
        self.to_process(osi)


class BbarBrick(ElementBase):
    op_type = 'bbarBrick'

    def __init__(self, osi, ele_nodes, mat, b1, b2, b3):
        self.ele_nodes = ele_nodes
        self.mat = mat.tag
        self.b1 = float(b1)
        self.b2 = float(b2)
        self.b3 = float(b3)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mat.tag, self.b1, self.b2, self.b3]
        self.to_process(osi)


class Brick20N(ElementBase):
    op_type = 'Brick20N'

    def __init__(self, osi, ele_nodes, mat, bf1, bf2, bf3, mass_den):
        self.ele_nodes = ele_nodes
        self.mat = mat.tag
        self.bf1 = float(bf1)
        self.bf2 = float(bf2)
        self.bf3 = float(bf3)
        self.mass_den = float(mass_den)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mat.tag, self.bf1, self.bf2, self.bf3, self.mass_den]
        self.to_process(osi)


class SSPbrick(ElementBase):
    op_type = 'SSPbrick'

    def __init__(self, osi, ele_nodes, mat, b1, b2, b3):
        self.ele_nodes = ele_nodes
        self.mat = mat.tag
        self.b1 = float(b1)
        self.b2 = float(b2)
        self.b3 = float(b3)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mat.tag, self.b1, self.b2, self.b3]
        self.to_process(osi)
