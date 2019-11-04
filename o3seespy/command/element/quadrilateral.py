from o3seespy.command.element.base_element import ElementBase


class Quad(ElementBase):

    def __init__(self, osi, ele_nodes, thick, type, mat, pressure, rho, b1, b2):
        self.ele_nodes = ele_nodes
        self.thick = float(thick)
        self.type = type
        self.mat = mat.tag
        self.pressure = float(pressure)
        self.rho = float(rho)
        self.b1 = float(b1)
        self.b2 = float(b2)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.thick, self.type, self.mat.tag, self.pressure, self.rho, self.b1, self.b2]
        self.to_process(osi)


class ShellMITC4(ElementBase):

    def __init__(self, osi, ele_nodes, sec):
        self.ele_nodes = ele_nodes
        self.sec = sec.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.sec.tag]
        self.to_process(osi)


class ShellDKGQ(ElementBase):

    def __init__(self, osi, ele_nodes, sec):
        self.ele_nodes = ele_nodes
        self.sec = sec.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.sec.tag]
        self.to_process(osi)


class ShellDKGT(ElementBase):

    def __init__(self, osi, ele_nodes, sec):
        self.ele_nodes = ele_nodes
        self.sec = sec.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.sec.tag]
        self.to_process(osi)


class ShellNLDKGQ(ElementBase):

    def __init__(self, osi, ele_nodes, sec):
        self.ele_nodes = ele_nodes
        self.sec = sec.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.sec.tag]
        self.to_process(osi)


class ShellNLDKGT(ElementBase):

    def __init__(self, osi, ele_nodes, sec):
        self.ele_nodes = ele_nodes
        self.sec = sec.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.sec.tag]
        self.to_process(osi)


class ShellNL(ElementBase):

    def __init__(self, osi, ele_nodes, sec):
        self.ele_nodes = ele_nodes
        self.sec = sec.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.sec.tag]
        self.to_process(osi)


class BbarQuad(ElementBase):

    def __init__(self, osi, ele_nodes, thick, mat):
        self.ele_nodes = ele_nodes
        self.thick = float(thick)
        self.mat = mat.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.thick, self.mat.tag]
        self.to_process(osi)


class EnhancedQuad(ElementBase):

    def __init__(self, osi, ele_nodes, thick, type, mat):
        self.ele_nodes = ele_nodes
        self.thick = float(thick)
        self.type = type
        self.mat = mat.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.thick, self.type, self.mat.tag]
        self.to_process(osi)


class SSPquad(ElementBase):

    def __init__(self, osi, ele_nodes, mat, type, thick, b1, b2):
        self.ele_nodes = ele_nodes
        self.mat = mat.tag
        self.type = type
        self.thick = float(thick)
        self.b1 = float(b1)
        self.b2 = float(b2)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mat.tag, self.type, self.thick, self.b1, self.b2]
        self.to_process(osi)
