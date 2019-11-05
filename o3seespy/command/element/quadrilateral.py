from o3seespy.command.element.base_element import ElementBase


class Quad(ElementBase):
    op_type = 'quad'

    def __init__(self, osi, ele_nodes, thick, type, mat, pressure, rho, b1, b2):
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.thick = float(thick)
        self.type = type
        self.mat = mat
        self.pressure = float(pressure)
        self.rho = float(rho)
        self.b1 = float(b1)
        self.b2 = float(b2)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.thick, self.type, self.mat.tag, self.pressure, self.rho, self.b1, self.b2]
        self.to_process(osi)


class ShellMITC4(ElementBase):
    op_type = 'ShellMITC4'

    def __init__(self, osi, ele_nodes, sec):
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.sec = sec
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.sec.tag]
        self.to_process(osi)


class ShellDKGQ(ElementBase):
    op_type = 'ShellDKGQ'

    def __init__(self, osi, ele_nodes, sec):
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.sec = sec
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.sec.tag]
        self.to_process(osi)


class ShellDKGT(ElementBase):
    op_type = 'ShellDKGT'

    def __init__(self, osi, ele_nodes, sec):
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.sec = sec
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.sec.tag]
        self.to_process(osi)


class ShellNLDKGQ(ElementBase):
    op_type = 'ShellNLDKGQ'

    def __init__(self, osi, ele_nodes, sec):
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.sec = sec
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.sec.tag]
        self.to_process(osi)


class ShellNLDKGT(ElementBase):
    op_type = 'ShellNLDKGT'

    def __init__(self, osi, ele_nodes, sec):
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.sec = sec
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.sec.tag]
        self.to_process(osi)


class ShellNL(ElementBase):
    op_type = 'ShellNL'

    def __init__(self, osi, ele_nodes, sec):
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.sec = sec
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.sec.tag]
        self.to_process(osi)


class BbarQuad(ElementBase):
    op_type = 'bbarQuad'

    def __init__(self, osi, ele_nodes, thick, mat):
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.thick = float(thick)
        self.mat = mat
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.thick, self.mat.tag]
        self.to_process(osi)


class EnhancedQuad(ElementBase):
    op_type = 'enhancedQuad'

    def __init__(self, osi, ele_nodes, thick, type, mat):
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.thick = float(thick)
        self.type = type
        self.mat = mat
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.thick, self.type, self.mat.tag]
        self.to_process(osi)


class SSPquad(ElementBase):
    op_type = 'SSPquad'

    def __init__(self, osi, ele_nodes, mat, type, thick, b1, b2):
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.mat = mat
        self.type = type
        self.thick = float(thick)
        self.b1 = float(b1)
        self.b2 = float(b2)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mat.tag, self.type, self.thick, self.b1, self.b2]
        self.to_process(osi)
