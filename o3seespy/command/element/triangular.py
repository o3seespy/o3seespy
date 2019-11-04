from o3seespy.command.element.base_element import ElementBase


class Tri31(ElementBase):
    op_type = 'Tri31'

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
