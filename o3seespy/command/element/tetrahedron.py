from o3seespy.command.element.base_element import ElementBase


class FourNodeTetrahedron(ElementBase):
    op_type = 'FourNodeTetrahedron'

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
