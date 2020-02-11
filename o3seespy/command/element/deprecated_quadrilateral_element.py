from o3seespy.command.element.base_element import ElementBase


class Quad(ElementBase):
    op_type = "quad"

    def __init__(self, osi, nodes, thick, mat_type, mat, pressure=0.0, rho=0.0, b1=0.0, b2=0.0):
        self.nodes = nodes
        self.thick = thick
        self.mat = mat
        self.mat_type = mat_type
        self.pressure = pressure
        self.rho = rho
        self.b1 = b1
        self.b2 = b2
        osi.n_ele += 1
        self._tag = osi.n_ele
        node_tags = [x.tag for x in self.nodes]
        self._parameters = [self.op_type, self._tag, *node_tags, self.thick, self.mat_type, self.mat.tag, self.pressure,
                            self.rho, self.b1, self.b2]
        # if self.pressure is not None:
        #     self._parameters.append(self.pressure)
        # if self.rho is not None:
        #     self._parameters.append(self.rho)
        # if self.b1 is not None:
        #     self._parameters.append(self.b1)
        # if self.b2 is not None:
        #     self._parameters.append(self.b2)

        self.to_process(osi)
