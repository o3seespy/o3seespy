from o3seespy.command.element.base_element import ElementBase


class PFEMElementBubble(ElementBase):
    op_type = 'PFEMElementBubble'

    def __init__(self, osi, nd1, nd2, nd3, nd4, rho, mu, b1, b2, b3, thickness, kappa):
        self.nd1 = int(nd1)
        self.nd2 = int(nd2)
        self.nd3 = int(nd3)
        self.nd4 = int(nd4)
        self.rho = float(rho)
        self.mu = float(mu)
        self.b1 = float(b1)
        self.b2 = float(b2)
        self.b3 = float(b3)
        self.thickness = float(thickness)
        self.kappa = float(kappa)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.nd1, self.nd2, self.nd3, self.nd4, self.rho, self.mu, self.b1, self.b2, self.b3, self.thickness, self.kappa]
        self.to_process(osi)


class PFEMElementCompressible(ElementBase):
    op_type = 'PFEMElementCompressible'

    def __init__(self, osi, nd1, nd2, nd3, nd4, rho, mu, b1, b2, thickness, kappa):
        self.nd1 = int(nd1)
        self.nd2 = int(nd2)
        self.nd3 = int(nd3)
        self.nd4 = int(nd4)
        self.rho = float(rho)
        self.mu = float(mu)
        self.b1 = float(b1)
        self.b2 = float(b2)
        self.thickness = float(thickness)
        self.kappa = float(kappa)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.nd1, self.nd2, self.nd3, self.nd4, self.rho, self.mu, self.b1, self.b2, self.thickness, self.kappa]
        self.to_process(osi)
