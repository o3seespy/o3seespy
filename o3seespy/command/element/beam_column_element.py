from o3seespy.command.element.base_element import ElementBase


class ForceBeamColumn(ElementBase):
    op_type = "forceBeamColumn"

    def __init__(self, osi, node_i, node_j, transf, integration, max_inter=10, tol=1e-12, mass=0.0):
        self.node_i = node_i
        self.node_j = node_j
        self.transf = transf
        self.integration = integration
        self.max_inter = max_inter  # TODO: this should be int, but docs say 'float'
        self.tol = float(tol)
        self.mass = float(mass)
        osi.n_eles += 1
        self._tag = osi.n_eles

        self._parameters = [self.op_type, self._tag, *[self.node_i.tag, self.node_j.tag], self.transf.tag,
                            self.integration.tag, self.max_inter, self.tol, self.mass]
        self.to_process(osi)
