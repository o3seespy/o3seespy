from o3seespy.command.element.base_element import ElementBase


class CatenaryCable(ElementBase):
    op_type = 'CatenaryCable'

    def __init__(self, osi, i_node, j_node, weight, big_e, big_a, l0, alpha, temperature_change, rho, error_tol, nsubsteps, mass_type):
        self.i_node = i_node
        self.j_node = j_node
        self.weight = float(weight)
        self.big_e = float(big_e)
        self.big_a = float(big_a)
        self.l0 = float(l0)
        self.alpha = float(alpha)
        self.temperature_change = float(temperature_change)
        self.rho = float(rho)
        self.error_tol = float(error_tol)
        self.nsubsteps = int(nsubsteps)
        self.mass_type = int(mass_type)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.i_node.tag, self.j_node.tag, self.weight, self.big_e, self.big_a, self.l0, self.alpha, self.temperature_change, self.rho, self.error_tol, self.nsubsteps, self.mass_type]
        self.to_process(osi)
