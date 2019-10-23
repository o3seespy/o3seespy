from o3seespy.base_model import OpenseesObject


class IntegratorBase(OpenseesObject):
    op_base_type = "integrator"


class Newmark(IntegratorBase):
    op_type = "Newmark"

    def __init__(self, osi, gamma=False, beta=False):
        self.gamma = gamma
        self.beta = beta
        self._parameters = [self.op_type, self.gamma, self.beta]
        self.to_process(osi)
