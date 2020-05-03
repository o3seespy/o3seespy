from o3seespy.base_model import OpenSeesObject


class AlgorithmBase(OpenSeesObject):
    op_base_type = "algorithm"


class Linear(AlgorithmBase):
    op_type = "Linear"

    def __init__(self, osi, secant=False, initial=False, factor_once=False):
        self.osi = osi
        self.secant = secant
        self.initial = initial
        self.factor_once = factor_once
        self._parameters = [self.op_type, self.secant, self.initial, self.factor_once]
        self.to_process(osi)


class Newton(AlgorithmBase):
    op_type = "Newton"

    def __init__(self, osi, secant=False, initial=False, initial_then_current=False):
        self.osi = osi
        self.secant = secant
        self.initial = initial
        self.initial_then_current = initial_then_current
        self._parameters = [self.op_type, self.secant, self.initial, self.initial_then_current]
        self.to_process(osi)
