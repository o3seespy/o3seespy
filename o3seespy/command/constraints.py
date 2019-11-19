from o3seespy.base_model import OpenseesObject


class ConstraintsBase(OpenseesObject):
    op_base_type = "constraints"


class Plain(ConstraintsBase):
    op_type = 'Plain'

    def __init__(self, osi):
        self._parameters = [self.op_type]
        self.to_process(osi)


class Lagrange(ConstraintsBase):
    op_type = 'Lagrange'

    def __init__(self, osi, alpha_m=1.0):
        self.alpha_m = float(alpha_m)
        self._parameters = [self.op_type, self.alpha_m]
        self.to_process(osi)


class Penalty(ConstraintsBase):
    op_type = 'Penalty'

    def __init__(self, osi, alpha_m=1.0):
        self.alpha_m = float(alpha_m)
        self._parameters = [self.op_type, self.alpha_m]
        self.to_process(osi)


class Transformation(ConstraintsBase):
    op_type = 'Transformation'

    def __init__(self, osi):
        self._parameters = [self.op_type]
        self.to_process(osi)
