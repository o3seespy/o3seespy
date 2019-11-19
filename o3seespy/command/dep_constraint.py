from o3seespy.base_model import OpenseesObject


class ConstraintsBase(OpenseesObject):
    op_base_type = "constraints"


class Transformation(ConstraintsBase):
    op_type = "Transformation"

    def __init__(self, osi):

        self._parameters = [self.op_type]
        self.to_process(osi)
