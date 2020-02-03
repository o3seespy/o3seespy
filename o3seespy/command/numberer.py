from o3seespy.base_model import OpenSeesObject


class NumbererBase(OpenSeesObject):
    op_base_type = "numberer"


class RCM(NumbererBase):
    op_type = "RCM"

    def __init__(self, osi):

        self._parameters = [self.op_type]
        self.to_process(osi)
