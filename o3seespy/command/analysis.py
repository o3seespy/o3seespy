from o3seespy.base_model import OpenseesObject


class AnalysisBase(OpenseesObject):
    op_base_type = "analysis"


class Transient(AnalysisBase):
    op_type = "Transient"

    def __init__(self, osi):

        self._parameters = [self.op_type]
        self.to_process(osi)
