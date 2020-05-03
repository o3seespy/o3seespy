from o3seespy.base_model import OpenSeesObject


class AnalysisBase(OpenSeesObject):
    op_base_type = "analysis"


class Static(AnalysisBase):
    op_type = "Static"

    def __init__(self, osi):
        self.osi = osi
        self._parameters = [self.op_type]
        self.to_process(osi)


class Transient(AnalysisBase):
    op_type = "Transient"

    def __init__(self, osi):
        self.osi = osi
        self._parameters = [self.op_type]
        self.to_process(osi)


class VariableTransient(AnalysisBase):
    op_type = "VariableTransient"

    def __init__(self, osi):
        self.osi = osi
        self._parameters = [self.op_type]
        self.to_process(osi)


class PFEM(AnalysisBase):
    op_type = "PFEM"

    def __init__(self, osi, dt_max, dt_min, gravity, ratio=0.5):
        self.osi = osi
        self.dt_max = dt_max
        self.dt_min = dt_min
        self.gravity = gravity
        self.ratio = ratio

        self._parameters = [self.op_type, self.dt_max, self.dt_min, self.gravity, self.ratio]
        self.to_process(osi)
