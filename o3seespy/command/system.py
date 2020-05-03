from o3seespy.base_model import OpenSeesObject


class SystemBase(OpenSeesObject):
    op_base_type = "system"


class ProfileSPD(SystemBase):
    op_type = "ProfileSPD"

    def __init__(self, osi):
        self.osi = osi
        self._parameters = [self.op_type]
        self.to_process(osi)


class SparseGeneral(SystemBase):
    op_type = "SparseGeneral"

    def __init__(self, osi):
        self.osi = osi
        self._parameters = [self.op_type]
        self.to_process(osi)


class BandGeneral(SystemBase):
    op_type = "BandGeneral"

    def __init__(self, osi):
        self.osi = osi
        self._parameters = [self.op_type]
        self.to_process(osi)


class FullGeneral(SystemBase):
    op_type = "FullGeneral"

    def __init__(self, osi):
        self.osi = osi
        self._parameters = [self.op_type]
        self.to_process(osi)
