from o3seespy.base_model import OpenseesObject


class SystemBase(OpenseesObject):
    op_base_type = "system"


class ProfileSPD(SystemBase):
    op_type = "ProfileSPD"

    def __init__(self, osi):

        self._parameters = [self.op_type]
        self.to_process(osi)


class SparseGeneral(SystemBase):
    op_type = "SparseGeneral"

    def __init__(self, osi):

        self._parameters = [self.op_type]
        self.to_process(osi)