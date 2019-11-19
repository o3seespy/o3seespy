from o3seespy.base_model import OpenseesObject


class BeamIntegrationBase(OpenseesObject):
    op_base_type = "beamIntegration"


class UserHinge(BeamIntegrationBase):
    op_type = 'UserHinge'

    def __init__(self, osi, sec_e, np_l, secs_l, locs_l, wts_l, np_r, secs_r, locs_r, wts_r):
        self.sec_e = sec_e
        self.np_l = int(np_l)
        self.secs_l = secs_l
        self.locs_l = locs_l
        self.wts_l = wts_l
        self.np_r = int(np_r)
        self.secs_r = secs_r
        self.locs_r = locs_r
        self.wts_r = wts_r
        osi.n_integ += 1
        self._tag = osi.n_integ
        self._parameters = [self.op_type, self._tag, self.sec_e.tag, self.np_l, *self.secs_l, *self.locs_l, *self.wts_l, self.np_r, *self.secs_r, *self.locs_r, *self.wts_r]
        self.to_process(osi)


class HingeMidpoint(BeamIntegrationBase):
    op_type = 'HingeMidpoint'

    def __init__(self, osi, sec_i, lp_i, sec_j, lp_j, sec_e):
        self.sec_i = sec_i
        self.lp_i = float(lp_i)
        self.sec_j = sec_j
        self.lp_j = float(lp_j)
        self.sec_e = sec_e
        osi.n_integ += 1
        self._tag = osi.n_integ
        self._parameters = [self.op_type, self._tag, self.sec_i.tag, self.lp_i, self.sec_j.tag, self.lp_j, self.sec_e.tag]
        self.to_process(osi)


class HingeRadau(BeamIntegrationBase):
    op_type = 'HingeRadau'

    def __init__(self, osi, sec_i, lp_i, sec_j, lp_j, sec_e):
        self.sec_i = sec_i
        self.lp_i = lp_i
        self.sec_j = sec_j
        self.lp_j = lp_j
        self.sec_e = sec_e
        osi.n_integ += 1
        self._tag = osi.n_integ
        self._parameters = [self.op_type, self._tag, self.sec_i.tag, self.lp_i, self.sec_j.tag, self.lp_j, self.sec_e.tag]
        self.to_process(osi)


class HingeRadauTwo(BeamIntegrationBase):
    op_type = 'HingeRadauTwo'

    def __init__(self, osi, sec_i, lp_i, sec_j, lp_j, sec_e):
        self.sec_i = sec_i
        self.lp_i = lp_i
        self.sec_j = sec_j
        self.lp_j = lp_j
        self.sec_e = sec_e
        osi.n_integ += 1
        self._tag = osi.n_integ
        self._parameters = [self.op_type, self._tag, self.sec_i.tag, self.lp_i, self.sec_j.tag, self.lp_j, self.sec_e.tag]
        self.to_process(osi)


class BeamhingeEndpoint(OpenseesObject):
    base_type = 'beamhingeEndpoint'

    def __init__(self, osi, lp_i, sec_j, lp_j, sec_e):
        self.lp_i = float(lp_i)
        self.sec_j = sec_j
        self.lp_j = float(lp_j)
        self.sec_e = sec_e
        osi.n_integ += 1
        self._tag = osi.n_integ
        self._parameters = [self._tag, self.lp_i, self.sec_j.tag, self.lp_j, self.sec_e.tag]
        self.to_process(osi)
