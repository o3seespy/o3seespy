from o3seespy.base_model import OpenSeesObject


class GroundMotionBase(OpenSeesObject):
    op_base_type = "groundMotion"


class Plain(GroundMotionBase):
    op_type = 'Plain'

    def __init__(self, osi, mtype, ts, fact=1.0):
        """
        Initial method for Plain groundMotion

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        mtype: str
            motion type: 'disp', 'vel', 'accel'

        Examples
        --------

        """
        self.osi = osi
        osi.n_gms += 1
        self._tag = osi.n_gms
        self.mtype = mtype
        self.ts = ts
        self.fact = fact
        self._parameters = [self._tag, self.op_type, f'-{self.mtype}', self.ts.tag, '-fact', self.fact]
        self.to_process(osi)


class PlainAccel(GroundMotionBase):
    op_type = 'Plain'

    def __init__(self, osi, ts, fact=1.0):
        super(PlainAccel, self).__init__(osi, 'accel', ts=ts, fact=fact)


class PlainVel(GroundMotionBase):
    op_type = 'Plain'

    def __init__(self, osi, ts, fact=1.0):
        super(PlainVel, self).__init__(osi, 'vel', ts=ts, fact=fact)


class PlainDisp(GroundMotionBase):
    op_type = 'Plain'

    def __init__(self, osi, ts, fact=1.0):
        super(PlainDisp, self).__init__(osi, 'disp', ts=ts, fact=fact)

