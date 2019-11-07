from o3seespy.base_model import OpenseesObject


class PatternBase(OpenseesObject):
    op_base_type = "pattern"


class Plain(PatternBase):
    op_type = 'Plain'

    def __init__(self, osi, ts, fact=None):
        self.ts = ts
        if fact is None:
            self.fact = None
        else:
            self.fact = float(fact)
        osi.n_pat += 1
        self._tag = osi.n_pat
        self._parameters = [self.op_type, self._tag, self.ts.tag]
        if getattr(self, 'fact') is not None:
            self._parameters += ['-fact', self.fact]
        self.to_process(osi)


class UniformExcitation(PatternBase):
    op_type = 'UniformExcitation'

    def __init__(self, osi, dir, disp_series=None, vel_series=None, accel_series=None, vel0=None, fact=None):
        self.dir = int(dir)
        self.disp_series = disp_series
        self.vel_series = vel_series
        self.accel_series = accel_series
        if vel0 is None:
            self.vel0 = None
        else:
            self.vel0 = float(vel0)
        if fact is None:
            self.fact = None
        else:
            self.fact = float(fact)
        osi.n_pat += 1
        self._tag = osi.n_pat
        self._parameters = [self.op_type, self._tag, self.dir]
        if getattr(self, 'disp_series') is not None:
            self._parameters += ['-disp', self.disp_series]
        if getattr(self, 'vel_series') is not None:
            self._parameters += ['-vel', self.vel_series]
        if getattr(self, 'accel_series') is not None:
            self._parameters += ['-accel', self.accel_series]
        if getattr(self, 'vel0') is not None:
            self._parameters += ['-vel0', self.vel0]
        if getattr(self, 'fact') is not None:
            self._parameters += ['-fact', self.fact]
        self.to_process(osi)


class MultipleSupport(PatternBase):
    op_type = 'MultipleSupport'

    def __init__(self, osi, ):
        osi.n_pat += 1
        self._tag = osi.n_pat
        self._parameters = [self.op_type, self._tag]
        self.to_process(osi)
