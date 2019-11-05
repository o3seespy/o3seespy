from o3seespy.base_model import OpenseesObject


class TimeSeriesBase(OpenseesObject):
    op_base_type = "timeSeries"


class Constant(TimeSeriesBase):
    op_type = 'Constant'

    def __init__(self, osi, factor=1.0):
        self.factor = float(factor)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, ]
        if getattr(self, 'factor') is not None:
            self._parameters += ['-factor', self.factor]
        self.to_process(osi)


class Linear(TimeSeriesBase):
    op_type = 'Linear'

    def __init__(self, osi, factor=1.0):
        self.factor = float(factor)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, ]
        if getattr(self, 'factor') is not None:
            self._parameters += ['-factor', self.factor]
        self.to_process(osi)


class Trig(TimeSeriesBase):
    op_type = 'Trig'

    def __init__(self, osi, t_start, t_end, period, factor=1.0, shift=0.0, zero_shift=0.0):
        self.t_start = float(t_start)
        self.t_end = float(t_end)
        self.period = float(period)
        self.factor = float(factor)
        self.shift = float(shift)
        self.zero_shift = float(zero_shift)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.t_start, self.t_end, self.period]
        if getattr(self, 'factor') is not None:
            self._parameters += ['-factor', self.factor]
        if getattr(self, 'shift') is not None:
            self._parameters += ['-shift', self.shift]
        if getattr(self, 'zero_shift') is not None:
            self._parameters += ['-zeroShift', self.zero_shift]
        self.to_process(osi)


class Triangle(TimeSeriesBase):
    op_type = 'Triangle'

    def __init__(self, osi, t_start, t_end, period, factor=1.0, shift=0.0, zero_shift=0.0):
        self.t_start = float(t_start)
        self.t_end = float(t_end)
        self.period = float(period)
        self.factor = float(factor)
        self.shift = float(shift)
        self.zero_shift = float(zero_shift)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.t_start, self.t_end, self.period]
        if getattr(self, 'factor') is not None:
            self._parameters += ['-factor', self.factor]
        if getattr(self, 'shift') is not None:
            self._parameters += ['-shift', self.shift]
        if getattr(self, 'zero_shift') is not None:
            self._parameters += ['-zeroShift', self.zero_shift]
        self.to_process(osi)


class Rectangular(TimeSeriesBase):
    op_type = 'Rectangular'

    def __init__(self, osi, t_start, t_end, factor=1.0):
        self.t_start = float(t_start)
        self.t_end = float(t_end)
        self.factor = float(factor)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.t_start, self.t_end]
        if getattr(self, 'factor') is not None:
            self._parameters += ['-factor', self.factor]
        self.to_process(osi)


class Pulse(TimeSeriesBase):
    op_type = 'Pulse'

    def __init__(self, osi, t_start, t_end, period, width=0.5, shift=0.0, factor=1.0, zero_shift=0.0):
        self.t_start = float(t_start)
        self.t_end = float(t_end)
        self.period = float(period)
        self.width = float(width)
        self.shift = float(shift)
        self.factor = float(factor)
        self.zero_shift = float(zero_shift)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.t_start, self.t_end, self.period]
        if getattr(self, 'width') is not None:
            self._parameters += ['-width', self.width]
        if getattr(self, 'shift') is not None:
            self._parameters += ['-shift', self.shift]
        if getattr(self, 'factor') is not None:
            self._parameters += ['-factor', self.factor]
        if getattr(self, 'zero_shift') is not None:
            self._parameters += ['-zeroShift', self.zero_shift]
        self.to_process(osi)


class Path(TimeSeriesBase):
    op_type = 'Path'

    def __init__(self, osi, dt=0.0, values=None, time=None, filepath='', file_time='', factor=1.0, start_time=0.0, use_last=False, prepend_zero=False):
        self.dt = float(dt)
        self.values = values
        self.time = time
        self.filepath = filepath
        self.file_time = file_time
        self.factor = float(factor)
        self.start_time = float(start_time)
        self.use_last = use_last
        self.prepend_zero = prepend_zero
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, ]
        if getattr(self, 'dt') is not None:
            self._parameters += ['-dt', self.dt]
        if getattr(self, 'values') is not None:
            self._parameters += ['-values', *self.values]
        if getattr(self, 'time') is not None:
            self._parameters += ['-time', *self.time]
        if getattr(self, 'filepath') is not None:
            self._parameters += ['-filepath', self.filepath]
        if getattr(self, 'file_time') is not None:
            self._parameters += ['-fileTime', self.file_time]
        if getattr(self, 'factor') is not None:
            self._parameters += ['-factor', self.factor]
        if getattr(self, 'start_time') is not None:
            self._parameters += ['-startTime', self.start_time]
        if getattr(self, 'use_last'):
            self._parameters += ['-useLast']
        if getattr(self, 'prepend_zero'):
            self._parameters += ['-prependZero']
        self.to_process(osi)
