from o3seespy.base_model import OpenseesObject


class TimeSeriesBase(OpenseesObject):
    op_base_type = "timeSeries"


class Constant(TimeSeriesBase):
    op_type = 'Constant'

    def __init__(self, osi, factor=None):
        if factor is None:
            self.factor = None
        else:
            self.factor = float(factor)
        osi.n_tseries += 1
        self._tag = osi.n_tseries
        self._parameters = [self.op_type, self._tag]
        if getattr(self, 'factor') is not None:
            self._parameters += ['-factor', self.factor]
        self.to_process(osi)


class Linear(TimeSeriesBase):
    op_type = 'Linear'

    def __init__(self, osi, factor=None):
        if factor is None:
            self.factor = None
        else:
            self.factor = float(factor)
        osi.n_tseries += 1
        self._tag = osi.n_tseries
        self._parameters = [self.op_type, self._tag]
        if getattr(self, 'factor') is not None:
            self._parameters += ['-factor', self.factor]
        self.to_process(osi)


class Trig(TimeSeriesBase):
    op_type = 'Trig'

    def __init__(self, osi, t_start, t_end, period, factor=None, shift=None, zero_shift=None):
        self.t_start = float(t_start)
        self.t_end = float(t_end)
        self.period = float(period)
        if factor is None:
            self.factor = None
        else:
            self.factor = float(factor)
        if shift is None:
            self.shift = None
        else:
            self.shift = float(shift)
        if zero_shift is None:
            self.zero_shift = None
        else:
            self.zero_shift = float(zero_shift)
        osi.n_tseries += 1
        self._tag = osi.n_tseries
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

    def __init__(self, osi, t_start, t_end, period, factor=None, shift=None, zero_shift=None):
        self.t_start = float(t_start)
        self.t_end = float(t_end)
        self.period = float(period)
        if factor is None:
            self.factor = None
        else:
            self.factor = float(factor)
        if shift is None:
            self.shift = None
        else:
            self.shift = float(shift)
        if zero_shift is None:
            self.zero_shift = None
        else:
            self.zero_shift = float(zero_shift)
        osi.n_tseries += 1
        self._tag = osi.n_tseries
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

    def __init__(self, osi, t_start, t_end, factor=None):
        self.t_start = float(t_start)
        self.t_end = float(t_end)
        if factor is None:
            self.factor = None
        else:
            self.factor = float(factor)
        osi.n_tseries += 1
        self._tag = osi.n_tseries
        self._parameters = [self.op_type, self._tag, self.t_start, self.t_end]
        if getattr(self, 'factor') is not None:
            self._parameters += ['-factor', self.factor]
        self.to_process(osi)


class Pulse(TimeSeriesBase):
    op_type = 'Pulse'

    def __init__(self, osi, t_start, t_end, period, width=None, shift=None, factor=None, zero_shift=None):
        self.t_start = float(t_start)
        self.t_end = float(t_end)
        self.period = float(period)
        if width is None:
            self.width = None
        else:
            self.width = float(width)
        if shift is None:
            self.shift = None
        else:
            self.shift = float(shift)
        if factor is None:
            self.factor = None
        else:
            self.factor = float(factor)
        if zero_shift is None:
            self.zero_shift = None
        else:
            self.zero_shift = float(zero_shift)
        osi.n_tseries += 1
        self._tag = osi.n_tseries
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

    def __init__(self, osi, dt=None, values=None, time=None, filepath=None, file_time=None, factor=None, start_time=None, use_last=False, prepend_zero=False):
        if dt is None:
            self.dt = None
        else:
            self.dt = float(dt)
        self.values = values
        self.time = time
        self.filepath = filepath
        self.file_time = file_time
        if factor is None:
            self.factor = None
        else:
            self.factor = float(factor)
        if start_time is None:
            self.start_time = None
        else:
            self.start_time = float(start_time)
        self.use_last = use_last
        self.prepend_zero = prepend_zero
        osi.n_tseries += 1
        self._tag = osi.n_tseries
        self._parameters = [self.op_type, self._tag]
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
