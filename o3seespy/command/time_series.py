from o3seespy.base_model import OpenSeesObject


class TimeSeriesBase(OpenSeesObject):
    op_base_type = "timeSeries"


class Constant(TimeSeriesBase):
    r"""
    The Constant TimeSeries Class
    
    This command is used to construct a TimeSeries object in which the load factor applied remains constant and is
    independent of the time in the domain, i.e. :math:`\lambda = f(t) = C`.
    """
    op_type = 'Constant'

    def __init__(self, osi, factor: float=None):
        """
        Initial method for Constant

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        factor: float, optional
            The load factor applied 

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.time_series.Constant(osi, factor=1.0)
        """
        self.osi = osi
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
    r"""
    The Linear TimeSeries Class
    
    This command is used to construct a TimeSeries object in which the load factor applied is linearly proportional to
    the time in the domain, i.e.:math:`\lambda = f(t) = cFactor * t`. 
    """
    op_type = 'Linear'

    def __init__(self, osi, factor: float=None):
        """
        Initial method for Linear

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        factor: float, optional
            Linear factor. 

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.time_series.Linear(osi, factor=1.0)
        """
        self.osi = osi
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
    r"""
    The Trig TimeSeries Class
    
    This command is used to construct a TimeSeries object in which the load factor is some trigonemtric function of the
    time in the domain

    .. math::

      \lambda = f(t) = 
      \begin{cases}
          cFactor * sin(\frac{2.0\pi(t-tStart)}{period}+\phi), &  tStart<=t<=tEnd\\
          0.0, & otherwise
      \end{cases}
      \phi = shift - \frac{period}{2.0\pi} * \arcsin(\frac{zeroShift}{cFactor})
    """
    op_type = 'Trig'

    def __init__(self, osi, t_start, t_end, period, factor: float=None, shift: float=None, zero_shift: float=None):
        """
        Initial method for Trig

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        t_start: float
            Starting time of non-zero load factor.
        t_end: float
            Ending time of non-zero load factor.
        period: float
            Characteristic period of sine wave.
        factor: float, optional
            Load factor. 
        shift: float, optional
            Phase shift in radians. 
        zero_shift: float, optional
            Zero shift. 

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.time_series.Trig(osi, t_start=1.0, t_end=1.0, period=1.0, factor=1.0, shift=0.0, zero_shift=0.0)
        """
        self.osi = osi
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
    r"""
    The Triangle TimeSeries Class
    
    This command is used to construct a TimeSeries object in which the load factor is some triangular function of the
    time in the domain.

    .. math::

      \lambda = f(t) = 
      \begin{cases}
          slope*k*period+zeroShift, & k < 0.25\\
      cFactor-slope*(k-0.25)*period+zeroShift, & k < 0.75\\
      -cFactor+slope*(k-0.75)*period+zeroShift, & k < 1.0\\
      0.0, & otherwise
      \end{cases}
    

    .. math::

      slope = \frac{cFactor}{period/4}
      k = \frac{t+\phi-tStart}{period}-floor(\frac{t+\phi-tStart}{period})
      \phi = shift - \frac{zeroShift}{slope}
    """
    op_type = 'Triangle'

    def __init__(self, osi, t_start, t_end, period, factor: float=None, shift: float=None, zero_shift: float=None):
        """
        Initial method for Triangle

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        t_start: float
            Starting time of non-zero load factor.
        t_end: float
            Ending time of non-zero load factor.
        period: float
            Characteristic period of sine wave.
        factor: float, optional
            Load factor. 
        shift: float, optional
            Phase shift in radians. 
        zero_shift: float, optional
            Zero shift. 

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.time_series.Triangle(osi, t_start=1.0, t_end=1.0, period=1.0, factor=1.0, shift=0.0, zero_shift=0.0)
        """
        self.osi = osi
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
    r"""
    The Rectangular TimeSeries Class
    
    This command is used to construct a TimeSeries object in which the load factor is constant for a specified period
    and 0 otherwise, i.e.

    .. math::

      \lambda = f(t) = 
      \begin{cases}
          cFactor, &  tStart<=t<=tEnd\\
      0.0, & otherwise
      \end{cases}
    """
    op_type = 'Rectangular'

    def __init__(self, osi, t_start, t_end, factor: float=None):
        """
        Initial method for Rectangular

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        t_start: float
            Starting time of non-zero load factor.
        t_end: float
            Ending time of non-zero load factor.
        factor: float, optional
            Load factor. 

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.time_series.Rectangular(osi, t_start=1.0, t_end=1.0, factor=1.0)
        """
        self.osi = osi
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
    r"""
    The Pulse TimeSeries Class
    
    This command is used to construct a TimeSeries object in which the load factor is some pulse function of the time in
    the domain.

    .. math::

      \lambda = f(t) = 
      \begin{cases}
          cFactor+zeroShift, &  k < width\\
      zeroshift, & k < 1\\
      0.0, & otherwise
      \end{cases}
    

    .. math::

      k = \frac{t+shift-tStart}{period}-floor(\frac{t+shift-tStart}{period})
    """
    op_type = 'Pulse'

    def __init__(self, osi, t_start, t_end, period, width: float=None, shift: float=None, factor: float=None, zero_shift: float=None):
        """
        Initial method for Pulse

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        t_start: float
            Starting time of non-zero load factor.
        t_end: float
            Ending time of non-zero load factor.
        period: float
            Characteristic period of pulse.
        width: float, optional
            Pulse width as a fraction of the period. (optinal)
        shift: float, optional
            Phase shift in seconds. 
        factor: float, optional
            Load factor. 
        zero_shift: float, optional
            Zero shift. 

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.time_series.Pulse(osi, t_start=1.0, t_end=1.0, period=1.0, width=0.5, shift=0.0, factor=1.0, zero_shift=0.0)
        """
        self.osi = osi
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
    """
    The Path TimeSeries Class
    
    The relationship between loadfactor and time is input by the user as a series of discrete points inthe 2d space
    (load factor, time). The input points can come from afile or from a list in the script. When the time specified
    does not matchany of the input points, linear interpolation is used between points.There are many ways to
    specify the load path, for example,
    """
    op_type = 'Path'

    def __init__(self, osi, dt: float=None, values: list=None, time: list=None, filepath: str=None, file_time: str=None, factor: float=None, start_time: float=None, use_last=False, prepend_zero=False):
        """
        Initial method for Path

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        dt: float, optional
            Time interval between specified points. 
        values: list, optional
            Load factor values in a |list|. 
        time: list, optional
            Time values in a |list|. 
        filepath: str, optional
            File containing the load factors values. 
        file_time: str, optional
            File containing the time values for corresponding load factors. 
        factor: float, optional
            A factor to multiply load factors by. 
        start_time: float, optional
            Provide a start time for provided load factors. 
        use_last: bool
            Use last value after the end of the series. 
        prepend_zero: bool
            Prepend a zero value to the series of load factors. 
        """
        self.osi = osi
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
