from o3seespy.base_model import OpenSeesObject


class PatternBase(OpenSeesObject):
    op_base_type = "pattern"


class Plain(PatternBase):
    """
    The Plain Pattern Class
    
    This commnand allows the user to construct a LoadPattern object. Each plain load pattern is associated with a
    TimeSeries object and can contain multiple NodalLoads, ElementalLoads and SP_Constraint objects. The command to
    generate LoadPattern object contains in { } the commands to generate all the loads and the single-point
    constraints in the pattern. To construct a load pattern and populate it, the following command is used:
    """
    op_type = 'Plain'

    def __init__(self, osi, ts, fact: float=None):
        """
        Initial method for Plain

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ts: obj
            The object of the time series to be used in the load pattern
        fact: float, optional
            Constant factor. 

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> ts = o3.time_series.Linear(osi, factor=1.0)
        >>> o3.pattern.Plain(osi, ts=ts, fact=1.0)
        """
        self.osi = osi
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
    """
    The UniformExcitation Pattern Class
    
    The UniformExcitation pattern allows the user to apply a uniform excitation to a model acting in a certain
    direction. The command is as follows:
    """
    op_type = 'UniformExcitation'

    def __init__(self, osi, dir, disp_series=None, vel_series=None, accel_series=None, vel0: float=None, fact: float=None):
        """
        Initial method for UniformExcitation

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        dir: int
            Direction in which ground motion acts #. corresponds to translation along the global x axis #. corresponds
            to translation along the global y axis #. corresponds to translation along the global z axis #. corresponds to
            rotation about the global x axis #. corresponds to rotation about the global y axis #. corresponds to
            rotation about the global z axis
        disp_series: obj, optional
            Object of the timeseries series defining the displacement history. 
        vel_series: obj, optional
            Object of the timeseries series defining the velocity history. 
        accel_series: obj, optional
            Object of the timeseries series defining the acceleration history. 
        vel0: float, optional
            The initial velocity (optional, default=0.0)
        fact: float, optional
            Constant factor (optional, default=1.0)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> ts = o3.time_series.Linear(osi, factor=1.0)
        >>> o3.pattern.UniformExcitation(osi, dir=1, accel_series=ts, vel0=1.0, fact=1.0)
        """
        self.osi = osi
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
            self._parameters += ['-disp', self.disp_series.tag]
        if getattr(self, 'vel_series') is not None:
            self._parameters += ['-vel', self.vel_series.tag]
        if getattr(self, 'accel_series') is not None:
            self._parameters += ['-accel', self.accel_series.tag]
        if getattr(self, 'vel0') is not None:
            self._parameters += ['-vel0', self.vel0]
        if getattr(self, 'fact') is not None:
            self._parameters += ['-fact', self.fact]
        self.to_process(osi)


class MultipleSupport(PatternBase):
    """
    The MultipleSupport Pattern Class
    
    The Multi-Support pattern allows similar or different prescribed ground motions to be input at various supports in
    the structure. In OpenSees, the prescribed motion is applied using single-point constraints, the single-point
    constraints taking their constraint value from user created ground motions.
    """
    op_type = 'MultipleSupport'

    def __init__(self, osi):
        """
        Initial method for MultipleSupport

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.pattern.MultipleSupport(osi)
        """
        self.osi = osi
        osi.n_pat += 1
        self._tag = osi.n_pat
        self._parameters = [self.op_type, self._tag]
        self.to_process(osi)
