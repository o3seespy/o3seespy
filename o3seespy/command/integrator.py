from o3seespy.base_model import OpenSeesObject


class IntegratorBase(OpenSeesObject):
    op_base_type = "integrator"


class LoadControl(IntegratorBase):
    """
    The LoadControl Integrator Class
    
    Create a OpenSees LoadControl integrator object.
    """
    op_type = 'LoadControl'

    def __init__(self, osi, incr, num_iter=1, min_incr: float=None, max_incr: float=None):
        r"""
        Initial method for LoadControl

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        incr: float
            Load factor increment :math:`\lambda`.
        num_iter: int, optional
            Number of iterations the user would like to occur in the solution algorithm. 
        min_incr: float (default=True), optional
            Min stepsize the user will allow :math:`\lambda_{min}`. 
        max_incr: float (default=True), optional
            Max stepsize the user will allow :math:`\lambda_{max}`. 

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.integrator.LoadControl(osi, incr=1.0, num_iter=1, min_incr=None, max_incr=None)
        """
        self.osi = osi
        self.incr = float(incr)
        self.num_iter = int(num_iter)
        if min_incr is None:
            self.min_incr = None
        else:
            self.min_incr = float(min_incr)
        if max_incr is None:
            self.max_incr = None
        else:
            self.max_incr = float(max_incr)
        self._parameters = [self.op_type, self.incr, self.num_iter]
        special_pms = ['min_incr', 'max_incr']
        packets = [False, False]
        for i, pm in enumerate(special_pms):
            if getattr(self, pm) is not None:
                if packets[i]:
                    self._parameters += [*getattr(self, pm)]
                else:
                    self._parameters += [getattr(self, pm)]
            else:
                break
        self.to_process(osi)


class DisplacementControl(IntegratorBase):
    """
    The DisplacementControl Integrator Class
    
    Create a DisplacementControl integrator.  In an analysis step with Displacement Control we seek to determine the
    time step that will result in a displacement increment for a particular degree-of-freedom at a node to be a
    prescribed value.
    """
    op_type = 'DisplacementControl'

    def __init__(self, osi, node, dof, incr, num_iter=1, d_umin: float=None, d_umax: float=None):
        r"""
        Initial method for DisplacementControl

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        node: obj
            Object of node whose response controls solution
        dof: int
            Degree of freedom at the node, 1 through ndf.
        incr: float
            First displacement increment :math:`\delta u_{dof}`.
        num_iter: int, optional
            Number of iterations the user would like to occur in the solution algorithm. 
        d_umin: None (default=True), optional
            
        d_umax: None (default=True), optional
            

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> node = o3.node.Node(osi, 0.0, 0.0)
        >>> o3.integrator.DisplacementControl(osi, node, dof=1, incr=1.0, num_iter=1, d_umin=None, d_umax=None)
        """
        self.osi = osi
        self.node = node
        self.dof = int(dof)
        self.incr = float(incr)
        self.num_iter = int(num_iter)
        self.d_umin = d_umin
        self.d_umax = d_umax
        self._parameters = [self.op_type, self.node.tag, self.dof, self.incr, self.num_iter]
        special_pms = ['d_umin', 'd_umax']
        packets = [False, False]
        for i, pm in enumerate(special_pms):
            if getattr(self, pm) is not None:
                if packets[i]:
                    self._parameters += [*getattr(self, pm)]
                else:
                    self._parameters += [getattr(self, pm)]
            else:
                break
        self.to_process(osi)


class ParallelDisplacementControl(IntegratorBase):
    """
    The ParallelDisplacementControl Integrator Class
    
    Create a Parallel version of DisplacementControl integrator.  In an analysis step with Displacement Control we seek
    to determine the time step that will result in a displacement increment for a particular degree-of-freedom at a node
    to be a prescribed value.
    """
    op_type = 'ParallelDisplacementControl'

    def __init__(self, osi, node, dof, incr, num_iter=1, d_umin: float=None, d_umax: float=None):
        r"""
        Initial method for ParallelDisplacementControl

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        node: obj
            Object of node whose response controls solution
        dof: int
            Degree of freedom at the node, 1 through ndf.
        incr: float
            First displacement increment :math:`\delta u_{dof}`.
        num_iter: int, optional
            Number of iterations the user would like to occur in the solution algorithm. 
        d_umin: None (default=True), optional
            
        d_umax: None (default=True), optional
            

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> node = o3.node.Node(osi, 0.0, 0.0)
        >>> o3.integrator.ParallelDisplacementControl(osi, node, dof=1, incr=1.0, num_iter=1, d_umin=None, d_umax=None)
        """
        self.osi = osi
        self.node = node
        self.dof = int(dof)
        self.incr = float(incr)
        self.num_iter = int(num_iter)
        self.d_umin = d_umin
        self.d_umax = d_umax
        self._parameters = [self.op_type, self.node.tag, self.dof, self.incr, self.num_iter]
        special_pms = ['d_umin', 'd_umax']
        packets = [False, False]
        for i, pm in enumerate(special_pms):
            if getattr(self, pm) is not None:
                if packets[i]:
                    self._parameters += [*getattr(self, pm)]
                else:
                    self._parameters += [getattr(self, pm)]
            else:
                break
        self.to_process(osi)


class MinUnbalDispNorm(IntegratorBase):
    """
    The MinUnbalDispNorm Integrator Class
    
    Create a MinUnbalDispNorm integrator.
    """
    op_type = 'MinUnbalDispNorm'

    def __init__(self, osi, dlambda1, jd=1, min_lambda: float=None, max_lambda: float=None, det: float=None):
        """
        Initial method for MinUnbalDispNorm

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        dlambda1: float
            First load increment (pseudo-time step) at the first iteration in the next invocation of the analysis
            command.
        jd: int, optional
            Factor relating first load increment at subsequent time steps. 
        min_lambda: float (default=True), optional
            Min load increment. 
        max_lambda: float (default=True), optional
            Max load increment. 
        det: None (default=True), optional
            

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.integrator.MinUnbalDispNorm(osi, dlambda1=1.0, jd=1, min_lambda=None, max_lambda=None, det=None)
        """
        self.osi = osi
        self.dlambda1 = float(dlambda1)
        self.jd = int(jd)
        if min_lambda is None:
            self.min_lambda = None
        else:
            self.min_lambda = float(min_lambda)
        if max_lambda is None:
            self.max_lambda = None
        else:
            self.max_lambda = float(max_lambda)
        self.det = det
        self._parameters = [self.op_type, self.dlambda1, self.jd]
        special_pms = ['min_lambda', 'max_lambda', 'det']
        packets = [False, False, False]
        for i, pm in enumerate(special_pms):
            if getattr(self, pm) is not None:
                if packets[i]:
                    self._parameters += [*getattr(self, pm)]
                else:
                    self._parameters += [getattr(self, pm)]
            else:
                break
        self.to_process(osi)


class ArcLength(IntegratorBase):
    """
    The ArcLength Integrator Class
    
    Create a ArcLength integrator. In an analysis step with ArcLength we seek to determine the time step that will
    result in our constraint equation being satisfied.
    """
    op_type = 'ArcLength'

    def __init__(self, osi, s, alpha):
        r"""
        Initial method for ArcLength

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        s: float
            The arclength.
        alpha: float
            :math:`\alpha` a scaling factor on the reference loads.

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.integrator.ArcLength(osi, s=1.0, alpha=1.0)
        """
        self.osi = osi
        self.s = float(s)
        self.alpha = float(alpha)
        self._parameters = [self.op_type, self.s, self.alpha]
        self.to_process(osi)


class CentralDifference(IntegratorBase):
    r"""
    The CentralDifference Integrator Class
    
    Create a centralDifference integrator.#. The calculation of :math:`U_t + \Delta t`, is based on using the
    equilibrium equation at time t. For this reason the method is called an explicit integration method.#. If
    there is no rayleigh damping and the C matrix is 0, for a diagonal mass matrix a diagonal solver may and
    should be used.#. For stability, :math:`\frac{\Delta t}{T_n} < \frac{1}{\pi}`
    """
    op_type = 'CentralDifference'

    def __init__(self, osi):
        """
        Initial method for CentralDifference

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.integrator.CentralDifference(osi)
        """
        self.osi = osi
        self._parameters = [self.op_type]
        self.to_process(osi)


class Newmark(IntegratorBase):
    """
    The Newmark Integrator Class
    
    Create a Newmark integrator.
    """
    op_type = 'Newmark'

    def __init__(self, osi, gamma, beta, form: str=None):
        r"""
        Initial method for Newmark

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        gamma: float
            :math:`\gamma` factor.
        beta: float
            :math:`\beta` factor.
        form: str, optional
            Flag to indicate which variable to be used as primary variable  * ``'d'`` -- displacement (default) *
            ``'v'`` -- velocity * ``'a'`` -- acceleration

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.integrator.Newmark(osi, gamma=1.0, beta=1.0, form=1)
        """
        self.osi = osi
        self.gamma = float(gamma)
        self.beta = float(beta)
        self.form = form
        self._parameters = [self.op_type, self.gamma, self.beta]
        if getattr(self, 'form') is not None:
            self._parameters += ['-formD', self.form]
        self.to_process(osi)


class HHT(IntegratorBase):
    """
    The HHT Integrator Class
    
    Create a Hilber-Hughes-Taylor (HHT) integrator. This is an implicit method that allows for energy dissipation and
    second order accuracy (which is not possible with the regular Newmark object). Depending on choices of input
    parameters, the method can be unconditionally stable.
    """
    op_type = 'HHT'

    def __init__(self, osi, alpha, gamma: float=None, beta: float=None):
        r"""
        Initial method for HHT

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        alpha: float
            :math:`\alpha` factor.
        gamma: float (default=True), optional
            :math:`\gamma` factor. 
        beta: float (default=True), optional
            :math:`\beta` factor. 

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.integrator.HHT(osi, alpha=1.0, gamma=1.0, beta=1.0)
        """
        self.osi = osi
        self.alpha = float(alpha)
        if gamma is None:
            self.gamma = None
        else:
            self.gamma = float(gamma)
        if beta is None:
            self.beta = None
        else:
            self.beta = float(beta)
        self._parameters = [self.op_type, self.alpha]
        special_pms = ['gamma', 'beta']
        packets = [False, False]
        for i, pm in enumerate(special_pms):
            if getattr(self, pm) is not None:
                if packets[i]:
                    self._parameters += [*getattr(self, pm)]
                else:
                    self._parameters += [getattr(self, pm)]
            else:
                break
        self.to_process(osi)


class GeneralizedAlpha(IntegratorBase):
    r"""
    The GeneralizedAlpha Integrator Class
    
    Create a GeneralizedAlpha integrator. This is an implicit method that like the HHT method allows for high frequency
    energy dissipation and second order accuracy, i.e. :math:`\Delta t^2`. Depending on choices of input parameters, the
    method can be unconditionally stable.
    """
    op_type = 'GeneralizedAlpha'

    def __init__(self, osi, alpha_m, alpha_f, gamma: float=None, beta: float=None):
        r"""
        Initial method for GeneralizedAlpha

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        alpha_m: float
            :math:`\alpha_m` factor.
        alpha_f: float
            :math:`\alpha_f` factor.
        gamma: float (default=True), optional
            :math:`\gamma` factor. 
        beta: float (default=True), optional
            :math:`\beta` factor. 

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.integrator.GeneralizedAlpha(osi, alpha_m=1.0, alpha_f=1.0, gamma=1.0, beta=1.0)
        """
        self.osi = osi
        self.alpha_m = float(alpha_m)
        self.alpha_f = float(alpha_f)
        if gamma is None:
            self.gamma = None
        else:
            self.gamma = float(gamma)
        if beta is None:
            self.beta = None
        else:
            self.beta = float(beta)
        self._parameters = [self.op_type, self.alpha_m, self.alpha_f]
        special_pms = ['gamma', 'beta']
        packets = [False, False]
        for i, pm in enumerate(special_pms):
            if getattr(self, pm) is not None:
                if packets[i]:
                    self._parameters += [*getattr(self, pm)]
                else:
                    self._parameters += [getattr(self, pm)]
            else:
                break
        self.to_process(osi)


class TRBDF2(IntegratorBase):
    """
    The TRBDF2 Integrator Class
    
    Create a TRBDF2 integrator. The TRBDF2 integrator is a composite scheme that alternates between the Trapezoidal
    scheme and a 3 point backward Euler scheme. It does this in an attempt to conserve energy and momentum, something
    Newmark does not always do.As opposed to dividing the time-step in 2 as outlined in the `Bathe2007`_, we just
    switch alternate between the 2 integration strategies,i.e. the time step in our implementation is double
    that described in the `Bathe2007`_.
    """
    op_type = 'TRBDF2'

    def __init__(self, osi):
        """
        Initial method for TRBDF2

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.integrator.TRBDF2(osi)
        """
        self.osi = osi
        self._parameters = [self.op_type]
        self.to_process(osi)


class ExplicitDifference(IntegratorBase):
    r"""
    The ExplicitDifference Integrator Class
    
    Create a ExplicitDifference integrator.#. When using Rayleigh damping, the damping ratio of high vibration modes is
    overrated, and the critical time step size will be much smaller. Hence Modal damping is more suitable for this
    method.#. There should be no zero element on the diagonal of the mass matrix when using this method.#.
    Diagonal solver should be used when lumped mass matrix is used because the equations are uncoupled.#.
    For stability, :math:`\Delta t \leq \left(\sqrt{\zeta^2+1}-\zeta\right)\frac{2}{\omega}`
    """
    op_type = 'ExplicitDifference'

    def __init__(self, osi):
        """
        Initial method for ExplicitDifference

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.integrator.ExplicitDifference(osi)
        """
        self.osi = osi
        self._parameters = [self.op_type]
        self.to_process(osi)
