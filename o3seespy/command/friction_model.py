from o3seespy.base_model import OpenSeesObject


class FrictionModelBase(OpenSeesObject):
    op_base_type = "frictionModel"


class Coulomb(FrictionModelBase):
    """
    The Coulomb FrictionModel Class
    
    This command is used to construct a `Coulomb friction <http://en.wikipedia.org/wiki/Friction>`_ model object.
    Coulomb's Law of Friction states that kinetic friction is independent of the sliding velocity.
    """
    op_type = 'Coulomb'

    def __init__(self, osi, mu):
        """
        Initial method for Coulomb

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        mu: float
            Coefficient of friction

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.friction_model.Coulomb(osi, mu=1.0)
        """
        self.osi = osi
        self.mu = float(mu)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.mu]
        self.to_process(osi)


class VelDependent(FrictionModelBase):
    """
    The VelDependent FrictionModel Class
    
    This command is used to construct a VelDependent friction model object. It is useful for modeling the behavior of
    `PTFE <http://en.wikipedia.org/wiki/Polytetrafluoroethylene>`_ or PTFE-like materials sliding on a stainless steel
    surface. For a detailed presentation on the velocity dependence of such interfaces please refer to Constantinou
    et al. (1999).
    """
    op_type = 'VelDependent'

    def __init__(self, osi, mu_slow, mu_fast, trans_rate):
        """
        Initial method for VelDependent

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        mu_slow: float
            Coefficient of friction at low velocity
        mu_fast: float
            Coefficient of friction at high velocity
        trans_rate: float
            Transition rate from low to high velocity

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.friction_model.VelDependent(osi, mu_slow=1.0, mu_fast=1.0, trans_rate=1.0)
        """
        self.osi = osi
        self.mu_slow = float(mu_slow)
        self.mu_fast = float(mu_fast)
        self.trans_rate = float(trans_rate)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.mu_slow, self.mu_fast, self.trans_rate]
        self.to_process(osi)


class VelNormalFrcDep(FrictionModelBase):
    """
    The VelNormalFrcDep FrictionModel Class
    
    This command is used to construct a VelNormalFrcDep friction model object.
    """
    op_type = 'VelNormalFrcDep'

    def __init__(self, osi, a_slow, n_slow, a_fast, n_fast, alpha0, alpha1, alpha2, max_mu_fact):
        r"""
        Initial method for VelNormalFrcDep

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        a_slow: float
            Constant for coefficient of friction at low velocity
        n_slow: float
            Exponent for coefficient of friction at low velocity
        a_fast: float
            Constant for coefficient of friction at high velocity
        n_fast: float
            Exponent for coefficient of friction at high velocity
        alpha0: float
            Constant rate parameter coefficient
        alpha1: float
            Linear rate parameter coefficient
        alpha2: float
            Quadratic rate parameter coefficient
        max_mu_fact: float
            Factor for determining the maximum coefficient of friction. this value prevents the friction coefficient
            from exceeding an unrealistic maximum value when the normal force becomes very small. the maximum friction
            coefficient is determined from μfast, for example :math:`\mu \leq maxmufac*μfast`.

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.friction_model.VelNormalFrcDep(osi, a_slow=1.0, n_slow=1.0, a_fast=1.0, n_fast=1.0, alpha0=1.0, alpha1=1.0, alpha2=1.0, max_mu_fact=1.0)
        """
        self.osi = osi
        self.a_slow = float(a_slow)
        self.n_slow = float(n_slow)
        self.a_fast = float(a_fast)
        self.n_fast = float(n_fast)
        self.alpha0 = float(alpha0)
        self.alpha1 = float(alpha1)
        self.alpha2 = float(alpha2)
        self.max_mu_fact = float(max_mu_fact)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.a_slow, self.n_slow, self.a_fast, self.n_fast, self.alpha0, self.alpha1, self.alpha2, self.max_mu_fact]
        self.to_process(osi)


class VelPressureDep(FrictionModelBase):
    """
    The VelPressureDep FrictionModel Class
    
    This command is used to construct a VelPressureDep friction model object.
    """
    op_type = 'VelPressureDep'

    def __init__(self, osi, mu_slow, mu_fast0, big_a, delta_mu, alpha, trans_rate):
        """
        Initial method for VelPressureDep

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        mu_slow: float
            Coefficient of friction at low velocity
        mu_fast0: float
            Initial coefficient of friction at high velocity
        big_a: float
            Nominal contact area
        delta_mu: float
            Pressure parameter calibrated from experimental data
        alpha: float
            Pressure parameter calibrated from experimental data
        trans_rate: float
            Transition rate from low to high velocity

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.friction_model.VelPressureDep(osi, mu_slow=1.0, mu_fast0=1.0, big_a=1.0, delta_mu=1.0, alpha=1.0, trans_rate=1.0)
        """
        self.osi = osi
        self.mu_slow = float(mu_slow)
        self.mu_fast0 = float(mu_fast0)
        self.big_a = float(big_a)
        self.delta_mu = float(delta_mu)
        self.alpha = float(alpha)
        self.trans_rate = float(trans_rate)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.mu_slow, self.mu_fast0, self.big_a, self.delta_mu, self.alpha, self.trans_rate]
        self.to_process(osi)


class VelDepMultiLinear(FrictionModelBase):
    """
    The VelDepMultiLinear FrictionModel Class
    
    This command is used to construct a VelDepMultiLinear friction model object. The friction-velocity relationship is
    given by a multi-linear curve that is define by a set of points. The slope given by the last two specified points on
    the positive velocity axis is extrapolated to infinite positive velocities. Velocity and friction points need to be
    equal or larger than zero (no negative values should be defined). The number of provided velocity points needs to
    be equal to the number of provided friction points.
    """
    op_type = 'VelDepMultiLinear'

    def __init__(self, osi, vel_points: list=None, frn_points: list=None):
        """
        Initial method for VelDepMultiLinear

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        vel_points: list, optional
            List of velocity points along friction-velocity curve
        frn_points: list, optional
            List of friction points along friction-velocity curve

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> vel_points = [0.0, 1.0]
        >>> frn_points = [1.0, 1.0]
        >>> o3.friction_model.VelDepMultiLinear(osi, vel_points=vel_points, frn_points=frn_points)
        """
        self.osi = osi
        self.vel_points = vel_points
        self.frn_points = frn_points
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag]
        if getattr(self, 'vel_points') is not None:
            self._parameters += ['-vel', *self.vel_points]
        if getattr(self, 'frn_points') is not None:
            self._parameters += ['-frn', *self.frn_points]
        self.to_process(osi)
