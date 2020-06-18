from o3seespy.base_model import OpenSeesObject


class BeamIntegrationBase(OpenSeesObject):
    op_base_type = "beamIntegration"


class Lobatto(BeamIntegrationBase):
    """
    The Lobatto BeamIntegration Class
    
    Create a Gauss-Lobatto beamIntegration object.Gauss-Lobatto integration is the most common approach for evaluating
    the response of`forceBeamColumn-Element` (`Neuenhofer and Filippou 1997`_) because it places an integration point at
    each end of the element, where bending moments are largest in the absence of interior element loads.
    """
    op_type = 'Lobatto'

    def __init__(self, osi, sec, big_n):
        """
        Initial method for Lobatto

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        sec: obj
            A previous-defined section object.
        big_n: int
            Number of integration points along the element.

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> sec = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
        >>> o3.beam_integration.Lobatto(osi, sec=sec, big_n=1)
        """
        self.osi = osi
        self.sec = sec
        self.big_n = int(big_n)
        osi.n_integ += 1
        self._tag = osi.n_integ
        self._parameters = [self.op_type, self._tag, self.sec.tag, self.big_n]
        self.to_process(osi)


class Legendre(BeamIntegrationBase):
    """
    The Legendre BeamIntegration Class
    
    Create a Gauss-Legendre beamIntegration object.Gauss-Legendre integration is more accurate than Gauss-Lobatto;
    however, it is not commonin force-based elements because there are no integration points at the element ends.of
    each integration point are tabulated in references on numerical analysis.The force deformation response at
    each integration point is defined by the section.The order of accuracy for Gauss-Legendre integration is
    2N-1.Arguments and examples see `Lobatto-BeamIntegration`.
    """
    op_type = 'Legendre'

    def __init__(self, osi, sec, big_n):
        """
        Initial method for Legendre

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        sec: obj
            
        big_n: None
            

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> sec = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
        >>> o3.beam_integration.Legendre(osi, sec=sec, big_n=1)
        """
        self.osi = osi
        self.sec = sec
        self.big_n = big_n
        osi.n_integ += 1
        self._tag = osi.n_integ
        self._parameters = [self.op_type, self._tag, self.sec.tag, self.big_n]
        self.to_process(osi)


class NewtonCotes(BeamIntegrationBase):
    """
    The NewtonCotes BeamIntegration Class
    
    Create a Newton-Cotes beamIntegration object.Newton-Cotes places integration points uniformly along the element,
    including a point ateach end of the element.spaced integration points are tabulated in references on numerical
    analysis. The force deformationresponse at each integration point is defined by the section.The order of
    accuracy for Gauss-Radau integration is N-1.Arguments and examples see `Lobatto-BeamIntegration`.
    """
    op_type = 'NewtonCotes'

    def __init__(self, osi, sec, big_n):
        """
        Initial method for NewtonCotes

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        sec: obj
            
        big_n: None
            

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> sec = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
        >>> o3.beam_integration.NewtonCotes(osi, sec=sec, big_n=1)
        """
        self.osi = osi
        self.sec = sec
        self.big_n = big_n
        osi.n_integ += 1
        self._tag = osi.n_integ
        self._parameters = [self.op_type, self._tag, self.sec.tag, self.big_n]
        self.to_process(osi)


class Radau(BeamIntegrationBase):
    """
    The Radau BeamIntegration Class
    
    Create a Gauss-Radau beamIntegration object.Gauss-Radau integration is not common in force-based elements because it
    places an integration point at only one end of the element; however, it forms the basis for optimal plastichinge
    integration methods.numerical analysis. The force-deformation response at each integration point is definedby
    the section. The order of accuracy for Gauss-Radau integration is 2N-2.Arguments and examples see
    `Lobatto-BeamIntegration`.
    """
    op_type = 'Radau'

    def __init__(self, osi, sec, big_n):
        """
        Initial method for Radau

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        sec: obj
            
        big_n: None
            

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> sec = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
        >>> o3.beam_integration.Radau(osi, sec=sec, big_n=1)
        """
        self.osi = osi
        self.sec = sec
        self.big_n = big_n
        osi.n_integ += 1
        self._tag = osi.n_integ
        self._parameters = [self.op_type, self._tag, self.sec.tag, self.big_n]
        self.to_process(osi)


class Trapezoidal(BeamIntegrationBase):
    """
    The Trapezoidal BeamIntegration Class
    
    Create a Trapezoidal beamIntegration object.Arguments and examples see `Lobatto-BeamIntegration`.
    """
    op_type = 'Trapezoidal'

    def __init__(self, osi, sec, big_n):
        """
        Initial method for Trapezoidal

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        sec: obj
            
        big_n: None
            

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> sec = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
        >>> o3.beam_integration.Trapezoidal(osi, sec=sec, big_n=1)
        """
        self.osi = osi
        self.sec = sec
        self.big_n = big_n
        osi.n_integ += 1
        self._tag = osi.n_integ
        self._parameters = [self.op_type, self._tag, self.sec.tag, self.big_n]
        self.to_process(osi)


class CompositeSimpson(BeamIntegrationBase):
    """
    The CompositeSimpson BeamIntegration Class
    
    Create a CompositeSimpson beamIntegration object.Arguments and examples see `Lobatto-BeamIntegration`.
    """
    op_type = 'CompositeSimpson'

    def __init__(self, osi, sec, big_n):
        """
        Initial method for CompositeSimpson

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        sec: obj
            
        big_n: None
            

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> sec = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
        >>> o3.beam_integration.CompositeSimpson(osi, sec=sec, big_n=1)
        """
        self.osi = osi
        self.sec = sec
        self.big_n = big_n
        osi.n_integ += 1
        self._tag = osi.n_integ
        self._parameters = [self.op_type, self._tag, self.sec.tag, self.big_n]
        self.to_process(osi)


class UserDefined(BeamIntegrationBase):
    """
    The UserDefined BeamIntegration Class
    
    Create a UserDefined beamIntegration object.This option allows user-specified locations and weights of the
    integration points.
    """
    op_type = 'UserDefined'

    def __init__(self, osi, big_n, secs, locs, wts):
        """
        Initial method for UserDefined

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        big_n: int
            Number of integration points along the element.
        secs: list
            A list previous-defined section objects.
        locs: list
            Locations of integration points along the element.
        wts: list
            Weights of integration points.

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> secs = [o3.section.Elastic2D(osi, e_mod=1.0, area=1.0, iz=1.0),
        >>>          o3.section.Elastic2D(osi, e_mod=1.0, area=1.0, iz=1.0)]
        >>> o3.beam_integration.UserDefined(osi, big_n=2, secs=secs, locs=[0.2, 0.9], wts=[0.5, 0.5])
        """
        self.osi = osi
        self.big_n = int(big_n)
        self.secs = [x.tag for x in secs]
        self.locs = locs
        self.wts = wts
        osi.n_integ += 1
        self._tag = osi.n_integ
        self._parameters = [self.op_type, self._tag, self.big_n, *self.secs, *self.locs, *self.wts]
        self.to_process(osi)


class FixedLocation(BeamIntegrationBase):
    r"""
    The FixedLocation BeamIntegration Class
    
    Create a FixedLocation beamIntegration object.This option allows user-specified locations of the integration points.
    The associated integrationweights are computed by the method of undetermined coefficients (Vandermondesystem)

    .. math::

      \sum^N_{i=1}x_i^{j-1}w_i = \int_0^1x^{j-1}dx = \frac{1}{j},\qquad (j=1,...,N)


    Note that `NewtonCotes-BeamIntegration` integration is recovered when the integration point locations are equally
    spaced.
    """
    op_type = 'FixedLocation'

    def __init__(self, osi, big_n, secs, locs):
        """
        Initial method for FixedLocation

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        big_n: int
            Number of integration points along the element.
        secs: list
            A list previous-defined section objects.
        locs: list
            Locations of integration points along the element.

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> secs = [o3.section.Elastic2D(osi, e_mod=1.0, area=1.0, iz=1.0),
        >>>         o3.section.Elastic2D(osi, e_mod=1.0, area=1.0, iz=1.0)]
        >>> o3.beam_integration.FixedLocation(osi, big_n=2, secs=secs, locs=[0.2, 0.9])
        """
        self.osi = osi
        self.big_n = int(big_n)
        self.secs = [x.tag for x in secs]
        self.locs = locs
        osi.n_integ += 1
        self._tag = osi.n_integ
        self._parameters = [self.op_type, self._tag, self.big_n, *self.secs, *self.locs]
        self.to_process(osi)


class LowOrder(BeamIntegrationBase):
    r"""
    The LowOrder BeamIntegration Class
    
    Create a LowOrder beamIntegration object.This option is a generalization of the `FixedLocation-BeamIntegration` and
    `UserDefined-BeamIntegration` integration approaches and is useful for moving load analysis (`Kidarsa, Scott and
    Higgins 2008`_). The locations of the integration points are user defined,while a selected number of weights
    are specified and the remaining weights arecomputed by the method of undetermined coefficients.

    .. math::

      \sum_{i=1}^{N_f}x_{fi}^{j-1}w_{fi}=\frac{1}{j}-\sum_{i=1}^{N_c}x_{ci}^{j-1}w_{ci}
    """
    op_type = 'LowOrder'

    def __init__(self, osi, big_n, secs, locs, wts):
        """
        Initial method for LowOrder

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        big_n: int
            Number of integration points along the element.
        secs: list
            A list previous-defined section objects.
        locs: list
            Locations of integration points along the element.
        wts: list
            Weights of integration points.

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> secs = [o3.section.Elastic2D(osi, e_mod=1.0, area=1.0, iz=1.0),
        >>>         o3.section.Elastic2D(osi, e_mod=1.0, area=1.0, iz=1.0)]
        >>> o3.beam_integration.LowOrder(osi, big_n=2, secs=secs, locs=[0.2, 0.9], wts=[0.5, 0.5])
        """
        self.osi = osi
        self.big_n = int(big_n)
        self.secs = [x.tag for x in secs]
        self.locs = locs
        self.wts = wts
        osi.n_integ += 1
        self._tag = osi.n_integ
        self._parameters = [self.op_type, self._tag, self.big_n, *self.secs, *self.locs, *self.wts]
        self.to_process(osi)


class MidDistance(BeamIntegrationBase):
    r"""
    The MidDistance BeamIntegration Class
    
    Create a MidDistance beamIntegration object.This option allows user-specified locations of the integration points.
    The associated integration weights are determined from the midpoints between adjacent integration point
    locations.:math:`w_i=(x_{i+1}-x_{i-1})/2` for :math:`i=2...N-1`, :math:`w_1=(x_1+x_2)/2`, and
    :math:`w_N=1-(x_{N-1}+x_N)/2`.
    """
    op_type = 'MidDistance'

    def __init__(self, osi, big_n, secs, locs):
        """
        Initial method for MidDistance

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        big_n: int
            Number of integration points along the element.
        secs: list
            A list previous-defined section objects.
        locs: list
            Locations of integration points along the element.

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> secs = [o3.section.Elastic2D(osi, e_mod=1.0, area=1.0, iz=1.0),
        >>>         o3.section.Elastic2D(osi, e_mod=1.0, area=1.0, iz=1.0)]
        >>> o3.beam_integration.MidDistance(osi, big_n=2, secs=secs, locs=[0.2, 0.9])
        """
        self.osi = osi
        self.big_n = int(big_n)
        self.secs = [x.tag for x in secs]
        self.locs = locs
        osi.n_integ += 1
        self._tag = osi.n_integ
        self._parameters = [self.op_type, self._tag, self.big_n, *self.secs, *self.locs]
        self.to_process(osi)


class UserHinge(BeamIntegrationBase):
    """
    The UserHinge BeamIntegration Class
    
    Create a UserHinge beamIntegration object.
    """
    op_type = 'UserHinge'

    def __init__(self, osi, sec_e, np_l, secs_ls, locs_l, wts_l, np_r, secs_rs, locs_r, wts_r):
        """
        Initial method for UserHinge

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        sec_e: obj
            A previous-defined section objects for non-hinge area.
        np_l: int
            Number of integration points along the left hinge.
        secs_ls: list
            A list of previous-defined section objects for left hinge area.
        locs_l: list
            A list of locations of integration points for left hinge area.
        wts_l: list
            A list of weights of integration points for left hinge area.
        np_r: int
            Number of integration points along the right hinge.
        secs_rs: list
            A list of previous-defined section objects for right hinge area.
        locs_r: list
            A list of locations of integration points for right hinge area.
        wts_r: list
            A list of weights of integration points for right hinge area.

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> sec_e = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
        >>> secs_l = [o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)]
        >>> secs_r = [o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)]
        >>> o3.beam_integration.UserHinge(osi, sec_e=sec_e, np_l=1, secs_ls=secs_l, locs_l=[1], wts_l=[1], np_r=1,
        >>>                               secs_rs=secs_r, locs_r=[1], wts_r=[1])
        """
        self.osi = osi
        self.sec_e = sec_e
        self.np_l = int(np_l)
        self.secs_ls = [x.tag for x in secs_ls]
        self.locs_l = locs_l
        self.wts_l = wts_l
        self.np_r = int(np_r)
        self.secs_rs = [x.tag for x in secs_rs]
        self.locs_r = locs_r
        self.wts_r = wts_r
        osi.n_integ += 1
        self._tag = osi.n_integ
        self._parameters = [self.op_type, self._tag, self.sec_e.tag, self.np_l, *self.secs_ls, *self.locs_l, *self.wts_l, self.np_r, *self.secs_rs, *self.locs_r, *self.wts_r]
        self.to_process(osi)


class HingeMidpoint(BeamIntegrationBase):
    """
    The HingeMidpoint BeamIntegration Class
    
    Create a HingeMidpoint beamIntegration object.Midpoint integration over each hinge region is the most accurate
    one-point integration rule;however, it does not place integration points at the element ends and there is a small
    integrationerror for linear curvature distributions along the element.
    """
    op_type = 'HingeMidpoint'

    def __init__(self, osi, sec_i, lp_i, sec_j, lp_j, sec_e):
        """
        Initial method for HingeMidpoint

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        sec_i: obj
            A previous-defined section object for hinge at i.
        lp_i: float
            The plastic hinge length at i.
        sec_j: obj
            A previous-defined section object for hinge at j.
        lp_j: float
            The plastic hinge length at j.
        sec_e: obj
            A previous-defined section object for the element interior.

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> sec_i = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
        >>> sec_j = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
        >>> sec_e = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
        >>> o3.beam_integration.HingeMidpoint(osi, sec_i=sec_i, lp_i=1.0, sec_j=sec_j, lp_j=1.0, sec_e=sec_e)
        """
        self.osi = osi
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
    """
    The HingeRadau BeamIntegration Class
    
    Create a HingeRadau beamIntegration object.Modified two-point Gauss-Radau integration over each hinge region places
    an integration point atthe element ends and at 8/3 the hinge length inside the element. This approach representslinear
    curvature distributions exactly and the characteristic length for softening plastic hinges is equal to the assumed
    palstic hinge length.Arguments and examples see `HingeMidPoint-BeamIntegration`.
    """
    op_type = 'HingeRadau'

    def __init__(self, osi, sec_i, lp_i, sec_j, lp_j, sec_e):
        """
        Initial method for HingeRadau

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        sec_i: obj
            
        lp_i: None
            
        sec_j: obj
            
        lp_j: None
            
        sec_e: obj
            

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> sec_i = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
        >>> sec_j = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
        >>> sec_e = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
        >>> o3.beam_integration.HingeRadau(osi, sec_i=sec_i, lp_i=1.0, sec_j=sec_j, lp_j=1.0, sec_e=sec_e)
        """
        self.osi = osi
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
    """
    The HingeRadauTwo BeamIntegration Class
    
    Create a HingeRadauTwo beamIntegration object.Two-point Gauss-Radau integration over each hinge region places an
    integrationpoint at the element ends and at 2/3 the hinge length inside the element. This approachrepresents linear
    curvature distributions exactly; however, the characteristic length for softeningplastic hinges is not equal to
    the assumed plastic hinge length (equals 1/4 of the plastic hinge length).Arguments and examples see
    `HingeMidPoint-BeamIntegration`.
    """
    op_type = 'HingeRadauTwo'

    def __init__(self, osi, sec_i, lp_i, sec_j, lp_j, sec_e):
        """
        Initial method for HingeRadauTwo

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        sec_i: obj
            
        lp_i: None
            
        sec_j: obj
            
        lp_j: None
            
        sec_e: obj
            

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> sec_i = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
        >>> sec_j = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
        >>> sec_e = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
        >>> o3.beam_integration.HingeRadauTwo(osi, sec_i=sec_i, lp_i=1.0, sec_j=sec_j, lp_j=1.0, sec_e=sec_e)
        """
        self.osi = osi
        self.sec_i = sec_i
        self.lp_i = lp_i
        self.sec_j = sec_j
        self.lp_j = lp_j
        self.sec_e = sec_e
        osi.n_integ += 1
        self._tag = osi.n_integ
        self._parameters = [self.op_type, self._tag, self.sec_i.tag, self.lp_i, self.sec_j.tag, self.lp_j, self.sec_e.tag]
        self.to_process(osi)


class BeamhingeEndpoint(OpenSeesObject):
    """
    The beamhingeEndpoint BeamhingeEndpoint Class
    
    Create a HingeEndpoint beamIntegration object.Endpoint integration over each hinge region moves the integration
    points to the element ends;however, there is a large integration error for linear curvature distributions along
    the element.
    """
    op_base_type = 'beamhingeEndpoint'

    def __init__(self, osi, lp_i, sec_j, lp_j, sec_e):
        """
        Initial method for beamhingeEndpoint

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        lp_i: float
            The plastic hinge length at i.
        sec_j: obj
            A previous-defined section object for hinge at j.
        lp_j: float
            The plastic hinge length at j.
        sec_e: obj
            A previous-defined section object for the element interior.
        """
        self.osi = osi
        self.lp_i = float(lp_i)
        self.sec_j = sec_j
        self.lp_j = float(lp_j)
        self.sec_e = sec_e
        osi.n_integ += 1
        self._tag = osi.n_integ
        self._parameters = [self._tag, self.lp_i, self.sec_j.tag, self.lp_j, self.sec_e.tag]
        self.to_process(osi)
