from o3seespy.command.element.base_element import ElementBase


class ElastomericBearingPlasticity2D(ElementBase):
    """
    The ElastomericBearingPlasticity2D Element Class
    
    This command is used to construct an elastomericBearing element object, which is defined by two nodes. The element
    can have zero length or the appropriate bearing height. The bearing has unidirectional (2D) or coupled (3D)
    plasticity properties for the shear deformations, and force-deformation behaviors defined by
    UniaxialMaterials in the remaining two (2D) or four (3D) directions. By default (sDratio =
    0.5) P-Delta moments are equally distributed to the two end-nodes. To avoid the
    introduction of artificial viscous damping in the isolation system (sometimes
    referred to as "damping leakage in the isolation system"), the bearing
    element does not contribute to the Rayleigh damping by default. If
    the element has non-zero length, the local x-axis is determined
    from the nodal geometry unless the optional x-axis vector is
    specified in which case the nodal geometry is ignored and
    the user-defined orientation is utilized.

    For a two-dimensional problem
    """
    op_type = 'elastomericBearingPlasticity'

    def __init__(self, osi, ele_nodes, k_init, qd, alpha1, alpha2, mu, p_mat=None, mz_mat=None, do_rayleigh=False, orient=None, mass: float=None, shear_dist: float=None):
        """
        Initial method for ElastomericBearingPlasticity2D

        Parameters
        ----------
        ele_nodes: listi
            A list of two element nodes
        k_init: float
            Initial elastic stiffness in local shear direction
        qd: float
            Characteristic strength
        alpha1: float
            Post yield stiffness ratio of linear hardening component
        alpha2: float
            Post yield stiffness ratio of non-linear hardening component
        mu: float
            Exponent of non-linear hardening component
        p_mat: obj
            Tag associated with previously-defined uniaxialmaterial in axial direction
        mz_mat: obj
            Tag associated with previously-defined uniaxialmaterial in moment direction around local z-axis
        do_rayleigh: str
            To include rayleigh damping from the bearing (optional, default = no rayleigh damping contribution)
        orient: None
            
        mass: float
            Element mass (optional, default = 0.0)
        shear_dist: float
            Shear distance from inode as a fraction of the element length (optional, default = 0.5)
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.k_init = float(k_init)
        self.qd = float(qd)
        self.alpha1 = float(alpha1)
        self.alpha2 = float(alpha2)
        self.mu = float(mu)
        self.p_mat = p_mat
        self.mz_mat = mz_mat
        self.do_rayleigh = do_rayleigh
        self.orient = orient
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        if shear_dist is None:
            self.shear_dist = None
        else:
            self.shear_dist = float(shear_dist)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.k_init, self.qd, self.alpha1, self.alpha2, self.mu]
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat.tag]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat.tag]
        if getattr(self, 'do_rayleigh'):
            self._parameters += ['-doRayleigh']
        if getattr(self, 'orient') is not None:
            self._parameters += ['-orient', *self.orient]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        if getattr(self, 'shear_dist') is not None:
            self._parameters += ['-shearDist', self.shear_dist]
        self.to_process(osi)


class ElastomericBearingPlasticity3D(ElementBase):
    """
    The ElastomericBearingPlasticity3D Element Class
    
    This command is used to construct an elastomericBearing element object, which is defined by two nodes. The element
    can have zero length or the appropriate bearing height. The bearing has unidirectional (2D) or coupled (3D)
    plasticity properties for the shear deformations, and force-deformation behaviors defined by
    UniaxialMaterials in the remaining two (2D) or four (3D) directions. By default (sDratio =
    0.5) P-Delta moments are equally distributed to the two end-nodes. To avoid the
    introduction of artificial viscous damping in the isolation system (sometimes
    referred to as "damping leakage in the isolation system"), the bearing
    element does not contribute to the Rayleigh damping by default. If
    the element has non-zero length, the local x-axis is determined
    from the nodal geometry unless the optional x-axis vector is
    specified in which case the nodal geometry is ignored and
    the user-defined orientation is utilized.

    For a three-dimensional problem
    """
    op_type = 'elastomericBearingPlasticity'

    def __init__(self, osi, ele_nodes, k_init, qd, alpha1, alpha2, mu, p_mat=None, t_mat=None, my_mat=None, mz_mat=None, do_rayleigh=False, orient=None, mass: float=None, shear_dist: float=None):
        """
        Initial method for ElastomericBearingPlasticity3D

        Parameters
        ----------
        ele_nodes: listi
            A list of two element nodes
        k_init: float
            Initial elastic stiffness in local shear direction
        qd: float
            Characteristic strength
        alpha1: float
            Post yield stiffness ratio of linear hardening component
        alpha2: float
            Post yield stiffness ratio of non-linear hardening component
        mu: float
            Exponent of non-linear hardening component
        p_mat: obj
            Tag associated with previously-defined uniaxialmaterial in axial direction
        t_mat: obj
            Tag associated with previously-defined uniaxialmaterial in torsional direction
        my_mat: obj
            Tag associated with previously-defined uniaxialmaterial in moment direction around local y-axis
        mz_mat: obj
            Tag associated with previously-defined uniaxialmaterial in moment direction around local z-axis
        do_rayleigh: str
            To include rayleigh damping from the bearing (optional, default = no rayleigh damping contribution)
        orient: None
            
        mass: float
            Element mass (optional, default = 0.0)
        shear_dist: float
            Shear distance from inode as a fraction of the element length (optional, default = 0.5)
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.k_init = float(k_init)
        self.qd = float(qd)
        self.alpha1 = float(alpha1)
        self.alpha2 = float(alpha2)
        self.mu = float(mu)
        self.p_mat = p_mat
        self.t_mat = t_mat
        self.my_mat = my_mat
        self.mz_mat = mz_mat
        self.do_rayleigh = do_rayleigh
        self.orient = orient
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        if shear_dist is None:
            self.shear_dist = None
        else:
            self.shear_dist = float(shear_dist)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.k_init, self.qd, self.alpha1, self.alpha2, self.mu]
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat.tag]
        if getattr(self, 't_mat') is not None:
            self._parameters += ['-T', self.t_mat.tag]
        if getattr(self, 'my_mat') is not None:
            self._parameters += ['-My', self.my_mat.tag]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat.tag]
        if getattr(self, 'do_rayleigh'):
            self._parameters += ['-doRayleigh']
        if getattr(self, 'orient') is not None:
            self._parameters += ['-orient', *self.orient]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        if getattr(self, 'shear_dist') is not None:
            self._parameters += ['-shearDist', self.shear_dist]
        self.to_process(osi)


class ElastomericBearingBoucWen2D(ElementBase):
    """
    The ElastomericBearingBoucWen2D Element Class
    
    This command is used to construct an elastomericBearing element object, which is defined by two nodes. The element
    can have zero length or the appropriate bearing height. The bearing has unidirectional (2D) or coupled (3D)
    plasticity properties for the shear deformations, and force-deformation behaviors defined by
    UniaxialMaterials in the remaining two (2D) or four (3D) directions. By default (sDratio =
    0.5) P-Delta moments are equally distributed to the two end-nodes. To avoid the
    introduction of artificial viscous damping in the isolation system (sometimes
    referred to as "damping leakage in the isolation system"), the bearing
    element does not contribute to the Rayleigh damping by default. If
    the element has non-zero length, the local x-axis is determined
    from the nodal geometry unless the optional x-axis vector is
    specified in which case the nodal geometry is ignored and
    the user-defined orientation is utilized.

    For a two-dimensional problem
    """
    op_type = 'elastomericBearingBoucWen'

    def __init__(self, osi, ele_nodes, k_init, qd, alpha1, alpha2, mu, eta, beta, gamma, p_mat=None, mz_mat=None, orient_vals=None, shear_dist: float=None, do_rayleigh=False, mass: float=None):
        """
        Initial method for ElastomericBearingBoucWen2D

        Parameters
        ----------
        ele_nodes: listi
            A list of two element nodes
        k_init: float
            Initial elastic stiffness in local shear direction
        qd: float
            Characteristic strength
        alpha1: float
            Post yield stiffness ratio of linear hardening component
        alpha2: float
            Post yield stiffness ratio of non-linear hardening component
        mu: float
            Exponent of non-linear hardening component
        eta: float
            Yielding exponent (sharpness of hysteresis loop corners) (default = 1.0)
        beta: float
            First hysteretic shape parameter (default = 0.5)
        gamma: float
            Second hysteretic shape parameter (default = 0.5)
        p_mat: obj
            Tag associated with previously-defined uniaxialmaterial in axial direction
        mz_mat: obj
            Tag associated with previously-defined uniaxialmaterial in moment direction around local z-axis
        orient_vals: listi
            Vector components in global coordinates defining local x-axis (optional), vector components in global
            coordinates defining local y-axis (optional)
        shear_dist: float
            Shear distance from inode as a fraction of the element length (optional, default = 0.5)
        do_rayleigh: str
            To include rayleigh damping from the bearing (optional, default = no rayleigh damping contribution)
        mass: float
            Element mass (optional, default = 0.0)
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.k_init = float(k_init)
        self.qd = float(qd)
        self.alpha1 = float(alpha1)
        self.alpha2 = float(alpha2)
        self.mu = float(mu)
        self.eta = float(eta)
        self.beta = float(beta)
        self.gamma = float(gamma)
        self.p_mat = p_mat
        self.mz_mat = mz_mat
        self.orient_vals = orient_vals
        if shear_dist is None:
            self.shear_dist = None
        else:
            self.shear_dist = float(shear_dist)
        self.do_rayleigh = do_rayleigh
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.k_init, self.qd, self.alpha1, self.alpha2, self.mu, self.eta, self.beta, self.gamma]
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat.tag]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat.tag]
        if getattr(self, 'orient_vals') is not None:
            self._parameters += ['-orient', *self.orient_vals]
        if getattr(self, 'shear_dist') is not None:
            self._parameters += ['-shearDist', self.shear_dist]
        if getattr(self, 'do_rayleigh'):
            self._parameters += ['-doRayleigh']
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        self.to_process(osi)


class ElastomericBearingBoucWen3D(ElementBase):
    """
    The ElastomericBearingBoucWen3D Element Class
    
    This command is used to construct an elastomericBearing element object, which is defined by two nodes. The element
    can have zero length or the appropriate bearing height. The bearing has unidirectional (2D) or coupled (3D)
    plasticity properties for the shear deformations, and force-deformation behaviors defined by
    UniaxialMaterials in the remaining two (2D) or four (3D) directions. By default (sDratio =
    0.5) P-Delta moments are equally distributed to the two end-nodes. To avoid the
    introduction of artificial viscous damping in the isolation system (sometimes
    referred to as "damping leakage in the isolation system"), the bearing
    element does not contribute to the Rayleigh damping by default. If
    the element has non-zero length, the local x-axis is determined
    from the nodal geometry unless the optional x-axis vector is
    specified in which case the nodal geometry is ignored and
    the user-defined orientation is utilized.

    For a three-dimensional problem
    """
    op_type = 'elastomericBearingBoucWen'

    def __init__(self, osi, ele_nodes, k_init, qd, alpha1, alpha2, mu, eta, beta, gamma, p_mat=None, t_mat=None, my_mat=None, mz_mat=None, orient_vals=None, shear_dist: float=None, do_rayleigh=False, mass: float=None):
        """
        Initial method for ElastomericBearingBoucWen3D

        Parameters
        ----------
        ele_nodes: listi
            A list of two element nodes
        k_init: float
            Initial elastic stiffness in local shear direction
        qd: float
            Characteristic strength
        alpha1: float
            Post yield stiffness ratio of linear hardening component
        alpha2: float
            Post yield stiffness ratio of non-linear hardening component
        mu: float
            Exponent of non-linear hardening component
        eta: float
            Yielding exponent (sharpness of hysteresis loop corners) (default = 1.0)
        beta: float
            First hysteretic shape parameter (default = 0.5)
        gamma: float
            Second hysteretic shape parameter (default = 0.5)
        p_mat: obj
            Tag associated with previously-defined uniaxialmaterial in axial direction
        t_mat: obj
            Tag associated with previously-defined uniaxialmaterial in torsional direction
        my_mat: obj
            Tag associated with previously-defined uniaxialmaterial in moment direction around local y-axis
        mz_mat: obj
            Tag associated with previously-defined uniaxialmaterial in moment direction around local z-axis
        orient_vals: listi
            Vector components in global coordinates defining local x-axis (optional), vector components in global
            coordinates defining local y-axis (optional)
        shear_dist: float
            Shear distance from inode as a fraction of the element length (optional, default = 0.5)
        do_rayleigh: str
            To include rayleigh damping from the bearing (optional, default = no rayleigh damping contribution)
        mass: float
            Element mass (optional, default = 0.0)
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.k_init = float(k_init)
        self.qd = float(qd)
        self.alpha1 = float(alpha1)
        self.alpha2 = float(alpha2)
        self.mu = float(mu)
        self.eta = float(eta)
        self.beta = float(beta)
        self.gamma = float(gamma)
        self.p_mat = p_mat
        self.t_mat = t_mat
        self.my_mat = my_mat
        self.mz_mat = mz_mat
        self.orient_vals = orient_vals
        if shear_dist is None:
            self.shear_dist = None
        else:
            self.shear_dist = float(shear_dist)
        self.do_rayleigh = do_rayleigh
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.k_init, self.qd, self.alpha1, self.alpha2, self.mu, self.eta, self.beta, self.gamma]
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat.tag]
        if getattr(self, 't_mat') is not None:
            self._parameters += ['-T', self.t_mat.tag]
        if getattr(self, 'my_mat') is not None:
            self._parameters += ['-My', self.my_mat.tag]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat.tag]
        if getattr(self, 'orient_vals') is not None:
            self._parameters += ['-orient', *self.orient_vals]
        if getattr(self, 'shear_dist') is not None:
            self._parameters += ['-shearDist', self.shear_dist]
        if getattr(self, 'do_rayleigh'):
            self._parameters += ['-doRayleigh']
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        self.to_process(osi)


class FlatSliderBearingmaxIter2D(ElementBase):
    """
    The FlatSliderBearingmaxIter2D Element Class
    
    This command is used to construct a flatSliderBearing element object, which is defined by two nodes. The iNode
    represents the flat sliding surface and the jNode represents the slider. The element can have zero length or the
    appropriate bearing height. The bearing has unidirectional (2D) or coupled (3D) friction properties for the
    shear deformations, and force-deformation behaviors defined by UniaxialMaterials in the remaining two (2D)
    or four (3D) directions. To capture the uplift behavior of the bearing, the user-specified
    UniaxialMaterial in the axial direction is modified for no-tension behavior. By default
    (sDratio = 0.0) P-Delta moments are entirely transferred to the flat sliding surface
    (iNode). It is important to note that rotations of the flat sliding surface
    (rotations at the iNode) affect the shear behavior of the bearing. To
    avoid the introduction of artificial viscous damping in the
    isolation system (sometimes referred to as "damping
    leakage in the isolation system"), the bearing
    element does not contribute to the Rayleigh
    damping by default. If the element has
    non-zero length, the local x-axis is
    determined from the nodal geometry
    unless the optional x-axis vector
    is specified in which case the
    nodal geometry is ignored and the user-defined orientation is utilized.

    For a two-dimensional problem
    """
    op_type = 'flatSliderBearing'

    def __init__(self, osi, ele_nodes, frn_mdl, k_init, iter, tol, p_mat=None, mz_mat=None, do_rayleigh=False, orient=None, mass: float=None, shear_dist: float=None):
        """
        Initial method for FlatSliderBearingmaxIter2D

        Parameters
        ----------
        ele_nodes: listi
            A list of two element nodes
        frn_mdl: obj
            Tag associated with previously-defined frictionmodel
        k_init: float
            Initial elastic stiffness in local shear direction
        iter: int
            Maximum number of iterations to undertake to satisfy element equilibrium (optional, default = 20)
        tol: float
            Convergence tolerance to satisfy element equilibrium (optional, default = 1e-8)
        p_mat: obj
            Tag associated with previously-defined uniaxialmaterial in axial direction
        mz_mat: obj
            Tag associated with previously-defined uniaxialmaterial in moment direction around local z-axis
        do_rayleigh: str
            To include rayleigh damping from the bearing (optional, default = no rayleigh damping contribution)
        orient: None
            
        mass: float
            Element mass (optional, default = 0.0)
        shear_dist: float
            Shear distance from inode as a fraction of the element length (optional, default = 0.0)
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.frn_mdl = frn_mdl
        self.k_init = float(k_init)
        self.p_mat = p_mat
        self.mz_mat = mz_mat
        self.do_rayleigh = do_rayleigh
        self.iter = int(iter)
        self.tol = float(tol)
        self.orient = orient
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        if shear_dist is None:
            self.shear_dist = None
        else:
            self.shear_dist = float(shear_dist)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.frn_mdl.tag, self.k_init, '-maxIter', self.iter, self.tol]
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat.tag]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat.tag]
        if getattr(self, 'do_rayleigh'):
            self._parameters += ['-doRayleigh']
        if getattr(self, 'orient') is not None:
            self._parameters += ['-orient', *self.orient]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        if getattr(self, 'shear_dist') is not None:
            self._parameters += ['-shearDist', self.shear_dist]
        self.to_process(osi)


class FlatSliderBearing3D(ElementBase):
    """
    The FlatSliderBearing3D Element Class
    
    This command is used to construct a flatSliderBearing element object, which is defined by two nodes. The iNode
    represents the flat sliding surface and the jNode represents the slider. The element can have zero length or the
    appropriate bearing height. The bearing has unidirectional (2D) or coupled (3D) friction properties for the
    shear deformations, and force-deformation behaviors defined by UniaxialMaterials in the remaining two (2D)
    or four (3D) directions. To capture the uplift behavior of the bearing, the user-specified
    UniaxialMaterial in the axial direction is modified for no-tension behavior. By default
    (sDratio = 0.0) P-Delta moments are entirely transferred to the flat sliding surface
    (iNode). It is important to note that rotations of the flat sliding surface
    (rotations at the iNode) affect the shear behavior of the bearing. To
    avoid the introduction of artificial viscous damping in the
    isolation system (sometimes referred to as "damping
    leakage in the isolation system"), the bearing
    element does not contribute to the Rayleigh
    damping by default. If the element has
    non-zero length, the local x-axis is
    determined from the nodal geometry
    unless the optional x-axis vector
    is specified in which case the
    nodal geometry is ignored and the user-defined orientation is utilized.

    For a three-dimensional problem
    """
    op_type = 'flatSliderBearing'

    def __init__(self, osi, ele_nodes, frn_mdl, k_init, p_mat=None, t_mat=None, my_mat=None, mz_mat=None, do_rayleigh=False, max_iter=None, tol: float=None, orient=None, mass: float=None, shear_dist: float=None):
        """
        Initial method for FlatSliderBearing3D

        Parameters
        ----------
        ele_nodes: listi
            A list of two element nodes
        frn_mdl: obj
            Tag associated with previously-defined frictionmodel
        k_init: float
            Initial elastic stiffness in local shear direction
        p_mat: obj
            Tag associated with previously-defined uniaxialmaterial in axial direction
        t_mat: obj
            Tag associated with previously-defined uniaxialmaterial in torsional direction
        my_mat: obj
            Tag associated with previously-defined uniaxialmaterial in moment direction around local y-axis
        mz_mat: obj
            Tag associated with previously-defined uniaxialmaterial in moment direction around local z-axis
        do_rayleigh: str
            To include rayleigh damping from the bearing (optional, default = no rayleigh damping contribution)
        max_iter: None
            
        tol: float
            Convergence tolerance to satisfy element equilibrium (optional, default = 1e-8)
        orient: None
            
        mass: float
            Element mass (optional, default = 0.0)
        shear_dist: float
            Shear distance from inode as a fraction of the element length (optional, default = 0.0)
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.frn_mdl = frn_mdl
        self.k_init = float(k_init)
        self.p_mat = p_mat
        self.t_mat = t_mat
        self.my_mat = my_mat
        self.mz_mat = mz_mat
        self.do_rayleigh = do_rayleigh
        self.max_iter = max_iter
        if tol is None:
            self.tol = None
        else:
            self.tol = float(tol)
        self.orient = orient
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        if shear_dist is None:
            self.shear_dist = None
        else:
            self.shear_dist = float(shear_dist)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.frn_mdl.tag, self.k_init]
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat.tag]
        if getattr(self, 't_mat') is not None:
            self._parameters += ['-T', self.t_mat.tag]
        if getattr(self, 'my_mat') is not None:
            self._parameters += ['-My', self.my_mat.tag]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat.tag]
        if getattr(self, 'do_rayleigh'):
            self._parameters += ['-doRayleigh']
        if getattr(self, 'max_iter') is not None:
            self._parameters += ['-iter', self.max_iter]
        if getattr(self, 'tol') is not None:
            if getattr(self, 'max_iter') is None:
                raise ValueError('Cannot set: tol and not: max_iter')
            self._parameters += [self.tol]
        if getattr(self, 'orient') is not None:
            self._parameters += ['-orient', *self.orient]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        if getattr(self, 'shear_dist') is not None:
            self._parameters += ['-shearDist', self.shear_dist]
        self.to_process(osi)


class SingleFPBearing2D(ElementBase):
    """
    The SingleFPBearing2D Element Class
    
    This command is used to construct a singleFPBearing element object, which is defined by two nodes. The iNode
    represents the concave sliding surface and the jNode represents the articulated slider. The element can have
    zero length or the appropriate bearing height. The bearing has unidirectional (2D) or coupled (3D) friction
    properties (with post-yield stiffening due to the concave sliding surface) for the shear deformations, and
    force-deformation behaviors defined by UniaxialMaterials in the remaining two (2D) or four (3D)
    directions. To capture the uplift behavior of the bearing, the user-specified UniaxialMaterial
    in the axial direction is modified for no-tension behavior. By default (sDratio = 0.0)
    P-Delta moments are entirely transferred to the concave sliding surface (iNode). It
    is important to note that rotations of the concave sliding surface (rotations at
    the iNode) affect the shear behavior of the bearing. To avoid the introduction
    of artificial viscous damping in the isolation system (sometimes referred to
    as "damping leakage in the isolation system"), the bearing element does not
    contribute to the Rayleigh damping by default. If the element has non-zero
    length, the local x-axis is determined from the nodal geometry unless the
    optional x-axis vector is specified in which case the nodal geometry is
    ignored and the user-defined orientation is utilized.

    For a two-dimensional problem
    """
    op_type = 'singleFPBearing'

    def __init__(self, osi, ele_nodes, frn_mdl, reff, k_init, p_mat=None, mz_mat=None, do_rayleigh=False, max_iter: int=None, tol: float=None, orient=None, mass: float=None, shear_dist: float=None):
        """
        Initial method for SingleFPBearing2D

        Parameters
        ----------
        ele_nodes: listi
            A list of two element nodes
        frn_mdl: obj
            Tag associated with previously-defined frictionmodel
        reff: float
            Effective radius of concave sliding surface
        k_init: float
            Initial elastic stiffness in local shear direction
        p_mat: obj
            Tag associated with previously-defined uniaxialmaterial in axial direction
        mz_mat: obj
            Tag associated with previously-defined uniaxialmaterial in moment direction around local z-axis
        do_rayleigh: str
            To include rayleigh damping from the bearing (optional, default = no rayleigh damping contribution)
        max_iter: int
            Maximum number of iterations to undertake to satisfy element equilibrium (optional, default = 20)
        tol: float
            Convergence tolerance to satisfy element equilibrium (optional, default = 1e-8)
        orient: None
            
        mass: float
            Element mass (optional, default = 0.0)
        shear_dist: float
            Shear distance from inode as a fraction of the element length (optional, default = 0.0)
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.frn_mdl = frn_mdl
        self.reff = float(reff)
        self.k_init = float(k_init)
        self.p_mat = p_mat
        self.mz_mat = mz_mat
        self.do_rayleigh = do_rayleigh
        if max_iter is None:
            self.max_iter = None
        else:
            self.max_iter = int(max_iter)
        if tol is None:
            self.tol = None
        else:
            self.tol = float(tol)
        self.orient = orient
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        if shear_dist is None:
            self.shear_dist = None
        else:
            self.shear_dist = float(shear_dist)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.frn_mdl.tag, self.reff, self.k_init]
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat.tag]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat.tag]
        if getattr(self, 'do_rayleigh'):
            self._parameters += ['-doRayleigh']
        if getattr(self, 'max_iter') is not None:
            self._parameters += ['-iter', self.max_iter]
        if getattr(self, 'tol') is not None:
            if getattr(self, 'max_iter') is None:
                raise ValueError('Cannot set: tol and not: max_iter')
            self._parameters += [self.tol]
        if getattr(self, 'orient') is not None:
            self._parameters += ['-orient', *self.orient]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        if getattr(self, 'shear_dist') is not None:
            self._parameters += ['-shearDist', self.shear_dist]
        self.to_process(osi)


class SingleFPBearing3D(ElementBase):
    """
    The SingleFPBearing3D Element Class
    
    This command is used to construct a singleFPBearing element object, which is defined by two nodes. The iNode
    represents the concave sliding surface and the jNode represents the articulated slider. The element can have
    zero length or the appropriate bearing height. The bearing has unidirectional (2D) or coupled (3D) friction
    properties (with post-yield stiffening due to the concave sliding surface) for the shear deformations, and
    force-deformation behaviors defined by UniaxialMaterials in the remaining two (2D) or four (3D)
    directions. To capture the uplift behavior of the bearing, the user-specified UniaxialMaterial
    in the axial direction is modified for no-tension behavior. By default (sDratio = 0.0)
    P-Delta moments are entirely transferred to the concave sliding surface (iNode). It
    is important to note that rotations of the concave sliding surface (rotations at
    the iNode) affect the shear behavior of the bearing. To avoid the introduction
    of artificial viscous damping in the isolation system (sometimes referred to
    as "damping leakage in the isolation system"), the bearing element does not
    contribute to the Rayleigh damping by default. If the element has non-zero
    length, the local x-axis is determined from the nodal geometry unless the
    optional x-axis vector is specified in which case the nodal geometry is
    ignored and the user-defined orientation is utilized.

    For a three-dimensional problem
    """
    op_type = 'singleFPBearing'

    def __init__(self, osi, ele_nodes, frn_mdl, reff, k_init, p_mat=None, t_mat=None, my_mat=None, mz_mat=None, do_rayleigh=False, max_iter: int=None, tol: float=None, orient=None, mass: float=None, shear_dist: float=None):
        """
        Initial method for SingleFPBearing3D

        Parameters
        ----------
        ele_nodes: listi
            A list of two element nodes
        frn_mdl: obj
            Tag associated with previously-defined frictionmodel
        reff: float
            Effective radius of concave sliding surface
        k_init: float
            Initial elastic stiffness in local shear direction
        p_mat: obj
            Tag associated with previously-defined uniaxialmaterial in axial direction
        t_mat: obj
            Tag associated with previously-defined uniaxialmaterial in torsional direction
        my_mat: obj
            Tag associated with previously-defined uniaxialmaterial in moment direction around local y axis
        mz_mat: obj
            Tag associated with previously-defined uniaxialmaterial in moment direction around local z-axis
        do_rayleigh: str
            To include rayleigh damping from the bearing (optional, default = no rayleigh damping contribution)
        max_iter: int
            Maximum number of iterations to undertake to satisfy element equilibrium (optional, default = 20)
        tol: float
            Convergence tolerance to satisfy element equilibrium (optional, default = 1e-8)
        orient: None
            
        mass: float
            Element mass (optional, default = 0.0)
        shear_dist: float
            Shear distance from inode as a fraction of the element length (optional, default = 0.0)
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.frn_mdl = frn_mdl
        self.reff = float(reff)
        self.k_init = float(k_init)
        self.p_mat = p_mat
        self.t_mat = t_mat
        self.my_mat = my_mat
        self.mz_mat = mz_mat
        self.do_rayleigh = do_rayleigh
        if max_iter is None:
            self.max_iter = None
        else:
            self.max_iter = int(max_iter)
        if tol is None:
            self.tol = None
        else:
            self.tol = float(tol)
        self.orient = orient
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        if shear_dist is None:
            self.shear_dist = None
        else:
            self.shear_dist = float(shear_dist)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.frn_mdl.tag, self.reff, self.k_init]
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat.tag]
        if getattr(self, 't_mat') is not None:
            self._parameters += ['-T', self.t_mat.tag]
        if getattr(self, 'my_mat') is not None:
            self._parameters += ['-My', self.my_mat.tag]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat.tag]
        if getattr(self, 'do_rayleigh'):
            self._parameters += ['-doRayleigh']
        if getattr(self, 'max_iter') is not None:
            self._parameters += ['-iter', self.max_iter]
        if getattr(self, 'tol') is not None:
            if getattr(self, 'max_iter') is None:
                raise ValueError('Cannot set: tol and not: max_iter')
            self._parameters += [self.tol]
        if getattr(self, 'orient') is not None:
            self._parameters += ['-orient', *self.orient]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        if getattr(self, 'shear_dist') is not None:
            self._parameters += ['-shearDist', self.shear_dist]
        self.to_process(osi)


class TFP(ElementBase):
    """
    The TFP Element Class
    
    This command is used to construct a Triple Friction Pendulum Bearing element object, which is defined by two nodes.
    The element can have zero length or the appropriate bearing height. The bearing has unidirectional (2D) or coupled
    (3D) friction properties (with post-yield stiffening due to the concave sliding surface) for the shear
    deformations, and force-deformation behaviors defined by UniaxialMaterials in the remaining two (2D)
    or four (3D) directions. To capture the uplift behavior of the bearing, the user-specified
    UniaxialMaterial in the axial direction is modified for no-tension behavior. P-Delta
    moments are entirely transferred to the concave sliding surface (iNode). It is
    important to note that rotations of the concave sliding surface (rotations at
    the iNode) affect the shear behavior of the bearing. If the element has
    non-zero length, the local x-axis is determined from the nodal
    geometry unless the optional x-axis vector is specified in
    which case the nodal geometry is ignored and the user-defined orientation is utilized.

    
    """
    op_type = 'TFP'

    def __init__(self, osi, ele_nodes, r1, r2, r3, r4, d1, d2, d3, d4, d1, d2, d3, d4, mu1, mu2, mu3, mu4, h1, h2, h3, h4, h0, col_load, big_k):
        """
        Initial method for TFP

        Parameters
        ----------
        ele_nodes: listi
            A list of two element nodes
        r1: float
            Radius of inner bottom sliding surface
        r2: float
            Radius of inner top sliding surface
        r3: float
            Radius of outer bottom sliding surface
        r4: float
            Radius of outer top sliding surface
        d1: float
            Diameter of inner bottom sliding surface
        d2: float
            Diameter of inner top sliding surface
        d3: float
            Diameter of outer bottom sliding surface
        d4: float
            Diameter of outer top sliding surface
        d1: float
            Diameter of inner slider
        d2: float
            Diameter of inner slider
        d3: float
            Diameter of outer bottom slider
        d4: float
            Diameter of outer top slider
        mu1: float
            Friction coefficient of inner bottom sliding surface
        mu2: float
            Friction coefficient of inner top sliding surface
        mu3: float
            Friction coefficient of outer bottom sliding surface
        mu4: float
            Friction coefficient of outer top sliding surface
        h1: float
            Height from inner bottom sliding surface to center of bearing
        h2: float
            Height from inner top sliding surface to center of bearing
        h3: float
            Height from outer bottom sliding surface to center of bearing
        h4: float
            Height from inner top sliding surface to center of bearing
        h0: float
            Total height of bearing
        col_load: float
            Initial axial load on bearing (only used for first time step then load come from model)
        big_k: float
            Optional, stiffness of spring in vertical dirn (dof 2 if ndm= 2, dof 3 if ndm = 3) (default=1.0e15)
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.r1 = float(r1)
        self.r2 = float(r2)
        self.r3 = float(r3)
        self.r4 = float(r4)
        self.d1 = float(d1)
        self.d2 = float(d2)
        self.d3 = float(d3)
        self.d4 = float(d4)
        self.d1 = float(d1)
        self.d2 = float(d2)
        self.d3 = float(d3)
        self.d4 = float(d4)
        self.mu1 = float(mu1)
        self.mu2 = float(mu2)
        self.mu3 = float(mu3)
        self.mu4 = float(mu4)
        self.h1 = float(h1)
        self.h2 = float(h2)
        self.h3 = float(h3)
        self.h4 = float(h4)
        self.h0 = float(h0)
        self.col_load = float(col_load)
        self.big_k = float(big_k)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.r1, self.r2, self.r3, self.r4, self.d1, self.d2, self.d3, self.d4, self.d1, self.d2, self.d3, self.d4, self.mu1, self.mu2, self.mu3, self.mu4, self.h1, self.h2, self.h3, self.h4, self.h0, self.col_load, self.big_k]
        self.to_process(osi)


class TripleFrictionPendulum(ElementBase):
    """
    The TripleFrictionPendulum Element Class
    
    
    """
    op_type = 'TripleFrictionPendulum'

    def __init__(self, osi, ele_nodes, frn_tag1, frn_tag2, frn_tag3, vert_mat, rot_z_mat, rot_x_mat, rot_y_mat, l1, l2, l3, d1, d2, d3, big_w, uy, kvt, min_fv, tol):
        """
        Initial method for TripleFrictionPendulum

        Parameters
        ----------
        ele_nodes: listi
            A list of two element nodes
        frn_tag1: int
            = tags associated with previously-defined frictionmodels at the three sliding interfaces
        frn_tag2: int
            = tags associated with previously-defined frictionmodels at the three sliding interfaces
        frn_tag3: int
            = tags associated with previously-defined frictionmodels at the three sliding interfaces
        vert_mat: obj
            = pre-defined material tag for compression behavior of the bearing
        rot_z_mat: obj
            = pre-defined material tags for rotational behavior about 3-axis, 1-axis and 2-axis, respectively.
        rot_x_mat: obj
            = pre-defined material tags for rotational behavior about 3-axis, 1-axis and 2-axis, respectively.
        rot_y_mat: obj
            = pre-defined material tags for rotational behavior about 3-axis, 1-axis and 2-axis, respectively.
        l1: float
            = effective radii. li = r_i - h_i (see figure 1)
        l2: float
            = effective radii. li = r_i - h_i (see figure 1)
        l3: float
            = effective radii. li = r_i - h_i (see figure 1)
        d1: float
            = displacement limits of pendulums (figure 1). displacement limit of the bearing is 2   ``d1`` +   ``d2`` + 
             ``d3`` +   ``l1``.   ``d3``/   ``l3`` -   ``l1``.   ``d2``/   ``l2``
        d2: float
            = displacement limits of pendulums (figure 1). displacement limit of the bearing is 2   ``d1`` +   ``d2`` + 
             ``d3`` +   ``l1``.   ``d3``/   ``l3`` -   ``l1``.   ``d2``/   ``l2``
        d3: float
            = displacement limits of pendulums (figure 1). displacement limit of the bearing is 2   ``d1`` +   ``d2`` + 
             ``d3`` +   ``l1``.   ``d3``/   ``l3`` -   ``l1``.   ``d2``/   ``l2``
        big_w: float
            = axial force used for the first trial of the first analysis step.
        uy: float
            = lateral displacement where sliding of the bearing starts. recommended value = 0.25 to 1 mm. a smaller
            value may cause convergence problem.
        kvt: float
            = tension stiffness k_vt of the bearing.
        min_fv: None
            
        tol: float
            = relative tolerance for checking the convergence of the element. recommended value = 1.e-10 to 1.e-3.
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.frn_tag1 = int(frn_tag1)
        self.frn_tag2 = int(frn_tag2)
        self.frn_tag3 = int(frn_tag3)
        self.vert_mat = vert_mat
        self.rot_z_mat = rot_z_mat
        self.rot_x_mat = rot_x_mat
        self.rot_y_mat = rot_y_mat
        self.l1 = float(l1)
        self.l2 = float(l2)
        self.l3 = float(l3)
        self.d1 = float(d1)
        self.d2 = float(d2)
        self.d3 = float(d3)
        self.big_w = float(big_w)
        self.uy = float(uy)
        self.kvt = float(kvt)
        self.min_fv = min_fv
        self.tol = float(tol)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.frn_tag1, self.frn_tag2, self.frn_tag3, self.vert_mat.tag, self.rot_z_mat.tag, self.rot_x_mat.tag, self.rot_y_mat.tag, self.l1, self.l2, self.l3, self.d1, self.d2, self.d3, self.big_w, self.uy, self.kvt, self.min_fv, self.tol]
        self.to_process(osi)


class MultipleShearSpring(ElementBase):
    """
    The MultipleShearSpring Element Class
    
    This command is used to construct a multipleShearSpring (MSS) element object, which is defined by two nodes. This
    element consists of a series of identical shear springs arranged radially to represent the isotropic behavior in the
    local y-z plane.

    
    """
    op_type = 'multipleShearSpring'

    def __init__(self, osi, ele_nodes, n_spring, mat=None, lim: float=None, mass: float=None, orient=None):
        """
        Initial method for MultipleShearSpring

        Parameters
        ----------
        ele_nodes: listi
            A list of two element nodes
        n_spring: int
            Number of springs
        mat: obj
            Tag associated with previously-defined uniaxialmaterial object
        lim: float
            Minimum deformation to calculate equivalent coefficient (see note 1)
        mass: float
            Element mass
        orient: None
            
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.n_spring = int(n_spring)
        self.mat = mat
        if lim is None:
            self.lim = None
        else:
            self.lim = float(lim)
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        self.orient = orient
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.n_spring]
        if getattr(self, 'mat') is not None:
            self._parameters += ['-mat', self.mat.tag]
        if getattr(self, 'lim') is not None:
            self._parameters += ['-lim', self.lim]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        if getattr(self, 'orient') is not None:
            self._parameters += ['-orient', *self.orient]
        self.to_process(osi)


class KikuchiBearingadjustPDOutput(ElementBase):
    """
    The KikuchiBearingadjustPDOutput Element Class
    
    This command is used to construct a KikuchiBearing element object, which is defined by two nodes. This element
    consists of multiple shear spring model (MSS) and multiple normal spring model (MNS).

    
    """
    op_type = 'KikuchiBearing'

    def __init__(self, osi, ele_nodes, total_rubber, ci, cj, shape: float=None, size: float=None, total_height: float=None, n_mss: int=None, mat_mss=None, lim_disp: float=None, n_mns: int=None, mat_mns=None, lamb: float=None, no_pd_input=False, no_tilt=False, orient=None, mass: float=None):
        """
        Initial method for KikuchiBearingadjustPDOutput

        Parameters
        ----------
        ele_nodes: listi
            A list of two element nodes
        total_rubber: float
            Total rubber thickness
        ci: float
            P-delta moment adjustment for reaction force (default:    ``ci`` =0.5,    ``cj`` =0.5)
        cj: float
            P-delta moment adjustment for reaction force (default:    ``ci`` =0.5,    ``cj`` =0.5)
        shape: float
            Following shapes are available: round, square
        size: float
            Diameter (round shape), length of edge (square shape)
        total_height: float
            Total height of the bearing (defaulut: distance between inode and jnode)
        n_mss: int
            Number of springs in mss = nmss
        mat_mss: obj
            Mattag for mss
        lim_disp: float
            Minimum deformation to calculate equivalent coefficient of mss (see note 1)
        n_mns: int
            Number of springs in mns = nmns*nmns (for round and square shape)
        mat_mns: obj
            Mattag for mns
        lamb: float
            Parameter to calculate compression modulus distribution on mns (see note 2)
        no_pd_input: str
            Not consider p-delta moment
        no_tilt: str
            Not consider tilt of rigid link
        orient: None
            
        mass: float
            Element mass
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        if shape is None:
            self.shape = None
        else:
            self.shape = float(shape)
        if size is None:
            self.size = None
        else:
            self.size = float(size)
        self.total_rubber = float(total_rubber)
        if total_height is None:
            self.total_height = None
        else:
            self.total_height = float(total_height)
        if n_mss is None:
            self.n_mss = None
        else:
            self.n_mss = int(n_mss)
        self.mat_mss = mat_mss
        if lim_disp is None:
            self.lim_disp = None
        else:
            self.lim_disp = float(lim_disp)
        if n_mns is None:
            self.n_mns = None
        else:
            self.n_mns = int(n_mns)
        self.mat_mns = mat_mns
        if lamb is None:
            self.lamb = None
        else:
            self.lamb = float(lamb)
        self.no_pd_input = no_pd_input
        self.no_tilt = no_tilt
        self.ci = float(ci)
        self.cj = float(cj)
        self.orient = orient
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.total_rubber, '-adjustPDOutput', self.ci, self.cj]
        if getattr(self, 'shape') is not None:
            self._parameters += ['-shape', self.shape]
        if getattr(self, 'size') is not None:
            self._parameters += ['-size', self.size]
        if getattr(self, 'total_height') is not None:
            self._parameters += ['-totalHeight', self.total_height]
        if getattr(self, 'n_mss') is not None:
            self._parameters += ['-nMSS', self.n_mss]
        if getattr(self, 'mat_mss') is not None:
            self._parameters += ['-matMSS', self.mat_mss.tag]
        if getattr(self, 'lim_disp') is not None:
            self._parameters += ['-limDisp', self.lim_disp]
        if getattr(self, 'n_mns') is not None:
            self._parameters += ['-nMNS', self.n_mns]
        if getattr(self, 'mat_mns') is not None:
            self._parameters += ['-matMNS', self.mat_mns.tag]
        if getattr(self, 'lamb') is not None:
            self._parameters += ['-lambda', self.lamb]
        if getattr(self, 'no_pd_input'):
            self._parameters += ['-noPDInput']
        if getattr(self, 'no_tilt'):
            self._parameters += ['-noTilt']
        if getattr(self, 'orient') is not None:
            self._parameters += ['-orient', *self.orient]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        self.to_process(osi)

class KikuchiBearingdoBalance(ElementBase):
    """
    The KikuchiBearingdoBalance Element Class
    
    This command is used to construct a KikuchiBearing element object, which is defined by two nodes. This element
    consists of multiple shear spring model (MSS) and multiple normal spring model (MNS).

    
    """
    op_type = 'KikuchiBearing'

    def __init__(self, osi, ele_nodes, total_rubber, lim_fo, lim_fi, n_iter, shape: float=None, size: float=None, total_height: float=None, n_mss: int=None, mat_mss=None, lim_disp: float=None, n_mns: int=None, mat_mns=None, lamb: float=None, no_pd_input=False, no_tilt=False, orient=None, mass: float=None):
        """
        Initial method for KikuchiBearingdoBalance

        Parameters
        ----------
        ele_nodes: listi
            A list of two element nodes
        total_rubber: float
            Total rubber thickness
        lim_fo: float
            Tolerance of external unbalanced force (   ``limfo``), tolorance of internal unbalanced force (  
            ``limfi``), number of iterations to get rid of internal unbalanced force (   ``niter``)
        lim_fi: float
            Tolerance of external unbalanced force (   ``limfo``), tolorance of internal unbalanced force (  
            ``limfi``), number of iterations to get rid of internal unbalanced force (   ``niter``)
        n_iter: float
            Tolerance of external unbalanced force (   ``limfo``), tolorance of internal unbalanced force (  
            ``limfi``), number of iterations to get rid of internal unbalanced force (   ``niter``)
        shape: float
            Following shapes are available: round, square
        size: float
            Diameter (round shape), length of edge (square shape)
        total_height: float
            Total height of the bearing (defaulut: distance between inode and jnode)
        n_mss: int
            Number of springs in mss = nmss
        mat_mss: obj
            Mattag for mss
        lim_disp: float
            Minimum deformation to calculate equivalent coefficient of mss (see note 1)
        n_mns: int
            Number of springs in mns = nmns*nmns (for round and square shape)
        mat_mns: obj
            Mattag for mns
        lamb: float
            Parameter to calculate compression modulus distribution on mns (see note 2)
        no_pd_input: str
            Not consider p-delta moment
        no_tilt: str
            Not consider tilt of rigid link
        orient: None
            
        mass: float
            Element mass
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        if shape is None:
            self.shape = None
        else:
            self.shape = float(shape)
        if size is None:
            self.size = None
        else:
            self.size = float(size)
        self.total_rubber = float(total_rubber)
        if total_height is None:
            self.total_height = None
        else:
            self.total_height = float(total_height)
        if n_mss is None:
            self.n_mss = None
        else:
            self.n_mss = int(n_mss)
        self.mat_mss = mat_mss
        if lim_disp is None:
            self.lim_disp = None
        else:
            self.lim_disp = float(lim_disp)
        if n_mns is None:
            self.n_mns = None
        else:
            self.n_mns = int(n_mns)
        self.mat_mns = mat_mns
        if lamb is None:
            self.lamb = None
        else:
            self.lamb = float(lamb)
        self.no_pd_input = no_pd_input
        self.no_tilt = no_tilt
        self.lim_fo = float(lim_fo)
        self.lim_fi = float(lim_fi)
        self.n_iter = float(n_iter)
        self.orient = orient
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.total_rubber, '-doBalance', self.lim_fo, self.lim_fi, self.n_iter]
        if getattr(self, 'shape') is not None:
            self._parameters += ['-shape', self.shape]
        if getattr(self, 'size') is not None:
            self._parameters += ['-size', self.size]
        if getattr(self, 'total_height') is not None:
            self._parameters += ['-totalHeight', self.total_height]
        if getattr(self, 'n_mss') is not None:
            self._parameters += ['-nMSS', self.n_mss]
        if getattr(self, 'mat_mss') is not None:
            self._parameters += ['-matMSS', self.mat_mss.tag]
        if getattr(self, 'lim_disp') is not None:
            self._parameters += ['-limDisp', self.lim_disp]
        if getattr(self, 'n_mns') is not None:
            self._parameters += ['-nMNS', self.n_mns]
        if getattr(self, 'mat_mns') is not None:
            self._parameters += ['-matMNS', self.mat_mns.tag]
        if getattr(self, 'lamb') is not None:
            self._parameters += ['-lambda', self.lamb]
        if getattr(self, 'no_pd_input'):
            self._parameters += ['-noPDInput']
        if getattr(self, 'no_tilt'):
            self._parameters += ['-noTilt']
        if getattr(self, 'orient') is not None:
            self._parameters += ['-orient', *self.orient]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        self.to_process(osi)


class YamamotoBiaxialHDRcoRS(ElementBase):
    """
    The YamamotoBiaxialHDRcoRS Element Class
    
    This command is used to construct a YamamotoBiaxialHDR element object, which is defined by two nodes. This element
    can be used to represent the isotropic behavior of high-damping rubber bearing in the local y-z plane.

    
    """
    op_type = 'YamamotoBiaxialHDR'

    def __init__(self, osi, ele_nodes, tp, d_do, d_di, hr, cr, cs, orient=None, mass: float=None):
        """
        Initial method for YamamotoBiaxialHDRcoRS

        Parameters
        ----------
        ele_nodes: listi
            A list of two element nodes
        tp: int
            Compound type = 1 : x0.6r manufactured by bridgestone corporation.
        d_do: float
            Outer diameter [m]
        d_di: float
            Bore diameter [m]
        hr: float
            Total thickness of rubber layer [m] optional data
        cr: float
            Coefficients for shear stress components of τr and τs
        cs: float
            Coefficients for shear stress components of τr and τs
        orient: None
            
        mass: float
            Element mass [kg]
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.tp = int(tp)
        self.d_do = float(d_do)
        self.d_di = float(d_di)
        self.hr = float(hr)
        self.cr = float(cr)
        self.cs = float(cs)
        self.orient = orient
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.tp, self.d_do, self.d_di, self.hr, '-coRS', self.cr, self.cs]
        if getattr(self, 'orient') is not None:
            self._parameters += ['--orient', *self.orient]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        self.to_process(osi)


class ElastomericX(ElementBase):
    """
    The ElastomericX Element Class
    
    This command is used to construct an ElastomericX bearing element object in three-dimension. The 3D continuum
    geometry of an elastomeric bearing is modeled as a 2-node, 12 DOF discrete element. This elements extends the
    formulation of Elastomeric_Bearing_(Bouc-Wen)_Element element. However, instead of the user providing
    material models as input arguments, it only requires geometric and material properties of an
    elastomeric bearing as arguments. The material models in six direction are formulated
    within the element from input arguments. The time-dependent values of mechanical
    properties (e.g., shear stiffness, buckling load capacity) can also be recorded
    using the "parameters" recorder.

    For 3D problem
    """
    op_type = 'ElastomericX'

    def __init__(self, osi, ele_nodes, fy, alpha, gr, kbulk, d1, d2, ts, tr, n, x1, x2, x3, y1, y2, y3, kc, phi_m, ac, s_dratio, m, cd, tc, tag1, tag2, tag3, tag4):
        """
        Initial method for ElastomericX

        Parameters
        ----------
        ele_nodes: listi
            A list of two element nodes
        fy: float
            Yield strength
        alpha: float
            Post-yield stiffness ratio
        gr: float
            Shear modulus of elastomeric bearing
        kbulk: float
            Bulk modulus of rubber
        d1: float
            Internal diameter
        d2: float
            Outer diameter (excluding cover thickness)
        ts: float
            Single steel shim layer thickness
        tr: float
            Single  rubber layer thickness
        n: int
            Number of rubber layers
        x1: float
            Vector components in global coordinates defining local x-axis (optional)
        x2: float
            Vector components in global coordinates defining local x-axis (optional)
        x3: float
            Vector components in global coordinates defining local x-axis (optional)
        y1: float
            Vector components in global coordinates defining local y-axis (optional)
        y2: float
            Vector components in global coordinates defining local y-axis (optional)
        y3: float
            Vector components in global coordinates defining local y-axis (optional)
        kc: float
            Cavitation parameter (optional, default = 10.0)
        phi_m: float
            Damage parameter (optional, default = 0.5)
        ac: float
            Strength reduction parameter (optional, default = 1.0)
        s_dratio: float
            Shear distance from inode as a fraction of the element length (optional, default = 0.5)
        m: float
            Element mass (optional, default = 0.0)
        cd: float
            Viscous damping parameter (optional, default = 0.0)
        tc: float
            Cover thickness (optional, default = 0.0)
        tag1: float
            Tag to include cavitation and post-cavitation (optional, default = 0)
        tag2: float
            Tag to include buckling load variation (optional, default = 0)
        tag3: float
            Tag to include horizontal stiffness variation (optional, default = 0)
        tag4: float
            Tag to include vertical stiffness variation (optional, default = 0)
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.fy = float(fy)
        self.alpha = float(alpha)
        self.gr = float(gr)
        self.kbulk = float(kbulk)
        self.d1 = float(d1)
        self.d2 = float(d2)
        self.ts = float(ts)
        self.tr = float(tr)
        self.n = int(n)
        self.x1 = float(x1)
        self.x2 = float(x2)
        self.x3 = float(x3)
        self.y1 = float(y1)
        self.y2 = float(y2)
        self.y3 = float(y3)
        self.kc = float(kc)
        self.phi_m = float(phi_m)
        self.ac = float(ac)
        self.s_dratio = float(s_dratio)
        self.m = float(m)
        self.cd = float(cd)
        self.tc = float(tc)
        self.tag1 = float(tag1)
        self.tag2 = float(tag2)
        self.tag3 = float(tag3)
        self.tag4 = float(tag4)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.fy, self.alpha, self.gr, self.kbulk, self.d1, self.d2, self.ts, self.tr, self.n, self.x1, self.x2, self.x3, self.y1, self.y2, self.y3, self.kc, self.phi_m, self.ac, self.s_dratio, self.m, self.cd, self.tc, self.tag1, self.tag2, self.tag3, self.tag4]
        self.to_process(osi)


class LeadRubberX(ElementBase):
    """
    The LeadRubberX Element Class
    
    This command is used to construct a LeadRubberX bearing element object in three-dimension. The 3D continuum geometry
    of a lead rubber bearing is modeled as a 2-node, 12 DOF discrete element. It extends the formulation of ElastomericX by
    including strength degradation in lead rubber bearing due to heating of the lead-core. The LeadRubberX element
    requires only the geometric and material properties of an elastomeric bearing as arguments. The material
    models in six direction are formulated within the element from input arguments. The time-dependent
    values of mechanical properties (e.g., shear stiffness, buckling load capacity, temperature in
    the lead-core, yield strength) can also be recorded using the "parameters" recorder.

    
    """
    op_type = 'LeadRubberX'

    def __init__(self, osi, ele_nodes, fy, alpha, gr, kbulk, d1, d2, ts, tr, n, x1, x2, x3, y1, y2, y3, kc, phi_m, ac, s_dratio, m, cd, tc, q_l, c_l, k_s, a_s, tag1, tag2, tag3, tag4, tag5):
        """
        Initial method for LeadRubberX

        Parameters
        ----------
        ele_nodes: listi
            A list of two element nodes
        fy: float
            Yield strength
        alpha: float
            Post-yield stiffness ratio
        gr: float
            Shear modulus of elastomeric bearing
        kbulk: float
            Bulk modulus of rubber
        d1: float
            Internal diameter
        d2: float
            Outer diameter (excluding cover thickness)
        ts: float
            Single steel shim layer thickness
        tr: float
            Single rubber layer thickness
        n: int
            Number of rubber layers
        x1: float
            Vector components in global coordinates defining local x-axis (optional)
        x2: float
            Vector components in global coordinates defining local x-axis (optional)
        x3: float
            Vector components in global coordinates defining local x-axis (optional)
        y1: float
            Vector components in global coordinates defining local y-axis (optional)
        y2: float
            Vector components in global coordinates defining local y-axis (optional)
        y3: float
            Vector components in global coordinates defining local y-axis (optional)
        kc: float
            Cavitation parameter (optional, default = 10.0)
        phi_m: float
            Damage parameter (optional, default = 0.5)
        ac: float
            Strength reduction parameter (optional, default = 1.0)
        s_dratio: float
            Shear distance from inode as a fraction of the element length (optional, default = 0.5)
        m: float
            Element mass (optional, default = 0.0)
        cd: float
            Viscous damping parameter (optional, default = 0.0)
        tc: float
            Cover thickness (optional, default = 0.0)
        q_l: float
            Density of lead (optional, default = 11200 kg/m3)
        c_l: float
            Specific heat of lead (optional, default = 130 n-m/kg oc)
        k_s: float
            Thermal conductivity of steel (optional, default = 50 w/m oc)
        a_s: float
            Thermal diffusivity of steel (optional, default = 1.41e-05 m2/s)
        tag1: int
            Tag to include cavitation and post-cavitation (optional, default = 0)
        tag2: int
            Tag to include buckling load variation (optional, default = 0)
        tag3: int
            Tag to include horizontal stiffness variation (optional, default = 0)
        tag4: int
            Tag to include vertical stiffness variation (optional, default = 0)
        tag5: int
            Tag to include strength degradation in shear due to heating of lead core (optional, default = 0)
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.fy = float(fy)
        self.alpha = float(alpha)
        self.gr = float(gr)
        self.kbulk = float(kbulk)
        self.d1 = float(d1)
        self.d2 = float(d2)
        self.ts = float(ts)
        self.tr = float(tr)
        self.n = int(n)
        self.x1 = float(x1)
        self.x2 = float(x2)
        self.x3 = float(x3)
        self.y1 = float(y1)
        self.y2 = float(y2)
        self.y3 = float(y3)
        self.kc = float(kc)
        self.phi_m = float(phi_m)
        self.ac = float(ac)
        self.s_dratio = float(s_dratio)
        self.m = float(m)
        self.cd = float(cd)
        self.tc = float(tc)
        self.q_l = float(q_l)
        self.c_l = float(c_l)
        self.k_s = float(k_s)
        self.a_s = float(a_s)
        self.tag1 = int(tag1)
        self.tag2 = int(tag2)
        self.tag3 = int(tag3)
        self.tag4 = int(tag4)
        self.tag5 = int(tag5)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.fy, self.alpha, self.gr, self.kbulk, self.d1, self.d2, self.ts, self.tr, self.n, self.x1, self.x2, self.x3, self.y1, self.y2, self.y3, self.kc, self.phi_m, self.ac, self.s_dratio, self.m, self.cd, self.tc, self.q_l, self.c_l, self.k_s, self.a_s, self.tag1, self.tag2, self.tag3, self.tag4, self.tag5]
        self.to_process(osi)


class HDR(ElementBase):
    """
    The HDR Element Class
    
    For 3D problem
    """
    op_type = 'HDR'

    def __init__(self, osi, ele_nodes, gr, kbulk, d1, d2, ts, tr, n, a1, a2, a3, b1, b2, b3, c1, c2, c3, c4, x1, x2, x3, y1, y2, y3, kc, phi_m, ac, s_dratio, m, tc):
        """
        Initial method for HDR

        Parameters
        ----------
        ele_nodes: listi
            A list of two element nodes
        gr: float
            Shear modulus of elastomeric bearing
        kbulk: float
            Bulk modulus of rubber
        d1: float
            Internal diameter
        d2: float
            Outer diameter (excluding cover thickness)
        ts: float
            Single steel shim layer thickness
        tr: float
            Single rubber layer thickness
        n: int
            Number of rubber layers
        a1: float
            Parameters of the grant model
        a2: float
            Parameters of the grant model
        a3: float
            Parameters of the grant model
        b1: float
            Parameters of the grant model
        b2: float
            Parameters of the grant model
        b3: float
            Parameters of the grant model
        c1: float
            Parameters of the grant model
        c2: float
            Parameters of the grant model
        c3: float
            Parameters of the grant model
        c4: float
            Parameters of the grant model
        x1: float
            Vector components in global coordinates defining local x-axis (optional)
        x2: float
            Vector components in global coordinates defining local x-axis (optional)
        x3: float
            Vector components in global coordinates defining local x-axis (optional)
        y1: float
            Vector components in global coordinates defining local y-axis (optional)
        y2: float
            Vector components in global coordinates defining local y-axis (optional)
        y3: float
            Vector components in global coordinates defining local y-axis (optional)
        kc: float
            Cavitation parameter (optional, default = 10.0)
        phi_m: float
            Damage parameter (optional, default = 0.5)
        ac: float
            Strength reduction parameter (optional, default = 1.0)
        s_dratio: float
            Shear distance from inode as a fraction of the element length (optional, default = 0.5)
        m: float
            Element mass (optional, default = 0.0)
        tc: float
            Cover thickness (optional, default = 0.0)
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.gr = float(gr)
        self.kbulk = float(kbulk)
        self.d1 = float(d1)
        self.d2 = float(d2)
        self.ts = float(ts)
        self.tr = float(tr)
        self.n = int(n)
        self.a1 = float(a1)
        self.a2 = float(a2)
        self.a3 = float(a3)
        self.b1 = float(b1)
        self.b2 = float(b2)
        self.b3 = float(b3)
        self.c1 = float(c1)
        self.c2 = float(c2)
        self.c3 = float(c3)
        self.c4 = float(c4)
        self.x1 = float(x1)
        self.x2 = float(x2)
        self.x3 = float(x3)
        self.y1 = float(y1)
        self.y2 = float(y2)
        self.y3 = float(y3)
        self.kc = float(kc)
        self.phi_m = float(phi_m)
        self.ac = float(ac)
        self.s_dratio = float(s_dratio)
        self.m = float(m)
        self.tc = float(tc)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.gr, self.kbulk, self.d1, self.d2, self.ts, self.tr, self.n, self.a1, self.a2, self.a3, self.b1, self.b2, self.b3, self.c1, self.c2, self.c3, self.c4, self.x1, self.x2, self.x3, self.y1, self.y2, self.y3, self.kc, self.phi_m, self.ac, self.s_dratio, self.m, self.tc]
        self.to_process(osi)


class FPBearingPTV(ElementBase):
    """
    The FPBearingPTV Element Class
    
    The FPBearingPTV command creates a single Friction Pendulum bearing element, which is capable of accounting for the
    changes in the coefficient of friction at the sliding surface with instantaneous values of the sliding velocity, axial
    pressure and temperature at the sliding surface. The constitutive modelling is similar to the existing
    singleFPBearing element, otherwise. The FPBearingPTV element has been verified and validated in
    accordance with the ASME guidelines, details of which are presented in Chapter 4 of Kumar et al. (2015a).

    
    """
    op_type = 'FPBearingPTV'

    def __init__(self, osi, ele_nodes, mu_ref, is_pressure_dependent, p_ref, is_temperature_dependent, diffusivity, conductivity, is_velocity_dependent, rate_parameter, reffective_fp, radius__contact, k_initial, the_material_a, the_material_b, the_material_c, the_material_d, x1, x2, x3, y1, y2, y3, shear_dist, do_rayleigh, mass, iter, tol, unit):
        """
        Initial method for FPBearingPTV

        Parameters
        ----------
        ele_nodes: listi
            A list of two element nodes
        mu_ref: float
            Reference coefficient of friction
        is_pressure_dependent: int
            1 if the coefficient of friction is a function of instantaneous axial pressure
        p_ref: float
            Reference axial pressure (the bearing pressure under static loads)
        is_temperature_dependent: int
            1 if the coefficient of friction is a function of instantaneous temperature at the sliding surface
        diffusivity: float
            Thermal diffusivity of steel
        conductivity: float
            Thermal conductivity of steel
        is_velocity_dependent: int
            1 if the coefficient of friction is a function of instantaneous velocity at the sliding surface
        rate_parameter: float
            The exponent that determines the shape of the coefficient of friction vs. sliding velocity curve
        reffective_fp: float
            Effective radius of curvature of the sliding surface of the fpbearing
        radius__contact: float
            Radius of contact area at the sliding surface
        k_initial: float
            Lateral  stiffness of the sliding bearing before sliding begins
        the_material_a: int
            Tag for the uniaxial material in the axial direction
        the_material_b: int
            Tag for the uniaxial material in the torsional direction
        the_material_c: int
            Tag for the uniaxial material for rocking about local y axis
        the_material_d: int
            Tag for the uniaxial material for rocking about local z axis
        x1: float
            Vector components to define local x axis
        x2: float
            Vector components to define local x axis
        x3: float
            Vector components to define local x axis
        y1: float
            Vector components to define local y axis
        y2: float
            Vector components to define local y axis
        y3: float
            Vector components to define local y axis
        shear_dist: float
            Shear distance from inode as a fraction of the length of the element
        do_rayleigh: int
            To include rayleigh damping from the bearing
        mass: float
            Element mass
        iter: int
            Maximum number of iterations to satisfy the equilibrium of element
        tol: float
            Convergence tolerance to satisfy the equilibrium of the element
        unit: int
            Tag to identify the unit from the list below. * ``1``: n, m, s, c * ``2``: kn, m, s, c * ``3``: n, mm, s, c
            * ``4``: kn, mm, s, c * ``5``: lb, in, s, c * ``6``: kip, in, s, c * ``7``: lb, ft, s, c * ``8``: kip, ft, s, c
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.mu_ref = float(mu_ref)
        self.is_pressure_dependent = int(is_pressure_dependent)
        self.p_ref = float(p_ref)
        self.is_temperature_dependent = int(is_temperature_dependent)
        self.diffusivity = float(diffusivity)
        self.conductivity = float(conductivity)
        self.is_velocity_dependent = int(is_velocity_dependent)
        self.rate_parameter = float(rate_parameter)
        self.reffective_fp = float(reffective_fp)
        self.radius__contact = float(radius__contact)
        self.k_initial = float(k_initial)
        self.the_material_a = int(the_material_a)
        self.the_material_b = int(the_material_b)
        self.the_material_c = int(the_material_c)
        self.the_material_d = int(the_material_d)
        self.x1 = float(x1)
        self.x2 = float(x2)
        self.x3 = float(x3)
        self.y1 = float(y1)
        self.y2 = float(y2)
        self.y3 = float(y3)
        self.shear_dist = float(shear_dist)
        self.do_rayleigh = int(do_rayleigh)
        self.mass = float(mass)
        self.iter = int(iter)
        self.tol = float(tol)
        self.unit = int(unit)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mu_ref, self.is_pressure_dependent, self.p_ref, self.is_temperature_dependent, self.diffusivity, self.conductivity, self.is_velocity_dependent, self.rate_parameter, self.reffective_fp, self.radius__contact, self.k_initial, self.the_material_a, self.the_material_b, self.the_material_c, self.the_material_d, self.x1, self.x2, self.x3, self.y1, self.y2, self.y3, self.shear_dist, self.do_rayleigh, self.mass, self.iter, self.tol, self.unit]
        self.to_process(osi)
