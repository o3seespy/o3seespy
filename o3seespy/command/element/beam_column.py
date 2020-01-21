from o3seespy.command.element.base_element import ElementBase


class ElasticBeamColumn2D(ElementBase):
    """
    The ElasticBeamColumn2D Element Class
    
    This command is used to construct an elasticBeamColumn element object. The arguments for the construction of an
    elastic beam-column element depend on the dimension of the problem, (ndm)

    For a two-dimensional problem
    """
    op_type = 'elasticBeamColumn'

    def __init__(self, osi, ele_nodes, big_a, big_e, iz, transf, mass=None, c_mass=False):
        """
        Initial method for ElasticBeamColumn2D

        Parameters
        ----------
        ele_nodes: listi
            A list of two element nodes
        big_a: float
            Cross-sectional area of element
        big_e: float
            Young's modulus
        iz: float
            Second moment of area about the local z-axis
        transf: obj
            Identifier for previously-defined coordinate-transformation (crdtransf) object adfg afe asfasffg asffgrgrg
            szfrgr3gr asgrr
        mass: float
            Element mass per unit length (optional, default = 0.0)
        c_mass: str
            To form consistent mass matrix (optional, default = lumped mass matrix)
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.big_a = float(big_a)
        self.big_e = float(big_e)
        self.iz = float(iz)
        self.transf = transf
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        self.c_mass = c_mass
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.big_a, self.big_e, self.iz, self.transf.tag]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        if getattr(self, 'c_mass'):
            self._parameters += ['-cMass']
        self.to_process(osi)


class ElasticBeamColumn3D(ElementBase):
    """
    The ElasticBeamColumn3D Element Class
    
    This command is used to construct an elasticBeamColumn element object. The arguments for the construction of an
    elastic beam-column element depend on the dimension of the problem, (ndm)

    For a three-dimensional problem
    """
    op_type = 'elasticBeamColumn'

    def __init__(self, osi, ele_nodes, big_a, big_e, big_g, big_j, iy, iz, transf, mass=None, c_mass=False):
        """
        Initial method for ElasticBeamColumn3D

        Parameters
        ----------
        ele_nodes: listi
            A list of two element nodes
        big_a: float
            Cross-sectional area of element
        big_e: float
            Young's modulus
        big_g: float
            Shear modulus
        big_j: float
            Torsional moment of inertia of cross section
        iy: float
            Second moment of area about the local y-axis
        iz: float
            Second moment of area about the local z-axis
        transf: obj
            Identifier for previously-defined coordinate-transformation (crdtransf) object adfg afe asfasffg asffgrgrg
            szfrgr3gr asgrr
        mass: float
            Element mass per unit length (optional, default = 0.0)
        c_mass: str
            To form consistent mass matrix (optional, default = lumped mass matrix)
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.big_a = float(big_a)
        self.big_e = float(big_e)
        self.big_g = float(big_g)
        self.big_j = float(big_j)
        self.iy = float(iy)
        self.iz = float(iz)
        self.transf = transf
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        self.c_mass = c_mass
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.big_a, self.big_e, self.big_g, self.big_j, self.iy, self.iz, self.transf.tag]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        if getattr(self, 'c_mass'):
            self._parameters += ['-cMass']
        self.to_process(osi)


class ModElasticBeam2Dmass(ElementBase):
    """
    The ModElasticBeam2Dmass Element Class
    
    This command is used to construct a ModElasticBeam2d element object. The arguments for the construction of an
    elastic beam-column element with stiffness modifiers is applicable for 2-D problems. This element should be used
    for modelling of a structural element with an equivalent combination of one elastic element with
    stiffness-proportional damping, and two springs at its two ends with no stiffness proportional
    damping to represent a prismatic section. The modelling technique is based on a number of
    analytical studies discussed in Zareian and Medina (2010) and Zareian and Krawinkler
    (2009) and is utilized in order to solve problems related to numerical damping in
    dynamic analysis of frame structures with concentrated plasticity springs.

    
    """
    op_type = 'ModElasticBeam2d'

    def __init__(self, osi, ele_nodes, big_a, big_e, iz, k11, k33, k44, transf, mass_dens, c_mass=False):
        """
        Initial method for ModElasticBeam2Dmass

        Parameters
        ----------
        ele_nodes: listi
            A list of two element nodes
        big_a: float
            Cross-sectional area of element
        big_e: float
            Young's modulus
        iz: float
            Second moment of area about the local z-axis
        k11: float
            Stiffness modifier for translation
        k33: float
            Stiffness modifier for translation
        k44: float
            Stiffness modifier for rotation
        transf: obj
            Identifier for previously-defined coordinate-transformation (crdtransf) object
        mass_dens: float
            Element mass per unit length (optional, default = 0.0)
        c_mass: str
            To form consistent mass matrix (optional, default = lumped mass matrix)
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.big_a = float(big_a)
        self.big_e = float(big_e)
        self.iz = float(iz)
        self.k11 = float(k11)
        self.k33 = float(k33)
        self.k44 = float(k44)
        self.transf = transf
        self.mass_dens = float(mass_dens)
        self.c_mass = c_mass
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.big_a, self.big_e, self.iz, self.k11, self.k33, self.k44, self.transf.tag, '-mass', self.mass_dens]
        if getattr(self, 'c_mass'):
            self._parameters += ['-cMass']
        self.to_process(osi)


class ElasticTimoshenkoBeam2D(ElementBase):
    """
    The ElasticTimoshenkoBeam2D Element Class
    
    This command is used to construct an ElasticTimoshenkoBeam element object. A Timoshenko beam is a frame member that
    accounts for shear deformations. The arguments for the construction of an elastic Timoshenko beam element depend on
    the dimension of the problem, ndm:

    For a two-dimensional problem:
    """
    op_type = 'ElasticTimoshenkoBeam'

    def __init__(self, osi, ele_nodes, big_e, big_g, big_a, iz, avy, transf, mass=None, c_mass=False):
        """
        Initial method for ElasticTimoshenkoBeam2D

        Parameters
        ----------
        ele_nodes: listi
            A list of two element nodes
        big_e: float
            Young's modulus
        big_g: float
            Shear modulus
        big_a: float
            Cross-sectional area of element
        iz: float
            Second moment of area about the local z-axis
        avy: float
            Shear area for the local y-axis
        transf: obj
            Identifier for previously-defined coordinate-transformation (crdtransf) object
        mass: float
            Element mass per unit length (optional, default = 0.0)
        c_mass: str
            To form consistent mass matrix (optional, default = lumped mass matrix)
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.big_e = float(big_e)
        self.big_g = float(big_g)
        self.big_a = float(big_a)
        self.iz = float(iz)
        self.avy = float(avy)
        self.transf = transf
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        self.c_mass = c_mass
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.big_e, self.big_g, self.big_a, self.iz, self.avy, self.transf.tag]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        if getattr(self, 'c_mass'):
            self._parameters += ['-cMass']
        self.to_process(osi)


class ElasticTimoshenkoBeam3D(ElementBase):
    """
    The ElasticTimoshenkoBeam3D Element Class
    
    This command is used to construct an ElasticTimoshenkoBeam element object. A Timoshenko beam is a frame member that
    accounts for shear deformations. The arguments for the construction of an elastic Timoshenko beam element depend on
    the dimension of the problem, ndm:

    For a three-dimensional problem:
    """
    op_type = 'ElasticTimoshenkoBeam'

    def __init__(self, osi, ele_nodes, big_e, big_g, big_a, iz, jx, iy, iz_2, avy, avz, transf, mass=None, c_mass=False):
        """
        Initial method for ElasticTimoshenkoBeam3D

        Parameters
        ----------
        ele_nodes: listi
            A list of two element nodes
        big_e: float
            Young's modulus
        big_g: float
            Shear modulus
        big_a: float
            Cross-sectional area of element
        iz: float
            Second moment of area about the local z-axis
        jx: float
            Torsional moment of inertia of cross section
        iy: float
            Second moment of area about the local y-axis
        iz_2: None
            
        avy: float
            Shear area for the local y-axis
        avz: float
            Shear area for the local z-axis
        transf: obj
            Identifier for previously-defined coordinate-transformation (crdtransf) object
        mass: float
            Element mass per unit length (optional, default = 0.0)
        c_mass: str
            To form consistent mass matrix (optional, default = lumped mass matrix)
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.big_e = float(big_e)
        self.big_g = float(big_g)
        self.big_a = float(big_a)
        self.iz = float(iz)
        self.jx = float(jx)
        self.iy = float(iy)
        self.iz_2 = iz_2
        self.avy = float(avy)
        self.avz = float(avz)
        self.transf = transf
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        self.c_mass = c_mass
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.big_e, self.big_g, self.big_a, self.iz, self.jx, self.iy, self.iz_2, self.avy, self.avz, self.transf.tag]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        if getattr(self, 'c_mass'):
            self._parameters += ['-cMass']
        self.to_process(osi)



class DispBeamColumn(ElementBase):
    """
    The DispBeamColumn Element Class
    
    Create a ForceBeamColumn element.
    """
    op_type = 'dispBeamColumn'

    def __init__(self, osi, ele_nodes, transf, integration, c_mass=False, mass=None):
        """
        Initial method for DispBeamColumn

        Parameters
        ----------
        ele_nodes: listi
            List of two node tags
        transf: obj
            Tag of transformation
        integration: obj
            Tag of :func:`beamintegration`
        c_mass: None
            
        mass: float
            Element mass density (per unit length), from which a lumped-mass matrix is formed (optional)
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.transf = transf
        self.integration = integration
        self.c_mass = c_mass
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.transf.tag, self.integration.tag]
        if getattr(self, 'c_mass'):
            self._parameters += ['-cMass']
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        self.to_process(osi)


class ForceBeamColumn(ElementBase):
    """
    The ForceBeamColumn Element Class
    
    Create a ForceBeamColumn element.
    """
    op_type = 'forceBeamColumn'

    def __init__(self, osi, ele_nodes, transf, integration, max_iter=None, tol=None, mass=None):
        """
        Initial method for ForceBeamColumn

        Parameters
        ----------
        ele_nodes: listi
            A list of two element nodes
        transf: obj
            Tag of transformation
        integration: obj
            Tag of :func:`beamintegration`
        max_iter: int
            Maximum number of iterations to undertake to satisfy element compatibility (optional)
        tol: float
            Tolerance for satisfaction of element compatibility (optional)
        mass: float
            Element mass density (per unit length), from which a lumped-mass matrix is formed (optional)
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.transf = transf
        self.integration = integration
        if max_iter is None:
            self.max_iter = None
        else:
            self.max_iter = int(max_iter)
        if tol is None:
            self.tol = None
        else:
            self.tol = float(tol)
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.transf.tag, self.integration.tag]
        if getattr(self, 'max_iter') is not None:
            self._parameters += ['-iter', self.max_iter]
        if getattr(self, 'tol') is not None:
            if getattr(self, 'max_iter') is None:
                raise ValueError('Cannot set: tol and not: max_iter')
            self._parameters += [self.tol]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        self.to_process(osi)


class NonlinearBeamColumn(ElementBase):
    """
    The NonlinearBeamColumn Element Class
    
    Create a nonlinearBeamColumn element. This element is for backward compatability.
    """
    op_type = 'nonlinearBeamColumn'

    def __init__(self, osi, ele_nodes, num_intgr_pts, sec, transf, max_iter=None, tol=None, mass=None, int_type=None):
        """
        Initial method for NonlinearBeamColumn

        Parameters
        ----------
        ele_nodes: listi
            A list of two element nodes
        num_intgr_pts: int
            Number of integration points.
        sec: obj
            Tag of section
        transf: obj
            Tag of transformation
        max_iter: int
            Maximum number of iterations to undertake to satisfy element compatibility (optional)
        tol: float
            Tolerance for satisfaction of element compatibility (optional)
        mass: float
            Element mass density (per unit length), from which a lumped-mass matrix is formed (optional)
        int_type: str
            Integration type (optional, default is ``'lobatto'``) * ``'lobatto'`` * ``'legendre'`` * ``'radau'`` *
            ``'newtoncotes'`` * ``'trapezoidal'``
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.num_intgr_pts = int(num_intgr_pts)
        self.sec = sec
        self.transf = transf
        if max_iter is None:
            self.max_iter = None
        else:
            self.max_iter = int(max_iter)
        if tol is None:
            self.tol = None
        else:
            self.tol = float(tol)
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        self.int_type = int_type
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.num_intgr_pts, self.sec.tag, self.transf.tag]
        if getattr(self, 'max_iter') is not None:
            self._parameters += ['-iter', self.max_iter]
        if getattr(self, 'tol') is not None:
            if getattr(self, 'max_iter') is None:
                raise ValueError('Cannot set: tol and not: max_iter')
            self._parameters += [self.tol]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        if getattr(self, 'int_type') is not None:
            self._parameters += ['-integration', self.int_type]
        self.to_process(osi)


class DispBeamColumnInt(ElementBase):
    """
    The DispBeamColumnInt Element Class
    
    This command is used to construct a dispBeamColumnInt element object, which is a distributed-plasticity,
    displacement-based beam-column element which includes interaction between flexural and shear components.

    
    """
    op_type = 'dispBeamColumnInt'

    def __init__(self, osi, ele_nodes, num_intgr_pts, sec, transf, c_rot, mass_dens=None):
        """
        Initial method for DispBeamColumnInt

        Parameters
        ----------
        ele_nodes: listi
            A list of two element nodes
        num_intgr_pts: int
            Number of integration points along the element.
        sec: obj
            Identifier for previously-defined section object
        transf: obj
            Identifier for previously-defined coordinate-transformation (crdtransf) object
        c_rot: float
            Identifier for element center of rotation (or center of curvature distribution). fraction of the height
            distance from bottom to the center of rotation (0 to 1)
        mass_dens: float
            Element mass density (per unit length), from which a lumped-mass matrix is formed (optional, default=0.0)
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.num_intgr_pts = int(num_intgr_pts)
        self.sec = sec
        self.transf = transf
        self.c_rot = float(c_rot)
        if mass_dens is None:
            self.mass_dens = None
        else:
            self.mass_dens = float(mass_dens)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.num_intgr_pts, self.sec.tag, self.transf.tag, self.c_rot]
        if getattr(self, 'mass_dens') is not None:
            self._parameters += ['-mass', self.mass_dens]
        self.to_process(osi)


class MVLEM(ElementBase):
    """
    The MVLEM Element Class
    
    The MVLEM element command is used to generate a two-dimensional Multiple-Vertical-Line-Element-Model (MVLEM; Vulcano
    et al., 1988; Orakcal et al., 2004, Kolozvari et al., 2015) for simulation of flexure-dominated RC wall behavior. A
    single model element incorporates six global degrees of freedom, three of each located at the center of rigid top
    and bottom beams, as illustrated in Figure 1a. The axial/flexural response of the MVLEM is simulated by a series
    of uniaxial elements (or macro-fibers) connected to the rigid beams at the top and bottom (e.g., floor) levels,
    whereas the shear response is described by a shear spring located at height ch from the bottom of the wall
    element (Figure 1a). Shear and flexural responses of the model element are uncoupled. The relative
    rotation between top and bottom faces of the wall element occurs about the point located on the
    central axis of the element at height ch (Figure 1b). Rotations and resulting transverse
    displacements are calculated based on the wall curvature, derived from section and
    material properties, corresponding to the bending moment at height ch of each
    element (Figure 1b). A value of c=0.4 was recommended by Vulcano et al.
    (1988) based on comparison of the model response with experimental results.

    
    """
    op_type = 'MVLEM'

    def __init__(self, osi, dens, ele_nodes, m, c, thick=None, widths=None, rho=None, mat_concrete_tags=None, mat_steel_tags=None, mat_shear=None):
        """
        Initial method for MVLEM

        Parameters
        ----------
        dens: float
            Wall density
        ele_nodes: listi
            A list of two element nodes
        m: int
            Number of element macro-fibers
        c: float
            Location of center of rotation from the inode, ``c`` = 0.4 (recommended)
        thick: listf
            A list of ``m`` macro-fiber thicknesses
        widths: listf
            A list of ``m`` macro-fiber widths
        rho: listf
            A list of m reinforcing ratios corresponding to macro-fibers; for each fiber: :math:`rho_i =
            a_{s,i}/a_{gross,i} (1 < i < m)`
        mat_concrete_tags: None
            A list of ``m`` uniaxialmaterial tags for concrete
        mat_steel_tags: None
            A list of ``m`` uniaxialmaterial tags for steel
        mat_shear: obj
            Tag of uniaxialmaterial for shear material
        """
        self.dens = float(dens)
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.m = int(m)
        self.c = float(c)
        self.thick = thick
        self.widths = widths
        self.rho = rho
        if mat_concrete_tags is None:
            self.mat_concrete_tags = None
        else:
            self.mat_concrete_tags = [x.tag for x in mat_concrete_tags]
        if mat_steel_tags is None:
            self.mat_steel_tags = None
        else:
            self.mat_steel_tags = [x.tag for x in mat_steel_tags]
        self.mat_shear = mat_shear
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.dens, *self.ele_nodes, self.m, self.c]
        if getattr(self, 'thick') is not None:
            self._parameters += ['-thick', *self.thick]
        if getattr(self, 'widths') is not None:
            self._parameters += ['-width', *self.widths]
        if getattr(self, 'rho') is not None:
            self._parameters += ['-rho', *self.rho]
        if getattr(self, 'mat_concrete_tags') is not None:
            self._parameters += ['-matConcrete', *self.mat_concrete_tags]
        if getattr(self, 'mat_steel_tags') is not None:
            self._parameters += ['-matSteel', *self.mat_steel_tags]
        if getattr(self, 'mat_shear') is not None:
            self._parameters += ['-matShear', self.mat_shear.tag]
        self.to_process(osi)


class SFIMVLEM(ElementBase):
    """
    The SFIMVLEM Element Class
    
    The SFI_MVLEM command is used to construct a Shear-Flexure Interaction Multiple-Vertical-Line-Element Model
    (SFI-MVLEM, Kolozvari et al., 2015a, b, c), which captures interaction between axial/flexural and shear
    behavior of RC structural walls and columns under cyclic loading. The SFI_MVLEM element (Figure 1)
    incorporates 2-D RC panel behavior described by the Fixed-Strut-Angle-Model (nDMaterial FSAM;
    Ulugtekin, 2010; Orakcal et al., 2012), into a 2-D macroscopic fiber-based model (MVLEM).
    The interaction between axial and shear behavior is captured at each RC panel
    (macro-fiber) level, which further incorporates interaction between shear
    and flexural behavior at the SFI_MVLEM element level.

    
    """
    op_type = 'SFI_MVLEM'

    def __init__(self, osi, ele_nodes, m, c, thick=None, widths=None, mat_tags=None):
        """
        Initial method for SFIMVLEM

        Parameters
        ----------
        ele_nodes: listi
            A list of two element nodes
        m: int
            Number of element macro-fibers
        c: float
            Location of center of rotation with from the inode, ``c`` = 0.4 (recommended)
        thick: None
            
        widths: None
            
        mat_tags: None
            
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.m = int(m)
        self.c = float(c)
        self.thick = thick
        self.widths = widths
        self.mat_tags = mat_tags
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.m, self.c]
        if getattr(self, 'thick') is not None:
            self._parameters += ['-thick', *self.thick]
        if getattr(self, 'widths') is not None:
            self._parameters += ['-width', *self.widths]
        if getattr(self, 'mat_tags') is not None:
            self._parameters += ['-mat', *self.mat_tags]
        self.to_process(osi)
