from o3seespy.command.element.base_element import ElementBase


class ElasticBeamColumn2D(ElementBase):
    """
    The ElasticBeamColumn2D Element Class
    
    This command is used to construct an elasticBeamColumn element object. The arguments for the construction of an
    elastic beam-column element depend on the dimension of the problem, (ndm)

    For a two-dimensional problem
    """
    op_type = 'elasticBeamColumn'

    def __init__(self, osi, ele_nodes, area, e_mod, iz, transf, mass: float=None, c_mass=False, release_code=None):
        """
        Initial method for ElasticBeamColumn2D

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of two element nodes
        area: float
            Cross-sectional area of element
        e_mod: float
            Young's modulus
        iz: float
            Second moment of area about the local z-axis
        transf: obj
            Identifier for previously-defined coordinate-transformation (crdtransf) object
        mass: float, optional
            Element mass per unit length (optional, default = 0.0)
        c_mass: bool
            To form consistent mass matrix (optional, default = lumped mass matrix)
        release_code: None, optional
            

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
        >>> transf = o3.geom_transf.Linear2D(osi, [])
        >>> o3.element.ElasticBeamColumn2D(osi, ele_nodes=ele_nodes, area=1.0, e_mod=1.0, iz=1.0, transf=transf, mass=1.0)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.area = float(area)
        self.e_mod = float(e_mod)
        self.iz = float(iz)
        self.transf = transf
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        self.c_mass = c_mass
        self.release_code = release_code
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.area, self.e_mod, self.iz, self.transf.tag]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        if getattr(self, 'c_mass'):
            self._parameters += ['-cMass']
        if getattr(self, 'release_code') is not None:
            self._parameters += ['-release', self.release_code]
        self.to_process(osi)


class ElasticBeamColumn3D(ElementBase):
    """
    The ElasticBeamColumn3D Element Class
    
    This command is used to construct an elasticBeamColumn element object. The arguments for the construction of an
    elastic beam-column element depend on the dimension of the problem, (ndm)

    For a three-dimensional problem
    """
    op_type = 'elasticBeamColumn'

    def __init__(self, osi, ele_nodes, area, e_mod, g_mod, jxx, iy, iz, transf, mass: float=None, c_mass=False):
        """
        Initial method for ElasticBeamColumn3D

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of two element nodes
        area: float
            Cross-sectional area of element
        e_mod: float
            Young's modulus
        g_mod: float
            Shear modulus
        jxx: float
            Torsional moment of inertia of cross section
        iy: float
            Second moment of area about the local y-axis
        iz: float
            Second moment of area about the local z-axis
        transf: obj
            Identifier for previously-defined coordinate-transformation (crdtransf) object
        mass: float, optional
            Element mass per unit length (optional, default = 0.0)
        c_mass: bool
            To form consistent mass matrix (optional, default = lumped mass matrix)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
        >>> transf = o3.geom_transf.Linear2D(osi, [])
        >>> o3.element.ElasticBeamColumn3D(osi, ele_nodes=ele_nodes, area=1.0, e_mod=1.0, g_mod=1.0, jxx=1.0, iy=1.0, iz=1.0, transf=transf, mass=1.0)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.area = float(area)
        self.e_mod = float(e_mod)
        self.g_mod = float(g_mod)
        self.jxx = float(jxx)
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
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.area, self.e_mod, self.g_mod, self.jxx, self.iy, self.iz, self.transf.tag]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        if getattr(self, 'c_mass'):
            self._parameters += ['-cMass']
        self.to_process(osi)


class ModElasticBeam2D(ElementBase):
    """
    The ModElasticBeam2D Element Class
    
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

    def __init__(self, osi, ele_nodes, area, e_mod, iz, k11, k33, k44, transf, c_mass=False, mass: float=None):
        """
        Initial method for ModElasticBeam2D

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of two element nodes
        area: float
            Cross-sectional area of element
        e_mod: float
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
        c_mass: bool
            To form consistent mass matrix (optional, default = lumped mass matrix)
        mass: float, optional
            Element mass per unit length (optional, default = 0.0)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.area = float(area)
        self.e_mod = float(e_mod)
        self.iz = float(iz)
        self.k11 = float(k11)
        self.k33 = float(k33)
        self.k44 = float(k44)
        self.transf = transf
        self.c_mass = c_mass
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.area, self.e_mod, self.iz, self.k11, self.k33, self.k44, self.transf.tag]
        if getattr(self, 'c_mass'):
            self._parameters += ['-cMass']
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
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

    def __init__(self, osi, ele_nodes, e_mod, g_mod, area, iz, avy, transf, c_mass=False, mass: float=None):
        """
        Initial method for ElasticTimoshenkoBeam2D

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of two element nodes
        e_mod: float
            Young's modulus
        g_mod: float
            Shear modulus
        area: float
            Cross-sectional area of element
        iz: float
            Second moment of area about the local z-axis
        avy: float
            Shear area for the local y-axis
        transf: obj
            Identifier for previously-defined coordinate-transformation (crdtransf) object
        c_mass: bool
            To form consistent mass matrix (optional, default = lumped mass matrix)
        mass: float, optional
            Element mass per unit length (optional, default = 0.0)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
        >>> transf = o3.geom_transf.Linear2D(osi, [])
        >>> o3.element.ElasticTimoshenkoBeam2D(osi, ele_nodes=ele_nodes, e_mod=1.0, g_mod=1.0, area=1.0, iz=1.0, avy=1.0, transf=transf, mass=1.0)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.e_mod = float(e_mod)
        self.g_mod = float(g_mod)
        self.area = float(area)
        self.iz = float(iz)
        self.avy = float(avy)
        self.transf = transf
        self.c_mass = c_mass
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.e_mod, self.g_mod, self.area, self.iz, self.avy, self.transf.tag]
        if getattr(self, 'c_mass'):
            self._parameters += ['-cMass']
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
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

    def __init__(self, osi, ele_nodes, e_mod, g_mod, area, iz, jxx, iy, iz_2, avy, avz, transf, c_mass=False, mass: float=None):
        """
        Initial method for ElasticTimoshenkoBeam3D

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of two element nodes
        e_mod: float
            Young's modulus
        g_mod: float
            Shear modulus
        area: float
            Cross-sectional area of element
        iz: float
            Second moment of area about the local z-axis
        jxx: float
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
        c_mass: bool
            To form consistent mass matrix (optional, default = lumped mass matrix)
        mass: float, optional
            Element mass per unit length (optional, default = 0.0)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
        >>> transf = o3.geom_transf.Linear2D(osi, [])
        >>> o3.element.ElasticTimoshenkoBeam3D(osi, ele_nodes=ele_nodes, e_mod=1.0, g_mod=1.0, area=1.0, iz=1.0, jxx=1.0, iy=1.0, iz_2=1, avy=1.0, avz=1.0, transf=transf, mass=1.0, c_mass=True)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.e_mod = float(e_mod)
        self.g_mod = float(g_mod)
        self.area = float(area)
        self.iz = float(iz)
        self.jxx = float(jxx)
        self.iy = float(iy)
        self.iz_2 = iz_2
        self.avy = float(avy)
        self.avz = float(avz)
        self.transf = transf
        self.c_mass = c_mass
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.e_mod, self.g_mod, self.area, self.iz, self.jxx, self.iy, self.iz_2, self.avy, self.avz, self.transf.tag]
        if getattr(self, 'c_mass'):
            self._parameters += ['-cMass']
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        self.to_process(osi)



class DispBeamColumn(ElementBase):
    """
    The DispBeamColumn Element Class
    
    Create a dispBeamColumn element.
    """
    op_type = 'dispBeamColumn'

    def __init__(self, osi, ele_nodes, transf, integration, c_mass=False, mass: float=None):
        """
        Initial method for DispBeamColumn

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            List of two node objects
        transf: obj
            Object of transformation
        integration: obj
            Object of :func:`beamintegration`
        c_mass: None
            
        mass: float, optional
            Element mass density (per unit length), from which a lumped-mass matrix is formed 

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> i_node = o3.node.Node(osi, 0.0, 0.0)
        >>> j_node = o3.node.Node(osi, 0.0, 1.0)
        >>> transf = o3.geom_transf.Linear2D(osi, [])
        >>> sec = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
        >>> integration = o3.beam_integration.Lobatto(osi, sec, 5)
        >>> o3.element.DispBeamColumn(osi, ele_nodes=[i_node, j_node], transf=transf, integration=integration, c_mass=1, mass=0.0)
        """
        self.osi = osi
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

    def __init__(self, osi, ele_nodes, transf, integration, max_iter: int=None, tol: float=None, mass: float=None):
        """
        Initial method for ForceBeamColumn

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of two element nodes
        transf: obj
            Object of transformation
        integration: obj
            Object of :func:`beamintegration`
        max_iter: int, optional
            Maximum number of iterations to undertake to satisfy element compatibility 
        tol: float, optional
            Tolerance for satisfaction of element compatibility 
        mass: float, optional
            Element mass density (per unit length), from which a lumped-mass matrix is formed 

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> i_node = o3.node.Node(osi, 0.0, 0.0)
        >>> j_node = o3.node.Node(osi, 0.0, 1.0)
        >>> transf = o3.geom_transf.Linear2D(osi, [])
        >>> sec = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
        >>> integration = o3.beam_integration.Lobatto(osi, sec, 5)
        >>> o3.element.ForceBeamColumn(osi, ele_nodes=[i_node, j_node], transf=transf, integration=integration, max_iter=10, tol=1e-12, mass=0.0)
        """
        self.osi = osi
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

    def __init__(self, osi, ele_nodes, num_intgr_pts, sec, transf, max_iter: int=None, tol: float=None, mass: float=None, int_type: str=None):
        """
        Initial method for NonlinearBeamColumn

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of two element nodes
        num_intgr_pts: int
            Number of integration points.
        sec: obj
            Object of section
        transf: obj
            Object of transformation
        max_iter: int, optional
            Maximum number of iterations to undertake to satisfy element compatibility 
        tol: float, optional
            Tolerance for satisfaction of element compatibility 
        mass: float, optional
            Element mass density (per unit length), from which a lumped-mass matrix is formed 
        int_type: str, optional
            Integration type (optional, default is ``'lobatto'``) * ``'lobatto'`` * ``'legendre'`` * ``'radau'`` *
            ``'newtoncotes'`` * ``'trapezoidal'``

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> i_node = o3.node.Node(osi, 0.0, 0.0)
        >>> j_node = o3.node.Node(osi, 0.0, 1.0)
        >>> sec = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
        >>> transf = o3.geom_transf.Linear2D(osi, [])
        >>> o3.element.NonlinearBeamColumn(osi, ele_nodes=[i_node, j_node], num_intgr_pts=1, sec=sec, transf=transf,
        >>>                                max_iter=10, tol=1e-12, mass=0.0, int_type="radau")
        """
        self.osi = osi
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

    def __init__(self, osi, ele_nodes, num_intgr_pts, sec, transf, c_rot, mass: float=None):
        """
        Initial method for DispBeamColumnInt

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
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
        mass: float, optional
            Element mass density (per unit length), from which a lumped-mass matrix is formed (optional, default=0.0)

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(2)]
        >>> sec = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
        >>> transf = o3.geom_transf.Linear2D(osi, [])
        >>> o3.element.DispBeamColumnInt(osi, ele_nodes=ele_nodes, num_intgr_pts=4, sec=sec, transf=transf, c_rot=1.0, mass=1.0)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.num_intgr_pts = int(num_intgr_pts)
        self.sec = sec
        self.transf = transf
        self.c_rot = float(c_rot)
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.num_intgr_pts, self.sec.tag, self.transf.tag, self.c_rot]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
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

    def __init__(self, osi, dens, ele_nodes, m, c, thick: list=None, widths: list=None, rho: list=None, mat_concretes: list=None, mat_steels: list=None, mat_shear=None):
        r"""
        Initial method for MVLEM

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        dens: float
            Wall density
        ele_nodes: list
            A list of two element nodes
        m: int
            Number of element macro-fibers
        c: float
            Location of center of rotation from the inode, ``c`` = 0.4 (recommended)
        thick: list, optional
            A list of ``m`` macro-fiber thicknesses
        widths: list, optional
            A list of ``m`` macro-fiber widths
        rho: list, optional
            A list of m reinforcing ratios corresponding to macro-fibers; for each fiber: :math:`rho_i =
            a_{s,i}/a_{gross,i} (1 < i < m)`
        mat_concretes: list, optional
            A list of ``m`` uniaxial_material objects for concrete
        mat_steels: list, optional
            A list of ``m`` uniaxial_material objects for steel
        mat_shear: obj, optional
            Object of uniaxial_material for shear material

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
        >>> mat_conc = [o3.uniaxial_material.Concrete01(osi, fpc=1.0, epsc0=1.0, fpcu=1.0, eps_u=1.0),
        >>>             o3.uniaxial_material.Concrete01(osi, fpc=1.0, epsc0=1.0, fpcu=1.0, eps_u=1.0)]
        >>> mat_steel = [o3.uniaxial_material.Steel02(osi, fy=1.0, e0=1.0, b=1.0, params=[15, 0.925, 0.15])]
        >>> mat_shear = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
        >>> o3.element.MVLEM(osi, dens=1.0, ele_nodes=ele_nodes, m=1, c=1.0, thick=[1.0, 1.0], widths=[1, 1], rho=[1., 1.],
        >>>                  mat_concretes=mat_conc, mat_steels=mat_steel, mat_shear=mat_shear)
        """
        self.osi = osi
        self.dens = float(dens)
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.m = int(m)
        self.c = float(c)
        self.thick = thick
        self.widths = widths
        self.rho = rho
        if mat_concretes is None:
            self.mat_concretes = None
        else:
            self.mat_concretes = [x.tag for x in mat_concretes]
        if mat_steels is None:
            self.mat_steels = None
        else:
            self.mat_steels = [x.tag for x in mat_steels]
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
        if getattr(self, 'mat_concretes') is not None:
            self._parameters += ['-matConcrete', *self.mat_concretes]
        if getattr(self, 'mat_steels') is not None:
            self._parameters += ['-matSteel', *self.mat_steels]
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

    def __init__(self, osi, ele_nodes, m, c, thick=None, widths=None, mats=None):
        """
        Initial method for SFIMVLEM

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of two element nodes
        m: int
            Number of element macro-fibers
        c: float
            Location of center of rotation with from the inode, ``c`` = 0.4 (recommended)
        thick: None, optional
            
        widths: None, optional
            
        mats: None, optional
            

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
        >>> mats = [o3.uniaxial_material.Elastic(osi, 1.0, 1.0), o3.uniaxial_material.Elastic(osi, 1.0, 1.0)]
        >>> mat_tags = [x.tag for x in mats]  # TODO: should pass in mats not mat tags
        >>> o3.element.SFIMVLEM(osi, ele_nodes=ele_nodes, m=1, c=1.0, thick=[1.0, 1.0], widths=[1.0, 1.0], mat_tags=mat_tags)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.m = int(m)
        self.c = float(c)
        self.thick = thick
        self.widths = widths
        self.mats = mats
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.m, self.c]
        if getattr(self, 'thick') is not None:
            self._parameters += ['-thick', *self.thick]
        if getattr(self, 'widths') is not None:
            self._parameters += ['-width', *self.widths]
        if getattr(self, 'mats') is not None:
            self._parameters += ['-mat', *self.mats]
        self.to_process(osi)
