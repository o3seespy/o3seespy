from o3seespy.command.element.base_element import ElementBase


class ZeroLength(ElementBase):
    """
    The ZeroLength Element Class
    
    This command is used to construct a zeroLength element object, which is defined by two nodes at the same location.
    The nodes are connected by multiple UniaxialMaterial objects to represent the force-deformation relationship for the
    element.
    """
    op_type = 'zeroLength'

    def __init__(self, osi, ele_nodes, mats: list=None, dirs: list=None, r_flag: float=None, orient: list=None):
        """
        Initial method for ZeroLength

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of two element nodes
        mats: list, optional
            A list of objects associated with previously-defined uniaxial_materials
        dirs: list, optional
            A list of material directions: * 1,2,3 - translation along local x,y,z axes, respectively; * 4,5,6 -
            rotation about local x,y,z axes, respectively
        r_flag: float, optional
            Optional, default = 0 * ``rflag`` = 0 no rayleigh damping (default)
        orient: list, optional
            

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [0, 0]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(2)]
        >>> bilinear_mat = o3.uniaxial_material.Steel01(osi, fy=300., e0=200.0e3, b=0.01)
        >>> o3.element.ZeroLength(osi, ele_nodes, mats=[bilinear_mat], dirs=[o3.cc.DOF2D_X], r_flag=1)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        if mats is None:
            self.mats = None
        else:
            self.mats = [x.tag for x in mats]
        self.dirs = dirs
        if r_flag is None:
            self.r_flag = None
        else:
            self.r_flag = float(r_flag)
        self.orient = orient
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes]
        if getattr(self, 'mats') is not None:
            self._parameters += ['-mat', *self.mats]
        if getattr(self, 'dirs') is not None:
            self._parameters += ['-dir', *self.dirs]
        if getattr(self, 'r_flag') is not None:
            self._parameters += ['-doRayleigh', self.r_flag]
        if getattr(self, 'orient') is not None:
            self._parameters += ['--orient', *self.orient]
        self.to_process(osi)


class ZeroLengthND(ElementBase):
    """
    The ZeroLengthND Element Class
    
    This command is used to construct a zeroLengthND element object, which is defined by two nodes at the same location.
    The nodes are connected by a single NDMaterial object to represent the force-deformation relationship for the element.
    """
    op_type = 'zeroLengthND'

    def __init__(self, osi, ele_nodes, mat, uni, orient: list=None):
        """
        Initial method for ZeroLengthND

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of two element nodes
        mat: obj
            Object associated with previously-defined ndmaterial object
        uni: obj
            Object associated with previously-defined uniaxial_material object which may be used to represent uncoupled
            behavior orthogonal to the plane of the ndmaterial response. see notes 2 and 3.
        orient: list, optional
            

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [0, 0]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(2)]
        >>> mat = o3.nd_material.ElasticIsotropic(osi, 1.0, 0.3)
        >>> uni = o3.uniaxial_material.Elastic(osi, 1.0)
        >>> o3.element.ZeroLengthND(osi, ele_nodes=ele_nodes, mat=mat, uni=uni, orient=[1, 2, 3, 4, 5, 6])
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.mat = mat
        self.uni = uni
        self.orient = orient
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mat.tag, self.uni.tag]
        if getattr(self, 'orient') is not None:
            self._parameters += ['--orient', *self.orient]
        self.to_process(osi)


class ZeroLengthSection(ElementBase):
    """
    The ZeroLengthSection Element Class
    
    This command is used to construct a zero length element object, which is defined by two nodes at the same location.
    The nodes are connected by a single section object to represent the force-deformation relationship for the element.
    """
    op_type = 'zeroLengthSection'

    def __init__(self, osi, ele_nodes, sec, r_flag: float=None, orient: list=None):
        """
        Initial method for ZeroLengthSection

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of two element nodes
        sec: obj
            Object associated with previously-defined section object
        r_flag: float, optional
            Optional, default = 0 * ``rflag`` = 0 no rayleigh damping (default) * ``rflag`` = 1 include rayleigh damping
        orient: list, optional
            

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [0, 0]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(2)]
        >>> sec = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
        >>> o3.element.ZeroLengthSection(osi, ele_nodes=ele_nodes, sec=sec, r_flag=1.0, orient=[1])
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.sec = sec
        if r_flag is None:
            self.r_flag = None
        else:
            self.r_flag = float(r_flag)
        self.orient = orient
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.sec.tag]
        if getattr(self, 'r_flag') is not None:
            self._parameters += ['-doRayleigh', self.r_flag]
        if getattr(self, 'orient') is not None:
            self._parameters += ['--orient', *self.orient]
        self.to_process(osi)


class CoupledZeroLength(ElementBase):
    """
    The CoupledZeroLength Element Class
    
    
    """
    op_type = 'CoupledZeroLength'

    def __init__(self, osi, ele_nodes, dirn1, dirn2, mat, r_flag=1):
        """
        Initial method for CoupledZeroLength

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of two element nodes
        dirn1: int
            The two directions, 1 through ndof.
        dirn2: int
            The two directions, 1 through ndof.
        mat: obj
            Objects associated with previously-defined uniaxial_material
        r_flag: float, optional
            Optional, default = 0 * ``rflag`` = 0 no rayleigh damping (default) * ``rflag`` = 1 include rayleigh damping

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [0, 0]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(2)]
        >>> mat = o3.uniaxial_material.Elastic(osi, 1.0)
        >>> o3.element.CoupledZeroLength(osi, ele_nodes=ele_nodes, dirn1=1, dirn2=1, mat=mat, r_flag=1)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.dirn1 = int(dirn1)
        self.dirn2 = int(dirn2)
        self.mat = mat
        self.r_flag = float(r_flag)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.dirn1, self.dirn2, self.mat.tag, self.r_flag]
        self.to_process(osi)


class ZeroLengthContact2Dnormal(ElementBase):
    """
    The ZeroLengthContact2Dnormal Element Class
    
    This command is used to construct a zeroLengthContact2D element, which is Node-to-node frictional contact element
    used in two dimensional analysis and three dimensional analysis:
    """
    op_type = 'zeroLengthContact2D'

    def __init__(self, osi, ele_nodes, kn, kt, mu, nx, ny):
        """
        Initial method for ZeroLengthContact2Dnormal

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of a constrained and a retained nodes
        kn: float
            Penalty in normal direction
        kt: float
            Penalty in tangential direction
        mu: float
            Friction coefficient
        nx: None
            
        ny: None
            

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [0, 0]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(2)]
        >>> o3.element.ZeroLengthContact2Dnormal(osi, ele_nodes=ele_nodes, kn=1.0, kt=1.0, mu=1.0, nx=1, ny=0)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.kn = float(kn)
        self.kt = float(kt)
        self.mu = float(mu)
        self.nx = nx
        self.ny = ny
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.kn, self.kt, self.mu, '-normal', self.nx, self.ny]
        self.to_process(osi)

class ZeroLengthContact3D(ElementBase):
    """
    The ZeroLengthContact3D Element Class
    
    This command is used to construct a zeroLengthContact3D element, which is Node-to-node frictional contact element
    used in two dimensional analysis and three dimensional analysis:
    """
    op_type = 'zeroLengthContact3D'

    def __init__(self, osi, ele_nodes, kn, kt, mu, c, dir):
        """
        Initial method for ZeroLengthContact3D

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of a constrained and a retained nodes
        kn: float
            Penalty in normal direction
        kt: float
            Penalty in tangential direction
        mu: float
            Friction coefficient
        c: float
            Cohesion (not available in 2d)
        dir: int
            Direction flag of the contact plane (3d), it can be: * 1 out normal of the master plane pointing to +x
            direction * 2 out normal of the master plane pointing to +y direction * 3 out normal of the master plane pointing
            to +z direction

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(2)]
        >>> o3.element.ZeroLengthContact3D(osi, ele_nodes=ele_nodes, kn=1.0, kt=1.0, mu=1.0, c=1.0, dir=1)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.kn = float(kn)
        self.kt = float(kt)
        self.mu = float(mu)
        self.c = float(c)
        self.dir = int(dir)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.kn, self.kt, self.mu, self.c, self.dir]
        self.to_process(osi)


class ZeroLengthContactNTS2D(ElementBase):
    """
    The ZeroLengthContactNTS2D Element Class
    
    
    """
    op_type = 'zeroLengthContactNTS2D'

    def __init__(self, osi, kn, kt, phi, s_nd_num: int=None, m_nd_num: int=None, nodes: list=None):
        """
        Initial method for ZeroLengthContactNTS2D

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        kn: float
            Penalty in normal direction
        kt: float
            Penalty in tangential direction
        phi: float
            Friction angle in degrees
        s_nd_num: int, optional
            Number of slave nodes
        m_nd_num: int, optional
            Number of master nodes
        nodes: list, optional
            Slave and master node objects respectively

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [0, 0]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(2)]
        >>> o3.element.ZeroLengthContactNTS2D(osi, s_nd_num=1, m_nd_num=1, nodes=ele_nodes, kn=1.0, kt=1.0, phi=1.0)
        """
        self.osi = osi
        if s_nd_num is None:
            self.s_nd_num = None
        else:
            self.s_nd_num = int(s_nd_num)
        if m_nd_num is None:
            self.m_nd_num = None
        else:
            self.m_nd_num = int(m_nd_num)
        if nodes is None:
            self.nodes = None
        else:
            self.nodes = [x.tag for x in nodes]
        self.kn = float(kn)
        self.kt = float(kt)
        self.phi = float(phi)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.kn, self.kt, self.phi]
        if getattr(self, 's_nd_num') is not None:
            self._parameters += ['-sNdNum', self.s_nd_num]
        if getattr(self, 'm_nd_num') is not None:
            self._parameters += ['-mNdNum', self.m_nd_num]
        if getattr(self, 'nodes') is not None:
            self._parameters += ['-Nodes', *self.nodes]
        self.to_process(osi)


class ZeroLengthInterface2Ddof(ElementBase):
    """
    The ZeroLengthInterface2Ddof Element Class
    
    
    """
    op_type = 'zeroLengthInterface2D'

    def __init__(self, osi, sdof, mdof, kn, kt, phi, s_nd_num: int=None, m_nd_num: int=None, nodes: list=None):
        """
        Initial method for ZeroLengthInterface2Ddof

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        sdof: int
            Slave and master degree of freedom
        mdof: int
            Slave and master degree of freedom
        kn: float
            Penalty in normal direction
        kt: float
            Penalty in tangential direction
        phi: float
            Friction angle in degrees
        s_nd_num: int, optional
            Number of slave nodes
        m_nd_num: int, optional
            Number of master nodes
        nodes: list, optional
            Slave and master node objects respectively

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> nodes = [1, 1]
        >>> o3.element.ZeroLengthInterface2Ddof(osi, s_nd_num=1, m_nd_num=1, sdof=1, mdof=1, nodes=nodes, kn=1.0, kt=1.0, phi=1.0)
        """
        self.osi = osi
        if s_nd_num is None:
            self.s_nd_num = None
        else:
            self.s_nd_num = int(s_nd_num)
        if m_nd_num is None:
            self.m_nd_num = None
        else:
            self.m_nd_num = int(m_nd_num)
        self.sdof = int(sdof)
        self.mdof = int(mdof)
        if nodes is None:
            self.nodes = None
        else:
            self.nodes = [x.tag for x in nodes]
        self.kn = float(kn)
        self.kt = float(kt)
        self.phi = float(phi)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, '-dof', self.sdof, self.mdof, self.kn, self.kt, self.phi]
        if getattr(self, 's_nd_num') is not None:
            self._parameters += ['-sNdNum', self.s_nd_num]
        if getattr(self, 'm_nd_num') is not None:
            self._parameters += ['-mNdNum', self.m_nd_num]
        if getattr(self, 'nodes') is not None:
            self._parameters += ['-Nodes', *self.nodes]
        self.to_process(osi)


class ZeroLengthImpact3D(ElementBase):
    """
    The ZeroLengthImpact3D Element Class
    
    This command constructs a node-to-node zero-length contact element in 3D space to simulate the impact/pounding and
    friction phenomena.
    """
    op_type = 'zeroLengthImpact3D'

    def __init__(self, osi, ele_nodes, direction, init_gap, friction_ratio, kt, kn, kn2, delta_y, cohesion):
        """
        Initial method for ZeroLengthImpact3D

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of a constrained and a retained nodes * ``1`` if out-normal vector of master plane points to +x
            direction * ``2`` if out-normal vector of master plane points to +y direction * ``3`` if out-normal vector of
            master plane points to +z direction
        direction: None
            
        init_gap: float
            Initial gap between master plane and slave plane
        friction_ratio: float
            Friction ratio in two tangential directions (parallel to master and slave planes)
        kt: float
            Penalty in two tangential directions
        kn: float
            Penalty in normal direction (normal to master and slave planes)
        kn2: float
            Penalty in normal direction after yielding based on hertz impact model
        delta_y: float
            Yield deformation based on hertz impact model
        cohesion: float
            Cohesion, if no cohesion, it is zero

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
        >>> o3.element.ZeroLengthImpact3D(osi, ele_nodes=ele_nodes, direction=1, init_gap=1.0, friction_ratio=1.0, kt=1.0, kn=1.0, kn2=1.0, delta_y=1.0, cohesion=1.0)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.direction = direction
        self.init_gap = float(init_gap)
        self.friction_ratio = float(friction_ratio)
        self.kt = float(kt)
        self.kn = float(kn)
        self.kn2 = float(kn2)
        self.delta_y = float(delta_y)
        self.cohesion = float(cohesion)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.direction, self.init_gap, self.friction_ratio, self.kt, self.kn, self.kn2, self.delta_y, self.cohesion]
        self.to_process(osi)
