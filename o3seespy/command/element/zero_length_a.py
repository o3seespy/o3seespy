from o3seespy.command.element.base_element import ElementBase


class ZeroLength(ElementBase):
    """
    The ZeroLength Element Class
    
    This command is used to construct a zeroLength element object, which is defined by two nodes at the same location.
    The nodes are connected by multiple UniaxialMaterial objects to represent the force-deformation relationship for the
    element.
    """
    op_type = 'zeroLength'

    def __init__(self, osi, ele_nodes, mats=None, dir_args=None, r_flag: float=None, orient=None):
        """
        Initial method for ZeroLength

        Parameters
        ----------
        ele_nodes: listi
            A list of two element nodes
        mats: None
            A list of tags associated with previously-defined uniaxialmaterials
        dir_args: listi
            A list of material directions: * 1,2,3 - translation along local x,y,z axes, respectively; * 4,5,6 -
            rotation about local x,y,z axes, respectively
        r_flag: float
            Optional, default = 0 * ``rflag`` = 0 no rayleigh damping (default)
        orient: None
            
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        if mats is None:
            self.mats = None
        else:
            self.mats = [x.tag for x in mats]
        self.dir_args = dir_args
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
        if getattr(self, 'dir_args') is not None:
            self._parameters += ['-dir', *self.dir_args]
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

    def __init__(self, osi, ele_nodes, mat, uni, orient=None):
        """
        Initial method for ZeroLengthND

        Parameters
        ----------
        ele_nodes: listi
            A list of two element nodes
        mat: obj
            Tag associated with previously-defined ndmaterial object
        uni: obj
            Tag associated with previously-defined uniaxialmaterial object which may be used to represent uncoupled
            behavior orthogonal to the plane of the ndmaterial response. see notes 2 and 3.
        orient: None
            
        """
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

    def __init__(self, osi, ele_nodes, sec, r_flag: float=None, orient=None):
        """
        Initial method for ZeroLengthSection

        Parameters
        ----------
        ele_nodes: listi
            A list of two element nodes
        sec: obj
            Tag associated with previously-defined section object
        r_flag: float
            Optional, default = 0 * ``rflag`` = 0 no rayleigh damping (default) * ``rflag`` = 1 include rayleigh damping
        orient: None
            
        """
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
        ele_nodes: listi
            A list of two element nodes
        dirn1: int
            The two directions, 1 through ndof.
        dirn2: int
            The two directions, 1 through ndof.
        mat: obj
            Tags associated with previously-defined uniaxialmaterial
        r_flag: float
            Optional, default = 0 * ``rflag`` = 0 no rayleigh damping (default) * ``rflag`` = 1 include rayleigh damping
        """
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
        ele_nodes: listi
            A list of a slave and a master nodes
        kn: float
            Penalty in normal direction
        kt: float
            Penalty in tangential direction
        mu: float
            Friction coefficient
        nx: None
            
        ny: None
            
        """
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
        ele_nodes: listi
            A list of a slave and a master nodes
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
        """
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

    def __init__(self, osi, kn, kt, phi, s_nd_num: int=None, m_nd_num: int=None, nodes=None):
        """
        Initial method for ZeroLengthContactNTS2D

        Parameters
        ----------
        kn: float
            Penalty in normal direction
        kt: float
            Penalty in tangential direction
        phi: float
            Friction angle in degrees
        s_nd_num: int
            Number of slave nodes
        m_nd_num: int
            Number of master nodes
        nodes: listi
            Slave and master node tags respectively
        """
        if s_nd_num is None:
            self.s_nd_num = None
        else:
            self.s_nd_num = int(s_nd_num)
        if m_nd_num is None:
            self.m_nd_num = None
        else:
            self.m_nd_num = int(m_nd_num)
        self.nodes = nodes
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

    def __init__(self, osi, sdof, mdof, kn, kt, phi, s_nd_num: int=None, m_nd_num: int=None, nodes=None):
        """
        Initial method for ZeroLengthInterface2Ddof

        Parameters
        ----------
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
        s_nd_num: int
            Number of slave nodes
        m_nd_num: int
            Number of master nodes
        nodes: listi
            Slave and master node tags respectively
        """
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
        self.nodes = nodes
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
        ele_nodes: listi
            A list of a slave and a master nodes * ``1`` if out-normal vector of master plane points to +x direction *
            ``2`` if out-normal vector of master plane points to +y direction * ``3`` if out-normal vector of master plane points
            to +z direction
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
        """
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
