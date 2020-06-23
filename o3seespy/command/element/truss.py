from o3seespy.command.element.base_element import ElementBase


class Truss(ElementBase):
    """
    The Truss Element Class
    
    This command is used to construct a truss element object. There are two ways to construct a truss element object:

 
      One way is to specify an area and a UniaxialMaterial identifier:
    """
    op_type = 'Truss'

    def __init__(self, osi, ele_nodes, big_a, mat, rho: float=None, c_flag: float=None, r_flag: float=None):
        """
        Initial method for Truss

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of two element nodes
        big_a: float
            Cross-sectional area of element
        mat: obj
            Object associated with previously-defined uniaxial_material
        rho: float, optional
            Mass per unit length, optional, default = 0.0
        c_flag: float, optional
            Consistent mass flag, optional, default = 0 * ``cflag`` = 0 lumped mass matrix (default) * ``cflag`` = 1
            consistent mass matrix
        r_flag: float, optional
            Rayleigh damping flag, optional, default = 0 * ``rflag`` = 0 no rayleigh damping (default) * ``rflag`` = 1
            include rayleigh damping

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [1, 0]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
        >>> mat = o3.uniaxial_material.Elastic(osi, 1.0)
        >>> o3.element.Truss(osi, ele_nodes=ele_nodes, big_a=1.0, mat=mat, rho=1.0, c_flag=1.0, r_flag=1.0)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.big_a = float(big_a)
        self.mat = mat
        if rho is None:
            self.rho = None
        else:
            self.rho = float(rho)
        if c_flag is None:
            self.c_flag = None
        else:
            self.c_flag = float(c_flag)
        if r_flag is None:
            self.r_flag = None
        else:
            self.r_flag = float(r_flag)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.big_a, self.mat.tag]
        if getattr(self, 'rho') is not None:
            self._parameters += ['-rho', self.rho]
        if getattr(self, 'c_flag') is not None:
            self._parameters += ['-cMass', self.c_flag]
        if getattr(self, 'r_flag') is not None:
            self._parameters += ['-doRayleigh', self.r_flag]
        self.to_process(osi)


class CorotTruss(ElementBase):
    """
    The CorotTruss Element Class
    
    This command is used to construct a corotational truss element object. There are two ways to construct a
    corotational truss element object:

    One way is to specify an area and a UniaxialMaterial identifier:
    """
    op_type = 'corotTruss'

    def __init__(self, osi, ele_nodes, big_a, mat, rho: float=None, c_flag: float=None, r_flag: float=None):
        """
        Initial method for CorotTruss

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of two element nodes
        big_a: float
            Cross-sectional area of element
        mat: obj
            Object associated with previously-defined uniaxial_material
        rho: float, optional
            Mass per unit length, optional, default = 0.0
        c_flag: float, optional
            Consistent mass flag, optional, default = 0 * ``cflag`` = 0 lumped mass matrix (default) * ``cflag`` = 1
            consistent mass matrix
        r_flag: float, optional
            Rayleigh damping flag, optional, default = 0 * ``rflag`` = 0 no rayleigh damping (default) * ``rflag`` = 1
            include rayleigh damping

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [1, 0]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
        >>> mat = o3.uniaxial_material.Elastic(osi, 1.0)
        >>> o3.element.CorotTruss(osi, ele_nodes=ele_nodes, big_a=1.0, mat=mat, rho=1.0, c_flag=1.0, r_flag=1.0)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.big_a = float(big_a)
        self.mat = mat
        if rho is None:
            self.rho = None
        else:
            self.rho = float(rho)
        if c_flag is None:
            self.c_flag = None
        else:
            self.c_flag = float(c_flag)
        if r_flag is None:
            self.r_flag = None
        else:
            self.r_flag = float(r_flag)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.big_a, self.mat.tag]
        if getattr(self, 'rho') is not None:
            self._parameters += ['-rho', self.rho]
        if getattr(self, 'c_flag') is not None:
            self._parameters += ['-cMass', self.c_flag]
        if getattr(self, 'r_flag') is not None:
            self._parameters += ['-doRayleigh', self.r_flag]
        self.to_process(osi)
