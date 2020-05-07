from o3seespy.command.element.base_element import ElementBase


class SSPquadUP(ElementBase):
    """
    The SSPquadUP Element Class
    
    This command is used to construct a SSPquadUP element object.

    
    """
    op_type = 'SSPquadUP'

    def __init__(self, osi, ele_nodes, mat, thick, f_bulk, f_den, k1, k2, void, alpha, b1=0.0, b2=0.0):
        """
        Initial method for SSPquadUP

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of four element nodes in counter-clockwise order
        mat: obj
            Unique integer object associated with previously-defined ndmaterial object
        thick: float
            Thickness of the element in out-of-plane direction
        f_bulk: float
            Bulk modulus of the pore fluid
        f_den: float
            Mass density of the pore fluid
        k1: float
            Permeability coefficients in global x- and y-directions, respectively
        k2: float
            Permeability coefficients in global x- and y-directions, respectively
        void: float
            Voids ratio
        alpha: float
            Spatial pressure field stabilization parameter (see discussion below for more information)
        b1: float, optional
            Constant body forces in global x- and y-directions, respectively (optional, default = 0.0) - see note 3
        b2: float, optional
            Constant body forces in global x- and y-directions, respectively (optional, default = 0.0) - see note 3

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> obj = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
        >>> ele_nodes = [o3.node.Node(osi, 0.0, x) for x in range(4)]
        >>> o3.element.SSPquadUP(osi, ele_nodes=ele_nodes, mat=obj, thick=1.0, f_bulk=1.0, f_den=1.0, k1=1.0, k2=1.0, void=1.0, alpha=1.0, b1=0.0, b2=0.0)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.mat = mat
        self.thick = float(thick)
        self.f_bulk = float(f_bulk)
        self.f_den = float(f_den)
        self.k1 = float(k1)
        self.k2 = float(k2)
        self.void = float(void)
        self.alpha = float(alpha)
        self.b1 = float(b1)
        self.b2 = float(b2)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mat.tag, self.thick, self.f_bulk, self.f_den, self.k1, self.k2, self.void, self.alpha, self.b1, self.b2]
        self.to_process(osi)


class SSPbrickUP(ElementBase):
    """
    The SSPbrickUP Element Class
    
    This command is used to construct a SSPbrickUP element object.

    
    """
    op_type = 'SSPbrickUP'

    def __init__(self, osi, ele_nodes, mat, f_bulk, f_den, k1, k2, k3, void, alpha, b1, b2, b3):
        """
        Initial method for SSPbrickUP

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of eight element nodes in counter-clockwise order
        mat: obj
            Unique integer object associated with previously-defined ndmaterial object
        f_bulk: float
            Bulk modulus of the pore fluid
        f_den: float
            Mass density of the pore fluid
        k1: float
            Permeability coefficients in global x-, y-, and z-directions, respectively
        k2: float
            Permeability coefficients in global x-, y-, and z-directions, respectively
        k3: float
            Permeability coefficients in global x-, y-, and z-directions, respectively
        void: float
            Voids ratio
        alpha: float
            Spatial pressure field stabilization parameter (see discussion below for more information)
        b1: float
            Constant body forces in global x-, y-, and z-directions, respectively (optional, default = 0.0) - see note 3
        b2: float
            Constant body forces in global x-, y-, and z-directions, respectively (optional, default = 0.0) - see note 3
        b3: float
            Constant body forces in global x-, y-, and z-directions, respectively (optional, default = 0.0) - see note 3

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> obj = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
        >>> ele_nodes = [o3.node.Node(osi, 0.0, x) for x in range(8)]
        >>> o3.element.SSPbrickUP(osi, ele_nodes=ele_nodes, mat=obj, f_bulk=1.0, f_den=1.0, k1=1.0, k2=1.0, k3=1.0, void=0.5, alpha=1.0e-5, b1=1.0, b2=1.0, b3=1.0)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.mat = mat
        self.f_bulk = float(f_bulk)
        self.f_den = float(f_den)
        self.k1 = float(k1)
        self.k2 = float(k2)
        self.k3 = float(k3)
        self.void = float(void)
        self.alpha = float(alpha)
        self.b1 = float(b1)
        self.b2 = float(b2)
        self.b3 = float(b3)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mat.tag, self.f_bulk, self.f_den, self.k1, self.k2, self.k3, self.void, self.alpha, self.b1, self.b2, self.b3]
        self.to_process(osi)
