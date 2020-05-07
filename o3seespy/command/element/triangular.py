from o3seespy.command.element.base_element import ElementBase


class Tri31(ElementBase):
    """
    The Tri31 Element Class
    
    This command is used to construct a constant strain triangular element (Tri31) which uses three nodes and one
    integration points.

    
    """
    op_type = 'Tri31'

    def __init__(self, osi, ele_nodes, thick, otype, mat, pressure, rho, b1, b2):
        """
        Initial method for Tri31

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of three element nodes in counter-clockwise order
        thick: float
            Element thickness
        otype: str
            String representing material behavior. the type parameter can be either ``'planestrain'`` or
            ``'planestress'``
        mat: obj
            Object of ndmaterial
        pressure: float
            Surface pressure (optional, default = 0.0)
        rho: float
            Element mass density (per unit volume) from which a lumped element mass matrix is computed (optional,
            default=0.0)
        b1: float
            Constant body forces defined in the domain (optional, default=0.0)
        b2: float
            Constant body forces defined in the domain (optional, default=0.0)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [1, 0], [1, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(3)]
        >>> mat = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
        >>> o3.element.Tri31(osi, ele_nodes=ele_nodes, thick=1.0, otype=o3.cc.PLANE_STRAIN, mat=mat, pressure=1.0, rho=1.0, b1=1.0, b2=1.0)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.thick = float(thick)
        self.otype = otype
        self.mat = mat
        self.pressure = float(pressure)
        self.rho = float(rho)
        self.b1 = float(b1)
        self.b2 = float(b2)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.thick, self.otype, self.mat.tag, self.pressure, self.rho, self.b1, self.b2]
        self.to_process(osi)
