from o3seespy.command.element.base_element import ElementBase


class SurfaceLoad(ElementBase):
    """
    The SurfaceLoad Element Class
    
    This command is used to construct a SurfaceLoad element object.

    
    """
    op_type = 'SurfaceLoad'

    def __init__(self, osi, ele_nodes, p):
        """
        Initial method for SurfaceLoad

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            The four nodes defining the element, input in counterclockwise order (-ndm 3 -ndf 3)
        p: float
            Applied pressure loading normal to the surface, outward is positive, inward is negative

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
        >>> o3.element.SurfaceLoad(osi, ele_nodes=ele_nodes, p=1.0)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.p = float(p)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.p]
        self.to_process(osi)


class VS3D4(ElementBase):
    """
    The VS3D4 Element Class
    
    This command is used to construct a four-node 3D viscous-spring boundary quad element object based on a bilinear
    isoparametric formulation.

    
    """
    op_type = 'VS3D4'

    def __init__(self, osi, ele_nodes, big_e, big_g, rho, big_r, alpha_n, alpha_t):
        """
        Initial method for VS3D4

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            4 end nodes
        big_e: float
            Young's modulus of element material
        big_g: float
            Shear modulus of element material
        rho: float
            Mass density of element material
        big_r: float
            Distance from the scattered wave source to the boundary
        alpha_n: float
            Correction parameter in the normal direction
        alpha_t: float
            Correction parameter in the tangential direction

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
        >>> o3.element.VS3D4(osi, ele_nodes=ele_nodes, big_e=1.0, big_g=1.0, rho=1.0, big_r=1.0, alpha_n=1.0, alpha_t=1.0)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.big_e = float(big_e)
        self.big_g = float(big_g)
        self.rho = float(rho)
        self.big_r = float(big_r)
        self.alpha_n = float(alpha_n)
        self.alpha_t = float(alpha_t)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.big_e, self.big_g, self.rho, self.big_r, self.alpha_n, self.alpha_t]
        self.to_process(osi)


class AC3D8(ElementBase):
    """
    The AC3D8 Element Class
    
    This command is used to construct an eight-node 3D brick acoustic element object based on a trilinear isoparametric
    formulation.

    
    """
    op_type = 'AC3D8'

    def __init__(self, osi, ele_nodes, mat):
        """
        Initial method for AC3D8

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            8 end nodes
        mat: obj
            Material object of previously defined nd material

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
        >>> mat = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
        >>> o3.element.AC3D8(osi, ele_nodes=ele_nodes, mat=mat)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.mat = mat
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mat.tag]
        self.to_process(osi)


class ASI3D8(ElementBase):
    """
    The ASI3D8 Element Class
    
    This command is used to construct an eight-node zero-thickness 3D brick acoustic-structure interface element object
    based on a bilinear isoparametric formulation. The nodes in the acoustic domain share the same coordinates with the
    nodes in the solid domain.

    
    """
    op_type = 'ASI3D8'

    def __init__(self, osi, ele_nodes1, ele_nodes2):
        """
        Initial method for ASI3D8

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes1: None
            
        ele_nodes2: None
            

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.element.ASI3D8(osi, ele_nodes1=1, ele_nodes2=1)
        """
        self.osi = osi
        self.ele_nodes1 = ele_nodes1
        self.ele_nodes2 = ele_nodes2
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes1, *self.ele_nodes2]
        self.to_process(osi)


class AV3D4(ElementBase):
    """
    The AV3D4 Element Class
    
    This command is used to construct a four-node 3D acoustic viscous boundary quad element object based on a bilinear
    isoparametric formulation.

    
    """
    op_type = 'AV3D4'

    def __init__(self, osi, ele_nodes, mat):
        """
        Initial method for AV3D4

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            4 end nodes
        mat: obj
            Material object of previously defined nd material

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
        >>> mat = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
        >>> o3.element.AV3D4(osi, ele_nodes=ele_nodes, mat=mat)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.mat = mat
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mat.tag]
        self.to_process(osi)
