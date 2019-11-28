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
        ele_nodes: listi
            The four nodes defining the element, input in counterclockwise order (-ndm 3 -ndf 3)
        p: float
            Applied pressure loading normal to the surface, outward is positive, inward is negative
        """
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
        ele_nodes: listi
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
        """
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
        ele_nodes: listi
            8 end nodes
        mat: obj
            Material tag of previously defined nd material
        """
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
        ele_nodes1: None
            
        ele_nodes2: None
            
        """
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
        ele_nodes: listi
            4 end nodes
        mat: obj
            Material tag of previously defined nd material
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.mat = mat
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mat.tag]
        self.to_process(osi)
