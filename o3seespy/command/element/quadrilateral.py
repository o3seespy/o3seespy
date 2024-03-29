from o3seespy.command.element.base_element import ElementBase


class Quad(ElementBase):
    """
    The Quad Element Class
    
    This command is used to construct a FourNodeQuad element object which uses a bilinear isoparametric formulation.

  
     
    """
    op_type = 'quad'

    def __init__(self, osi, ele_nodes, thick, otype, mat, pressure=0.0, rho=0.0, b1=0.0, b2=0.0):
        """
        Initial method for Quad

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of four element nodes in counter-clockwise order
        thick: float
            Element thickness
        otype: str
            String representing material behavior. the type parameter can be either ``'planestrain'`` or
            ``'planestress'``
        mat: obj
            Object of ndmaterial
        pressure: float, optional
            Surface pressure (optional, default = 0.0)
        rho: float, optional
            Element mass density (per unit volume) from which a lumped element mass matrix is computed (optional,
            default=0.0)
        b1: float, optional
            Constant body forces defined in the isoparametric domain (optional, default=0.0)
        b2: float, optional
            Constant body forces defined in the isoparametric domain (optional, default=0.0)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2, ndf=2)
        >>> obj = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
        >>> coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
        >>> o3.element.Quad(osi, ele_nodes=ele_nodes, thick=1.0, otype='PlaneStrain', mat=obj, pressure=1.0, rho=1.0, b1=0.0, b2=0.0)
        """
        self.osi = osi
        self.ele_node_tags = [x.tag for x in ele_nodes]
        self.ele_nodes = ele_nodes
        self.thick = float(thick)
        self.otype = otype
        self.mat = mat
        self.pressure = float(pressure)
        self.rho = float(rho)
        self.b1 = float(b1)
        self.b2 = float(b2)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_node_tags, self.thick, self.otype, self.mat.tag, self.pressure, self.rho, self.b1, self.b2]
        self.to_process(osi)


class ShellMITC4(ElementBase):
    """
    The ShellMITC4 Element Class
    
    This command is used to construct a ShellMITC4 element object, which uses a bilinear isoparametric formulation in
    combination with a modified shear interpolation to improve thin-plate bending performance.

    
    """
    op_type = 'ShellMITC4'

    def __init__(self, osi, ele_nodes, sec):
        """
        Initial method for ShellMITC4

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of four element nodes in counter-clockwise order
        sec: obj
            Object associated with previously-defined sectionforcedeformation object. currently must be either a
            ``'platefibersection'``, or ``'elasticmembraneplatesection'``

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
        >>> sec = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
        >>> o3.element.ShellMITC4(osi, ele_nodes=ele_nodes, sec=sec)
        """
        self.osi = osi
        self.ele_node_tags = [x.tag for x in ele_nodes]
        self.ele_nodes = ele_nodes
        self.sec = sec
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_node_tags, self.sec.tag]
        self.to_process(osi)


class ShellDKGQ(ElementBase):
    """
    The ShellDKGQ Element Class
    
    This command is used to construct a ShellDKGQ element object, which is a quadrilateral shell element based on the
    theory of generalized conforming element.

    
    """
    op_type = 'ShellDKGQ'

    def __init__(self, osi, ele_nodes, sec):
        """
        Initial method for ShellDKGQ

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of four element nodes in counter-clockwise order
        sec: obj
            Object associated with previously-defined sectionforcedeformation object. currently can be a
            ``'platefibersection'``, a ``'elasticmembraneplatesection'`` and a ``'layeredshell'`` section

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
        >>> sec = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
        >>> o3.element.ShellDKGQ(osi, ele_nodes=ele_nodes, sec=sec)
        """
        self.osi = osi
        self.ele_node_tags = [x.tag for x in ele_nodes]
        self.ele_nodes = ele_nodes
        self.sec = sec
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_node_tags, self.sec.tag]
        self.to_process(osi)


class ShellDKGT(ElementBase):
    """
    The ShellDKGT Element Class
    
    This command is used to construct a ShellDKGT element object, which is a triangular shell element based on the
    theory of generalized conforming element.

    
    """
    op_type = 'ShellDKGT'

    def __init__(self, osi, ele_nodes, sec):
        """
        Initial method for ShellDKGT

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of three element nodes in clockwise or counter-clockwise order
        sec: obj
            Object associated with previously-defined sectionforcedeformation object. currently can be a
            ``'platefibersection'``, a ``'elasticmembraneplatesection'`` and a ``'layeredshell'`` section

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [1, 0], [1, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
        >>> sec = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
        >>> o3.element.ShellDKGT(osi, ele_nodes=ele_nodes, sec=sec)
        """
        self.osi = osi
        self.ele_node_tags = [x.tag for x in ele_nodes]
        self.ele_nodes = ele_nodes
        self.sec = sec
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_node_tags, self.sec.tag]
        self.to_process(osi)


class ShellNLDKGQ(ElementBase):
    """
    The ShellNLDKGQ Element Class
    
    This command is used to construct a ShellNLDKGQ element object accounting for the geometric nonlinearity of large
    deformation using the updated Lagrangian formula, which is developed based on the ShellDKGQ element.

    
    """
    op_type = 'ShellNLDKGQ'

    def __init__(self, osi, ele_nodes, sec):
        """
        Initial method for ShellNLDKGQ

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of four element nodes in counter-clockwise order
        sec: obj
            Object associated with previously-defined sectionforcedeformation object. currently can be a
            ``'platefibersection'``, a ``'elasticmembraneplatesection'`` and a ``'layeredshell'`` section

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
        >>> sec = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
        >>> o3.element.ShellNLDKGQ(osi, ele_nodes=ele_nodes, sec=sec)
        """
        self.osi = osi
        self.ele_node_tags = [x.tag for x in ele_nodes]
        self.ele_nodes = ele_nodes
        self.sec = sec
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_node_tags, self.sec.tag]
        self.to_process(osi)


class ShellNLDKGT(ElementBase):
    """
    The ShellNLDKGT Element Class
    
    This command is used to construct a ShellNLDKGT element object accounting for the geometric nonlinearity of large
    deformation using the updated Lagrangian formula, which is developed based on the ShellDKGT element.

    
    """
    op_type = 'ShellNLDKGT'

    def __init__(self, osi, ele_nodes, sec):
        """
        Initial method for ShellNLDKGT

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of three element nodes in clockwise or counter-clockwise order around the element
        sec: obj
            Object associated with previously-defined sectionforcedeformation object. currently can be a
            ``'platefibersection'``, a ``'elasticmembraneplatesection'`` and a ``'layeredshell'`` section

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [1, 0], [1, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
        >>> sec = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
        >>> o3.element.ShellNLDKGT(osi, ele_nodes=ele_nodes, sec=sec)
        """
        self.osi = osi
        self.ele_node_tags = [x.tag for x in ele_nodes]
        self.ele_nodes = ele_nodes
        self.sec = sec
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_node_tags, self.sec.tag]
        self.to_process(osi)


class ShellNL(ElementBase):
    """
    The ShellNL Element Class
    
    
    """
    op_type = 'ShellNL'

    def __init__(self, osi, ele_nodes, sec):
        """
        Initial method for ShellNL

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of nine element nodes, input is the typical, firstly four corner nodes counter-clockwise, then
            mid-side nodes counter-clockwise and finally the central node
        sec: obj
            Object associated with previously-defined sectionforcedeformation object. currently can be a
            ``'platefibersection'``, a ``'elasticmembraneplatesection'`` and a ``'layeredshell'`` section

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [1, 0], [1, 1], [0, 1], [0.5, 0], [1, 0.5], [0.5, 1], [0, 0.5], [0.5, 0.5]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
        >>> sec = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
        >>> o3.element.ShellNL(osi, ele_nodes=ele_nodes, sec=sec)
        """
        self.osi = osi
        self.ele_node_tags = [x.tag for x in ele_nodes]
        self.ele_nodes = ele_nodes
        self.sec = sec
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_node_tags, self.sec.tag]
        self.to_process(osi)


class BbarQuad(ElementBase):
    """
    The BbarQuad Element Class
    
    This command is used to construct a four-node quadrilateral element object, which uses a bilinear isoparametric
    formulation along with a mixed volume/pressure B-bar assumption. This element is for plane strain problems only.

    
    """
    op_type = 'bbarQuad'

    def __init__(self, osi, ele_nodes, thick, mat):
        """
        Initial method for BbarQuad

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of four element nodes in counter-clockwise order
        thick: float
            Element thickness
        mat: obj
            Object of ndmaterial

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2, ndf=2)
        >>> mat = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
        >>> coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
        >>> o3.element.BbarQuad(osi, ele_nodes=ele_nodes, thick=1.0, mat=mat)
        """
        self.osi = osi
        self.ele_node_tags = [x.tag for x in ele_nodes]
        self.ele_nodes = ele_nodes
        self.thick = float(thick)
        self.mat = mat
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_node_tags, self.thick, self.mat.tag]
        self.to_process(osi)


class EnhancedQuad(ElementBase):
    """
    The EnhancedQuad Element Class
    
    This command is used to construct a four-node quadrilateral element, which uses a bilinear isoparametric formulation
    with enhanced strain modes.

    
    """
    op_type = 'enhancedQuad'

    def __init__(self, osi, ele_nodes, thick, otype, mat):
        """
        Initial method for EnhancedQuad

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of four element nodes in counter-clockwise order
        thick: float
            Element thickness
        otype: str
            String representing material behavior. valid options depend on the ndmaterial object and its available
            material formulations. the type parameter can be either ``'PlaneStrain'`` or ``'PlaneStress'``
        mat: obj
            Object of ndmaterial

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2, ndf=2)
        >>> obj = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
        >>> coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
        >>> o3.element.EnhancedQuad(osi, ele_nodes=ele_nodes, thick=1.0, otype='PlaneStress', mat=obj)
        """
        self.osi = osi
        self.ele_node_tags = [x.tag for x in ele_nodes]
        self.ele_nodes = ele_nodes
        self.thick = float(thick)
        self.otype = otype
        self.mat = mat
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_node_tags, self.thick, self.otype, self.mat.tag]
        self.to_process(osi)


class SSPquad(ElementBase):
    """
    The SSPquad Element Class
    
    This command is used to construct a SSPquad element object.

    
    """
    op_type = 'SSPquad'

    def __init__(self, osi, ele_nodes, mat, otype, thick, b1, b2):
        """
        Initial method for SSPquad

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of four element nodes in counter-clockwise order
        mat: obj
            Unique integer object associated with previously-defined ndmaterial object
        otype: str
            String to relay material behavior to the element, can be either ``'planestrain'`` or ``'planestress'``
        thick: float
            Thickness of the element in out-of-plane direction
        b1: float
            Constant body forces in global x- and y-directions, respectively (optional, default = 0.0)
        b2: float
            Constant body forces in global x- and y-directions, respectively (optional, default = 0.0)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> obj = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
        >>> coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
        >>> o3.element.SSPquad(osi, ele_nodes=ele_nodes, mat=obj, otype='PlaneStrain', thick=1.0, b1=0.0, b2=0.0)
        """
        self.osi = osi
        self.ele_node_tags = [x.tag for x in ele_nodes]
        self.ele_nodes = ele_nodes
        self.mat = mat
        self.otype = otype
        self.thick = float(thick)
        self.b1 = float(b1)
        self.b2 = float(b2)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_node_tags, self.mat.tag, self.otype, self.thick, self.b1, self.b2]
        self.to_process(osi)


class MVLEM3DCoR(ElementBase):
    """
    The MVLEM3DCoR Element Class
    
    | Developed and implemented by: | `Kristijan Kolozvari <mailto:kkolozvari@fullerton.edu>`_ (CSU Fullerton)| Kamiar
    Kalbasi (CSU Fullerton)| Kutay Orakcal (Bogazici University)| John Wallace (UCLA)The MVLEM_3D model (Figure 1a) is a
    three-dimensional four-node element with 24 DOFs for nonlinear analysis of flexure-controlled non-rectangular
    reinforced concrete walls subjected to multidirectional loading. The model is an extension of the
    two-dimensional, two-node Multiple-Vertical-Line-Element-Model (`MVLEM
    
    """
    op_type = 'MVLEM_3D'

    def __init__(self, osi, ele_nodes, m, c, thick: list=None, widths: list=None, rho: list=None, mat_concretes: list=None, mat_steels: list=None, mat_shear=None):
        r"""
        Initial method for MVLEM3DCoR

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of four element nodes defined in the counter-clockwise direction
        m: int
            Number of element uniaxial fibers
        c: float
            Location of center of rotation from the base (optional; default = 0.4 (recommended))
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
        """
        self.osi = osi
        self.ele_node_tags = [x.tag for x in ele_nodes]
        self.ele_nodes = ele_nodes
        self.m = int(m)
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
        self.c = float(c)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_node_tags, self.m, '-CoR', self.c]
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

class MVLEM3DThickMod(ElementBase):
    """
    The MVLEM3DThickMod Element Class
    
    | Developed and implemented by: | `Kristijan Kolozvari <mailto:kkolozvari@fullerton.edu>`_ (CSU Fullerton)| Kamiar
    Kalbasi (CSU Fullerton)| Kutay Orakcal (Bogazici University)| John Wallace (UCLA)The MVLEM_3D model (Figure 1a) is a
    three-dimensional four-node element with 24 DOFs for nonlinear analysis of flexure-controlled non-rectangular
    reinforced concrete walls subjected to multidirectional loading. The model is an extension of the
    two-dimensional, two-node Multiple-Vertical-Line-Element-Model (`MVLEM
    
    """
    op_type = 'MVLEM_3D'

    def __init__(self, osi, ele_nodes, m, t_mod, thick: list=None, widths: list=None, rho: list=None, mat_concretes: list=None, mat_steels: list=None, mat_shear=None):
        r"""
        Initial method for MVLEM3DThickMod

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of four element nodes defined in the counter-clockwise direction
        m: int
            Number of element uniaxial fibers
        t_mod: float
            Thickness multiplier (optional; default = 0.63 equivalent to 0.25ig for out-of-plane bending)
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
        """
        self.osi = osi
        self.ele_node_tags = [x.tag for x in ele_nodes]
        self.ele_nodes = ele_nodes
        self.m = int(m)
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
        self.t_mod = float(t_mod)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_node_tags, self.m, '-ThickMod', self.t_mod]
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

class MVLEM3DPoisson(ElementBase):
    """
    The MVLEM3DPoisson Element Class
    
    | Developed and implemented by: | `Kristijan Kolozvari <mailto:kkolozvari@fullerton.edu>`_ (CSU Fullerton)| Kamiar
    Kalbasi (CSU Fullerton)| Kutay Orakcal (Bogazici University)| John Wallace (UCLA)The MVLEM_3D model (Figure 1a) is a
    three-dimensional four-node element with 24 DOFs for nonlinear analysis of flexure-controlled non-rectangular
    reinforced concrete walls subjected to multidirectional loading. The model is an extension of the
    two-dimensional, two-node Multiple-Vertical-Line-Element-Model (`MVLEM
    
    """
    op_type = 'MVLEM_3D'

    def __init__(self, osi, ele_nodes, m, nu, thick: list=None, widths: list=None, rho: list=None, mat_concretes: list=None, mat_steels: list=None, mat_shear=None):
        r"""
        Initial method for MVLEM3DPoisson

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of four element nodes defined in the counter-clockwise direction
        m: int
            Number of element uniaxial fibers
        nu: float
            Poisson ratio for out-of-plane bending (optional; default = 0.25)
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
        """
        self.osi = osi
        self.ele_node_tags = [x.tag for x in ele_nodes]
        self.ele_nodes = ele_nodes
        self.m = int(m)
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
        self.nu = float(nu)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_node_tags, self.m, '-Poisson', self.nu]
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

class MVLEM3DDensity(ElementBase):
    """
    The MVLEM3DDensity Element Class
    
    | Developed and implemented by: | `Kristijan Kolozvari <mailto:kkolozvari@fullerton.edu>`_ (CSU Fullerton)| Kamiar
    Kalbasi (CSU Fullerton)| Kutay Orakcal (Bogazici University)| John Wallace (UCLA)The MVLEM_3D model (Figure 1a) is a
    three-dimensional four-node element with 24 DOFs for nonlinear analysis of flexure-controlled non-rectangular
    reinforced concrete walls subjected to multidirectional loading. The model is an extension of the
    two-dimensional, two-node Multiple-Vertical-Line-Element-Model (`MVLEM
    
    """
    op_type = 'MVLEM_3D'

    def __init__(self, osi, ele_nodes, m, dens, thick: list=None, widths: list=None, rho: list=None, mat_concretes: list=None, mat_steels: list=None, mat_shear=None):
        r"""
        Initial method for MVLEM3DDensity

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of four element nodes defined in the counter-clockwise direction
        m: int
            Number of element uniaxial fibers
        dens: float
            Density (optional; default = 0.0)
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
        """
        self.osi = osi
        self.ele_node_tags = [x.tag for x in ele_nodes]
        self.ele_nodes = ele_nodes
        self.m = int(m)
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
        self.dens = float(dens)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_node_tags, self.m, '-Density', self.dens]
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


class SFIMVLEM3DCoR(ElementBase):
    """
    The SFIMVLEM3DCoR Element Class
    
    | Developed and implemented by: | `Kristijan Kolozvari <mailto:kkolozvari@fullerton.edu>`_ (CSU Fullerton)| Kamiar
    Kalbasi (CSU Fullerton)| Kutay Orakcal (Bogazici University)| John Wallace (UCLA)The SFI-MVLEM-3D model (Figure 1a)
    is a three-dimensional four-node element with 24 DOFs that incorporates axial-flexural-shear interaction and can
    be used for nonlinear analysis of non-rectangular reinforced concrete walls subjected to multidirectional
    loading. The SFI-MVLEM-3D model is an extension of the two-dimensional, two-node
    Shear-Flexure-Interaction Multiple-Vertical-Line-Element-Model (`SFI-MVLEM
    
    """
    op_type = 'SFI_MVLEM_3D'

    def __init__(self, osi, ele_nodes, m, c, thicks: list=None, widths: list=None, mats: list=None):
        """
        Initial method for SFIMVLEM3DCoR

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of four element nodes defined in the counter-clockwise direction
        m: int
            Number of element uniaxial fibers
        c: float
            Location of center of rotation from the base (optional; default = 0.4 (recommended))
        thicks: list, optional
            A list of ``m`` macro-fiber thicknesses
        widths: list, optional
            A list of ``m`` macro-fiber widths
        mats: list, optional
            A list of ``m`` material objects corresponding to ndmaterial fsam
        """
        self.osi = osi
        self.ele_node_tags = [x.tag for x in ele_nodes]
        self.ele_nodes = ele_nodes
        self.m = int(m)
        self.thicks = thicks
        self.widths = widths
        if mats is None:
            self.mats = None
        else:
            self.mats = [x.tag for x in mats]
        self.c = float(c)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_node_tags, self.m, '-CoR', self.c]
        if getattr(self, 'thicks') is not None:
            self._parameters += ['-thick', *self.thicks]
        if getattr(self, 'widths') is not None:
            self._parameters += ['-width', *self.widths]
        if getattr(self, 'mats') is not None:
            self._parameters += ['-mat', *self.mats]
        self.to_process(osi)

class SFIMVLEM3DThickMod(ElementBase):
    """
    The SFIMVLEM3DThickMod Element Class
    
    | Developed and implemented by: | `Kristijan Kolozvari <mailto:kkolozvari@fullerton.edu>`_ (CSU Fullerton)| Kamiar
    Kalbasi (CSU Fullerton)| Kutay Orakcal (Bogazici University)| John Wallace (UCLA)The SFI-MVLEM-3D model (Figure 1a)
    is a three-dimensional four-node element with 24 DOFs that incorporates axial-flexural-shear interaction and can
    be used for nonlinear analysis of non-rectangular reinforced concrete walls subjected to multidirectional
    loading. The SFI-MVLEM-3D model is an extension of the two-dimensional, two-node
    Shear-Flexure-Interaction Multiple-Vertical-Line-Element-Model (`SFI-MVLEM
    
    """
    op_type = 'SFI_MVLEM_3D'

    def __init__(self, osi, ele_nodes, m, t_mod, thicks: list=None, widths: list=None, mats: list=None):
        """
        Initial method for SFIMVLEM3DThickMod

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of four element nodes defined in the counter-clockwise direction
        m: int
            Number of element uniaxial fibers
        t_mod: float
            Thickness multiplier (optional; default = 0.63 equivalent to 0.25ig for out-of-plane bending)
        thicks: list, optional
            A list of ``m`` macro-fiber thicknesses
        widths: list, optional
            A list of ``m`` macro-fiber widths
        mats: list, optional
            A list of ``m`` material objects corresponding to ndmaterial fsam
        """
        self.osi = osi
        self.ele_node_tags = [x.tag for x in ele_nodes]
        self.ele_nodes = ele_nodes
        self.m = int(m)
        self.thicks = thicks
        self.widths = widths
        if mats is None:
            self.mats = None
        else:
            self.mats = [x.tag for x in mats]
        self.t_mod = float(t_mod)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_node_tags, self.m, '-ThickMod', self.t_mod]
        if getattr(self, 'thicks') is not None:
            self._parameters += ['-thick', *self.thicks]
        if getattr(self, 'widths') is not None:
            self._parameters += ['-width', *self.widths]
        if getattr(self, 'mats') is not None:
            self._parameters += ['-mat', *self.mats]
        self.to_process(osi)

class SFIMVLEM3DPoisson(ElementBase):
    """
    The SFIMVLEM3DPoisson Element Class
    
    | Developed and implemented by: | `Kristijan Kolozvari <mailto:kkolozvari@fullerton.edu>`_ (CSU Fullerton)| Kamiar
    Kalbasi (CSU Fullerton)| Kutay Orakcal (Bogazici University)| John Wallace (UCLA)The SFI-MVLEM-3D model (Figure 1a)
    is a three-dimensional four-node element with 24 DOFs that incorporates axial-flexural-shear interaction and can
    be used for nonlinear analysis of non-rectangular reinforced concrete walls subjected to multidirectional
    loading. The SFI-MVLEM-3D model is an extension of the two-dimensional, two-node
    Shear-Flexure-Interaction Multiple-Vertical-Line-Element-Model (`SFI-MVLEM
    
    """
    op_type = 'SFI_MVLEM_3D'

    def __init__(self, osi, ele_nodes, m, nu, thicks: list=None, widths: list=None, mats: list=None):
        """
        Initial method for SFIMVLEM3DPoisson

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of four element nodes defined in the counter-clockwise direction
        m: int
            Number of element uniaxial fibers
        nu: float
            Poisson ratio for out-of-plane bending (optional; default = 0.25)
        thicks: list, optional
            A list of ``m`` macro-fiber thicknesses
        widths: list, optional
            A list of ``m`` macro-fiber widths
        mats: list, optional
            A list of ``m`` material objects corresponding to ndmaterial fsam
        """
        self.osi = osi
        self.ele_node_tags = [x.tag for x in ele_nodes]
        self.ele_nodes = ele_nodes
        self.m = int(m)
        self.thicks = thicks
        self.widths = widths
        if mats is None:
            self.mats = None
        else:
            self.mats = [x.tag for x in mats]
        self.nu = float(nu)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_node_tags, self.m, '-Poisson', self.nu]
        if getattr(self, 'thicks') is not None:
            self._parameters += ['-thick', *self.thicks]
        if getattr(self, 'widths') is not None:
            self._parameters += ['-width', *self.widths]
        if getattr(self, 'mats') is not None:
            self._parameters += ['-mat', *self.mats]
        self.to_process(osi)

class SFIMVLEM3DDensity(ElementBase):
    """
    The SFIMVLEM3DDensity Element Class
    
    | Developed and implemented by: | `Kristijan Kolozvari <mailto:kkolozvari@fullerton.edu>`_ (CSU Fullerton)| Kamiar
    Kalbasi (CSU Fullerton)| Kutay Orakcal (Bogazici University)| John Wallace (UCLA)The SFI-MVLEM-3D model (Figure 1a)
    is a three-dimensional four-node element with 24 DOFs that incorporates axial-flexural-shear interaction and can
    be used for nonlinear analysis of non-rectangular reinforced concrete walls subjected to multidirectional
    loading. The SFI-MVLEM-3D model is an extension of the two-dimensional, two-node
    Shear-Flexure-Interaction Multiple-Vertical-Line-Element-Model (`SFI-MVLEM
    
    """
    op_type = 'SFI_MVLEM_3D'

    def __init__(self, osi, ele_nodes, m, dens, thicks: list=None, widths: list=None, mats: list=None):
        """
        Initial method for SFIMVLEM3DDensity

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of four element nodes defined in the counter-clockwise direction
        m: int
            Number of element uniaxial fibers
        dens: float
            Density (optional; default = 0.0)
        thicks: list, optional
            A list of ``m`` macro-fiber thicknesses
        widths: list, optional
            A list of ``m`` macro-fiber widths
        mats: list, optional
            A list of ``m`` material objects corresponding to ndmaterial fsam
        """
        self.osi = osi
        self.ele_node_tags = [x.tag for x in ele_nodes]
        self.ele_nodes = ele_nodes
        self.m = int(m)
        self.thicks = thicks
        self.widths = widths
        if mats is None:
            self.mats = None
        else:
            self.mats = [x.tag for x in mats]
        self.dens = float(dens)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_node_tags, self.m, '-Density', self.dens]
        if getattr(self, 'thicks') is not None:
            self._parameters += ['-thick', *self.thicks]
        if getattr(self, 'widths') is not None:
            self._parameters += ['-width', *self.widths]
        if getattr(self, 'mats') is not None:
            self._parameters += ['-mat', *self.mats]
        self.to_process(osi)
