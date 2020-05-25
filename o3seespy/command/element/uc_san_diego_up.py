from o3seespy.command.element.base_element import ElementBase


class QuadUP(ElementBase):
    """
    The QuadUP Element Class
    
    FourNodeQuadUP is a four-node plane-strain element using bilinear isoparametric formulation. This element is
    implemented for simulating dynamic response of solid-fluid fully coupled material, based on Biot's theory of
    porous medium. Each element node has 3 degrees-of-freedom (DOF): DOF 1 and 2 for solid displacement (u) and
    DOF 3 for fluid pressure (p).

    
    """
    op_type = 'quadUP'

    def __init__(self, osi, ele_nodes, thick, mat, bulk, fmass, h_perm, v_perm, b1=0, b2=0, t=0):
        r"""
        Initial method for QuadUP

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of four element nodes in counter-clockwise order
        thick: float
            Element thickness
        mat: obj
            Object of an ndmaterial object (previously defined) of which the element is composed
        bulk: float
            Combined undrained bulk modulus bc relating changes in pore pressure and volumetric strain, may be
            approximated by: :math:`b_c \approx b_f/n` where :math:`b_f` is the bulk modulus of fluid phase
            (:math:`2.2\times 10^6` kpa (or :math:`3.191\times 10^5` psi) for water), and n the initial porosity.
        fmass: float
            Fluid mass density
        h_perm: float
            Permeability coefficient in horizontal and vertical directions respectively.
        v_perm: float
            Permeability coefficient in horizontal and vertical directions respectively.
        b1: float, optional
            Optional gravity acceleration components in horizontal and vertical directions respectively (defaults are
            0.0)
        b2: float, optional
            Optional gravity acceleration components in horizontal and vertical directions respectively (defaults are
            0.0)
        t: float, optional
            Optional uniform element normal traction, positive in tension (default is 0.0)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> mat = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
        >>> coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
        >>> o3.element.QuadUP(osi, ele_nodes=ele_nodes, thick=1.0, mat=mat, bulk=1.0, fmass=1.0, h_perm=1.0, v_perm=1.0, b1=0, b2=0, t=0)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.thick = float(thick)
        self.mat = mat
        self.bulk = float(bulk)
        self.fmass = float(fmass)
        self.h_perm = float(h_perm)
        self.v_perm = float(v_perm)
        self.b1 = float(b1)
        self.b2 = float(b2)
        self.t = float(t)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.thick, self.mat.tag, self.bulk, self.fmass, self.h_perm, self.v_perm, self.b1, self.b2, self.t]
        self.to_process(osi)


class BrickUP(ElementBase):
    """
    The BrickUP Element Class
    
    BrickUP is an 8-node hexahedral linear isoparametric element. Each node has 4 degrees-of-freedom (DOF): DOFs 1 to 3
    for solid displacement (u) and DOF 4 for fluid pressure (p). This element is implemented for simulating dynamic
    response of solid-fluid fully coupled material, based on Biot's theory of porous medium.

    
    """
    op_type = 'brickUP'

    def __init__(self, osi, ele_nodes, mat, bulk, fmass, perm_x, perm_y, perm_z, b_x=0, b_y=0, b_z=0):
        r"""
        Initial method for BrickUP

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of eight element nodes
        mat: obj
            Object of an ndmaterial object (previously defined) of which the element is composed
        bulk: float
            Combined undrained bulk modulus bc relating changes in pore pressure and volumetric strain, may be
            approximated by: :math:`b_c \approx b_f/n` where :math:`b_f` is the bulk modulus of fluid phase
            (:math:`2.2\times 10^6` kpa (or :math:`3.191\times 10^5` psi) for water), and n the initial porosity.
        fmass: float
            Fluid mass density
        perm_x: float
            Permeability coefficients in x, y, and z directions respectively.
        perm_y: float
            Permeability coefficients in x, y, and z directions respectively.
        perm_z: float
            Permeability coefficients in x, y, and z directions respectively.
        b_x: float, optional
            Optional gravity acceleration components in x, y, and z directions directions respectively (defaults are
            0.0)
        b_y: float, optional
            Optional gravity acceleration components in x, y, and z directions directions respectively (defaults are
            0.0)
        b_z: float, optional
            Optional gravity acceleration components in x, y, and z directions directions respectively (defaults are
            0.0)

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=3, ndf=6)
        >>> coords = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
        >>> mat = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
        >>> o3.element.BrickUP(osi, ele_nodes=ele_nodes, mat=mat, bulk=1.0, fmass=1.0, perm_x=1.0, perm_y=1.0, perm_z=1.0, b_x=0, b_y=0, b_z=0)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.mat = mat
        self.bulk = float(bulk)
        self.fmass = float(fmass)
        self.perm_x = float(perm_x)
        self.perm_y = float(perm_y)
        self.perm_z = float(perm_z)
        self.b_x = float(b_x)
        self.b_y = float(b_y)
        self.b_z = float(b_z)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mat.tag, self.bulk, self.fmass, self.perm_x, self.perm_y, self.perm_z, self.b_x, self.b_y, self.b_z]
        self.to_process(osi)


class BbarQuadUP(ElementBase):
    """
    The BbarQuadUP Element Class
    
    bbarQuadUP is a four-node plane-strain mixed volume/pressure element, which uses a tri-linear isoparametric
    formulation. This element is implemented for simulating dynamic response of solid-fluid fully coupled
    material, based on Biot's theory of porous medium. Each element node has 3 degrees-of-freedom (DOF):
    DOF 1 and 2 for solid displacement (u) and DOF 3 for fluid pressure (p).

    
    """
    op_type = 'bbarQuadUP'

    def __init__(self, osi, ele_nodes, thick, mat, bulk, fmass, h_perm, v_perm, b1=0, b2=0, t=0):
        r"""
        Initial method for BbarQuadUP

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of four element nodes in counter-clockwise order
        thick: float
            Element thickness
        mat: obj
            Object of an ndmaterial object (previously defined) of which the element is composed
        bulk: float
            Combined undrained bulk modulus bc relating changes in pore pressure and volumetric strain, may be
            approximated by: :math:`b_c \approx b_f/n` where :math:`b_f` is the bulk modulus of fluid phase
            (:math:`2.2\times 10^6` kpa (or :math:`3.191\times 10^5` psi) for water), and n the initial porosity.
        fmass: float
            Fluid mass density
        h_perm: float
            Permeability coefficient in horizontal and vertical directions respectively.
        v_perm: float
            Permeability coefficient in horizontal and vertical directions respectively.
        b1: float, optional
            Optional gravity acceleration components in horizontal and vertical directions respectively (defaults are
            0.0)
        b2: float, optional
            Optional gravity acceleration components in horizontal and vertical directions respectively (defaults are
            0.0)
        t: float, optional
            Optional uniform element normal traction, positive in tension (default is 0.0)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> mat = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
        >>> coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
        >>> o3.element.BbarQuadUP(osi, ele_nodes=ele_nodes, thick=1.0, mat=mat, bulk=1.0, fmass=1.0, h_perm=1.0, v_perm=1.0, b1=0, b2=0, t=0)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.thick = float(thick)
        self.mat = mat
        self.bulk = float(bulk)
        self.fmass = float(fmass)
        self.h_perm = float(h_perm)
        self.v_perm = float(v_perm)
        self.b1 = float(b1)
        self.b2 = float(b2)
        self.t = float(t)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.thick, self.mat.tag, self.bulk, self.fmass, self.h_perm, self.v_perm, self.b1, self.b2, self.t]
        self.to_process(osi)


class BbarBrickUP(ElementBase):
    """
    The BbarBrickUP Element Class
    
    bbarBrickUP is a 8-node mixed volume/pressure element, which uses a tri-linear isoparametric formulation.Each node
    has 4 degrees-of-freedom (DOF): DOFs 1 to 3 for solid displacement (u) and DOF 4 for fluid pressure (p). This element
    is implemented for simulating dynamic response of solid-fluid fully coupled material, based on Biot's theory of
    porous medium.

    
    """
    op_type = 'bbarBrickUP'

    def __init__(self, osi, ele_nodes, mat, bulk, fmass, perm_x, perm_y, perm_z, b_x=0, b_y=0, b_z=0):
        r"""
        Initial method for BbarBrickUP

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of eight element nodes
        mat: obj
            Object of an ndmaterial object (previously defined) of which the element is composed
        bulk: float
            Combined undrained bulk modulus bc relating changes in pore pressure and volumetric strain, may be
            approximated by: :math:`b_c \approx b_f/n` where :math:`b_f` is the bulk modulus of fluid phase
            (:math:`2.2\times 10^6` kpa (or :math:`3.191\times 10^5` psi) for water), and n the initial porosity.
        fmass: float
            Fluid mass density
        perm_x: float
            Permeability coefficients in x, y, and z directions respectively.
        perm_y: float
            Permeability coefficients in x, y, and z directions respectively.
        perm_z: float
            Permeability coefficients in x, y, and z directions respectively.
        b_x: float, optional
            Optional gravity acceleration components in x, y, and z directions directions respectively (defaults are
            0.0)
        b_y: float, optional
            Optional gravity acceleration components in x, y, and z directions directions respectively (defaults are
            0.0)
        b_z: float, optional
            Optional gravity acceleration components in x, y, and z directions directions respectively (defaults are
            0.0)

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=3)
        >>> coords = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
        >>> mat = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
        >>> o3.element.BbarBrickUP(osi, ele_nodes=ele_nodes, mat=mat, bulk=1.0, fmass=1.0, perm_x=1.0, perm_y=1.0, perm_z=1.0, b_x=0, b_y=0, b_z=0)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.mat = mat
        self.bulk = float(bulk)
        self.fmass = float(fmass)
        self.perm_x = float(perm_x)
        self.perm_y = float(perm_y)
        self.perm_z = float(perm_z)
        self.b_x = float(b_x)
        self.b_y = float(b_y)
        self.b_z = float(b_z)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mat.tag, self.bulk, self.fmass, self.perm_x, self.perm_y, self.perm_z, self.b_x, self.b_y, self.b_z]
        self.to_process(osi)


class N94QuadUP(ElementBase):
    """
    The N94QuadUP Element Class
    
    Nine_Four_Node_QuadUP is a 9-node quadrilateral plane-strain element. The four corner nodes have 3
    degrees-of-freedom (DOF) each: DOF 1 and 2 for solid displacement (u) and DOF 3 for fluid pressure
    (p). The other five nodes have 2 DOFs each for solid displacement. This element is implemented
    for simulating dynamic response of solid-fluid fully coupled material, based on Biot's theory of porous medium.

    
    """
    op_type = '9_4_QuadUP'

    def __init__(self, osi, ele_nodes, thick, mat, bulk, fmass, h_perm, v_perm, b1=0, b2=0):
        r"""
        Initial method for N94QuadUP

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of nine element nodes
        thick: float
            Element thickness
        mat: obj
            Object of an ndmaterial object (previously defined) of which the element is composed
        bulk: float
            Combined undrained bulk modulus bc relating changes in pore pressure and volumetric strain, may be
            approximated by: :math:`b_c \approx b_f/n` where :math:`b_f` is the bulk modulus of fluid phase
            (:math:`2.2\times 10^6` kpa (or :math:`3.191\times 10^5` psi) for water), and n the initial porosity.
        fmass: float
            Fluid mass density
        h_perm: float
            Permeability coefficient in horizontal and vertical directions respectively.
        v_perm: float
            Permeability coefficient in horizontal and vertical directions respectively.
        b1: float, optional
            Optional gravity acceleration components in horizontal and vertical directions respectively (defaults are
            0.0)
        b2: float, optional
            Optional gravity acceleration components in horizontal and vertical directions respectively (defaults are
            0.0)

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1], [0.5, 0.5, 0.5]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
        >>> mat = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
        >>> o3.element.N94QuadUP(osi, ele_nodes=ele_nodes, thick=1.0, mat=mat, bulk=1.0, fmass=1.0, h_perm=1.0, v_perm=1.0, b1=0, b2=0)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.thick = float(thick)
        self.mat = mat
        self.bulk = float(bulk)
        self.fmass = float(fmass)
        self.h_perm = float(h_perm)
        self.v_perm = float(v_perm)
        self.b1 = float(b1)
        self.b2 = float(b2)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.thick, self.mat.tag, self.bulk, self.fmass, self.h_perm, self.v_perm, self.b1, self.b2]
        self.to_process(osi)


class N208BrickUP(ElementBase):
    """
    The N208BrickUP Element Class
    
    Twenty_Eight_Node_BrickUP is a 20-node hexahedral isoparametric element.The eight corner nodes have 4
    degrees-of-freedom (DOF) each: DOFs 1 to 3 for solid displacement (u) and DOF 4 for fluid pressure (p).
    The other nodes have 3 DOFs each for solid displacement. This element is implemented for simulating
    dynamic response of solid-fluid fully coupled material, based on Biot's theory of porous medium.

    
    """
    op_type = '20_8_BrickUP'

    def __init__(self, osi, ele_nodes, mat, bulk, fmass, perm_x, perm_y, perm_z, b_x=0, b_y=0, b_z=0):
        r"""
        Initial method for N208BrickUP

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of twenty element nodes
        mat: obj
            Object of an ndmaterial object (previously defined) of which the element is composed
        bulk: float
            Combined undrained bulk modulus bc relating changes in pore pressure and volumetric strain, may be
            approximated by: :math:`b_c \approx b_f/n` where :math:`b_f` is the bulk modulus of fluid phase
            (:math:`2.2\times 10^6` kpa (or :math:`3.191\times 10^5` psi) for water), and n the initial porosity.
        fmass: float
            Fluid mass density
        perm_x: float
            Permeability coefficients in x, y, and z directions respectively.
        perm_y: float
            Permeability coefficients in x, y, and z directions respectively.
        perm_z: float
            Permeability coefficients in x, y, and z directions respectively.
        b_x: float, optional
            Optional gravity acceleration components in x, y, and z directions directions respectively (defaults are
            0.0)
        b_y: float, optional
            Optional gravity acceleration components in x, y, and z directions directions respectively (defaults are
            0.0)
        b_z: float, optional
            Optional gravity acceleration components in x, y, and z directions directions respectively (defaults are
            0.0)

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=3)
        >>> coords = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
        >>> mat = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
        >>> o3.element.N208BrickUP(osi, ele_nodes=ele_nodes, mat=mat, bulk=1.0, fmass=1.0, perm_x=1.0, perm_y=1.0, perm_z=1.0, b_x=0, b_y=0, b_z=0)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.mat = mat
        self.bulk = float(bulk)
        self.fmass = float(fmass)
        self.perm_x = float(perm_x)
        self.perm_y = float(perm_y)
        self.perm_z = float(perm_z)
        self.b_x = float(b_x)
        self.b_y = float(b_y)
        self.b_z = float(b_z)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mat.tag, self.bulk, self.fmass, self.perm_x, self.perm_y, self.perm_z, self.b_x, self.b_y, self.b_z]
        self.to_process(osi)
