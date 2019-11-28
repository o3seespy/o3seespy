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
        """
        Initial method for QuadUP

        Parameters
        ----------
        ele_nodes: listi
            A list of four element nodes in counter-clockwise order
        thick: float
            Element thickness
        mat: obj
            Tag of an ndmaterial object (previously defined) of which the element is composed
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
        b1: float
            Optional gravity acceleration components in horizontal and vertical directions respectively (defaults are
            0.0)
        b2: float
            Optional gravity acceleration components in horizontal and vertical directions respectively (defaults are
            0.0)
        t: float
            Optional uniform element normal traction, positive in tension (default is 0.0)
        """
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
        """
        Initial method for BrickUP

        Parameters
        ----------
        ele_nodes: listi
            A list of eight element nodes
        mat: obj
            Tag of an ndmaterial object (previously defined) of which the element is composed
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
        b_x: float
            Optional gravity acceleration components in x, y, and z directions directions respectively (defaults are
            0.0)
        b_y: float
            Optional gravity acceleration components in x, y, and z directions directions respectively (defaults are
            0.0)
        b_z: float
            Optional gravity acceleration components in x, y, and z directions directions respectively (defaults are
            0.0)
        """
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
        """
        Initial method for BbarQuadUP

        Parameters
        ----------
        ele_nodes: listi
            A list of four element nodes in counter-clockwise order
        thick: float
            Element thickness
        mat: obj
            Tag of an ndmaterial object (previously defined) of which the element is composed
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
        b1: float
            Optional gravity acceleration components in horizontal and vertical directions respectively (defaults are
            0.0)
        b2: float
            Optional gravity acceleration components in horizontal and vertical directions respectively (defaults are
            0.0)
        t: float
            Optional uniform element normal traction, positive in tension (default is 0.0)
        """
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
        """
        Initial method for BbarBrickUP

        Parameters
        ----------
        ele_nodes: listi
            A list of eight element nodes
        mat: obj
            Tag of an ndmaterial object (previously defined) of which the element is composed
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
        b_x: float
            Optional gravity acceleration components in x, y, and z directions directions respectively (defaults are
            0.0)
        b_y: float
            Optional gravity acceleration components in x, y, and z directions directions respectively (defaults are
            0.0)
        b_z: float
            Optional gravity acceleration components in x, y, and z directions directions respectively (defaults are
            0.0)
        """
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
        """
        Initial method for N94QuadUP

        Parameters
        ----------
        ele_nodes: listi
            A list of nine element nodes
        thick: float
            Element thickness
        mat: obj
            Tag of an ndmaterial object (previously defined) of which the element is composed
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
        b1: float
            Optional gravity acceleration components in horizontal and vertical directions respectively (defaults are
            0.0)
        b2: float
            Optional gravity acceleration components in horizontal and vertical directions respectively (defaults are
            0.0)
        """
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
        """
        Initial method for N208BrickUP

        Parameters
        ----------
        ele_nodes: listi
            A list of twenty element nodes
        mat: obj
            Tag of an ndmaterial object (previously defined) of which the element is composed
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
        b_x: float
            Optional gravity acceleration components in x, y, and z directions directions respectively (defaults are
            0.0)
        b_y: float
            Optional gravity acceleration components in x, y, and z directions directions respectively (defaults are
            0.0)
        b_z: float
            Optional gravity acceleration components in x, y, and z directions directions respectively (defaults are
            0.0)
        """
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
