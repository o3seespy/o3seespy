from o3seespy.command.element.base_element import ElementBase


class StdBrick(ElementBase):
    """
    The StdBrick Element Class
    
    This element is used to construct an eight-node brick element object, which uses a trilinear isoparametric
    formulation.

    
    """
    op_type = 'stdBrick'

    def __init__(self, osi, ele_nodes, mat, b1, b2, b3):
        """
        Initial method for StdBrick

        Parameters
        ----------
        ele_nodes: listi
            A list of eight element nodes in bottom and top faces and in counter-clockwise order
        mat: obj
            Tag of ndmaterial
        b1: float
            Body forces in global x,y,z directions
        b2: float
            Body forces in global x,y,z directions
        b3: float
            Body forces in global x,y,z directions
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.mat = mat
        self.b1 = float(b1)
        self.b2 = float(b2)
        self.b3 = float(b3)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mat.tag, self.b1, self.b2, self.b3]
        self.to_process(osi)


class BbarBrick(ElementBase):
    """
    The BbarBrick Element Class
    
    This command is used to construct an eight-node mixed volume/pressure brick element object, which uses a trilinear
    isoparametric formulation.

    
    """
    op_type = 'bbarBrick'

    def __init__(self, osi, ele_nodes, mat, b1, b2, b3):
        """
        Initial method for BbarBrick

        Parameters
        ----------
        ele_nodes: listi
            A list of eight element nodes in bottom and top faces and in counter-clockwise order
        mat: obj
            Tag of ndmaterial
        b1: float
            Body forces in global x,y,z directions
        b2: float
            Body forces in global x,y,z directions
        b3: float
            Body forces in global x,y,z directions
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.mat = mat
        self.b1 = float(b1)
        self.b2 = float(b2)
        self.b3 = float(b3)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mat.tag, self.b1, self.b2, self.b3]
        self.to_process(osi)


class Brick20N(ElementBase):
    """
    The Brick20N Element Class
    
    The element is used to construct a twenty-node three dimensional element object

    
    """
    op_type = 'Brick20N'

    def __init__(self, osi, ele_nodes, mat, bf1, bf2, bf3, mass_den):
        """
        Initial method for Brick20N

        Parameters
        ----------
        ele_nodes: listi
            A list of twenty element nodes, input order is shown in notes below
        mat: obj
            Material tag associated with previsouly-defined ndmaterial object
        bf1: float
            Body force in the direction of global coordinates x, y and z
        bf2: float
            Body force in the direction of global coordinates x, y and z
        bf3: float
            Body force in the direction of global coordinates x, y and z
        mass_den: float
            Mass density (mass/volume)
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.mat = mat
        self.bf1 = float(bf1)
        self.bf2 = float(bf2)
        self.bf3 = float(bf3)
        self.mass_den = float(mass_den)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mat.tag, self.bf1, self.bf2, self.bf3, self.mass_den]
        self.to_process(osi)


class SSPbrick(ElementBase):
    """
    The SSPbrick Element Class
    
    This command is used to construct a SSPbrick element object.

    
    """
    op_type = 'SSPbrick'

    def __init__(self, osi, ele_nodes, mat, b1, b2, b3):
        """
        Initial method for SSPbrick

        Parameters
        ----------
        ele_nodes: listi
            A list of eight element nodes in bottom and top faces and in counter-clockwise order
        mat: obj
            Unique integer tag associated with previously-defined ndmaterial object
        b1: float
            Constant body forces in global x-, y-, and z-directions, respectively (optional, default = 0.0)
        b2: float
            Constant body forces in global x-, y-, and z-directions, respectively (optional, default = 0.0)
        b3: float
            Constant body forces in global x-, y-, and z-directions, respectively (optional, default = 0.0)
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.mat = mat
        self.b1 = float(b1)
        self.b2 = float(b2)
        self.b3 = float(b3)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mat.tag, self.b1, self.b2, self.b3]
        self.to_process(osi)
