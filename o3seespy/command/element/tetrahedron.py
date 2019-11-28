from o3seespy.command.element.base_element import ElementBase


class FourNodeTetrahedron(ElementBase):
    """
    The FourNodeTetrahedron Element Class
    
    This command is used to construct a standard four-node tetrahedron element objec with one-point Gauss integration.


       
    """
    op_type = 'FourNodeTetrahedron'

    def __init__(self, osi, ele_nodes, mat, b1, b2, b3):
        """
        Initial method for FourNodeTetrahedron

        Parameters
        ----------
        ele_nodes: listi
            A list of four element nodes
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
