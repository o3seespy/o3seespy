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
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of four element nodes
        mat: obj
            Object of ndmaterial
        b1: float
            Body forces in global x,y,z directions
        b2: float
            Body forces in global x,y,z directions
        b3: float
            Body forces in global x,y,z directions

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
        >>> mat = o3.nd_material.ElasticIsotropic(osi, 1, 0.45)
        >>> o3.element.FourNodeTetrahedron(osi, ele_nodes=ele_nodes, mat=mat, b1=1.0, b2=1.0, b3=1.0)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.mat = mat
        self.b1 = float(b1)
        self.b2 = float(b2)
        self.b3 = float(b3)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mat.tag, self.b1, self.b2, self.b3]
        self.to_process(osi)
