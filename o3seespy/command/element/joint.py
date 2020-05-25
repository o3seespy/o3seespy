from o3seespy.command.element.base_element import ElementBase


class BeamColumnJoint(ElementBase):
    """
    The BeamColumnJoint Element Class
    
    This command is used to construct a two-dimensional beam-column-joint element object. The element may be used with
    both two-dimensional and three-dimensional structures; however, load is transferred only in the plane of the
    element.

    
    """
    op_type = 'beamColumnJoint'

    def __init__(self, osi, ele_nodes, mat1, mat2, mat3, mat4, mat5, mat6, mat7, mat8, mat9, mat10, mat11, mat12, mat13, ele_height_fac=1.0, ele_width_fac=1.0):
        """
        Initial method for BeamColumnJoint

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of four element nodes
        mat1: obj
            Uniaxial material object for left bar-slip spring at node 1
        mat2: obj
            Uniaxial material object for right bar-slip spring at node 1
        mat3: obj
            Uniaxial material object for interface-shear spring at node 1
        mat4: obj
            Uniaxial material object for lower bar-slip spring at node 2
        mat5: obj
            Uniaxial material object for upper bar-slip spring at node 2
        mat6: obj
            Uniaxial material object for interface-shear spring at node 2
        mat7: obj
            Uniaxial material object for left bar-slip spring at node 3
        mat8: obj
            Uniaxial material object for right bar-slip spring at node 3
        mat9: obj
            Uniaxial material object for interface-shear spring at node 3
        mat10: obj
            Uniaxial material object for lower bar-slip spring at node 4
        mat11: obj
            Uniaxial material object for upper bar-slip spring at node 4
        mat12: obj
            Uniaxial material object for interface-shear spring at node 4
        mat13: obj
            Uniaxial material object for shear-panel
        ele_height_fac: float, optional
            Floating point value (as a ratio to the total height of the element) to be considered for determination of
            the distance in between the tension-compression couples (optional, default: 1.0)
        ele_width_fac: float, optional
            Floating point value (as a ratio to the total width of the element) to be considered for determination of
            the distance in between the tension-compression couples (optional, default: 1.0)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
        >>> mats = [o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None) for x in range(13)]
        >>> o3.element.BeamColumnJoint(osi, ele_nodes, *mats, ele_height_fac=1.0, ele_width_fac=1.0)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.mat1 = mat1
        self.mat2 = mat2
        self.mat3 = mat3
        self.mat4 = mat4
        self.mat5 = mat5
        self.mat6 = mat6
        self.mat7 = mat7
        self.mat8 = mat8
        self.mat9 = mat9
        self.mat10 = mat10
        self.mat11 = mat11
        self.mat12 = mat12
        self.mat13 = mat13
        self.ele_height_fac = float(ele_height_fac)
        self.ele_width_fac = float(ele_width_fac)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mat1.tag, self.mat2.tag, self.mat3.tag, self.mat4.tag, self.mat5.tag, self.mat6.tag, self.mat7.tag, self.mat8.tag, self.mat9.tag, self.mat10.tag, self.mat11.tag, self.mat12.tag, self.mat13.tag, self.ele_height_fac, self.ele_width_fac]
        self.to_process(osi)


class ElasticTubularJoint(ElementBase):
    """
    The ElasticTubularJoint Element Class
    
    This command is used to construct an ElasticTubularJoint element object, which models joint flexibility of tubular
    joints in two dimensional analysis of any structure having tubular joints.

    
    """
    op_type = 'ElasticTubularJoint'

    def __init__(self, osi, ele_nodes, brace__diameter, brace__angle, big_e, chord__diameter, chord__thickness, chord__angle):
        """
        Initial method for ElasticTubularJoint

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of two element nodes
        brace__diameter: float
            Outer diameter of brace
        brace__angle: float
            Angle between brace and chord axis 0 < brace_angle < 90
        big_e: float
            Young's modulus
        chord__diameter: float
            Outer diameter of chord
        chord__thickness: float
            Thickness of chord
        chord__angle: float
            Angle between chord axis and global x-axis 0 < chord_angle < 180

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
        >>> o3.element.ElasticTubularJoint(osi, ele_nodes=ele_nodes, brace__diameter=1.0, brace__angle=1.0, big_e=1.0, chord__diameter=1.0, chord__thickness=1.0, chord__angle=1.0)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.brace__diameter = float(brace__diameter)
        self.brace__angle = float(brace__angle)
        self.big_e = float(big_e)
        self.chord__diameter = float(chord__diameter)
        self.chord__thickness = float(chord__thickness)
        self.chord__angle = float(chord__angle)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.brace__diameter, self.brace__angle, self.big_e, self.chord__diameter, self.chord__thickness, self.chord__angle]
        self.to_process(osi)


class Joint2D(ElementBase):
    """
    The Joint2D Element Class
    
    This command is used to construct a two-dimensional beam-column-joint element object. The two dimensional
    beam-column joint is idealized as a parallelogram shaped shear panel with adjacent elements connected to its
    mid-points. The midpoints of the parallelogram are referred to as external nodes. These nodes are the only
    analysis components that connect the joint element to the surrounding structure.

    
    """
    op_type = 'Joint2D'

    def __init__(self, osi, ele_nodes, mat1, mat2, mat3, mat4, mat_c, lrg_dsp, dmg, dmg1dmg2dmg3dmg4dmg_c=None):
        """
        Initial method for Joint2D

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of five element nodes = ``[nd1,nd2,nd3,nd4,ndc]``. ``ndc`` is the central node of beam-column joint.
            (the object ``ndc`` is used to generate the internal node, thus, the node should not exist in the domain or be used by
            any other node)
        mat1: int
            Uniaxial material object for interface rotational spring at node 1. use a zero object to indicate the case
            that a beam-column element is rigidly framed to the joint. 
        mat2: int
            Uniaxial material object for interface rotational spring at node 2. use a zero object to indicate the case
            that a beam-column element is rigidly framed to the joint. 
        mat3: int
            Uniaxial material object for interface rotational spring at node 3. use a zero object to indicate the case
            that a beam-column element is rigidly framed to the joint. 
        mat4: int
            Uniaxial material object for interface rotational spring at node 4. use a zero object to indicate the case
            that a beam-column element is rigidly framed to the joint. 
        mat_c: int
            Uniaxial material object for rotational spring of the central node that describes shear panel behavior
        lrg_dsp: obj
            An integer indicating the flag for considering large deformations: * ``0`` - for small deformations and
            constant geometry * ``1`` - for large deformations and time varying geometry * ``2`` - for large deformations
            ,time varying geometry and length correction
        dmg: obj
            Damage model object
        dmg1dmg2dmg3dmg4dmg_c: None, optional
            

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [1, 0], [1, 1], [0, 1], [0.5, 0.5]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
        >>> o3.element.Joint2D(osi, ele_nodes=ele_nodes, mat1=1, mat2=1, mat3=1, mat4=1, mat_c=1, lrg_dsp='', dmg='', dmg1dmg2dmg3dmg4dmg_c=1)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.mat1 = int(mat1)
        self.mat2 = int(mat2)
        self.mat3 = int(mat3)
        self.mat4 = int(mat4)
        self.mat_c = int(mat_c)
        self.lrg_dsp = lrg_dsp
        self.dmg = dmg
        self.dmg1dmg2dmg3dmg4dmg_c = dmg1dmg2dmg3dmg4dmg_c
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mat1, self.mat2, self.mat3, self.mat4, self.mat_c, self.lrg_dsp.tag, self.dmg.tag]
        if getattr(self, 'dmg1dmg2dmg3dmg4dmg_c') is not None:
            self._parameters += ['-damage', self.dmg1dmg2dmg3dmg4dmg_c]
        self.to_process(osi)
