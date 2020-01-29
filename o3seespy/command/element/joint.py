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
        ele_nodes: listi
            A list of four element nodes
        mat1: int
            Uniaxial material tag for left bar-slip spring at node 1
        mat2: int
            Uniaxial material tag for right bar-slip spring at node 1
        mat3: int
            Uniaxial material tag for interface-shear spring at node 1
        mat4: int
            Uniaxial material tag for lower bar-slip spring at node 2
        mat5: int
            Uniaxial material tag for upper bar-slip spring at node 2
        mat6: int
            Uniaxial material tag for interface-shear spring at node 2
        mat7: int
            Uniaxial material tag for left bar-slip spring at node 3
        mat8: int
            Uniaxial material tag for right bar-slip spring at node 3
        mat9: int
            Uniaxial material tag for interface-shear spring at node 3
        mat10: int
            Uniaxial material tag for lower bar-slip spring at node 4
        mat11: int
            Uniaxial material tag for upper bar-slip spring at node 4
        mat12: int
            Uniaxial material tag for interface-shear spring at node 4
        mat13: int
            Uniaxial material tag for shear-panel
        ele_height_fac: float
            Floating point value (as a ratio to the total height of the element) to be considered for determination of
            the distance in between the tension-compression couples (optional, default: 1.0)
        ele_width_fac: float
            Floating point value (as a ratio to the total width of the element) to be considered for determination of
            the distance in between the tension-compression couples (optional, default: 1.0)
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.mat1 = int(mat1)
        self.mat2 = int(mat2)
        self.mat3 = int(mat3)
        self.mat4 = int(mat4)
        self.mat5 = int(mat5)
        self.mat6 = int(mat6)
        self.mat7 = int(mat7)
        self.mat8 = int(mat8)
        self.mat9 = int(mat9)
        self.mat10 = int(mat10)
        self.mat11 = int(mat11)
        self.mat12 = int(mat12)
        self.mat13 = int(mat13)
        self.ele_height_fac = float(ele_height_fac)
        self.ele_width_fac = float(ele_width_fac)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mat1, self.mat2, self.mat3, self.mat4, self.mat5, self.mat6, self.mat7, self.mat8, self.mat9, self.mat10, self.mat11, self.mat12, self.mat13, self.ele_height_fac, self.ele_width_fac]
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
        ele_nodes: listi
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
        """
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
        ele_nodes: listi
            A list of five element nodes = ``[nd1,nd2,nd3,nd4,ndc]``. ``ndc`` is the central node of beam-column joint.
            (the tag ``ndc`` is used to generate the internal node, thus, the node should not exist in the domain or be used by
            any other node)
        mat1: int
            Uniaxial material tag for interface rotational spring at node 1. use a zero tag to indicate the case that a
            beam-column element is rigidly framed to the joint. (optional)
        mat2: int
            Uniaxial material tag for interface rotational spring at node 2. use a zero tag to indicate the case that a
            beam-column element is rigidly framed to the joint. (optional)
        mat3: int
            Uniaxial material tag for interface rotational spring at node 3. use a zero tag to indicate the case that a
            beam-column element is rigidly framed to the joint. (optional)
        mat4: int
            Uniaxial material tag for interface rotational spring at node 4. use a zero tag to indicate the case that a
            beam-column element is rigidly framed to the joint. (optional)
        mat_c: int
            Uniaxial material tag for rotational spring of the central node that describes shear panel behavior
        lrg_dsp: obj
            An integer indicating the flag for considering large deformations: * ``0`` - for small deformations and
            constant geometry * ``1`` - for large deformations and time varying geometry * ``2`` - for large deformations
            ,time varying geometry and length correction
        dmg: obj
            Damage model tag
        dmg1dmg2dmg3dmg4dmg_c: None
            
        """
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
