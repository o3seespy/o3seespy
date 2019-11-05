from o3seespy.command.element.base_element import ElementBase


class BeamColumnJoint(ElementBase):
    op_type = 'beamColumnJoint'

    def __init__(self, osi, ele_nodes, mat1, mat2, mat3, mat4, mat5, mat6, mat7, mat8, mat9, mat10, mat11, mat12, mat13, ele_height_fac=1.0, ele_width_fac=1.0):
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
    op_type = 'ElasticTubularJoint'

    def __init__(self, osi, ele_nodes, brace__diameter, brace__angle, big_e, chord__diameter, chord__thickness, chord__angle):
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
    op_type = 'Joint2D'

    def __init__(self, osi, ele_nodes, mat1, mat2, mat3, mat4, mat_c, lrg_dsp):
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.mat1 = int(mat1)
        self.mat2 = int(mat2)
        self.mat3 = int(mat3)
        self.mat4 = int(mat4)
        self.mat_c = int(mat_c)
        self.lrg_dsp = lrg_dsp
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mat1, self.mat2, self.mat3, self.mat4, self.mat_c, self.lrg_dsp.tag]
        self.to_process(osi)
