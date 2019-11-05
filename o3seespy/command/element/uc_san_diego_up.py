from o3seespy.command.element.base_element import ElementBase


class QuadUP(ElementBase):
    op_type = 'quadUP'

    def __init__(self, osi, ele_nodes, thick, mat, bulk, fmass, h_perm, v_perm, b1=0, b2=0, t=0):
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
    op_type = 'brickUP'

    def __init__(self, osi, ele_nodes, mat, bulk, fmass, perm_x, perm_y, perm_z, b_x=0, b_y=0, b_z=0):
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
    op_type = 'bbarQuadUP'

    def __init__(self, osi, ele_nodes, thick, mat, bulk, fmass, h_perm, v_perm, b1=0, b2=0, t=0):
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
    op_type = 'bbarBrickUP'

    def __init__(self, osi, ele_nodes, mat, bulk, fmass, perm_x, perm_y, perm_z, b_x=0, b_y=0, b_z=0):
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
    op_type = '9_4_QuadUP'

    def __init__(self, osi, ele_nodes, thick, mat, bulk, fmass, h_perm, v_perm, b1=0, b2=0):
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
    op_type = '20_8_BrickUP'

    def __init__(self, osi, ele_nodes, mat, bulk, fmass, perm_x, perm_y, perm_z, b_x=0, b_y=0, b_z=0):
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
