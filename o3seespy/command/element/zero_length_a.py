from o3seespy.command.element.base_element import ElementBase


class ZeroLengthdoRayleigh(ElementBase):
    op_type = 'zeroLength'

    def __init__(self, osi, ele_nodes, mat_tags=None, dir_args=None, r_flag=0):
        self.ele_nodes = ele_nodes
        self.mat_tags = mat_tags
        self.dir_args = dir_args
        self.r_flag = float(r_flag)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, '-doRayleigh', self.r_flag]
        if getattr(self, 'mat_tags') is not None:
            self._parameters += ['-mat', *self.mat_tags]
        if getattr(self, 'dir_args') is not None:
            self._parameters += ['-dir', *self.dir_args]
        self.to_process(osi)

class ZeroLengthorient(ElementBase):
    op_type = 'zeroLength'

    def __init__(self, osi, ele_nodes, mat_tags=None, dir_args=None, vecx, vecyp):
        self.ele_nodes = ele_nodes
        self.mat_tags = mat_tags
        self.dir_args = dir_args
        self.vecx = vecx
        self.vecyp = vecyp
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, '-orient', *self.vecx, *self.vecyp]
        if getattr(self, 'mat_tags') is not None:
            self._parameters += ['-mat', *self.mat_tags]
        if getattr(self, 'dir_args') is not None:
            self._parameters += ['-dir', *self.dir_args]
        self.to_process(osi)


class ZeroLengthNDorient(ElementBase):
    op_type = 'zeroLengthND'

    def __init__(self, osi, ele_nodes, mat, uni, vecx, vecyp):
        self.ele_nodes = ele_nodes
        self.mat = mat.tag
        self.uni = uni.tag
        self.vecx = vecx
        self.vecyp = vecyp
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mat.tag, self.uni.tag, '-orient', *self.vecx, self.vecyp]
        self.to_process(osi)


class ZeroLengthSectionorient(ElementBase):
    op_type = 'zeroLengthSection'

    def __init__(self, osi, ele_nodes, sec, vecx, vecyp):
        self.ele_nodes = ele_nodes
        self.sec = sec.tag
        self.vecx = vecx
        self.vecyp = vecyp
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.sec.tag, '-orient', *self.vecx, *self.vecyp]
        self.to_process(osi)

class ZeroLengthSectiondoRayleigh(ElementBase):
    op_type = 'zeroLengthSection'

    def __init__(self, osi, ele_nodes, sec, r_flag):
        self.ele_nodes = ele_nodes
        self.sec = sec.tag
        self.r_flag = float(r_flag)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.sec.tag, '-doRayleigh', self.r_flag]
        self.to_process(osi)


class CoupledZeroLength(ElementBase):
    op_type = 'CoupledZeroLength'

    def __init__(self, osi, ele_nodes, dirn1, dirn2, mat, r_flag=1):
        self.ele_nodes = ele_nodes
        self.dirn1 = int(dirn1)
        self.dirn2 = int(dirn2)
        self.mat = mat.tag
        self.r_flag = float(r_flag)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.dirn1, self.dirn2, self.mat.tag, self.r_flag]
        self.to_process(osi)


class ZeroLengthContact2Dnormal(ElementBase):
    op_type = 'zeroLengthContact2D'

    def __init__(self, osi, ele_nodes, kn, kt, mu, nx, ny):
        self.ele_nodes = ele_nodes
        self.kn = float(kn)
        self.kt = float(kt)
        self.mu = float(mu)
        self.nx = nx
        self.ny = ny
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.kn, self.kt, self.mu, '-normal', self.nx, self.ny]
        self.to_process(osi)

class ZeroLengthContact3D(ElementBase):
    op_type = 'zeroLengthContact3D'

    def __init__(self, osi, ele_nodes, kn, kt, mu, c, dir):
        self.ele_nodes = ele_nodes
        self.kn = float(kn)
        self.kt = float(kt)
        self.mu = float(mu)
        self.c = float(c)
        self.dir = int(dir)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.kn, self.kt, self.mu, self.c, self.dir]
        self.to_process(osi)


class ZeroLengthContactNTS2D(ElementBase):
    op_type = 'zeroLengthContactNTS2D'

    def __init__(self, osi, s_nd_num=None, m_nd_num=None, nodes=None, kn, kt, phi):
        self.s_nd_num = int(s_nd_num)
        self.m_nd_num = int(m_nd_num)
        self.nodes = nodes
        self.kn = float(kn)
        self.kt = float(kt)
        self.phi = float(phi)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.kn, self.kt, self.phi]
        if getattr(self, 's_nd_num') is not None:
            self._parameters += ['-sNdNum', self.s_nd_num]
        if getattr(self, 'm_nd_num') is not None:
            self._parameters += ['-mNdNum', self.m_nd_num]
        if getattr(self, 'nodes') is not None:
            self._parameters += ['-Nodes', *self.nodes]
        self.to_process(osi)


class ZeroLengthInterface2Ddof(ElementBase):
    op_type = 'zeroLengthInterface2D'

    def __init__(self, osi, s_nd_num=None, m_nd_num=None, sdof, mdof, nodes=None, kn, kt, phi):
        self.s_nd_num = int(s_nd_num)
        self.m_nd_num = int(m_nd_num)
        self.sdof = int(sdof)
        self.mdof = int(mdof)
        self.nodes = nodes
        self.kn = float(kn)
        self.kt = float(kt)
        self.phi = float(phi)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, '-dof', self.sdof, self.mdof, self.kn, self.kt, self.phi]
        if getattr(self, 's_nd_num') is not None:
            self._parameters += ['-sNdNum', self.s_nd_num]
        if getattr(self, 'm_nd_num') is not None:
            self._parameters += ['-mNdNum', self.m_nd_num]
        if getattr(self, 'nodes') is not None:
            self._parameters += ['-Nodes', *self.nodes]
        self.to_process(osi)


class ZeroLengthImpact3D(ElementBase):
    op_type = 'zeroLengthImpact3D'

    def __init__(self, osi, ele_nodes, direction, init_gap, friction_ratio, kt, kn, kn2, delta_y, cohesion):
        self.ele_nodes = ele_nodes
        self.direction = direction
        self.init_gap = float(init_gap)
        self.friction_ratio = float(friction_ratio)
        self.kt = float(kt)
        self.kn = float(kn)
        self.kn2 = float(kn2)
        self.delta_y = float(delta_y)
        self.cohesion = float(cohesion)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.direction, self.init_gap, self.friction_ratio, self.kt, self.kn, self.kn2, self.delta_y, self.cohesion]
        self.to_process(osi)
