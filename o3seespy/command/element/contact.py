from o3seespy.command.element.base_element import ElementBase


class SimpleContact2D(ElementBase):
    op_type = 'SimpleContact2D'

    def __init__(self, osi, i_node, j_node, s_node, l_node, mat, g_tol, f_tol):
        self.i_node = int(i_node)
        self.j_node = int(j_node)
        self.s_node = int(s_node)
        self.l_node = int(l_node)
        self.mat = mat
        self.g_tol = float(g_tol)
        self.f_tol = float(f_tol)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.i_node, self.j_node, self.s_node, self.l_node, self.mat.tag, self.g_tol, self.f_tol]
        self.to_process(osi)


class SimpleContact3D(ElementBase):
    op_type = 'SimpleContact3D'

    def __init__(self, osi, i_node, j_node, k_node, l_node, s_node, lagr_node, mat, g_tol, f_tol):
        self.i_node = int(i_node)
        self.j_node = int(j_node)
        self.k_node = int(k_node)
        self.l_node = int(l_node)
        self.s_node = int(s_node)
        self.lagr_node = int(lagr_node)
        self.mat = mat
        self.g_tol = float(g_tol)
        self.f_tol = float(f_tol)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.i_node, self.j_node, self.k_node, self.l_node, self.s_node, self.lagr_node, self.mat.tag, self.g_tol, self.f_tol]
        self.to_process(osi)


class BeamContact2D(ElementBase):
    op_type = 'BeamContact2D'

    def __init__(self, osi, i_node, j_node, s_node, l_node, mat, width, g_tol, f_tol, c_flag):
        self.i_node = int(i_node)
        self.j_node = int(j_node)
        self.s_node = int(s_node)
        self.l_node = int(l_node)
        self.mat = mat
        self.width = float(width)
        self.g_tol = float(g_tol)
        self.f_tol = float(f_tol)
        self.c_flag = int(c_flag)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.i_node, self.j_node, self.s_node, self.l_node, self.mat.tag, self.width, self.g_tol, self.f_tol, self.c_flag]
        self.to_process(osi)


class BeamContact3D(ElementBase):
    op_type = 'BeamContact3D'

    def __init__(self, osi, i_node, j_node, s_node, l_node, radius, crd_transf, mat, g_tol, f_tol, c_flag):
        self.i_node = int(i_node)
        self.j_node = int(j_node)
        self.s_node = int(s_node)
        self.l_node = int(l_node)
        self.radius = float(radius)
        self.crd_transf = int(crd_transf)
        self.mat = mat
        self.g_tol = float(g_tol)
        self.f_tol = float(f_tol)
        self.c_flag = int(c_flag)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.i_node, self.j_node, self.s_node, self.l_node, self.radius, self.crd_transf, self.mat.tag, self.g_tol, self.f_tol, self.c_flag]
        self.to_process(osi)


class BeamEndContact3D(ElementBase):
    op_type = 'BeamEndContact3D'

    def __init__(self, osi, i_node, j_node, s_node, l_node, radius, g_tol, f_tol, c_flag):
        self.i_node = int(i_node)
        self.j_node = int(j_node)
        self.s_node = int(s_node)
        self.l_node = int(l_node)
        self.radius = float(radius)
        self.g_tol = float(g_tol)
        self.f_tol = float(f_tol)
        self.c_flag = float(c_flag)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.i_node, self.j_node, self.s_node, self.l_node, self.radius, self.g_tol, self.f_tol, self.c_flag]
        self.to_process(osi)
