from o3seespy.command.element.base_element import ElementBase


class ElastomericBearingPlasticityP(ElementBase):

    def __init__(self, osi, ele_nodes, k_init, qd, alpha1, alpha2, mu, mat, p_mat=None, t_mat=None, my_mat=None, mz_mat=None):
        self.ele_nodes = ele_nodes
        self.k_init = float(k_init)
        self.qd = float(qd)
        self.alpha1 = float(alpha1)
        self.alpha2 = float(alpha2)
        self.mu = float(mu)
        self.mat = mat.tag
        self.p_mat = p_mat.tag
        self.t_mat = t_mat.tag
        self.my_mat = my_mat.tag
        self.mz_mat = mz_mat.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.k_init, self.qd, self.alpha1, self.alpha2, self.mu, '-P', self.mat.tag]
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat]
        if getattr(self, 't_mat') is not None:
            self._parameters += ['-T', self.t_mat]
        if getattr(self, 'my_mat') is not None:
            self._parameters += ['-My', self.my_mat]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat]
        self.to_process(osi)

class ElastomericBearingPlasticityMz(ElementBase):

    def __init__(self, osi, ele_nodes, k_init, qd, alpha1, alpha2, mu, mat_tag_2, p_mat=None, t_mat=None, my_mat=None, mz_mat=None):
        self.ele_nodes = ele_nodes
        self.k_init = float(k_init)
        self.qd = float(qd)
        self.alpha1 = float(alpha1)
        self.alpha2 = float(alpha2)
        self.mu = float(mu)
        self.mat_tag_2 = mat_tag_2
        self.p_mat = p_mat.tag
        self.t_mat = t_mat.tag
        self.my_mat = my_mat.tag
        self.mz_mat = mz_mat.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.k_init, self.qd, self.alpha1, self.alpha2, self.mu, '-Mz', self.mat_tag_2]
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat]
        if getattr(self, 't_mat') is not None:
            self._parameters += ['-T', self.t_mat]
        if getattr(self, 'my_mat') is not None:
            self._parameters += ['-My', self.my_mat]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat]
        self.to_process(osi)

class ElastomericBearingPlasticityorient(ElementBase):

    def __init__(self, osi, ele_nodes, k_init, qd, alpha1, alpha2, mu, x1, x2, x3, y1, y2, y3, p_mat=None, t_mat=None, my_mat=None, mz_mat=None):
        self.ele_nodes = ele_nodes
        self.k_init = float(k_init)
        self.qd = float(qd)
        self.alpha1 = float(alpha1)
        self.alpha2 = float(alpha2)
        self.mu = float(mu)
        self.x1 = float(x1)
        self.x2 = float(x2)
        self.x3 = float(x3)
        self.y1 = float(y1)
        self.y2 = float(y2)
        self.y3 = float(y3)
        self.p_mat = p_mat.tag
        self.t_mat = t_mat.tag
        self.my_mat = my_mat.tag
        self.mz_mat = mz_mat.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.k_init, self.qd, self.alpha1, self.alpha2, self.mu, '-orient', self.x1, self.x2, self.x3, self.y1, self.y2, self.y3]
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat]
        if getattr(self, 't_mat') is not None:
            self._parameters += ['-T', self.t_mat]
        if getattr(self, 'my_mat') is not None:
            self._parameters += ['-My', self.my_mat]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat]
        self.to_process(osi)

class ElastomericBearingPlasticityshearDist(ElementBase):

    def __init__(self, osi, ele_nodes, k_init, qd, alpha1, alpha2, mu, s_dratio, doRayleigh=False, p_mat=None, t_mat=None, my_mat=None, mz_mat=None):
        self.ele_nodes = ele_nodes
        self.k_init = float(k_init)
        self.qd = float(qd)
        self.alpha1 = float(alpha1)
        self.alpha2 = float(alpha2)
        self.mu = float(mu)
        self.s_dratio = float(s_dratio)
        self.do_rayleigh = do_rayleigh
        self.p_mat = p_mat.tag
        self.t_mat = t_mat.tag
        self.my_mat = my_mat.tag
        self.mz_mat = mz_mat.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.k_init, self.qd, self.alpha1, self.alpha2, self.mu, '-shearDist', self.s_dratio]
        if getattr(self, 'do_rayleigh') is not None:
            self._parameters += ['-do_rayleigh']
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat]
        if getattr(self, 't_mat') is not None:
            self._parameters += ['-T', self.t_mat]
        if getattr(self, 'my_mat') is not None:
            self._parameters += ['-My', self.my_mat]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat]
        self.to_process(osi)

class ElastomericBearingPlasticitymass(ElementBase):

    def __init__(self, osi, ele_nodes, k_init, qd, alpha1, alpha2, mu, m, p_mat=None, t_mat=None, my_mat=None, mz_mat=None):
        self.ele_nodes = ele_nodes
        self.k_init = float(k_init)
        self.qd = float(qd)
        self.alpha1 = float(alpha1)
        self.alpha2 = float(alpha2)
        self.mu = float(mu)
        self.m = float(m)
        self.p_mat = p_mat.tag
        self.t_mat = t_mat.tag
        self.my_mat = my_mat.tag
        self.mz_mat = mz_mat.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.k_init, self.qd, self.alpha1, self.alpha2, self.mu, '-mass', self.m]
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat]
        if getattr(self, 't_mat') is not None:
            self._parameters += ['-T', self.t_mat]
        if getattr(self, 'my_mat') is not None:
            self._parameters += ['-My', self.my_mat]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat]
        self.to_process(osi)


class ElastomericBearingBoucWentTagMz(ElementBase):

    def __init__(self, osi, ele_nodes, k_init, qd, alpha1, alpha2, mu, eta, beta, mma'-P=False, mat, beat, gamma, p_mat=None, t_mat=None, my_mat=None, mz_mat=None):
        self.ele_nodes = ele_nodes
        self.k_init = float(k_init)
        self.qd = float(qd)
        self.alpha1 = float(alpha1)
        self.alpha2 = float(alpha2)
        self.mu = float(mu)
        self.eta = float(eta)
        self.beta = float(beta)
        self.mma'-p = mma'-p
        self.mat = mat.tag
        self.beat = beat
        self.gamma = float(gamma)
        self.p_mat = p_mat.tag
        self.t_mat = t_mat.tag
        self.my_mat = my_mat.tag
        self.mz_mat = mz_mat.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.k_init, self.qd, self.alpha1, self.alpha2, self.mu, self.eta, self.beta, '-tTag'-Mz', self.mat.tag, self.beat, self.gamma]
        if getattr(self, 'mma'-p') is not None:
            self._parameters += ['-mma'-p']
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat]
        if getattr(self, 't_mat') is not None:
            self._parameters += ['-T', self.t_mat]
        if getattr(self, 'my_mat') is not None:
            self._parameters += ['-My', self.my_mat]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat]
        self.to_process(osi)

class ElastomericBearingBoucWenorient(ElementBase):

    def __init__(self, osi, ele_nodes, k_init, qd, alpha1, alpha2, mu, eta, beta, mma'-P=False, x1, x2, x3, y1, y2, y3, beat, gamma, p_mat=None, t_mat=None, my_mat=None, mz_mat=None):
        self.ele_nodes = ele_nodes
        self.k_init = float(k_init)
        self.qd = float(qd)
        self.alpha1 = float(alpha1)
        self.alpha2 = float(alpha2)
        self.mu = float(mu)
        self.eta = float(eta)
        self.beta = float(beta)
        self.mma'-p = mma'-p
        self.x1 = float(x1)
        self.x2 = float(x2)
        self.x3 = float(x3)
        self.y1 = float(y1)
        self.y2 = float(y2)
        self.y3 = float(y3)
        self.beat = beat
        self.gamma = float(gamma)
        self.p_mat = p_mat.tag
        self.t_mat = t_mat.tag
        self.my_mat = my_mat.tag
        self.mz_mat = mz_mat.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.k_init, self.qd, self.alpha1, self.alpha2, self.mu, self.eta, self.beta, '-orient', self.x1, self.x2, self.x3, self.y1, self.y2, self.y3, self.beat, self.gamma]
        if getattr(self, 'mma'-p') is not None:
            self._parameters += ['-mma'-p']
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat]
        if getattr(self, 't_mat') is not None:
            self._parameters += ['-T', self.t_mat]
        if getattr(self, 'my_mat') is not None:
            self._parameters += ['-My', self.my_mat]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat]
        self.to_process(osi)

class ElastomericBearingBoucWenshearDist(ElementBase):

    def __init__(self, osi, ele_nodes, k_init, qd, alpha1, alpha2, mu, eta, beta, mma'-P=False, s_dratio, doRayleigh=False, beat, gamma, p_mat=None, t_mat=None, my_mat=None, mz_mat=None):
        self.ele_nodes = ele_nodes
        self.k_init = float(k_init)
        self.qd = float(qd)
        self.alpha1 = float(alpha1)
        self.alpha2 = float(alpha2)
        self.mu = float(mu)
        self.eta = float(eta)
        self.beta = float(beta)
        self.mma'-p = mma'-p
        self.s_dratio = float(s_dratio)
        self.do_rayleigh = do_rayleigh
        self.beat = beat
        self.gamma = float(gamma)
        self.p_mat = p_mat.tag
        self.t_mat = t_mat.tag
        self.my_mat = my_mat.tag
        self.mz_mat = mz_mat.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.k_init, self.qd, self.alpha1, self.alpha2, self.mu, self.eta, self.beta, '-shearDist', self.s_dratio, self.beat, self.gamma]
        if getattr(self, 'mma'-p') is not None:
            self._parameters += ['-mma'-p']
        if getattr(self, 'do_rayleigh') is not None:
            self._parameters += ['-do_rayleigh']
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat]
        if getattr(self, 't_mat') is not None:
            self._parameters += ['-T', self.t_mat]
        if getattr(self, 'my_mat') is not None:
            self._parameters += ['-My', self.my_mat]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat]
        self.to_process(osi)

class ElastomericBearingBoucWenmass(ElementBase):

    def __init__(self, osi, ele_nodes, k_init, qd, alpha1, alpha2, mu, eta, beta, mma'-P=False, m, beat, gamma, p_mat=None, t_mat=None, my_mat=None, mz_mat=None):
        self.ele_nodes = ele_nodes
        self.k_init = float(k_init)
        self.qd = float(qd)
        self.alpha1 = float(alpha1)
        self.alpha2 = float(alpha2)
        self.mu = float(mu)
        self.eta = float(eta)
        self.beta = float(beta)
        self.mma'-p = mma'-p
        self.m = float(m)
        self.beat = beat
        self.gamma = float(gamma)
        self.p_mat = p_mat.tag
        self.t_mat = t_mat.tag
        self.my_mat = my_mat.tag
        self.mz_mat = mz_mat.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.k_init, self.qd, self.alpha1, self.alpha2, self.mu, self.eta, self.beta, '-mass', self.m, self.beat, self.gamma]
        if getattr(self, 'mma'-p') is not None:
            self._parameters += ['-mma'-p']
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat]
        if getattr(self, 't_mat') is not None:
            self._parameters += ['-T', self.t_mat]
        if getattr(self, 'my_mat') is not None:
            self._parameters += ['-My', self.my_mat]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat]
        self.to_process(osi)

class ElastomericBearingBoucWenmas(ElementBase):

    def __init__(self, osi, ele_nodes, k_init, qd, alpha1, alpha2, mu, eta, beta, mma'-P=False, m, beat, gamma, p_mat=None, t_mat=None, my_mat=None, mz_mat=None):
        self.ele_nodes = ele_nodes
        self.k_init = float(k_init)
        self.qd = float(qd)
        self.alpha1 = float(alpha1)
        self.alpha2 = float(alpha2)
        self.mu = float(mu)
        self.eta = float(eta)
        self.beta = float(beta)
        self.mma'-p = mma'-p
        self.m = float(m)
        self.beat = beat
        self.gamma = float(gamma)
        self.p_mat = p_mat.tag
        self.t_mat = t_mat.tag
        self.my_mat = my_mat.tag
        self.mz_mat = mz_mat.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.k_init, self.qd, self.alpha1, self.alpha2, self.mu, self.eta, self.beta, '-mas', self.m, self.beat, self.gamma]
        if getattr(self, 'mma'-p') is not None:
            self._parameters += ['-mma'-p']
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat]
        if getattr(self, 't_mat') is not None:
            self._parameters += ['-T', self.t_mat]
        if getattr(self, 'my_mat') is not None:
            self._parameters += ['-My', self.my_mat]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat]
        self.to_process(osi)


class FlatSliderBearingP(ElementBase):

    def __init__(self, osi, ele_nodes, frn_mdl, k_init, mat, p_mat=None, t_mat=None, my_mat=None, mz_mat=None):
        self.ele_nodes = ele_nodes
        self.frn_mdl = frn_mdl.tag
        self.k_init = float(k_init)
        self.mat = mat.tag
        self.p_mat = p_mat.tag
        self.t_mat = t_mat.tag
        self.my_mat = my_mat.tag
        self.mz_mat = mz_mat.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.frn_mdl.tag, self.k_init, '-P', self.mat.tag]
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat]
        if getattr(self, 't_mat') is not None:
            self._parameters += ['-T', self.t_mat]
        if getattr(self, 'my_mat') is not None:
            self._parameters += ['-My', self.my_mat]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat]
        self.to_process(osi)

class FlatSliderBearingMz(ElementBase):

    def __init__(self, osi, ele_nodes, frn_mdl, k_init, mat_tag_2, p_mat=None, t_mat=None, my_mat=None, mz_mat=None):
        self.ele_nodes = ele_nodes
        self.frn_mdl = frn_mdl.tag
        self.k_init = float(k_init)
        self.mat_tag_2 = mat_tag_2
        self.p_mat = p_mat.tag
        self.t_mat = t_mat.tag
        self.my_mat = my_mat.tag
        self.mz_mat = mz_mat.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.frn_mdl.tag, self.k_init, '-Mz', self.mat_tag_2]
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat]
        if getattr(self, 't_mat') is not None:
            self._parameters += ['-T', self.t_mat]
        if getattr(self, 'my_mat') is not None:
            self._parameters += ['-My', self.my_mat]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat]
        self.to_process(osi)

class FlatSliderBearingorient(ElementBase):

    def __init__(self, osi, ele_nodes, frn_mdl, k_init, x1, x2, x3, y1, y2, y3, p_mat=None, t_mat=None, my_mat=None, mz_mat=None):
        self.ele_nodes = ele_nodes
        self.frn_mdl = frn_mdl.tag
        self.k_init = float(k_init)
        self.x1 = float(x1)
        self.x2 = float(x2)
        self.x3 = float(x3)
        self.y1 = float(y1)
        self.y2 = float(y2)
        self.y3 = float(y3)
        self.p_mat = p_mat.tag
        self.t_mat = t_mat.tag
        self.my_mat = my_mat.tag
        self.mz_mat = mz_mat.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.frn_mdl.tag, self.k_init, '-orient', self.x1, self.x2, self.x3, self.y1, self.y2, self.y3]
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat]
        if getattr(self, 't_mat') is not None:
            self._parameters += ['-T', self.t_mat]
        if getattr(self, 'my_mat') is not None:
            self._parameters += ['-My', self.my_mat]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat]
        self.to_process(osi)

class FlatSliderBearingshearDist(ElementBase):

    def __init__(self, osi, ele_nodes, frn_mdl, k_init, s_dratio, doRayleigh=False, p_mat=None, t_mat=None, my_mat=None, mz_mat=None):
        self.ele_nodes = ele_nodes
        self.frn_mdl = frn_mdl.tag
        self.k_init = float(k_init)
        self.s_dratio = float(s_dratio)
        self.do_rayleigh = do_rayleigh
        self.p_mat = p_mat.tag
        self.t_mat = t_mat.tag
        self.my_mat = my_mat.tag
        self.mz_mat = mz_mat.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.frn_mdl.tag, self.k_init, '-shearDist', self.s_dratio]
        if getattr(self, 'do_rayleigh') is not None:
            self._parameters += ['-do_rayleigh']
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat]
        if getattr(self, 't_mat') is not None:
            self._parameters += ['-T', self.t_mat]
        if getattr(self, 'my_mat') is not None:
            self._parameters += ['-My', self.my_mat]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat]
        self.to_process(osi)

class FlatSliderBearingmass(ElementBase):

    def __init__(self, osi, ele_nodes, frn_mdl, k_init, m, p_mat=None, t_mat=None, my_mat=None, mz_mat=None):
        self.ele_nodes = ele_nodes
        self.frn_mdl = frn_mdl.tag
        self.k_init = float(k_init)
        self.m = float(m)
        self.p_mat = p_mat.tag
        self.t_mat = t_mat.tag
        self.my_mat = my_mat.tag
        self.mz_mat = mz_mat.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.frn_mdl.tag, self.k_init, '-mass', self.m]
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat]
        if getattr(self, 't_mat') is not None:
            self._parameters += ['-T', self.t_mat]
        if getattr(self, 'my_mat') is not None:
            self._parameters += ['-My', self.my_mat]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat]
        self.to_process(osi)

class FlatSliderBearingiter(ElementBase):

    def __init__(self, osi, ele_nodes, frn_mdl, k_init, max_iter, tol, p_mat=None, t_mat=None, my_mat=None, mz_mat=None):
        self.ele_nodes = ele_nodes
        self.frn_mdl = frn_mdl.tag
        self.k_init = float(k_init)
        self.max_iter = int(max_iter)
        self.tol = float(tol)
        self.p_mat = p_mat.tag
        self.t_mat = t_mat.tag
        self.my_mat = my_mat.tag
        self.mz_mat = mz_mat.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.frn_mdl.tag, self.k_init, '-iter', self.max_iter, self.tol]
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat]
        if getattr(self, 't_mat') is not None:
            self._parameters += ['-T', self.t_mat]
        if getattr(self, 'my_mat') is not None:
            self._parameters += ['-My', self.my_mat]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat]
        self.to_process(osi)


class SingleFPBearingP(ElementBase):

    def __init__(self, osi, ele_nodes, frn_mdl, reff, k_init, mat, p_mat=None, t_mat=None, my_mat=None, mz_mat=None):
        self.ele_nodes = ele_nodes
        self.frn_mdl = frn_mdl.tag
        self.reff = float(reff)
        self.k_init = float(k_init)
        self.mat = mat.tag
        self.p_mat = p_mat.tag
        self.t_mat = t_mat.tag
        self.my_mat = my_mat.tag
        self.mz_mat = mz_mat.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.frn_mdl.tag, self.reff, self.k_init, '-P', self.mat.tag]
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat]
        if getattr(self, 't_mat') is not None:
            self._parameters += ['-T', self.t_mat]
        if getattr(self, 'my_mat') is not None:
            self._parameters += ['-My', self.my_mat]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat]
        self.to_process(osi)

class SingleFPBearingMz(ElementBase):

    def __init__(self, osi, ele_nodes, frn_mdl, reff, k_init, mat_tag_2, p_mat=None, t_mat=None, my_mat=None, mz_mat=None):
        self.ele_nodes = ele_nodes
        self.frn_mdl = frn_mdl.tag
        self.reff = float(reff)
        self.k_init = float(k_init)
        self.mat_tag_2 = mat_tag_2
        self.p_mat = p_mat.tag
        self.t_mat = t_mat.tag
        self.my_mat = my_mat.tag
        self.mz_mat = mz_mat.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.frn_mdl.tag, self.reff, self.k_init, '-Mz', self.mat_tag_2]
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat]
        if getattr(self, 't_mat') is not None:
            self._parameters += ['-T', self.t_mat]
        if getattr(self, 'my_mat') is not None:
            self._parameters += ['-My', self.my_mat]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat]
        self.to_process(osi)

class SingleFPBearingorient(ElementBase):

    def __init__(self, osi, ele_nodes, frn_mdl, reff, k_init, x1, x2, x3, y1, y2, y3, p_mat=None, t_mat=None, my_mat=None, mz_mat=None):
        self.ele_nodes = ele_nodes
        self.frn_mdl = frn_mdl.tag
        self.reff = float(reff)
        self.k_init = float(k_init)
        self.x1 = float(x1)
        self.x2 = float(x2)
        self.x3 = float(x3)
        self.y1 = float(y1)
        self.y2 = float(y2)
        self.y3 = float(y3)
        self.p_mat = p_mat.tag
        self.t_mat = t_mat.tag
        self.my_mat = my_mat.tag
        self.mz_mat = mz_mat.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.frn_mdl.tag, self.reff, self.k_init, '-orient', self.x1, self.x2, self.x3, self.y1, self.y2, self.y3]
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat]
        if getattr(self, 't_mat') is not None:
            self._parameters += ['-T', self.t_mat]
        if getattr(self, 'my_mat') is not None:
            self._parameters += ['-My', self.my_mat]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat]
        self.to_process(osi)

class SingleFPBearingshearDist(ElementBase):

    def __init__(self, osi, ele_nodes, frn_mdl, reff, k_init, s_dratio, doRayleigh=False, p_mat=None, t_mat=None, my_mat=None, mz_mat=None):
        self.ele_nodes = ele_nodes
        self.frn_mdl = frn_mdl.tag
        self.reff = float(reff)
        self.k_init = float(k_init)
        self.s_dratio = float(s_dratio)
        self.do_rayleigh = do_rayleigh
        self.p_mat = p_mat.tag
        self.t_mat = t_mat.tag
        self.my_mat = my_mat.tag
        self.mz_mat = mz_mat.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.frn_mdl.tag, self.reff, self.k_init, '-shearDist', self.s_dratio]
        if getattr(self, 'do_rayleigh') is not None:
            self._parameters += ['-do_rayleigh']
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat]
        if getattr(self, 't_mat') is not None:
            self._parameters += ['-T', self.t_mat]
        if getattr(self, 'my_mat') is not None:
            self._parameters += ['-My', self.my_mat]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat]
        self.to_process(osi)

class SingleFPBearingmass(ElementBase):

    def __init__(self, osi, ele_nodes, frn_mdl, reff, k_init, m, p_mat=None, t_mat=None, my_mat=None, mz_mat=None):
        self.ele_nodes = ele_nodes
        self.frn_mdl = frn_mdl.tag
        self.reff = float(reff)
        self.k_init = float(k_init)
        self.m = float(m)
        self.p_mat = p_mat.tag
        self.t_mat = t_mat.tag
        self.my_mat = my_mat.tag
        self.mz_mat = mz_mat.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.frn_mdl.tag, self.reff, self.k_init, '-mass', self.m]
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat]
        if getattr(self, 't_mat') is not None:
            self._parameters += ['-T', self.t_mat]
        if getattr(self, 'my_mat') is not None:
            self._parameters += ['-My', self.my_mat]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat]
        self.to_process(osi)

class SingleFPBearingiter(ElementBase):

    def __init__(self, osi, ele_nodes, frn_mdl, reff, k_init, max_iter, tol, p_mat=None, t_mat=None, my_mat=None, mz_mat=None):
        self.ele_nodes = ele_nodes
        self.frn_mdl = frn_mdl.tag
        self.reff = float(reff)
        self.k_init = float(k_init)
        self.max_iter = int(max_iter)
        self.tol = float(tol)
        self.p_mat = p_mat.tag
        self.t_mat = t_mat.tag
        self.my_mat = my_mat.tag
        self.mz_mat = mz_mat.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.frn_mdl.tag, self.reff, self.k_init, '-iter', self.max_iter, self.tol]
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat]
        if getattr(self, 't_mat') is not None:
            self._parameters += ['-T', self.t_mat]
        if getattr(self, 'my_mat') is not None:
            self._parameters += ['-My', self.my_mat]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat]
        self.to_process(osi)


class TFP(ElementBase):

    def __init__(self, osi, ele_nodes, r1, r2, r3, r4, d1, d2, d3, d4, d1, d2, d3, d4, mu1, mu2, mu3, mu4, h1, h2, h3, h4, h0, col_load, big_k):
        self.ele_nodes = ele_nodes
        self.r1 = float(r1)
        self.r2 = float(r2)
        self.r3 = float(r3)
        self.r4 = float(r4)
        self.d1 = float(d1)
        self.d2 = float(d2)
        self.d3 = float(d3)
        self.d4 = float(d4)
        self.d1 = float(d1)
        self.d2 = float(d2)
        self.d3 = float(d3)
        self.d4 = float(d4)
        self.mu1 = float(mu1)
        self.mu2 = float(mu2)
        self.mu3 = float(mu3)
        self.mu4 = float(mu4)
        self.h1 = float(h1)
        self.h2 = float(h2)
        self.h3 = float(h3)
        self.h4 = float(h4)
        self.h0 = float(h0)
        self.col_load = float(col_load)
        self.big_k = float(big_k)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.r1, self.r2, self.r3, self.r4, self.d1, self.d2, self.d3, self.d4, self.d1, self.d2, self.d3, self.d4, self.mu1, self.mu2, self.mu3, self.mu4, self.h1, self.h2, self.h3, self.h4, self.h0, self.col_load, self.big_k]
        self.to_process(osi)


class TripleFrictionPendulum(ElementBase):

    def __init__(self, osi, ele_nodes, frn_tag1, frn_tag2, frn_tag3, vert_mat, rot_z_mat, rot_x_mat, rot_y_mat, l1, l2, l3, d1, d2, d3, big_w, uy, kvt, min_fv, tol):
        self.ele_nodes = ele_nodes
        self.frn_tag1 = int(frn_tag1)
        self.frn_tag2 = int(frn_tag2)
        self.frn_tag3 = int(frn_tag3)
        self.vert_mat = vert_mat.tag
        self.rot_z_mat = rot_z_mat.tag
        self.rot_x_mat = rot_x_mat.tag
        self.rot_y_mat = rot_y_mat.tag
        self.l1 = float(l1)
        self.l2 = float(l2)
        self.l3 = float(l3)
        self.d1 = float(d1)
        self.d2 = float(d2)
        self.d3 = float(d3)
        self.big_w = float(big_w)
        self.uy = float(uy)
        self.kvt = float(kvt)
        self.min_fv = min_fv
        self.tol = float(tol)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.frn_tag1, self.frn_tag2, self.frn_tag3, self.vert_mat.tag, self.rot_z_mat.tag, self.rot_x_mat.tag, self.rot_y_mat.tag, self.l1, self.l2, self.l3, self.d1, self.d2, self.d3, self.big_w, self.uy, self.kvt, self.min_fv, self.tol]
        self.to_process(osi)


class MultipleShearSpringlim(ElementBase):

    def __init__(self, osi, ele_nodes, n_spring, mat=None, dsp):
        self.ele_nodes = ele_nodes
        self.n_spring = int(n_spring)
        self.mat = mat.tag
        self.dsp = float(dsp)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.n_spring, '-lim', self.dsp]
        if getattr(self, 'mat') is not None:
            self._parameters += ['-mat', self.mat]
        self.to_process(osi)

class MultipleShearSpringorient(ElementBase):

    def __init__(self, osi, ele_nodes, n_spring, mat=None, x1, x2, x3, yp1, yp2, yp3):
        self.ele_nodes = ele_nodes
        self.n_spring = int(n_spring)
        self.mat = mat.tag
        self.x1 = float(x1)
        self.x2 = float(x2)
        self.x3 = float(x3)
        self.yp1 = float(yp1)
        self.yp2 = float(yp2)
        self.yp3 = float(yp3)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.n_spring, '-orient', self.x1, self.x2, self.x3, self.yp1, self.yp2, self.yp3]
        if getattr(self, 'mat') is not None:
            self._parameters += ['-mat', self.mat]
        self.to_process(osi)

class MultipleShearSpringmass(ElementBase):

    def __init__(self, osi, ele_nodes, n_spring, mat=None, m):
        self.ele_nodes = ele_nodes
        self.n_spring = int(n_spring)
        self.mat = mat.tag
        self.m = float(m)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.n_spring, '-mass', self.m]
        if getattr(self, 'mat') is not None:
            self._parameters += ['-mat', self.mat]
        self.to_process(osi)


class KikuchiBearingorient(ElementBase):

    def __init__(self, osi, ele_nodes, shape=None, size=None, total_rubber, total_height=None, n_mss=None, mat_mss=None, lim_disp=None, n_mns=None, mat_mns=None, lamb=None, x1, x2, x3, yp1, yp2, yp3):
        self.ele_nodes = ele_nodes
        self.shape = float(shape)
        self.size = float(size)
        self.total_rubber = float(total_rubber)
        self.total_height = float(total_height)
        self.n_mss = int(n_mss)
        self.mat_mss = mat_mss.tag
        self.lim_disp = float(lim_disp)
        self.n_mns = int(n_mns)
        self.mat_mns = mat_mns.tag
        self.lamb = float(lamb)
        self.x1 = float(x1)
        self.x2 = float(x2)
        self.x3 = float(x3)
        self.yp1 = float(yp1)
        self.yp2 = float(yp2)
        self.yp3 = float(yp3)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.total_rubber, '-orient', self.x1, self.x2, self.x3, self.yp1, self.yp2, self.yp3]
        if getattr(self, 'shape') is not None:
            self._parameters += ['-shape', self.shape]
        if getattr(self, 'size') is not None:
            self._parameters += ['-size', self.size]
        if getattr(self, 'total_height') is not None:
            self._parameters += ['-totalHeight', self.total_height]
        if getattr(self, 'n_mss') is not None:
            self._parameters += ['-nMSS', self.n_mss]
        if getattr(self, 'mat_mss') is not None:
            self._parameters += ['-matMSS', self.mat_mss]
        if getattr(self, 'lim_disp') is not None:
            self._parameters += ['-limDisp', self.lim_disp]
        if getattr(self, 'n_mns') is not None:
            self._parameters += ['-nMNS', self.n_mns]
        if getattr(self, 'mat_mns') is not None:
            self._parameters += ['-matMNS', self.mat_mns]
        if getattr(self, 'lamb') is not None:
            self._parameters += ['-lambda', self.lamb]
        self.to_process(osi)

class KikuchiBearingmass(ElementBase):

    def __init__(self, osi, ele_nodes, shape=None, size=None, total_rubber, total_height=None, n_mss=None, mat_mss=None, lim_disp=None, n_mns=None, mat_mns=None, lamb=None, m, noPDInput=False, noTilt=False):
        self.ele_nodes = ele_nodes
        self.shape = float(shape)
        self.size = float(size)
        self.total_rubber = float(total_rubber)
        self.total_height = float(total_height)
        self.n_mss = int(n_mss)
        self.mat_mss = mat_mss.tag
        self.lim_disp = float(lim_disp)
        self.n_mns = int(n_mns)
        self.mat_mns = mat_mns.tag
        self.lamb = float(lamb)
        self.m = float(m)
        self.no_pd_input = no_pd_input
        self.no_tilt = no_tilt
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.total_rubber, '-mass', self.m]
        if getattr(self, 'shape') is not None:
            self._parameters += ['-shape', self.shape]
        if getattr(self, 'size') is not None:
            self._parameters += ['-size', self.size]
        if getattr(self, 'total_height') is not None:
            self._parameters += ['-totalHeight', self.total_height]
        if getattr(self, 'n_mss') is not None:
            self._parameters += ['-nMSS', self.n_mss]
        if getattr(self, 'mat_mss') is not None:
            self._parameters += ['-matMSS', self.mat_mss]
        if getattr(self, 'lim_disp') is not None:
            self._parameters += ['-limDisp', self.lim_disp]
        if getattr(self, 'n_mns') is not None:
            self._parameters += ['-nMNS', self.n_mns]
        if getattr(self, 'mat_mns') is not None:
            self._parameters += ['-matMNS', self.mat_mns]
        if getattr(self, 'lamb') is not None:
            self._parameters += ['-lambda', self.lamb]
        if getattr(self, 'no_pd_input') is not None:
            self._parameters += ['-no_pd_input']
        if getattr(self, 'no_tilt') is not None:
            self._parameters += ['-no_tilt']
        self.to_process(osi)

class KikuchiBearingadjustPDOutput(ElementBase):

    def __init__(self, osi, ele_nodes, shape=None, size=None, total_rubber, total_height=None, n_mss=None, mat_mss=None, lim_disp=None, n_mns=None, mat_mns=None, lamb=None, ci, cj):
        self.ele_nodes = ele_nodes
        self.shape = float(shape)
        self.size = float(size)
        self.total_rubber = float(total_rubber)
        self.total_height = float(total_height)
        self.n_mss = int(n_mss)
        self.mat_mss = mat_mss.tag
        self.lim_disp = float(lim_disp)
        self.n_mns = int(n_mns)
        self.mat_mns = mat_mns.tag
        self.lamb = float(lamb)
        self.ci = float(ci)
        self.cj = float(cj)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.total_rubber, '-adjustPDOutput', self.ci, self.cj]
        if getattr(self, 'shape') is not None:
            self._parameters += ['-shape', self.shape]
        if getattr(self, 'size') is not None:
            self._parameters += ['-size', self.size]
        if getattr(self, 'total_height') is not None:
            self._parameters += ['-totalHeight', self.total_height]
        if getattr(self, 'n_mss') is not None:
            self._parameters += ['-nMSS', self.n_mss]
        if getattr(self, 'mat_mss') is not None:
            self._parameters += ['-matMSS', self.mat_mss]
        if getattr(self, 'lim_disp') is not None:
            self._parameters += ['-limDisp', self.lim_disp]
        if getattr(self, 'n_mns') is not None:
            self._parameters += ['-nMNS', self.n_mns]
        if getattr(self, 'mat_mns') is not None:
            self._parameters += ['-matMNS', self.mat_mns]
        if getattr(self, 'lamb') is not None:
            self._parameters += ['-lambda', self.lamb]
        self.to_process(osi)

class KikuchiBearingdoBalance(ElementBase):

    def __init__(self, osi, ele_nodes, shape=None, size=None, total_rubber, total_height=None, n_mss=None, mat_mss=None, lim_disp=None, n_mns=None, mat_mns=None, lamb=None, lim_fo, lim_fi, n_iter):
        self.ele_nodes = ele_nodes
        self.shape = float(shape)
        self.size = float(size)
        self.total_rubber = float(total_rubber)
        self.total_height = float(total_height)
        self.n_mss = int(n_mss)
        self.mat_mss = mat_mss.tag
        self.lim_disp = float(lim_disp)
        self.n_mns = int(n_mns)
        self.mat_mns = mat_mns.tag
        self.lamb = float(lamb)
        self.lim_fo = float(lim_fo)
        self.lim_fi = float(lim_fi)
        self.n_iter = float(n_iter)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.total_rubber, '-doBalance', self.lim_fo, self.lim_fi, self.n_iter]
        if getattr(self, 'shape') is not None:
            self._parameters += ['-shape', self.shape]
        if getattr(self, 'size') is not None:
            self._parameters += ['-size', self.size]
        if getattr(self, 'total_height') is not None:
            self._parameters += ['-totalHeight', self.total_height]
        if getattr(self, 'n_mss') is not None:
            self._parameters += ['-nMSS', self.n_mss]
        if getattr(self, 'mat_mss') is not None:
            self._parameters += ['-matMSS', self.mat_mss]
        if getattr(self, 'lim_disp') is not None:
            self._parameters += ['-limDisp', self.lim_disp]
        if getattr(self, 'n_mns') is not None:
            self._parameters += ['-nMNS', self.n_mns]
        if getattr(self, 'mat_mns') is not None:
            self._parameters += ['-matMNS', self.mat_mns]
        if getattr(self, 'lamb') is not None:
            self._parameters += ['-lambda', self.lamb]
        self.to_process(osi)


class YamamotoBiaxialHDRcoRS(ElementBase):

    def __init__(self, osi, ele_nodes, tp, d_do, d_di, hr, cr, cs):
        self.ele_nodes = ele_nodes
        self.tp = int(tp)
        self.d_do = float(d_do)
        self.d_di = float(d_di)
        self.hr = float(hr)
        self.cr = float(cr)
        self.cs = float(cs)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.tp, self.d_do, self.d_di, self.hr, '-coRS', self.cr, self.cs]
        self.to_process(osi)

class YamamotoBiaxialHDRorient(ElementBase):

    def __init__(self, osi, ele_nodes, tp, d_do, d_di, hr, vecx, vecyp):
        self.ele_nodes = ele_nodes
        self.tp = int(tp)
        self.d_do = float(d_do)
        self.d_di = float(d_di)
        self.hr = float(hr)
        self.vecx = vecx
        self.vecyp = vecyp
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.tp, self.d_do, self.d_di, self.hr, '-orient', *self.vecx, *self.vecyp]
        self.to_process(osi)

class YamamotoBiaxialHDRmass(ElementBase):

    def __init__(self, osi, ele_nodes, tp, d_do, d_di, hr, m):
        self.ele_nodes = ele_nodes
        self.tp = int(tp)
        self.d_do = float(d_do)
        self.d_di = float(d_di)
        self.hr = float(hr)
        self.m = float(m)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.tp, self.d_do, self.d_di, self.hr, '-mass', self.m]
        self.to_process(osi)


class ElastomericX(ElementBase):

    def __init__(self, osi, ele_nodes, fy, alpha, gr, kbulk, d1, d2, ts, tr, n, x1, x2, x3, y1, y2, y3, kc, phi_m, ac, s_dratio, m, cd, tc, tag1, tag2, tag3, tag4):
        self.ele_nodes = ele_nodes
        self.fy = float(fy)
        self.alpha = float(alpha)
        self.gr = float(gr)
        self.kbulk = float(kbulk)
        self.d1 = float(d1)
        self.d2 = float(d2)
        self.ts = float(ts)
        self.tr = float(tr)
        self.n = int(n)
        self.x1 = float(x1)
        self.x2 = float(x2)
        self.x3 = float(x3)
        self.y1 = float(y1)
        self.y2 = float(y2)
        self.y3 = float(y3)
        self.kc = float(kc)
        self.phi_m = float(phi_m)
        self.ac = float(ac)
        self.s_dratio = float(s_dratio)
        self.m = float(m)
        self.cd = float(cd)
        self.tc = float(tc)
        self.tag1 = float(tag1)
        self.tag2 = float(tag2)
        self.tag3 = float(tag3)
        self.tag4 = float(tag4)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.fy, self.alpha, self.gr, self.kbulk, self.d1, self.d2, self.ts, self.tr, self.n, self.x1, self.x2, self.x3, self.y1, self.y2, self.y3, self.kc, self.phi_m, self.ac, self.s_dratio, self.m, self.cd, self.tc, self.tag1, self.tag2, self.tag3, self.tag4]
        self.to_process(osi)


class LeadRubberX(ElementBase):

    def __init__(self, osi, ele_nodes, fy, alpha, gr, kbulk, d1, d2, ts, tr, n, x1, x2, x3, y1, y2, y3, kc, phi_m, ac, s_dratio, m, cd, tc, q_l, c_l, k_s, a_s, tag1, tag2, tag3, tag4, tag5):
        self.ele_nodes = ele_nodes
        self.fy = float(fy)
        self.alpha = float(alpha)
        self.gr = float(gr)
        self.kbulk = float(kbulk)
        self.d1 = float(d1)
        self.d2 = float(d2)
        self.ts = float(ts)
        self.tr = float(tr)
        self.n = int(n)
        self.x1 = float(x1)
        self.x2 = float(x2)
        self.x3 = float(x3)
        self.y1 = float(y1)
        self.y2 = float(y2)
        self.y3 = float(y3)
        self.kc = float(kc)
        self.phi_m = float(phi_m)
        self.ac = float(ac)
        self.s_dratio = float(s_dratio)
        self.m = float(m)
        self.cd = float(cd)
        self.tc = float(tc)
        self.q_l = float(q_l)
        self.c_l = float(c_l)
        self.k_s = float(k_s)
        self.a_s = float(a_s)
        self.tag1 = int(tag1)
        self.tag2 = int(tag2)
        self.tag3 = int(tag3)
        self.tag4 = int(tag4)
        self.tag5 = int(tag5)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.fy, self.alpha, self.gr, self.kbulk, self.d1, self.d2, self.ts, self.tr, self.n, self.x1, self.x2, self.x3, self.y1, self.y2, self.y3, self.kc, self.phi_m, self.ac, self.s_dratio, self.m, self.cd, self.tc, self.q_l, self.c_l, self.k_s, self.a_s, self.tag1, self.tag2, self.tag3, self.tag4, self.tag5]
        self.to_process(osi)


class HDR(ElementBase):

    def __init__(self, osi, ele_nodes, gr, kbulk, d1, d2, ts, tr, n, a1, a2, a3, b1, b2, b3, c1, c2, c3, c4, x1, x2, x3, y1, y2, y3, kc, phi_m, ac, s_dratio, m, tc):
        self.ele_nodes = ele_nodes
        self.gr = float(gr)
        self.kbulk = float(kbulk)
        self.d1 = float(d1)
        self.d2 = float(d2)
        self.ts = float(ts)
        self.tr = float(tr)
        self.n = int(n)
        self.a1 = float(a1)
        self.a2 = float(a2)
        self.a3 = float(a3)
        self.b1 = float(b1)
        self.b2 = float(b2)
        self.b3 = float(b3)
        self.c1 = float(c1)
        self.c2 = float(c2)
        self.c3 = float(c3)
        self.c4 = float(c4)
        self.x1 = float(x1)
        self.x2 = float(x2)
        self.x3 = float(x3)
        self.y1 = float(y1)
        self.y2 = float(y2)
        self.y3 = float(y3)
        self.kc = float(kc)
        self.phi_m = float(phi_m)
        self.ac = float(ac)
        self.s_dratio = float(s_dratio)
        self.m = float(m)
        self.tc = float(tc)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.gr, self.kbulk, self.d1, self.d2, self.ts, self.tr, self.n, self.a1, self.a2, self.a3, self.b1, self.b2, self.b3, self.c1, self.c2, self.c3, self.c4, self.x1, self.x2, self.x3, self.y1, self.y2, self.y3, self.kc, self.phi_m, self.ac, self.s_dratio, self.m, self.tc]
        self.to_process(osi)


class FPBearingPTV(ElementBase):

    def __init__(self, osi, ele_nodes, mu_ref, is_pressure_dependent, p_ref, is_temperature_dependent, diffusivity, conductivity, is_velocity_dependent, rate_parameter, reffective_fp, radius__contact, k_initial, the_material_a, the_material_b, the_material_c, the_material_d, x1, x2, x3, y1, y2, y3, shear_dist, do_rayleigh, mass, iter, tol, unit):
        self.ele_nodes = ele_nodes
        self.mu_ref = float(mu_ref)
        self.is_pressure_dependent = int(is_pressure_dependent)
        self.p_ref = float(p_ref)
        self.is_temperature_dependent = int(is_temperature_dependent)
        self.diffusivity = float(diffusivity)
        self.conductivity = float(conductivity)
        self.is_velocity_dependent = int(is_velocity_dependent)
        self.rate_parameter = float(rate_parameter)
        self.reffective_fp = float(reffective_fp)
        self.radius__contact = float(radius__contact)
        self.k_initial = float(k_initial)
        self.the_material_a = int(the_material_a)
        self.the_material_b = int(the_material_b)
        self.the_material_c = int(the_material_c)
        self.the_material_d = int(the_material_d)
        self.x1 = float(x1)
        self.x2 = float(x2)
        self.x3 = float(x3)
        self.y1 = float(y1)
        self.y2 = float(y2)
        self.y3 = float(y3)
        self.shear_dist = float(shear_dist)
        self.do_rayleigh = int(do_rayleigh)
        self.mass = float(mass)
        self.iter = int(iter)
        self.tol = float(tol)
        self.unit = int(unit)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.mu_ref, self.is_pressure_dependent, self.p_ref, self.is_temperature_dependent, self.diffusivity, self.conductivity, self.is_velocity_dependent, self.rate_parameter, self.reffective_fp, self.radius__contact, self.k_initial, self.the_material_a, self.the_material_b, self.the_material_c, self.the_material_d, self.x1, self.x2, self.x3, self.y1, self.y2, self.y3, self.shear_dist, self.do_rayleigh, self.mass, self.iter, self.tol, self.unit]
        self.to_process(osi)
