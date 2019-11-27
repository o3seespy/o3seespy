from o3seespy.command.uniaxial_material.base_material import UniaxialMaterialBase



class Elastic(UniaxialMaterialBase):
    op_type = 'Elastic'

    def __init__(self, osi, big_e, eta=0.0, eneg=None):
        self.big_e = float(big_e)
        self.eta = float(eta)
        if eneg is None:
            self.eneg = None
        else:
            self.eneg = float(eneg)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.big_e, self.eta]
        special_pms = ['eneg']
        packets = [False]
        for i, pm in enumerate(special_pms):
            if getattr(self, pm) is not None:
                if packets[i]:
                    self._parameters += [*getattr(self, pm)]
                else:
                    self._parameters += [getattr(self, pm)]
        self.to_process(osi)


class ElasticPP(UniaxialMaterialBase):
    op_type = 'ElasticPP'

    def __init__(self, osi, big_e, epsy_p, epsy_n=None, eps0=0.0):
        self.big_e = float(big_e)
        self.epsy_p = float(epsy_p)
        if epsy_n is None:
            self.epsy_n = None
        else:
            self.epsy_n = float(epsy_n)
        self.eps0 = float(eps0)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.big_e, self.epsy_p]
        special_pms = ['epsy_n', 'eps0']
        packets = [False, False]
        for i, pm in enumerate(special_pms):
            if getattr(self, pm) is not None:
                if packets[i]:
                    self._parameters += [*getattr(self, pm)]
                else:
                    self._parameters += [getattr(self, pm)]
        self.to_process(osi)


class ElasticPPGap(UniaxialMaterialBase):
    op_type = 'ElasticPPGap'

    def __init__(self, osi, big_e, fy, gap, eta=0.0, damage='noDamage'):
        self.big_e = float(big_e)
        self.fy = float(fy)
        self.gap = float(gap)
        self.eta = float(eta)
        self.damage = damage
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.big_e, self.fy, self.gap, self.eta, self.damage]
        self.to_process(osi)


class ENT(UniaxialMaterialBase):
    op_type = 'ENT'

    def __init__(self, osi, big_e):
        self.big_e = float(big_e)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.big_e]
        self.to_process(osi)


class Parallel(UniaxialMaterialBase):
    op_type = 'Parallel'

    def __init__(self, osi, tags, factor_args=None):
        self.tags = tags
        self.factor_args = factor_args
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, *self.tags]
        if getattr(self, 'factor_args') is not None:
            self._parameters += ['-factor', *self.factor_args]
        self.to_process(osi)


class Series(UniaxialMaterialBase):
    op_type = 'Series'

    def __init__(self, osi, tags):
        self.tags = tags
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, *self.tags]
        self.to_process(osi)
