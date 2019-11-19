from o3seespy.base_model import OpenseesObject


class IntegratorBase(OpenseesObject):
    op_base_type = "integrator"


class CentralDifference(IntegratorBase):
    op_type = 'CentralDifference'

    def __init__(self, osi):
        self._parameters = [self.op_type]
        self.to_process(osi)


class Newmark(IntegratorBase):
    op_type = 'Newmark'

    def __init__(self, osi, gamma, beta, form=None):
        self.gamma = float(gamma)
        self.beta = float(beta)
        self.form = form
        self._parameters = [self.op_type, self.gamma, self.beta]
        if getattr(self, 'form') is not None:
            self._parameters += ['-formD', self.form]
        self.to_process(osi)


class HHT(IntegratorBase):
    op_type = 'HHT'

    def __init__(self, osi, alpha, gamma=None, beta=None):
        self.alpha = float(alpha)
        self.gamma = float(gamma)
        self.beta = float(beta)
        self._parameters = [self.op_type, self.alpha]
        special_pms = ['gamma', 'beta']
        packets = [False, False]
        for i, pm in enumerate(special_pms):
            if getattr(self, pm) is not None:
                if packets[i]:
                    self._parameters += [*getattr(self, pm)]
                else:
                    self._parameters += [getattr(self, pm)]
        self.to_process(osi)


class GeneralizedAlpha(IntegratorBase):
    op_type = 'GeneralizedAlpha'

    def __init__(self, osi, alpha_m, alpha_f, gamma=None, beta=None):
        self.alpha_m = float(alpha_m)
        self.alpha_f = float(alpha_f)
        self.gamma = float(gamma)
        self.beta = float(beta)
        self._parameters = [self.op_type, self.alpha_m, self.alpha_f]
        special_pms = ['gamma', 'beta']
        packets = [False, False]
        for i, pm in enumerate(special_pms):
            if getattr(self, pm) is not None:
                if packets[i]:
                    self._parameters += [*getattr(self, pm)]
                else:
                    self._parameters += [getattr(self, pm)]
        self.to_process(osi)


class TRBDF2(IntegratorBase):
    op_type = 'TRBDF2'

    def __init__(self, osi):
        self._parameters = [self.op_type]
        self.to_process(osi)


class ExplicitDifference(IntegratorBase):
    op_type = 'ExplicitDifference'

    def __init__(self, osi):
        self._parameters = [self.op_type]
        self.to_process(osi)
