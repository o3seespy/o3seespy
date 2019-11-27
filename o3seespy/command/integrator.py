from o3seespy.base_model import OpenseesObject


class IntegratorBase(OpenseesObject):
    op_base_type = "integrator"


class LoadControl(IntegratorBase):
    op_type = 'LoadControl'

    def __init__(self, osi, incr, num_iter=1, min_incr=None, max_incr=None):
        self.incr = float(incr)
        self.num_iter = int(num_iter)
        if min_incr is None:
            self.min_incr = None
        else:
            self.min_incr = float(min_incr)
        if max_incr is None:
            self.max_incr = None
        else:
            self.max_incr = float(max_incr)
        self._parameters = [self.op_type, self.incr, self.num_iter]
        special_pms = ['min_incr', 'max_incr']
        packets = [False, False]
        for i, pm in enumerate(special_pms):
            if getattr(self, pm) is not None:
                if packets[i]:
                    self._parameters += [*getattr(self, pm)]
                else:
                    self._parameters += [getattr(self, pm)]
        self.to_process(osi)


class DisplacementControl(IntegratorBase):
    op_type = 'DisplacementControl'

    def __init__(self, osi, nd, dof, incr, num_iter=1, d_umin=None, d_umax=None):
        self.nd = int(nd)
        self.dof = int(dof)
        self.incr = float(incr)
        self.num_iter = int(num_iter)
        self.d_umin = d_umin
        self.d_umax = d_umax
        self._parameters = [self.op_type, self.nd, self.dof, self.incr, self.num_iter]
        special_pms = ['d_umin', 'd_umax']
        packets = [False, False]
        for i, pm in enumerate(special_pms):
            if getattr(self, pm) is not None:
                if packets[i]:
                    self._parameters += [*getattr(self, pm)]
                else:
                    self._parameters += [getattr(self, pm)]
        self.to_process(osi)


class ParallelDisplacementControl(IntegratorBase):
    op_type = 'ParallelDisplacementControl'

    def __init__(self, osi, nd, dof, incr, num_iter=1, d_umin=None, d_umax=None):
        self.nd = int(nd)
        self.dof = int(dof)
        self.incr = float(incr)
        self.num_iter = int(num_iter)
        self.d_umin = d_umin
        self.d_umax = d_umax
        self._parameters = [self.op_type, self.nd, self.dof, self.incr, self.num_iter]
        special_pms = ['d_umin', 'd_umax']
        packets = [False, False]
        for i, pm in enumerate(special_pms):
            if getattr(self, pm) is not None:
                if packets[i]:
                    self._parameters += [*getattr(self, pm)]
                else:
                    self._parameters += [getattr(self, pm)]
        self.to_process(osi)


class MinUnbalDispNorm(IntegratorBase):
    op_type = 'MinUnbalDispNorm'

    def __init__(self, osi, dlambda1, jd=1, min_lambda=None, max_lambda=None, det=None):
        self.dlambda1 = float(dlambda1)
        self.jd = int(jd)
        if min_lambda is None:
            self.min_lambda = None
        else:
            self.min_lambda = float(min_lambda)
        if max_lambda is None:
            self.max_lambda = None
        else:
            self.max_lambda = float(max_lambda)
        self.det = det
        self._parameters = [self.op_type, self.dlambda1, self.jd]
        special_pms = ['min_lambda', 'max_lambda', 'det']
        packets = [False, False, False]
        for i, pm in enumerate(special_pms):
            if getattr(self, pm) is not None:
                if packets[i]:
                    self._parameters += [*getattr(self, pm)]
                else:
                    self._parameters += [getattr(self, pm)]
        self.to_process(osi)


class ArcLength(IntegratorBase):
    op_type = 'ArcLength'

    def __init__(self, osi, s, alpha):
        self.s = float(s)
        self.alpha = float(alpha)
        self._parameters = [self.op_type, self.s, self.alpha]
        self.to_process(osi)


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
        if gamma is None:
            self.gamma = None
        else:
            self.gamma = float(gamma)
        if beta is None:
            self.beta = None
        else:
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
        if gamma is None:
            self.gamma = None
        else:
            self.gamma = float(gamma)
        if beta is None:
            self.beta = None
        else:
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
