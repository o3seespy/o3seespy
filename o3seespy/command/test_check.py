from o3seespy.base_model import OpenseesObject


class TestBase(OpenseesObject):
    op_base_type = "test"


class EnergyIncr(TestBase):
    op_type = "EnergyIncr"

    def __init__(self, osi, tol, max_iter, p_flag=0, n_type=2):
        self.tol = float(tol)
        self.max_iter = int(max_iter)  # changed to avoid python function iter
        self.p_flag = int(p_flag)
        self.n_type = int(n_type)
        self._parameters = [self.op_type, self.tol, self.max_iter, self.p_flag, self.n_type]
        self.to_process(osi)


class FixedNumIter(TestBase):
    op_type = "FixedNumIter"

    def __init__(self, osi, max_iter, p_flag=0, n_type=2):
        self.max_iter = int(max_iter)  # changed to avoid python function iter
        self.p_flag = int(p_flag)
        self.n_type = int(n_type)
        self._parameters = [self.op_type, self.max_iter, self.p_flag, self.n_type]
        self.to_process(osi)


class NormDispIncr(TestBase):
    op_type = "NormDispIncr"

    def __init__(self, osi, tol, max_iter, p_flag=0, n_type=2):
        self.tol = float(tol)
        self.max_iter = int(max_iter)  # changed to avoid python function iter
        self.p_flag = int(p_flag)
        self.n_type = int(n_type)
        self._parameters = [self.op_type, self.tol, self.max_iter, self.p_flag, self.n_type]
        self.to_process(osi)


class RelativeNormDispIncr(TestBase):
    op_type = "RelativeNormDispIncr"

    def __init__(self, osi, tol, max_iter, p_flag=0, n_type=2):
        self.tol = float(tol)
        self.max_iter = int(max_iter)  # changed to avoid python function iter
        self.p_flag = int(p_flag)
        self.n_type = int(n_type)
        self._parameters = [self.op_type, self.tol, self.max_iter, self.p_flag, self.n_type]
        self.to_process(osi)


class NormUnbalance(TestBase):
    op_type = "NormUnbalance"

    def __init__(self, osi, tol, max_iter, p_flag=0, n_type=2, max_inc=-1):
        self.tol = float(tol)
        self.max_iter = int(max_iter)  # changed to avoid python function iter
        self.p_flag = int(p_flag)
        self.n_type = int(n_type)
        self.max_inc = int(max_inc)
        self._parameters = [self.op_type, self.tol, self.max_iter, self.p_flag, self.n_type, self.n_type, self.max_inc]
        self.to_process(osi)


class RelativeNormUnbalance(TestBase):
    op_type = "RelativeNormUnbalance"

    def __init__(self, osi, tol, max_iter, p_flag=0, n_type=2):
        self.tol = float(tol)
        self.max_iter = int(max_iter)  # changed to avoid python function iter
        self.p_flag = int(p_flag)
        self.n_type = int(n_type)
        self._parameters = [self.op_type, self.tol, self.max_iter, self.p_flag, self.n_type, self.n_type]
        self.to_process(osi)