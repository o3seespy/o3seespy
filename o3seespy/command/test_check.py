from o3seespy.base_model import OpenSeesObject


class TestBase(OpenSeesObject):
    op_base_type = "test"


class EnergyIncr(TestBase):
    op_type = "EnergyIncr"

    def __init__(self, osi, tol, max_iter, p_flag=0, n_type=2):
        self.osi = osi
        self.tol = float(tol)
        self.max_iter = int(max_iter)  # changed to avoid python function iter
        self.p_flag = int(p_flag)
        self.n_type = int(n_type)
        self._parameters = [self.op_type, self.tol, self.max_iter, self.p_flag, self.n_type]
        self.to_process(osi)


class FixedNumIter(TestBase):
    op_type = "FixedNumIter"

    def __init__(self, osi, max_iter, p_flag=0, n_type=2):
        self.osi = osi
        self.max_iter = int(max_iter)  # changed to avoid python function iter
        self.p_flag = int(p_flag)
        self.n_type = int(n_type)
        self._parameters = [self.op_type, self.max_iter, self.p_flag, self.n_type]
        self.to_process(osi)


class NormDispIncr(TestBase):
    op_type = "NormDispIncr"

    def __init__(self, osi, tol, max_iter, p_flag=0, n_type=2):
        self.osi = osi
        self.tol = float(tol)
        self.max_iter = int(max_iter)  # changed to avoid python function iter
        self.p_flag = int(p_flag)
        self.n_type = int(n_type)
        self._parameters = [self.op_type, self.tol, self.max_iter, self.p_flag, self.n_type]
        self.to_process(osi)


class RelativeNormDispIncr(TestBase):
    op_type = "RelativeNormDispIncr"

    def __init__(self, osi, tol, max_iter, p_flag=0, n_type=2):
        """
        Initial method for RelativeNormDispIncr

        Parameters
        ----------
        tol: float
            Tolerance criteria used to check for convergence.
        max_iter: int
            Max number of iterations to check
        p_flag: int
            Print flag (optional):
            * 0 print nothing.
            * 1 print information on norms each time ``test()`` is invoked.
            * 2 print information on norms and number of iterations at end of successful test.
            * 4 at each step it will print the norms and also the :math:`\\delta u` and :math:`R(u)` vectors.
            * 5 if it fails to converge at end of ``numiter``
                it will print an error message **but return a successfull test**.
        n_type: int
            Type of norm, (0 = max-norm, 1 = 1-norm, 2 = 2-norm). (optional)
        """
        self.osi = osi
        self.tol = float(tol)
        self.max_iter = int(max_iter)  # changed to avoid python function iter
        self.p_flag = int(p_flag)
        self.n_type = int(n_type)
        self._parameters = [self.op_type, self.tol, self.max_iter, self.p_flag, self.n_type]
        self.to_process(osi)


class NormUnbalance(TestBase):
    """
    The NormUnbalance Test Class

    Create a NormUnbalance test, which uses the norm of the right hand side of the matrix equation to determine if
    convergence has been reached.
    """
    op_type = 'NormUnbalance'

    def __init__(self, osi, tol, max_iter, p_flag=0, n_type=2, max_incr: int = None):
        """
        Initial method for NormUnbalance

        Parameters
        ----------
        tol: float
            Tolerance criteria used to check for convergence.
        max_iter: int
            Max number of iterations to check
        p_flag: int
            Print flag (optional): * 0 print nothing. * 1 print information on norms each time ``test()`` is invoked. *
            2 print information on norms and number of iterations at end of successful test. * 4 at each step it will print the
            norms and also the :math:`\\delta u` and :math:`r(u)` vectors. * 5 if it fails to converge at end of ``numiter``
            it will print an error message **but return a successfull test**.
        n_type: int
            Type of norm, (0 = max-norm, 1 = 1-norm, 2 = 2-norm). (optional)
        max_incr: int (default=True)
            Maximum times of error increasing. (optional)
        """
        self.osi = osi
        self.tol = float(tol)
        self.max_iter = int(max_iter)
        self.p_flag = int(p_flag)
        self.n_type = int(n_type)
        if max_incr is None:
            self.max_incr = None
        else:
            self.max_incr = int(max_incr)
        self._parameters = [self.op_type, self.tol, self.max_iter, self.p_flag, self.n_type]
        special_pms = ['max_incr']
        packets = [False]
        for i, pm in enumerate(special_pms):
            if getattr(self, pm) is not None:
                if packets[i]:
                    self._parameters += [*getattr(self, pm)]
                else:
                    self._parameters += [getattr(self, pm)]
            else:
                break
        self.to_process(osi)


class RelativeNormUnbalance(TestBase):
    op_type = "RelativeNormUnbalance"

    def __init__(self, osi, tol, max_iter, p_flag=0, n_type=2):
        self.osi = osi
        self.tol = float(tol)
        self.max_iter = int(max_iter)  # changed to avoid python function iter
        self.p_flag = int(p_flag)
        self.n_type = int(n_type)
        self._parameters = [self.op_type, self.tol, self.max_iter, self.p_flag, self.n_type, self.n_type]
        self.to_process(osi)
