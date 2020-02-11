from o3seespy.base_model import OpenSeesObject


class TestBase(OpenSeesObject):
    op_base_type = "test"


class NormUnbalance(TestBase):
    """
    The NormUnbalance Test Class
    
    Create a NormUnbalance test, which uses the norm of the right hand side of the matrix equation to determine if
    convergence has been reached.
    """
    op_type = 'NormUnbalance'

    def __init__(self, osi, tol, max_iter, p_flag=0, n_type=2, max_incr: float=None):
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


class NormDispIncr(TestBase):
    """
    The NormDispIncr Test Class
    
    Create a NormUnbalance test, which uses the norm of the left hand side solution vector of the matrix equation to
    determine if convergence has been reached.
    """
    op_type = 'NormDispIncr'

    def __init__(self, osi, tol, max_iter, p_flag=0, n_type=2):
        """
        Initial method for NormDispIncr

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
        """
        self.tol = float(tol)
        self.max_iter = int(max_iter)
        self.p_flag = int(p_flag)
        self.n_type = int(n_type)
        self._parameters = [self.op_type, self.tol, self.max_iter, self.p_flag, self.n_type]
        self.to_process(osi)


class EnergyIncr(TestBase):
    """
    The EnergyIncr Test Class
    
    Create a EnergyIncr test, which uses the dot product of the solution vector and norm of the right hand side of the
    matrix equation to determine if convergence has been reached.
    """
    op_type = 'EnergyIncr'

    def __init__(self, osi, tol, max_iter, p_flag=0, n_type=2):
        """
        Initial method for EnergyIncr

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
        """
        self.tol = float(tol)
        self.max_iter = int(max_iter)
        self.p_flag = int(p_flag)
        self.n_type = int(n_type)
        self._parameters = [self.op_type, self.tol, self.max_iter, self.p_flag, self.n_type]
        self.to_process(osi)


class RelativeNormUnbalance(TestBase):
    """
    The RelativeNormUnbalance Test Class
    
    Create a RelativeNormUnbalance test, which uses the relative norm of the right hand side of the matrix equation to
    determine if convergence has been reached.
    """
    op_type = 'RelativeNormUnbalance'

    def __init__(self, osi, tol, max_iter, p_flag=0, n_type=2):
        """
        Initial method for RelativeNormUnbalance

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
        """
        self.tol = float(tol)
        self.max_iter = int(max_iter)
        self.p_flag = int(p_flag)
        self.n_type = int(n_type)
        self._parameters = [self.op_type, self.tol, self.max_iter, self.p_flag, self.n_type]
        self.to_process(osi)


class RelativeNormDispIncr(TestBase):
    """
    The RelativeNormDispIncr Test Class
    
    Create a RelativeNormDispIncr test, which uses the relative of the solution vector of the matrix equation to
    determine if convergence has been reached.
    """
    op_type = 'RelativeNormDispIncr'

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
            Print flag (optional): * 0 print nothing. * 1 print information on norms each time ``test()`` is invoked. *
            2 print information on norms and number of iterations at end of successful test. * 4 at each step it will print the
            norms and also the :math:`\\delta u` and :math:`r(u)` vectors. * 5 if it fails to converge at end of ``numiter``
            it will print an error message **but return a successfull test**.
        n_type: int
            Type of norm, (0 = max-norm, 1 = 1-norm, 2 = 2-norm). (optional)
        """
        self.tol = float(tol)
        self.max_iter = int(max_iter)
        self.p_flag = int(p_flag)
        self.n_type = int(n_type)
        self._parameters = [self.op_type, self.tol, self.max_iter, self.p_flag, self.n_type]
        self.to_process(osi)


class RelativeTotalNormDispIncr(TestBase):
    """
    The RelativeTotalNormDispIncr Test Class
    
    Create a RelativeTotalNormDispIncr test, which uses the ratio of the current norm to the total norm (the sum of all
    the norms since last convergence) of the solution vector.
    """
    op_type = 'relativeTotalNormDispIncr'

    def __init__(self, osi, tol, max_iter, p_flag=0, n_type=2):
        """
        Initial method for RelativeTotalNormDispIncr

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
        """
        self.tol = float(tol)
        self.max_iter = int(max_iter)
        self.p_flag = int(p_flag)
        self.n_type = int(n_type)
        self._parameters = [self.op_type, self.tol, self.max_iter, self.p_flag, self.n_type]
        self.to_process(osi)


class RelativeEnergyIncr(TestBase):
    """
    The RelativeEnergyIncr Test Class
    
    Create a RelativeEnergyIncr test, which uses the relative dot product of the solution vector and norm of the right
    hand side of the matrix equation to determine if convergence has been reached.
    """
    op_type = 'RelativeEnergyIncr'

    def __init__(self, osi, tol, max_iter, p_flag=0, n_type=2):
        """
        Initial method for RelativeEnergyIncr

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
        """
        self.tol = float(tol)
        self.max_iter = int(max_iter)
        self.p_flag = int(p_flag)
        self.n_type = int(n_type)
        self._parameters = [self.op_type, self.tol, self.max_iter, self.p_flag, self.n_type]
        self.to_process(osi)


class FixedNumIter(TestBase):
    """
    The FixedNumIter Test Class
    
    Create a FixedNumIter test, that performs a fixed number of iterations without testing for convergence.
    """
    op_type = 'FixedNumIter'

    def __init__(self, osi, max_iter, p_flag=0, n_type=2):
        """
        Initial method for FixedNumIter

        Parameters
        ----------
        max_iter: int
            Max number of iterations to check
        p_flag: int
            Print flag (optional): * 0 print nothing. * 1 print information on norms each time ``test()`` is invoked. *
            2 print information on norms and number of iterations at end of successful test. * 4 at each step it will print the
            norms and also the :math:`\\delta u` and :math:`r(u)` vectors. * 5 if it fails to converge at end of ``numiter``
            it will print an error message **but return a successfull test**.
        n_type: int
            Type of norm, (0 = max-norm, 1 = 1-norm, 2 = 2-norm). (optional)
        """
        self.max_iter = int(max_iter)
        self.p_flag = int(p_flag)
        self.n_type = int(n_type)
        self._parameters = [self.op_type, self.max_iter, self.p_flag, self.n_type]
        self.to_process(osi)


class NormDispAndUnbalance(TestBase):
    """
    The NormDispAndUnbalance Test Class
    
    Create a NormDispAndUnbalance test, which check if both
    """
    op_type = 'NormDispAndUnbalance'

    def __init__(self, osi, tol_incr, tol_r, max_iter, p_flag=0, n_type=2, maxincr=-1):
        """
        Initial method for NormDispAndUnbalance

        Parameters
        ----------
        tol_incr: float
            Tolerance for right hand residual
        tol_r: None
            
        max_iter: int
            Max number of iterations to check
        p_flag: int
            Print flag (optional): * 0 print nothing. * 1 print information on norms each time ``test()`` is invoked. *
            2 print information on norms and number of iterations at end of successful test. * 4 at each step it will print the
            norms and also the :math:`\\delta u` and :math:`r(u)` vectors. * 5 if it fails to converge at end of ``numiter``
            it will print an error message **but return a successfull test**.
        n_type: int
            Type of norm, (0 = max-norm, 1 = 1-norm, 2 = 2-norm). (optional)
        maxincr: int
            Maximum times of error increasing. (optional)
        """
        self.tol_incr = float(tol_incr)
        self.tol_r = tol_r
        self.max_iter = int(max_iter)
        self.p_flag = int(p_flag)
        self.n_type = int(n_type)
        self.maxincr = int(maxincr)
        self._parameters = [self.op_type, self.tol_incr, self.tol_r, self.max_iter, self.p_flag, self.n_type, self.maxincr]
        self.to_process(osi)


class NormDispOrUnbalance(TestBase):
    """
    The NormDispOrUnbalance Test Class
    
    Create a NormDispOrUnbalance test, which check if both
    """
    op_type = 'NormDispOrUnbalance'

    def __init__(self, osi, tol_incr, tol_r, max_iter, p_flag=0, n_type=2, maxincr=-1):
        """
        Initial method for NormDispOrUnbalance

        Parameters
        ----------
        tol_incr: float
            Tolerance for right hand residual
        tol_r: None
            
        max_iter: int
            Max number of iterations to check
        p_flag: int
            Print flag (optional): * 0 print nothing. * 1 print information on norms each time ``test()`` is invoked. *
            2 print information on norms and number of iterations at end of successful test. * 4 at each step it will print the
            norms and also the :math:`\\delta u` and :math:`r(u)` vectors. * 5 if it fails to converge at end of ``numiter``
            it will print an error message **but return a successfull test**.
        n_type: int
            Type of norm, (0 = max-norm, 1 = 1-norm, 2 = 2-norm). (optional)
        maxincr: int
            Maximum times of error increasing. (optional)
        """
        self.tol_incr = float(tol_incr)
        self.tol_r = tol_r
        self.max_iter = int(max_iter)
        self.p_flag = int(p_flag)
        self.n_type = int(n_type)
        self.maxincr = int(maxincr)
        self._parameters = [self.op_type, self.tol_incr, self.tol_r, self.max_iter, self.p_flag, self.n_type, self.maxincr]
        self.to_process(osi)
