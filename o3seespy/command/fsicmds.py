
def mesh(osi, tag, args):
    """
    Create a mesh object. See below for available mesh types... toctree:::maxdepth:
    2linemeshtrimeshquadmeshtetmeshpartmeshbgmesh

    Parameters
    ----------
    osi: o3seespy.OpenSeesInstance
    tag: None
            
    args: None
            
    """
    _parameters = [tag, *args]
    return osi.to_process("mesh", _parameters)


def remesh(osi):
    r"""
    * :math:`\alpha \ge 0` for updating moving mesh.* :math:`\alpha < 0` for updating background mesh.If there are nodes
    shared by different mesh in the domain,the principles to decide what element should be used fora triangle:#. If all 3
    nodes share the same mesh, use that mesh for the triangle.

    Parameters
    ----------
    osi: o3seespy.OpenSeesInstance
    """
    _parameters = []
    return osi.to_process("remesh", _parameters)


def integrator(osi):
    """
    Create a PFEM Integrator.

    Parameters
    ----------
    osi: o3seespy.OpenSeesInstance
    """
    _parameters = []
    return osi.to_process("integrator", _parameters)


def system(osi, compressible=False, mumps=False):
    """
    Create a incompressible PFEM system of equations using the Umfpack solver

    Parameters
    ----------
    osi: o3seespy.OpenSeesInstance
    compressible: None
            
    mumps: None
            
    """
    _parameters = []
    if compressible:
        _parameters += ['-compressible']
    if mumps:
        _parameters += ['-mumps']
    return osi.to_process("system", _parameters)


def test(osi, tolv, tolp, tolrv, tolrp, tolrelv, tolrelp, max_iter, maxincr, p_flag=0, n_type=2):
    """
    Create a PFEM test, which check both increments and residual forvelocities and pressures.

    Parameters
    ----------
    osi: o3seespy.OpenSeesInstance
    tolv: float
            Tolerance for velocity increments
    tolp: float
            Tolerance for pressure increments
    tolrv: float
            Tolerance for relative velocity increments
    tolrp: float
            Tolerance for relative pressure increments
    tolrelv: None
            
    tolrelp: None
            
    max_iter: int
            Max number of iterations to check
    maxincr: int
            Max times for error increasing
    p_flag: int, optional
            Print flag : * 0 print nothing. * 1 print information on norms each time ``test()`` is invoked. * 2 print
            information on norms and number of iterations at end of successful test. * 4 at each step it will print the norms
            and also the :math:`\delta u` and :math:`r(u)` vectors. * 5 if it fails to converge at end of ``numiter`` it
            will print an error message **but return a successfull test**.
    n_type: int, optional
            Type of norm, (0 = max-norm, 1 = 1-norm, 2 = 2-norm). 
    """
    tolv = float(tolv)
    tolp = float(tolp)
    tolrv = float(tolrv)
    tolrp = float(tolrp)
    max_iter = int(max_iter)
    maxincr = int(maxincr)
    p_flag = int(p_flag)
    n_type = int(n_type)
    _parameters = [tolv, tolp, tolrv, tolrp, tolrelv, tolrelp, max_iter, maxincr, p_flag, n_type]
    return osi.to_process("test", _parameters)


def analysis(osi, dtmax, dtmin, gravity, ratio=0.5):
    """
    Create a OpenSees PFEMAnalysis object. 

    Parameters
    ----------
    osi: o3seespy.OpenSeesInstance
    dtmax: float
            Maximum time steps.
    dtmin: float
            Mimimum time steps.
    gravity: float
            Gravity acceleration used to move isolated particles.
    ratio: float, optional
            The ratio to reduce time steps if it was not converged. 
    """
    dtmax = float(dtmax)
    dtmin = float(dtmin)
    gravity = float(gravity)
    ratio = float(ratio)
    _parameters = [dtmax, dtmin, gravity, ratio]
    return osi.to_process("analysis", _parameters)
