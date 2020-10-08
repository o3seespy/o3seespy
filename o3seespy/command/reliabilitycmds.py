
def random_variable(osi, dist, mean: float=None, stdv: float=None, start_point: float=None, params: list=None):
    """
    Create a random variable with user specified distribution

    Parameters
    ----------
    osi: o3seespy.OpenSeesInstance
    dist: str
            Random variable distribution * ``'normal'`` * ``'lognormal'`` * ``'gamma'`` * ``'shiftedexponential'`` *
            ``'shiftedrayleigh'`` * ``'exponential'`` * ``'rayleigh'`` * ``'uniform'`` * ``'beta'`` * ``'type1largestvalue'`` *
            ``'type1smallestvalue'`` * ``'type2largestvalue'`` * ``'type3smallestvalue'`` * ``'chisquare'`` * ``'gumbel'`` *
            ``'weibull'`` * ``'laplace'`` * ``'pareto'``
    mean: float, optional
            Mean value
    stdv: float, optional
            Standard deviation
    start_point: float, optional
            Starting point of the distribution
    params: list, optional
            A list of parameter objects
    """
    if mean is not None:
        mean = float(mean)
    if stdv is not None:
        stdv = float(stdv)
    if start_point is not None:
        start_point = float(start_point)
    _parameters = [dist]
    if mean is not None:
        _parameters += ['-mean', mean]
    if stdv is not None:
        _parameters += ['-stdv', stdv]
    if start_point is not None:
        _parameters += ['-startPoint', start_point]
    if params is not None:
        _parameters += ['-parameters', *params]
    return osi.to_process("randomVariable", _parameters)
