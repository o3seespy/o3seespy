from o3seespy.base_model import OpenSeesObject


class AlgorithmBase(OpenSeesObject):
    op_base_type = "algorithm"


class Linear(AlgorithmBase):
    op_type = "Linear"

    def __init__(self, osi, secant=False, initial=False, factor_once=False):
        self.osi = osi
        self.secant = secant
        self.initial = initial
        self.factor_once = factor_once
        self._parameters = [self.op_type, self.secant, self.initial, self.factor_once]
        self.to_process(osi)


class Newton(AlgorithmBase):
    op_type = "Newton"

    def __init__(self, osi, secant=False, initial=False, initial_then_current=False):
        self.osi = osi
        self.secant = secant
        self.initial = initial
        self.initial_then_current = initial_then_current
        self._parameters = [self.op_type, self.secant, self.initial, self.initial_then_current]
        self.to_process(osi)


class SecantNewton(AlgorithmBase):
    op_type = "SecantNewton"

    def __init__(self, osi, iterate='current', increment='current', maxDim=3):
        self.osi = osi
        self.iterate = iterate
        self.increment = increment
        self.maxDim = maxDim
        self._parameters = [self.op_type, self.iterate, self.increment, self.maxDim]
        self.to_process(osi)


class RaphsonNewton(AlgorithmBase):
    op_type = "RaphsonNewton"

    def __init__(self, osi, iterate='current', increment='current'):
        self.osi = osi
        self.iterate = iterate
        self.increment = increment
        self._parameters = [self.op_type, self.iterate, self.increment]
        self.to_process(osi)


class KrylovNewton(AlgorithmBase):
    op_type = "KrylovNewton"

    def __init__(self, osi, tang_inter='current', tang_incr='current', max_inter=3):
        """

        Parameters
        ----------
        osi
        tang_inter: str
            options are: 'current', 'initial', 'noTangent'
        tang_incr
        max_inter
        """
        self.osi = osi
        self.tang_inter = tang_inter
        self.tang_incr = tang_incr
        self.max_inter = max_inter
        self._parameters = [self.op_type, self.tang_inter, self.tang_incr, self.max_inter]
        self.to_process(osi)
