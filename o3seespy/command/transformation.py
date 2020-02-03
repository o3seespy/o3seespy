from o3seespy.base_model import OpenSeesObject


class TransformationBase(OpenSeesObject):
    op_base_type = "geomTransf"


class Linear(TransformationBase):
    op_type = 'Linear'

    def __init__(self, osi, vecxz, d_i=None, d_j=None):
        raise ValueError('Use geom_transf.Linear not transformation.Linear')
        self.vecxz = vecxz
        self.d_i = d_i
        self.d_j = d_j

        osi.n_transformation += 1
        self._tag = osi.n_transformation
        self._parameters = [self.op_type, self._tag, *self.vecxz]
        if d_i is not None:
            self._parameters += ['jntOffset', self.d_i, self.d_j]
        self.to_process(osi)

