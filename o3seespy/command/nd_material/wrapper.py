from o3seespy.command.nd_material.base_material import NDMaterialBase


class InitialStateAnalysisWrapper(NDMaterialBase):
    op_type = 'InitialStateAnalysisWrapper'

    def __init__(self, osi, n_d_mat, n_dim):
        self.n_d_mat = n_d_mat
        self.n_dim = int(n_dim)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.n_d_mat.tag, self.n_dim]
        self.to_process(osi)
