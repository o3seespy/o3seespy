from o3seespy.command.nd_material.base_material import NDMaterialBase


class InitialStateAnalysisWrapper(NDMaterialBase):
    """
    The InitialStateAnalysisWrapper NDMaterial Class
    
    The InitialStateAnalysisWrapper nDMaterial allows for the use of the InitialStateAnalysis command for setting
    initial conditions. The InitialStateAnalysisWrapper can be used with any nDMaterial. This material wrapper
    allows for the development of an initial stress field while maintaining the original geometry of the
    problem. An example analysis is provided below to demonstrate the use of this material wrapper object.
    """
    op_type = 'InitialStateAnalysisWrapper'

    def __init__(self, osi, n_d_mat, n_dim):
        """
        Initial method for InitialStateAnalysisWrapper

        Parameters
        ----------
        n_d_mat: obj
            The tag of the associated ndmaterial object
        n_dim: int
            Number of dimensions (2 for 2d, 3 for 3d)
        """
        self.n_d_mat = n_d_mat
        self.n_dim = int(n_dim)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.n_d_mat.tag, self.n_dim]
        self.to_process(osi)
