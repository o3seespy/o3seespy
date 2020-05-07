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
        osi: o3seespy.OpenSeesInstance
        n_d_mat: obj
            The object of the associated ndmaterial object
        n_dim: int
            Number of dimensions (2 for 2d, 3 for 3d)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> mat = o3.nd_material.ElasticIsotropic(osi, e_mod=1.0, nu=1.0, rho=0.0)
        >>> o3.nd_material.InitialStateAnalysisWrapper(osi, n_d_mat=mat, n_dim=1)
        """
        self.osi = osi
        self.n_d_mat = n_d_mat
        self.n_dim = int(n_dim)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.n_d_mat.tag, self.n_dim]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def set_update_material_stage(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'updateMaterialStage', value, ele, eles)

    def set_g_mod(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'shearModulus', value, ele, eles)

    def set_bulk_mod(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'bulkModulus', value, ele, eles)

    def set_friction_angle(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'frictionAngle', value, ele, eles)

    def set_cohesion(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'cohesion', value, ele, eles)


class InitStressNDMaterial(NDMaterialBase):
    """
    The InitStressNDMaterial NDMaterial Class
    
    This command is used to construct an Initial Stress material object.The stress-strain behaviour for this material is
    defined by another material.Initial Stress Material enables definition of initial stress for the material under
    consideration.The strain that corresponds to the initial stress will be calculated from the other material.
    """
    op_type = 'InitStressNDMaterial'

    def __init__(self, osi, other, init_stress, n_dim):
        """
        Initial method for InitStressNDMaterial

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        other: obj
            Object of the other material
        init_stress: float
            Initial stress
        n_dim: int
            Number of dimensions (e.g. if plane strain ndim=2)
        """
        self.osi = osi
        self.other = other
        self.init_stress = float(init_stress)
        self.n_dim = int(n_dim)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.other.tag, self.init_stress, self.n_dim]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class InitStrainNDMaterial(NDMaterialBase):
    """
    The InitStrainNDMaterial NDMaterial Class
    
    This command is used to construct an Initial Strain material object. The stress-strain behaviour for this material
    is defined by another material. Initial Strain Material enables definition of initial strains for the material under
    consideration. The stress that corresponds to the initial strain will be calculated from the other material.
    """
    op_type = 'InitStrainNDMaterial'

    def __init__(self, osi, other, init_strain, n_dim):
        """
        Initial method for InitStrainNDMaterial

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        other: obj
            Object of the other material
        init_strain: float
            Initial strain
        n_dim: float
            Number of dimensions
        """
        self.osi = osi
        self.other = other
        self.init_strain = float(init_strain)
        self.n_dim = float(n_dim)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.other.tag, self.init_strain, self.n_dim]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)
