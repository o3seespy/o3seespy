from o3seespy.command.nd_material.base_material import NDMaterialBase


class ElasticIsotropic(NDMaterialBase):
    op_type = "ElasticIsotropic"

    def __init__(self, osi, e_mod, nu, rho=None):
        """
        ElasticIsotropic material
        """
        self.e_mod = e_mod
        self.nu = nu
        self.rho = rho
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.e_mod, self.nu]
        if self.rho is not None:
            self._parameters.append(self.rho)
        self.to_process(osi)

