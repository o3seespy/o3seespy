from o3seespy.command.nd_material.base_material import NDMaterialBase


class ElasticIsotropic(NDMaterialBase):
    type = "ElasticIsotropic"

    def __init__(self, osi, e_mod, v, rho=None):
        """
        ElasticIsotropic material
        """
        self.e_mod = e_mod
        self.v = v
        self.rho = rho
        osi.n_mats += 1
        self._tag = osi.n_mats
        self._parameters = [self.op_type, self._tag, self.e_mod, self.v]
        if self.rho is not None:
            self._parameters.append(self.rho)
        self.to_process(osi)

