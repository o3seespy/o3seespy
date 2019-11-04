from o3seespy.command.uniaxial_material.base_material import UniaxialMaterialBase


class Concrete01(UniaxialMaterialBase):
    op_type = "Concrete01"
    
    def __init__(self, osi, fpc, epsc0, fpcu, eps_ult):
        """
        A uni-axial Kent-Scott-Park concrete material object

        Parameters
        ----------
        osi : opensees_pack.opensees_instance.OpenseesInstance object
            An instance of opensees
        fpc : float
            Concrete compressive strength at 28 days (compression is negative)
        epsc0 : float
            Concrete strain at maximum strength
        fpcu : float
            Concrete crushing strength
        eps_ult : float
            Concrete strain at crushing strength
        """
        self.fpc = fpc
        self.epsc0 = epsc0
        self.fpcu = fpcu
        self.eps_ult = eps_ult
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.fpc, self.epsc0, self.fpcu, self.eps_ult]
        self.to_process(osi)
