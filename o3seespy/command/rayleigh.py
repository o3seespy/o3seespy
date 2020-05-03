from o3seespy.base_model import OpenSeesObject


class Rayleigh(OpenSeesObject):
    op_base_type = "rayleigh"
    op_type = "rayleigh"

    def __init__(self, osi, alpha_m, beta_k, beta_k_init, beta_k_comm):
        """
        Assign Rayleigh damping to previously defined nodes

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        alpha_m : float
            factor applied to element or nodes mass matrix
        beta_k : float
            factor applied to element current stiffness matrix
        beta_k_init : float
            factor applied to element initial stiffness matrix
        beta_k_comm : float
            factor applied to element committed stiffness matrix
        """
        self.osi = osi
        self.alpha_m = float(alpha_m)
        self.beta_k = float(beta_k)
        self.beta_k_init = float(beta_k_init)
        self.beta_k_comm = float(beta_k_comm)
        self._parameters = [self.alpha_m, self.beta_k, self.beta_k_init, self.beta_k_comm]
        self.to_process(osi)

