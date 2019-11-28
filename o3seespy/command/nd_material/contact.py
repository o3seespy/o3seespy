from o3seespy.command.nd_material.base_material import NDMaterialBase


class ContactMaterial2D(NDMaterialBase):
    """
    The ContactMaterial2D NDMaterial Class
    
    This command is used to construct a ContactMaterial2D nDMaterial object.
    """
    op_type = 'ContactMaterial2D'

    def __init__(self, osi, mu, big_g, c, t):
        """
        Initial method for ContactMaterial2D

        Parameters
        ----------
        mu: float
            Interface frictional coefficient
        big_g: float
            Interface stiffness parameter
        c: float
            Interface cohesive intercept
        t: float
            Interface tensile strength
        """
        self.mu = float(mu)
        self.big_g = float(big_g)
        self.c = float(c)
        self.t = float(t)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.mu, self.big_g, self.c, self.t]
        self.to_process(osi)


class ContactMaterial3D(NDMaterialBase):
    """
    The ContactMaterial3D NDMaterial Class
    
    This command is used to construct a ContactMaterial3D nDMaterial object.
    """
    op_type = 'ContactMaterial3D'

    def __init__(self, osi, mu, big_g, c, t):
        """
        Initial method for ContactMaterial3D

        Parameters
        ----------
        mu: float
            Interface frictional coefficient
        big_g: float
            Interface stiffness parameter
        c: float
            Interface cohesive intercept
        t: float
            Interface tensile strength
        """
        self.mu = float(mu)
        self.big_g = float(big_g)
        self.c = float(c)
        self.t = float(t)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.mu, self.big_g, self.c, self.t]
        self.to_process(osi)
