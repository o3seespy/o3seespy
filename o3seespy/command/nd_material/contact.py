from o3seespy.command.nd_material.base_material import NDMaterialBase


class ContactMaterial2D(NDMaterialBase):
    """
    The ContactMaterial2D NDMaterial Class
    
    This command is used to construct a ContactMaterial2D nDMaterial object.
    """
    op_type = 'ContactMaterial2D'

    def __init__(self, osi, mu, g_mod, c, t):
        """
        Initial method for ContactMaterial2D

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        mu: float
            Interface frictional coefficient
        g_mod: float
            Interface stiffness parameter
        c: float
            Interface cohesive intercept
        t: float
            Interface tensile strength

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.nd_material.ContactMaterial2D(osi, mu=1.0, g_mod=1.0, c=1.0, t=1.0)
        """
        self.osi = osi
        self.mu = float(mu)
        self.g_mod = float(g_mod)
        self.c = float(c)
        self.t = float(t)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.mu, self.g_mod, self.c, self.t]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class ContactMaterial3D(NDMaterialBase):
    """
    The ContactMaterial3D NDMaterial Class
    
    This command is used to construct a ContactMaterial3D nDMaterial object.
    """
    op_type = 'ContactMaterial3D'

    def __init__(self, osi, mu, g_mod, c, t):
        """
        Initial method for ContactMaterial3D

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        mu: float
            Interface frictional coefficient
        g_mod: float
            Interface stiffness parameter
        c: float
            Interface cohesive intercept
        t: float
            Interface tensile strength

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=3)
        >>> o3.nd_material.ContactMaterial3D(osi, mu=1.0, g_mod=1.0, c=1.0, t=1.0)
        """
        self.osi = osi
        self.mu = float(mu)
        self.g_mod = float(g_mod)
        self.c = float(c)
        self.t = float(t)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.mu, self.g_mod, self.c, self.t]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)
