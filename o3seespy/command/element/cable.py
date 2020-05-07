from o3seespy.command.element.base_element import ElementBase


class CatenaryCable(ElementBase):
    """
    The CatenaryCable Element Class
    
    This command is used to construct a catenary cable element object.

    
    """
    op_type = 'CatenaryCable'

    def __init__(self, osi, i_node, j_node, weight, big_e, big_a, l0, alpha, temperature_change, rho, error_tol, nsubsteps, mass_type):
        """
        Initial method for CatenaryCable

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        i_node: obj
            End nodes (3 dof per node)
        j_node: obj
            End nodes (3 dof per node)
        weight: float
            Undefined
        big_e: float
            Elastic modulus of the cable material
        big_a: float
            Cross-sectional area of element
        l0: float
            Unstretched length of the cable
        alpha: float
            Coefficient of thermal expansion
        temperature_change: float
            Temperature change for the element
        rho: float
            Mass per unit length
        error_tol: float
            Allowed tolerance for within-element equilbrium (newton-rhapson iterations)
        nsubsteps: int
            Number of within-element substeps into which equilibrium iterations are subdivided (not number of steps to
            convergence)
        mass_type: int
            Mass matrix model to use (``masstype`` = 0 lumped mass matrix,    ``masstype`` = 1 rigid-body mass matrix
            (in development))

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> i_node = o3.node.Node(osi, 0.0, 0.0)
        >>> j_node = o3.node.Node(osi, 0.0, 1.0)
        >>> o3.element.CatenaryCable(osi, i_node=i_node, j_node=j_node, weight=1.0, big_e=1.0, big_a=1.0, l0=1.0, alpha=1.0, temperature_change=1.0, rho=1.0, error_tol=1.0, nsubsteps=1, mass_type=1)
        """
        self.osi = osi
        self.i_node = i_node
        self.j_node = j_node
        self.weight = float(weight)
        self.big_e = float(big_e)
        self.big_a = float(big_a)
        self.l0 = float(l0)
        self.alpha = float(alpha)
        self.temperature_change = float(temperature_change)
        self.rho = float(rho)
        self.error_tol = float(error_tol)
        self.nsubsteps = int(nsubsteps)
        self.mass_type = int(mass_type)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.i_node.tag, self.j_node.tag, self.weight, self.big_e, self.big_a, self.l0, self.alpha, self.temperature_change, self.rho, self.error_tol, self.nsubsteps, self.mass_type]
        self.to_process(osi)
