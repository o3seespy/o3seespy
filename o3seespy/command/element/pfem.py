from o3seespy.command.element.base_element import ElementBase


class PFEMElementBubble(ElementBase):
    """
    The PFEMElementBubble Element Class
    
    Create a PFEM Bubble element, which is a fluid element for FSI analysis.
    """
    op_type = 'PFEMElementBubble'

    def __init__(self, osi, ele_nodes, rho, mu, b1, b2, b3, thickness, kappa):
        """
        Initial method for PFEMElementBubble

        Parameters
        ----------
        ele_nodes: list
            A list of three or four element nodes, four are required for 3d
        rho: float
            Fluid density
        mu: float
            Fluid viscosity
        b1: float
            Body body acceleration in x direction
        b2: float
            Body body acceleration in y direction
        b3: float
            Body body acceleration in z direction (required for 3d)
        thickness: float
            Element thickness (required for 2d)
        kappa: float
            Fluid bulk modulus (optional)
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.rho = float(rho)
        self.mu = float(mu)
        self.b1 = float(b1)
        self.b2 = float(b2)
        self.b3 = float(b3)
        self.thickness = float(thickness)
        self.kappa = float(kappa)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.rho, self.mu, self.b1, self.b2, self.b3, self.thickness, self.kappa]
        self.to_process(osi)


class PFEMElementCompressible(ElementBase):
    """
    The PFEMElementCompressible Element Class
    
    Create a PFEM compressible element, which is a fluid element for FSI analysis.
    """
    op_type = 'PFEMElementCompressible'

    def __init__(self, osi, ele_nodes, rho, mu, b1, b2, thickness, kappa):
        """
        Initial method for PFEMElementCompressible

        Parameters
        ----------
        ele_nodes: list
            A list of four element nodes, last one is middle node
        rho: float
            Fluid density
        mu: float
            Fluid viscosity
        b1: float
            Body body acceleration in x direction
        b2: float
            Body body acceleration in y direction
        thickness: float
            Element thickness (optional)
        kappa: float
            Fluid bulk modulus (optional)
        """
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.rho = float(rho)
        self.mu = float(mu)
        self.b1 = float(b1)
        self.b2 = float(b2)
        self.thickness = float(thickness)
        self.kappa = float(kappa)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.rho, self.mu, self.b1, self.b2, self.thickness, self.kappa]
        self.to_process(osi)
