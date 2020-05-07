from o3seespy.command.nd_material.base_material import NDMaterialBase


class CycLiqCP(NDMaterialBase):
    """
    The CycLiqCP NDMaterial Class
    
    This command is used to construct a multi-dimensional material object that that follows the constitutive behavior of
    a cyclic elastoplasticity model for large post- liquefaction deformation.CycLiqCP material is a cyclic elastoplasticity
    model for large post-liquefaction deformation, and is implemented using a cutting plane algorithm. The model is
    capable of reproducing small to large deformation in the pre- to post-liquefaction regime. The elastic moduli
    of the model are pressure dependent. The plasticity in the model is developed within the framework of
    bounding surface plasticity, with special consideration to the formulation of reversible and
    irreversible dilatancy.The model does not take into consideration of the state of sand, and
    requires different parameters for sand under different densities and confining pressures.
    The surfaces (i.e. failure and maximum pre-stress) are considered as circles in the pi
    plane.The model has been validated against VELACS centrifuge model tests and has used
    on numerous simulations of liquefaction related problems.When this material is
    employed in regular solid elements (e.g., FourNodeQuad, Brick), it simulates
    drained soil response. When solid-fluid coupled elements (u-p elements and
    SSP u-p elements) are used, the model is able to simulate undrained and partially drained behavior of soil.
    """
    op_type = 'CycLiqCP'

    def __init__(self, osi, g0, kappa, h, mfc, dre1, mdc, dre2, rdr, alpha, dir, ein, rho):
        """
        Initial method for CycLiqCP

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        g0: float
            A constant related to elastic shear modulus
        kappa: float
            Bulk modulus
        h: float
            Model parameter for plastic modulus
        mfc: float
            Stress ratio at failure in triaxial compression
        dre1: float
            Coefficient for reversible dilatancy generation
        mdc: float
            Stress ratio at which the reversible dilatancy sign changes
        dre2: float
            Coefficient for reversible dilatancy release
        rdr: float
            Reference shear strain length
        alpha: float
            Parameter controlling the decrease rate of irreversible dilatancy
        dir: float
            Coefficient for irreversible dilatancy potential
        ein: float
            Initial void ratio
        rho: float
            Saturated mass density

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.nd_material.CycLiqCP(osi, g0=1.0, kappa=1.0, h=1.0, mfc=1.0, dre1=1.0, mdc=1.0, dre2=1.0, rdr=1.0, alpha=1.0, dir=1.0, ein=1.0, rho=1.0)
        """
        self.osi = osi
        self.g0 = float(g0)
        self.kappa = float(kappa)
        self.h = float(h)
        self.mfc = float(mfc)
        self.dre1 = float(dre1)
        self.mdc = float(mdc)
        self.dre2 = float(dre2)
        self.rdr = float(rdr)
        self.alpha = float(alpha)
        self.dir = float(dir)
        self.ein = float(ein)
        self.rho = float(rho)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.g0, self.kappa, self.h, self.mfc, self.dre1, self.mdc, self.dre2, self.rdr, self.alpha, self.dir, self.ein, self.rho]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def set_update_material_stage(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'updateMaterialStage', value, ele, eles)


class CycLiqCPSP(NDMaterialBase):
    """
    The CycLiqCPSP NDMaterial Class
    
    This command is used to construct a multi-dimensional material object that that follows the constitutive behavior of
    a cyclic elastoplasticity model for large post- liquefaction deformation.CycLiqCPSP material is a constitutive model
    for sand with special considerations for cyclic behaviour and accumulation of large post-liquefaction shear
    deformation, and is implemented using a cutting plane algorithm. The model: (1) achieves the simulation of
    post-liquefaction shear deformation based on its physics, allowing the unified description of pre- and
    post-liquefaction behavior of sand; (2) directly links the cyclic mobility of sand with reversible
    and irreversible dilatancy, enabling the unified description of monotonic and cyclic loading; (3)
    introduces critical state soil mechanics concepts to achieve unified modelling of sand under
    different states.The critical, maximum stress ratio and reversible dilatancy surfaces
    follow a rounded triangle in the pi plane similar to the Matsuoka-Nakai
    criterion.When this material is employed in regular solid elements
    (e.g., FourNodeQuad, Brick), it simulates drained soil response.
    When solid-fluid coupled elements (u-p elements and SSP u-p
    elements) are used, the model is able to simulate undrained and partially drained behavior of soil.
    """
    op_type = 'CycLiqCPSP'

    def __init__(self, osi, g0, kappa, h, big_m, dre1, dre2, rdr, alpha, dir, lambdac, ksi, e0, np, nd, ein, rho):
        """
        Initial method for CycLiqCPSP

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        g0: float
            A constant related to elastic shear modulus
        kappa: float
            Bulk modulus
        h: float
            Model parameter for plastic modulus
        big_m: float
            Critical state stress ratio
        dre1: float
            Coefficient for reversible dilatancy generation
        dre2: float
            Coefficient for reversible dilatancy release
        rdr: float
            Reference shear strain length
        alpha: float
            Parameter controlling the decrease rate of irreversible dilatancy
        dir: float
            Coefficient for irreversible dilatancy potential
        lambdac: float
            Critical state constant
        ksi: float
            Critical state constant
        e0: float
            Void ratio at pc=0
        np: float
            Material constant for peak mobilized stress ratio
        nd: float
            Material constant for reversible dilatancy generation stress ratio
        ein: float
            Initial void ratio
        rho: float
            Saturated mass density

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.nd_material.CycLiqCPSP(osi, g0=1.0, kappa=1.0, h=1.0, big_m=1.0, dre1=1.0, dre2=1.0, rdr=1.0, alpha=1.0, dir=1.0, lambdac=1.0, ksi=1.0, e0=1.0, np=1.0, nd=1.0, ein=1.0, rho=1.0)
        """
        self.osi = osi
        self.g0 = float(g0)
        self.kappa = float(kappa)
        self.h = float(h)
        self.big_m = float(big_m)
        self.dre1 = float(dre1)
        self.dre2 = float(dre2)
        self.rdr = float(rdr)
        self.alpha = float(alpha)
        self.dir = float(dir)
        self.lambdac = float(lambdac)
        self.ksi = float(ksi)
        self.e0 = float(e0)
        self.np = float(np)
        self.nd = float(nd)
        self.ein = float(ein)
        self.rho = float(rho)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.g0, self.kappa, self.h, self.big_m, self.dre1, self.dre2, self.rdr, self.alpha, self.dir, self.lambdac, self.ksi, self.e0, self.np, self.nd, self.ein, self.rho]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def set_update_material_stage(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'updateMaterialStage', value, ele, eles)
