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

    def __init__(self, osi, g0, kappa, h, mfc, dre1, mdc, dre2, gamma_dr, alpha, d_ir, e_init, den):
        """
        Initial method for CycLiqCP

        Notes:
         - Units are in kPa. E.g. atmospheric pressure is hard coded as 101 kPa.

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
        gamma_dr: float
            Reference shear strain length
        alpha: float
            Parameter controlling the decrease rate of irreversible dilatancy
        d_ir: float
            Coefficient for irreversible dilatancy potential
        e_init: float
            Initial void ratio
        den: float
            Saturated mass density

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.nd_material.CycLiqCP(osi, g0=1.0, kappa=1.0, h=1.0, mfc=1.0, dre1=1.0, mdc=1.0, dre2=1.0, gamma_dr=1.0, alpha=1.0, d_ir=1.0, e_init=1.0, den=1.0)
        """
        self.osi = osi
        self.g0 = float(g0)
        self.kappa = float(kappa)
        self.h = float(h)
        self.mfc = float(mfc)
        self.dre1 = float(dre1)
        self.mdc = float(mdc)
        self.dre2 = float(dre2)
        self.gamma_dr = float(gamma_dr)
        self.alpha = float(alpha)
        self.d_ir = float(d_ir)
        self.e_init = float(e_init)
        self.den = float(den)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.g0, self.kappa, self.h, self.mfc, self.dre1, self.mdc,
                            self.dre2, self.gamma_dr, self.alpha, self.d_ir, self.e_init, self.den]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def set_update_material_stage(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'updateMaterialStage', value, ele, eles)

    def update_to_nonlinear(self):
        from o3seespy import update_material_stage
        update_material_stage(self.osi, self, 1)

    def update_to_linear(self):
        from o3seespy import update_material_stage
        update_material_stage(self.osi, self, 2)

    def update_to_nonlinear_elastic(self):
        from o3seespy import update_material_stage
        update_material_stage(self.osi, self, 0)


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

    def __init__(self, osi, g0, kappa, h, big_m, dre1, dre2, gamma_dr, alpha, d_ir, lambda_c, ksi, e_0, n_p, n_d, e_init, den):
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
        gamma_dr: float
            Reference shear strain length (gamma_d,r)
        alpha: float
            Parameter controlling the decrease rate of irreversible dilatancy
        d_ir: float
            Coefficient for irreversible dilatancy potential
        lambda_c: float
            Critical state constant
        ksi: float
            Critical state constant
        e_0: float
            Void ratio at pc=0
        n_p: float
            Material constant for peak mobilized stress ratio
        n_d: float
            Material constant for reversible dilatancy generation stress ratio
        e_init: float
            Initial void ratio
        den: float
            Saturated mass density

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.nd_material.CycLiqCPSP(osi, g0=200, kappa=0.008, h=1.8, big_m=1.25, dre1=0.6, dre2=30, gamma_dr=0.05,
        >>> alpha=20, d_ir=1.4, lambda_c=0.019, ksi=0.7, e_0=0.934, n_p=1.1, n_d=7.8, e_init=0.87, den=1.6)
        """
        self.osi = osi
        self.g0 = float(g0)
        self.kappa = float(kappa)
        self.h = float(h)
        self.big_m = float(big_m)
        self.dre1 = float(dre1)
        self.dre2 = float(dre2)
        self.gamma_dr = float(gamma_dr)
        self.alpha = float(alpha)
        self.d_ir = float(d_ir)
        self.lambda_c = float(lambda_c)
        self.ksi = float(ksi)
        self.e_0 = float(e_0)
        self.n_p = float(n_p)
        self.n_d = float(n_d)
        self.e_init = float(e_init)
        self.den = float(den)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.g0, self.kappa, self.h, self.big_m, self.dre1, self.dre2,
                            self.gamma_dr, self.alpha, self.d_ir, self.lambda_c, self.ksi, self.e_0, self.n_p,
                            self.n_d, self.e_init, self.den]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def set_update_material_stage(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'updateMaterialStage', value, ele, eles)

    def update_to_linear(self):
        from o3seespy import update_material_stage
        update_material_stage(self.osi, self, 0)

    def update_to_nonlinear(self):
        from o3seespy import update_material_stage
        update_material_stage(self.osi, self, 1)
