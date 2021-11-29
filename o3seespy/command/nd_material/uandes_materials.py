from o3seespy.command.nd_material.base_material import NDMaterialBase


class SAniSandMS(NDMaterialBase):
    """
    The SAniSandMS NDMaterial Class

    This command is used to construct a multi-dimensional SAniSandMS material.
    """
    op_type = 'SAniSandMS'

    def __init__(self, osi, g0, nu, e_init, m_c, c_c, lambda_c, e_0, ksi, p_atm, m_yield, h_0, c_h, n_b, a_0, n_d,
                 zeta, mu0, beta, den, int_scheme, tan_type=0, tol_f=1.0e-7, tol_r=1.0e-7):
        r"""
        Initial method for SAniSandMS

        Note: When in elastic mode the shear modulus is equal to the shear modulus at the atmospheric pressure

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        g0: float
            Shear modulus constant
        nu: float
            Poisson ratio
        e_init: float
            Initial void ratio
        m_c: float
            Critical state stress ratio
        c_c: float
            Ratio of critical state stress ratio in extension and compression
        lambda_c: float
            Critical state line constant
        e_0: float
            Critical void ratio at p = 0
        ksi: float
            Critical state line constant
        p_atm: float
            Atmospheric pressure
        m_yield: float
            Yield surface constant (radius of yield surface in stress ratio space)
        h_0: float
            Constant parameter
        c_h: float
            Constant parameter
        n_b: float
            Bounding surface parameter, :math:`n_b \ge 0`
        a_0: float
            Dilatancy parameter
        n_d: float
            Dilatancy surface parameter :math:`n_d \ge 0`
        zeta: float
            Memory surface shrinkage parameter
        mu0: float
            Ratcheting parameter
        den: float
            Mass density of the material
        int_scheme: int, optional (default=1)
            Integration scheme type:
                * 1 = Modified Euler
                * 3 = Runge Kutta 4th order with error control (Reccommended)
        tan_type: int, optional (default=0)
            Tangent type: (appears to be unused - possibly only used in Implicit scheme)
                * 0: Elastic Tangent
                * 1: Contiuum ElastoPlastic Tangent
                * 2: Consistent ElastoPlastic Tangent
        jaco_type: int, optional (default=1)
            Jacobian type: (possibly only used in Implicit scheme)
                * 0: Finite Difference Jacobian
                * 1: Analytical Jacobian
        tol_f: float, optional (default=1.0e-7)
        tol_r: float, optional (default=1.0e-7)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.nd_material.SAniSandMS(osi, g0=1110.0, nu=0.05, e_init=0.72, m_c=1.27, c_c=0.712, lambda_c=0.049, e_0=0.845,
        >>> ksi=0.27, p_atm=101.3, m_yield=0.01, h_0=5.95, c_h=1.01, n_b=2.0, a_0=1.06, n_d=1.17, zeta=0.0005, mu0=260.0, beta=1, den=1.6)
        """
        self.osi = osi
        self.g0 = float(g0)
        self.nu = float(nu)
        self.e_init = float(e_init)
        self.m_c = float(m_c)
        self.c_c = float(c_c)
        self.lambda_c = float(lambda_c)
        self.e_0 = float(e_0)
        self.ksi = float(ksi)
        self.p_atm = float(p_atm)
        self.m_yield = float(m_yield)
        self.h_0 = float(h_0)
        self.c_h = float(c_h)
        self.n_b = float(n_b)
        self.a_0 = float(a_0)
        self.n_d = float(n_d)
        self.zeta = float(zeta)  # note DM04 has z_max, c_z
        self.mu0 = float(mu0)
        self.beta = beta
        self.den = float(den)
        self.int_scheme = int(int_scheme)
        self.tan_type = int(tan_type)
        dep_jaco_type = 1
        self.tol_f = float(tol_f)
        self.tol_r = float(tol_r)
        dep_fabric_flag = 1
        dep_flow_flag = 1
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.g0, self.nu, self.e_init, self.m_c, self.c_c, self.lambda_c,
                            self.e_0, self.ksi, self.p_atm, self.m_yield, self.h_0, self.c_h, self.n_b, self.a_0,
                            self.n_d, self.zeta, self.mu0, self.beta, self.den, dep_fabric_flag, dep_flow_flag,
                            self.int_scheme, self.tan_type, dep_jaco_type,
                            self.tol_f, self.tol_r]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def set_update_material_stage(self, stage):
        parameters = ['-material', self.tag, '-stage', stage]
        self.osi.to_process("updateMaterialStage", parameters)

    def set_material_state(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'materialState', value, ele, eles)

    def set_integration_scheme(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'IntegrationScheme', value, ele, eles)

    def set_ref_shear_modulus(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'refShearModulus', value, ele, eles)
    #
    # def set_g_mod(self, value, ele=None, eles=None):
    #     self.set_parameter(self.osi, 'ShearModulus', value, ele, eles)

    def set_nu(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'poissonRatio', value, ele, eles)

    def set_void_ratio(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'voidRatio', value, ele, eles)

    def update_to_nonlinear(self):
        from o3seespy import update_material_stage
        update_material_stage(self.osi, self, 1)

    def update_to_linear(self):
        from o3seespy import update_material_stage
        update_material_stage(self.osi, self, 0)
