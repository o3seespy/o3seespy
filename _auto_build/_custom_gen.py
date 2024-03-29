class NDMaterialBase(object):
    pass


class UniaxialMaterialBase(object):
    pass


class ElementBase(object):
    pass


class SectionBase(object):
    pass


class LayerBase(object):
    pass


class LayeredShell(SectionBase):
    """
    The LayeredShell Section Class

    This command will create the section of the multi-layer shell element, including the multi-dimensional concrete,
    reinforcement material and the corresponding thickness.
    """
    op_type = 'LayeredShell'

    def __init__(self, osi, mats):
        """
        Initial method for LayeredShell

        Parameters
        ----------
        mats: list
            A list of material objs and thickness, ``[[mat1,thk1], ..., [mat2,thk2]]``
        """
        self.mats = []
        for i, mat in enumerate(mats):
            # self.mats.append([mats[i][0].tag, mats[i][1]])
            self.mats.append(mats[i][0].tag)
            self.mats.append(mats[i][1])

        self.n_layers = int(len(self.mats))
        osi.n_sect += 1
        self._tag = osi.n_sect
        self._parameters = [self.op_type, self._tag, self.n_layers, *self.mats]
        self.to_process(osi)


class Aggregator(SectionBase):
    """
    The Aggregator Section Class

    This command is used to construct a SectionAggregator object which aggregates groups previously-defined
    UniaxialMaterial objects into a single section force-deformation model. Each UniaxialMaterial object
    represents the section force-deformation response for a particular section degree-of-freedom (dof).
    There is no interaction between responses in different dof directions. The aggregation can include
    one previously defined section.
    """
    op_type = 'Aggregator'

    def __init__(self, osi, mats, section=None):
        """
        Initial method for Aggregator

        Parameters
        ----------
        mats: list
            List of mat objs and dofs of previously-defined uniaxialmaterial objects, ``mats =
            [[mattag1,dof1],[mattag2,dof2],...]`` the force-deformation quantity to be modeled by this
            section object. one of the following section dof may be used: * ``'p'`` axial
            force-deformation * ``'mz'`` moment-curvature about section local z-axis *
            ``'vy'`` shear force-deformation along section local y-axis * ``'my'``
            moment-curvature about section local y-axis * ``'vz'`` shear
            force-deformation along section local z-axis * ``'t'`` torsion force-deformation
        section: obj
            Tag of previously-defined section object to which the uniaxialmaterial objects are aggregated as additional
            force-deformation relationships (optional)
        """
        self.mats = []
        for i, mat in enumerate(mats):
            self.mats.append(mats[i][0].tag)
            self.mats.append(mats[i][1])
        self.section = section
        osi.n_sect += 1
        self._tag = osi.n_sect
        self._parameters = [self.op_type, self._tag, *self.mats]
        if getattr(self, 'section') is not None:
            self._parameters += ['-section', self.section.tag]
        self.to_process(osi)


class ManzariDafalias(NDMaterialBase):
    """
    The ManzariDafalias NDMaterial Class

    This command is used to construct a multi-dimensional Manzari-Dafalias(2004) material.
    """
    op_type = 'ManzariDafalias'

    def __init__(self, osi, g0, nu, e_init, m_c, c_c, lambda_c, e_0, ksi, p_atm, m_yield, h_0, c_h, n_b, a_0, n_d,
                 z_max, c_z, den, int_scheme=1, tan_type=0, jaco_type=1, tol_f=1.0e-7, tol_r=1.0e-7):
        r"""
        Initial method for ManzariDafalias

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
        z_max: float
            Fabric-dilatancy tensor parameter
        c_z: float
            Fabric-dilatancy tensor parameter
        den: float
            Mass density of the material
        int_scheme: int, optional (default=1)
            Integration scheme type:
                * 0 = Modified Euler constraining maximum energy increment (See Note 6)
                * 1 = Modified Euler with error control
                * 2 = Backward Euler (Implicit)
                * 3 = Runge Kutta 4th order (inconsistent results to other methods - 3&6 produce same)
                * 4 = Forward Euler constraining maximum energy increment (See Note 6)
                * 5 = Forward Euler
                * 6 = Runge-Kutta 4-th order constraining maximum energy increment (inconsistent results to other methods - 3&6 produce same)
                * 7 = Not implemented. runs int_scheme=9. Modified Euler constraining maximum strain increment
                * 8 = Not implemented. runs int_scheme=9. Runge-Kutta 4-th order constraining maximum strain increment
                * 9 = Forward Euler constraining maximum strain increment
                * 45 = Runge Kutta 45 with error control after Sloan [very slow]
            Notes:
                1. To use implicit integration, `int_scheme` must be equal to 2
                2. IS=3,6 - RK4 methods. produce same results under small time steps (diff to others)
                3. IS=0,1 - Modified Euler methods. produce same results under small time steps and same as Backward Euler.
                4. IS=4,5,9 - Forward Euler methods. produce same results as 0, 1 and 2 under small time steps
                5. The maximum strain increment is hardcoded as 1e-5
                6. WARNING: The maximum energy increment is hardcoded as 1e-4 (not normalised, should use p_atm=101 (i.e. kPa))
                   also only does 2 substeps at half strain increment, if exceeds energy increment.
                9. IS=9 is the preferred method for global explicit methods since has no substepping in Euler step,
                   which is fine since time step typically small enough, however, has an additional strain check which
                   may result in substepping.
        tan_type: int, optional (default=0)
            Tangent type:
                * 0: Elastic Tangent
                * 1: Contiuum ElastoPlastic Tangent
                * 2: Consistent ElastoPlastic Tangent
        jaco_type: int, optional (default=1)
            Jacobian type:
                * 0: Finite Difference Jacobian
                * 1: Analytical Jacobian
        tol_f: float, optional (default=1.0e-7)
            Tolerance for evaluating whether stress state outside yield surface
        tol_r: float, optional (default=1.0e-7)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.nd_material.ManzariDafalias(osi, g0=1.0, nu=1.0, e_init=1.0, m_c=1.0, c_c=1.0, lambda_c=1.0, e_0=1.0,
        >>> ksi=1.0, p_atm=1.0, m_yield=1.0, h_0=1.0, c_h=1.0, n_b=1.0, a_0=1.0, n_d=1.0, z_max=1.0, c_z=1.0, den=1.0)
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
        self.z_max = float(z_max)
        self.c_z = float(c_z)
        self.den = float(den)
        self.int_scheme = int(int_scheme)
        self.tan_type = int(tan_type)
        self.jaco_type = int(jaco_type)
        self.tol_f = float(tol_f)
        self.tol_r = float(tol_r)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.g0, self.nu, self.e_init, self.m_c, self.c_c, self.lambda_c,
                            self.e_0, self.ksi, self.p_atm, self.m_yield, self.h_0, self.c_h, self.n_b, self.a_0,
                            self.n_d, self.z_max, self.c_z, self.den, self.int_scheme, self.tan_type, self.jaco_type,
                            self.tol_f, self.tol_r]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def set_update_material_stage(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'updateMaterialStage', value, ele, eles)

    def set_material_state(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'materialState', value, ele, eles)

    def set_integration_scheme(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'IntegrationScheme', value, ele, eles)

    def set_jacobian(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'Jacobian', value, ele, eles)

    def set_ref_shear_modulus(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'refShearModulus', value, ele, eles)
    #
    # def set_g_mod(self, value, ele=None, eles=None):
    #     self.set_parameter(self.osi, 'ShearModulus', value, ele, eles)

    def set_nu(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'poissonRatio', value, ele, eles)

    def set_void_ratio(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'voidRatio', value, ele, eles)

    def set_stress_correction(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'stressCorrection', value, ele, eles)

    def update_to_nonlinear(self):
        from o3seespy import update_material_stage
        update_material_stage(self.osi, self, 1)

    def update_to_linear(self):
        from o3seespy import update_material_stage
        update_material_stage(self.osi, self, 0)


class PressureIndependMultiYield(NDMaterialBase):
    op_type = "PressureIndependMultiYield"

    def __init__(self, osi,  nd, rho, g_mod_ref, bulk_mod_ref, cohesion, peak_strain, phi=0.,
                 p_ref=100., d=0., n_surf=20, strains=None, ratios=None):
        r"""
        PressureIndependMultiYield material

        Parameters
        ==========
        matTag: int
            integer tag identifying material
        nd: float
            Number of dimensions, 2 for plane-strain, and 3 for 3D analysis.
       rho: float
            Saturated soil mass density.
       g_mod_ref, float
            (:math:`G_0`) Reference low-strain shear modulus,
            specified at a reference mean effective confining pressure (`p_ref`).
       bulk_mod_ref: float
            (:math:`B_r`) Reference bulk modulus,
            specified at a reference mean effective confining pressure (`p_ref`)).
       cohesion: float
            (:math:`c`) Apparent cohesion at zero effective confinement.
       peak_strain: float
            (:math:`\gamma_{max}`) An octahedral shear strain at
            which the maximum shear strength is reached,
            specified at a reference mean effective confining
            pressure refPress of p'r (see below).
       phi: float
            (:math:`phi`) Friction angle at peak shear
            strength in degrees, optional (default is 0.0).
       p_ref: float
            (:math:`p'_ref`) Reference mean effective confining pressure at which
                          :math:`G_r`, :math:`B_r`, and :math:`\gamma_{max}`
                          are defined, optional (default is 100. kPa).
       d: float
            (:math:`d`) A positive constant defining variations
                        of :math:`G` and :math:`B` as a function of
                          instantaneous effective
                          confinement :math:`p'` (default is 0.0)

                          :math:`G=G_r(\frac{p'}{p'_ref})^d`

                          :math:`B=B_r(\frac{p'}{p'_ref})^d`

                          If :math:`\phi=0`, :math:`d` is reset to 0.0.

       n_surf: float, optional
            Number of yield surfaces, optional (must be less
            than 40, default is 20). The surfaces are generated
            based on the hyperbolic relation.
       strains: list
            Instead of automatic surfaces generation,
            you can define yield surfaces directly based on
            desired shear modulus reduction curve.
            provide a list of shear strains and corresponding shear modulus reduction ratios (`ratios`)
        ratios: list
            Instead of automatic surfaces generation,
            you can define yield surfaces directly based on
            desired shear modulus reduction curve.
            provide a list of shear strains and corresponding shear modulus reduction ratios (`ratios`)
        """
        self.osi = osi
        self.nd = nd
        self.rho = float(rho)
        self.g_mod_ref = float(g_mod_ref)
        self.bulk_mod_ref = float(bulk_mod_ref)
        self.cohesion = float(cohesion)
        self.peak_strain = float(peak_strain)
        self.phi = float(phi)
        self.p_ref = float(p_ref)
        self.d = float(d)

        if strains is not None:
            assert len(strains) == len(ratios)
            yield_surf = []
            for i in range(len(strains)):
                yield_surf.append(strains[i])
                yield_surf.append(ratios[i])
            self.yield_surf = yield_surf
            n_surf = -len(strains)  # from docs 'add a minus sign in front of noYieldSurf'
        else:
            self.yield_surf = None

        assert abs(n_surf) < 40
        self.n_surf = int(n_surf)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat

        self._parameters = [self.op_type, self._tag, self.nd, self.rho, self.g_mod_ref, self.bulk_mod_ref,
                            self.cohesion, self.peak_strain, self.phi, self.p_ref, self.d]

        if self.yield_surf is not None:
            self._parameters.append(self.n_surf)
            self._parameters += list(self.yield_surf)

        else:
            # self._keyword_args['noYieldSurf'] = self.no_yield_surf
            self._parameters.append(self.n_surf)
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def update_to_nonlinear(self):
        from o3seespy import update_material_stage
        update_material_stage(self.osi, self, 1)

    def update_to_linear(self):
        from o3seespy import update_material_stage
        update_material_stage(self.osi, self, 0)

    def set_nu(self, nu, ele=None, eles=None, adj_g_mod=False):
        if adj_g_mod:
            g_mod = 3 * self.bulk_mod_ref * (1 - 2 * nu) / (2 * (1 + nu))
            self.set_parameter(self.osi, 'shearModulus', g_mod, ele, eles)
        else:
            bulk_mod = 2 * self.g_mod_ref * (1 + nu) / (3 * (1 - 2 * nu))
            self.set_parameter(self.osi, 'bulkModulus', bulk_mod, ele, eles)

    @property
    def nu(self):
        return (3 * self.bulk_mod_ref - 2 * self.g_mod_ref) / (2 * (3 * self.bulk_mod_ref + self.g_mod_ref))


class PressureDependMultiYield(NDMaterialBase):
    """
    The PressureDependMultiYield NDMaterial Class

    PressureDependMultiYield material is an elastic-plastic material for simulating the essential response
    characteristics of pressure sensitive soil materials under general loading conditions. Such
    characteristics include dilatancy (shear-induced volume contraction or dilation) and
    non-flow liquefaction (cyclic mobility), typically exhibited in sands or silts during monotonic or cyclic loading.
    """
    op_type = 'PressureDependMultiYield'

    def __init__(self, osi, nd, rho, g_mod_ref, bulk_mod_ref, phi, peak_strain, p_ref,
                 d, pt_ang, con_rate, dil_rates, liquefac, n_surf=20.0, strains=None, ratios=None,
                 e_init=0.6, cs_params=None, c=0.3):
        r"""
        Initial method for PressureDependMultiYield

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        nd: int
            Number of dimensions, 2 for plane-strain, and 3 for 3d analysis.
        rho: float
            Saturated soil mass density.
        g_mod_ref, float
            (:math:`G_0`) Reference low-strain shear modulus,
            specified at a reference mean effective confining pressure (`p_ref`).
       bulk_mod_ref: float
            (:math:`B_r`) Reference bulk modulus,
            specified at a reference mean effective confining pressure (`p_ref`)).
        phi: float
            (:math:`phi`) friction angle at peak shear strength in degrees, optional (default is 0.0).
        peak_strain: float
            (:math:`\\gamma_{max}`) An octahedral shear strain at
            which the maximum shear strength is reached,
            specified at a reference mean effective confining
            pressure refPress of p'r (see below).
        p_ref: float
            (:math:`p'_r`) reference mean effective confining pressure at which :math:`g_r`, :math:`b_r`, and
            :math:`\gamma_{max}` are defined, optional (default is 100. kpa).
        d: float
            (:math:`d`) a positive constant defining variations of :math:`g` and :math:`b` as a function of
            instantaneous effective confinement :math:`p'` (default is 0.0) :math:`g=g_r(\frac{p'}{p'_r})^d`
            :math:`b=b_r(\frac{p'}{p'_r})^d` if :math:`\phi=0`, :math:`d` is reset to 0.0.
        pt_ang: float
            (:math:`\phi_{pt}`) phase transformation angle, in degrees.
        con_rate: float
            A non-negative constant defining the rate of shear-induced volume decrease (contraction) or pore pressure
            buildup. a larger value corresponds to faster contraction rate.
        dil_rates: list
            Non-negative constants defining the rate of shear-induced volume increase (dilation). larger values
            correspond to stronger dilation rate. ``dil_rates = [dilat1, dilat2]``.
        liquefac: list
            Parameters controlling the mechanism of liquefaction-induced perfectly plastic shear strain accumulation,
            i.e., cyclic mobility. set ``liquefac[0] = 0`` to deactivate this mechanism altogether. ``liquefac[0]`` defines the
            effective confining pressure (e.g., 10 kpa in si units or 1.45 psi in english units) below which the mechanism is
            in effect. smaller values should be assigned to denser sands. ``liquefac[1]`` defines the maximum amount of
            perfectly plastic shear strain developed at zero effective confinement during each loading phase. smaller
            values should be assigned to denser sands. ``liquefac[2]`` defines the maximum amount of biased
            perfectly plastic shear strain :math:`\gamma_b` accumulated at each loading phase under biased
            shear loading conditions, as :math:`\gamma_b=liquefac[1]\times liquefac[2]`. typically,
            :math:`liquefac[2]` takes a value between 0.0 and 3.0. smaller values should be
            assigned to denser sands. see the references listed at the end of this chapter for more information.
        n_surf: int, optional
            Number of yield surfaces, optional (must be less than 40, default is 20). the surfaces are generated based
            on the hyperbolic relation defined in note 2 below.
        strains: list
            Instead of automatic surfaces generation,
            you can define yield surfaces directly based on
            desired shear modulus reduction curve.
            provide a list of shear strains and corresponding shear modulus reduction ratios (`ratios`)
        ratios: list
            Instead of automatic surfaces generation,
            you can define yield surfaces directly based on
            desired shear modulus reduction curve.
            provide a list of shear strains and corresponding shear modulus reduction ratios (`ratios`)
        e_init: float, optional
            Initial void ratio, optional (default is 0.6).
        cs_params: list (default=[0.9, 0.02, 0.7, 101.0]), optional
            ``params=[cs1, cs2, cs3, pa]`` defining a straight critical-state line ec in e-p' space. if cs3=0, ec =
            cs1-cs2 log(p'/pa) else (li and wang, jgge, 124(12)), ec = cs1-cs2(p'/pa)cs3 where pa is atmospheric pressure for
            normalization (typically 101 kpa in si units, or 14.65 psi in english units). all four constants are optional
        c: float, optional
            Numerical constant (default value = 0.3 kPa)
        """
        self.osi = osi
        self.nd = int(nd)
        self.rho = float(rho)
        self.g_mod_ref = float(g_mod_ref)
        self.bulk_mod_ref = float(bulk_mod_ref)
        self.peak_strain = float(peak_strain)
        self.phi = float(phi)
        self.p_ref = float(p_ref)
        self.d = float(d)
        self.pt_ang = float(pt_ang)
        self.con_rate = float(con_rate)
        self.dil_rates = dil_rates
        self.liquefac = liquefac
        assert n_surf < 40
        self.n_surf = int(n_surf)
        if strains is not None:
            assert len(strains) == len(ratios)
            yield_surf = []
            for i in range(len(strains)):
                yield_surf.append(ratios[i])
                yield_surf.append(strains[i])
            self.yield_surf = yield_surf
        else:
            self.yield_surf = None
        self.e_init = float(e_init)
        self.cs_params = cs_params
        self.c = float(c)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.nd, self.rho, self.g_mod_ref, self.bulk_mod_ref,
                            self.phi, self.peak_strain, self.p_ref, self.d, self.pt_ang,
                            self.con_rate, *self.dil_rates, *self.liquefac, self.n_surf]
        if self.yield_surf is not None:
            self._parameters += list(self.yield_surf)
        special_pms = ['e_init', 'cs_params', 'c']
        packets = [False, True, False]
        for i, pm in enumerate(special_pms):
            if getattr(self, pm) is not None:
                if packets[i]:
                    self._parameters += [*getattr(self, pm)]
                else:
                    self._parameters += [getattr(self, pm)]
            else:
                break
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def update_to_nonlinear(self):
        from o3seespy import update_material_stage
        update_material_stage(self.osi, self, 1)

    def set_nu(self, nu, ele=None, eles=None, adj_g_mod=False):
        if adj_g_mod:
            g_mod = 3 * self.bulk_mod_ref * (1 - 2 * nu) / (2 * (1 + nu))
            self.set_parameter(self.osi, 'shearModulus', g_mod, ele, eles)
        else:
            bulk_mod = 2 * self.g_mod_ref * (1 + nu) / (3 * (1 - 2 * nu))
            self.set_parameter(self.osi, 'bulkModulus', bulk_mod, ele, eles)

    @property
    def nu(self):
        return (3 * self.bulk_mod_ref - 2 * self.g_mod_ref) / (2 * (3 * self.bulk_mod_ref + self.g_mod_ref))


class PressureDependMultiYield02(NDMaterialBase):
    r"""
    The PressureDependMultiYield02 NDMaterial Class

    PressureDependMultiYield02 material is modified from PressureDependMultiYield material, with:#. additional
    parameters (``contrac[2]`` and ``dilat[2]``) to account for :math:`K_{\sigma}` effect,#. a parameter to
    account for the influence of previous dilation history on subsequent contraction phase
    (``contrac[1]``), and#. modified logic related to permanent shear strain accumulation
    (``liquefac[0]`` and ``liquefac[1]``).
    """
    op_type = 'PressureDependMultiYield02'

    def __init__(self, osi, nd, rho, g_mod_ref, bulk_mod_ref, phi, peak_strain, p_ref,
                 d, pt_ang, con_rates, dil_rates, liquefac=(1., 0.), n_surf=20.0, strains=None, ratios=None,
                 e_init=0.6, cs_params=None, c=0.1):
        r"""
        Initial method for PressureDependMultiYield02

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        nd: int
            Number of dimensions, 2 for plane-strain, and 3 for 3d analysis.
        rho: float
            Saturated soil mass density.
        g_mod_ref, float
            (:math:`G_0`) Reference low-strain shear modulus,
            specified at a reference mean effective confining pressure (`p_ref`).
       bulk_mod_ref: float
            (:math:`B_r`) Reference bulk modulus,
            specified at a reference mean effective confining pressure (`p_ref`)).
        phi: float
            (:math:`phi`) friction angle at peak shear strength in degrees, optional (default is 0.0).
        peak_strain: float
            (:math:`\\gamma_{max}`) An octahedral shear strain at
            which the maximum shear strength is reached,
            specified at a reference mean effective confining
            pressure refPress of p'r (see below).
        p_ref: float
            (:math:`p'_r`) reference mean effective confining pressure at which :math:`g_r`, :math:`b_r`, and
            :math:`\gamma_{max}` are defined, optional (default is 100. kpa).
        d: float
            (:math:`d`) a positive constant defining variations of :math:`g` and :math:`b` as a function of
            instantaneous effective confinement :math:`p'` (default is 0.0) :math:`g=g_r(\frac{p'}{p'_r})^d`
            :math:`b=b_r(\frac{p'}{p'_r})^d` if :math:`\phi=0`, :math:`d` is reset to 0.0.
        pt_ang: float
            (:math:`\phi_{pt}`) phase transformation angle, in degrees.
        con_rates: list
            A list of constants defining contraction behaviour `[contrac1, contrac2, contrac3]`.
             `contrac1`: non-negative constant defining the rate of shear-induced volume decrease (contraction)
             or pore pressure buildup. a larger value corresponds to faster contraction rate.
             `contrac2`: A non-negative constant reflecting dilation history on contraction tendency.
             `contrac3`: A non-negative constant reflecting Kσ effect.
        dil_rates: list
            A list of constants defining contraction behaviour `[dilat1, dilat2, dilat3]`.
            `dilat1` and `dilat2`: Non-negative constants defining the rate of shear-induced volume
            increase (dilation). larger values correspond to stronger dilation rate.
            `dilat3`: A non-negative constant reflecting Kσ effect.
        liquefac: list
            Parameters controlling the mechanism of liquefaction-induced perfectly plastic shear strain accumulation,
            i.e., cyclic mobility. [`liquefac1`, `liquefac2`]. NOTE: Different from `PressureDependMultiYield`.
            `liquefac1`: Damage parameter to define accumulated permanent shear strain as a function of
            dilation history.
            `liquefac2`: Damage parameter to define biased accumulation of permanent shear strain as a
            function of load reversal history.
        n_surf: int, optional
            Number of yield surfaces, optional (must be less than 40, default is 20). the surfaces are generated based
            on the hyperbolic relation defined in note 2 below.
        strains: list
            Instead of automatic surfaces generation,
            you can define yield surfaces directly based on
            desired shear modulus reduction curve.
            provide a list of shear strains and corresponding shear modulus reduction ratios (`ratios`)
        ratios: list
            Instead of automatic surfaces generation,
            you can define yield surfaces directly based on
            desired shear modulus reduction curve.
            provide a list of shear strains and corresponding shear modulus reduction ratios (`ratios`)
        e_init: float, optional
            Initial void ratio, optional (default is 0.6).
        cs_params: list (default=[0.9, 0.02, 0.7, 101.0]), optional
            ``params=[cs1, cs2, cs3, pa]`` defining a straight critical-state line ec in e-p' space. if cs3=0, ec =
            cs1-cs2 log(p'/pa) else (li and wang, jgge, 124(12)), ec = cs1-cs2(p'/pa)cs3 where pa is atmospheric pressure for
            normalization (typically 101 kpa in si units, or 14.65 psi in english units). all four constants are optional
        c: float, optional
            Numerical constant (default value = 0.1 kPa)
        """
        self.osi = osi
        self.nd = int(nd)
        self.rho = float(rho)
        self.g_mod_ref = float(g_mod_ref)
        self.bulk_mod_ref = float(bulk_mod_ref)
        self.peak_strain = float(peak_strain)
        self.phi = float(phi)
        self.p_ref = float(p_ref)
        self.d = float(d)
        self.pt_ang = float(pt_ang)
        self.con_rates = con_rates
        contrac1 = con_rates[0]
        contrac2 = con_rates[1]
        contrac3 = con_rates[2]
        self.dil_rates = dil_rates
        dilat1 = dil_rates[0]
        dilat2 = dil_rates[1]
        dilat3 = dil_rates[2]
        self.liquefac = liquefac
        assert n_surf < 40
        self.n_surf = int(n_surf)
        if strains is not None:
            assert len(strains) == len(ratios)
            yield_surf = []
            for i in range(len(strains)):
                yield_surf.append(ratios[i])
                yield_surf.append(strains[i])
            self.yield_surf = yield_surf
        else:
            self.yield_surf = None
        self.e_init = float(e_init)
        self.cs_params = cs_params
        self.c = float(c)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.nd, self.rho, self.g_mod_ref, self.bulk_mod_ref,
                            self.phi, self.peak_strain, self.p_ref, self.d, self.pt_ang,
                            contrac1, contrac3, dilat1, dilat3, self.n_surf]
        if self.yield_surf is not None:
            self._parameters += list(self.yield_surf)

        self._parameters += [contrac2, dilat2]
        special_pms = ['liquefac', 'e_init', 'cs_params', 'c']
        packets = [True, False, True, False]
        for i, pm in enumerate(special_pms):
            if getattr(self, pm) is not None:
                if packets[i]:
                    self._parameters += [*getattr(self, pm)]
                else:
                    self._parameters += [getattr(self, pm)]
            else:
                break
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def update_to_nonlinear(self):
        from o3seespy import update_material_stage
        update_material_stage(self.osi, self, 1)

    def set_nu(self, nu, ele=None, eles=None, adj_g_mod=False):
        if adj_g_mod:
            g_mod = 3 * self.bulk_mod_ref * (1 - 2 * nu) / (2 * (1 + nu))
            self.set_parameter(self.osi, 'shearModulus', g_mod, ele, eles)
        else:
            bulk_mod = 2 * self.g_mod_ref * (1 + nu) / (3 * (1 - 2 * nu))
            self.set_parameter(self.osi, 'bulkModulus', bulk_mod, ele, eles)

    @property
    def nu(self):
        return (3 * self.bulk_mod_ref - 2 * self.g_mod_ref) / (2 * (3 * self.bulk_mod_ref + self.g_mod_ref))


class Steel01(UniaxialMaterialBase):
    """
    The Steel01 UniaxialMaterial Class

    This command is used to construct a uniaxial bilinear steel material object with kinematic hardening and optional
    isotropic hardening described by a non-linear evolution equation (REF: Fedeas).
    """
    op_type = "Steel01"

    def __init__(self, osi, fy: float, e0: float, b: float, a1=None, a2=None, a3=None, a4=None):
        """
        Initial method for Steel01

        Parameters
        ----------
        fy: float
            Yield strength
        e0: float
            Initial elastic tangent
        b: float
            Strain-hardening ratio (ratio between post-yield tangent and initial elastic tangent)
        a1: float
            Isotropic hardening parameter, increase of compression yield envelope as proportion of yield strength after
            a plastic strain of :math:`a_2*(f_y/e_0)` (optional)
        a2: float
            Isotropic hardening parameter
        a3: float
            Isotropic hardening parameter, increase of tension yield envelope as proportion of yield strength after a
            plastic strain of :math:`a_4*(f_y/e_0)`. (optional)
        a4: float
            Isotropic hardening parameter (see explanation
        """
        self.osi = osi
        self.fy = float(fy)
        self.e0 = float(e0)
        self.b = float(b)
        self.a_values = [a1, a2, a3, a4]
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.a4 = a4
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self.tag, self.fy, self.e0, self.b]
        for a in self.a_values:
            if a is None:
                break
            self._parameters.append(a)
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def set_fy(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'Fy', value, ele, eles)

    def set_e_mod(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'E', value, ele, eles)

    def set_b(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'b', value, ele, eles)

    def set_a1(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'a1', value, ele, eles)

    def set_a2(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'a2', value, ele, eles)

    def set_a3(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'a3', value, ele, eles)

    def set_a4(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'a4', value, ele, eles)


class Fiber(SectionBase):
    """
    The Fiber Section Class

    This command allows the user to construct a FiberSection object. Each FiberSection object is composed of Fibers,
    with each fiber containing a UniaxialMaterial, an area and a location (y,z). The dofs for 2D section are ``[P,
    Mz]``,for 3D are ``[P,Mz,My,T]``.
    """
    op_type = 'Fiber'

    def __init__(self, osi, gj: float = None, torsion_mat=None):
        """
        Initial method for Fiber

        Supports pre-building

        Parameters
        ----------
        gj: float
            Linear-elastic torsional stiffness assigned to the section
        torsion_mat: obj
            uniaxial_material object assigned to the section for torsional response (can be nonlinear)
        """
        self.osi = osi
        if gj is None:
            self.gj = None
        else:
            self.gj = float(gj)
        self.torsion_mat = torsion_mat
        osi.n_sect += 1
        self._tag = osi.n_sect
        self._parameters = [self.op_type, self._tag]
        if getattr(self, 'gj') is not None:
            self._parameters += ['-GJ', self.gj]
        if getattr(self, 'torsion_mat') is not None:
            self._parameters += ['-torsion', self.torsion_mat.tag]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class Circ(LayerBase):
    """
    The Circ Layer Class

    This command is used to construct a line of fibers along a circular arc
    """
    op_type = 'circ'

    def __init__(self, osi, mat, num_fiber, area_fiber, center, radius, ang=None):
        """
        Initial method for Circ

        Supports pre-building

        Parameters
        ----------
        mat: obj
            Material tag associated with this fiber (uniaxial_material for a fiber section and nd_material for use
            in an nd fiber section).
        num_fiber: int
            Number of fibers along line
        area_fiber: float
            Area of each fiber
        center: listf
            Y & z-coordinates of center of circular arc
        radius: float
            Radius of circlular arc
        ang: listf
            Starting and ending angle (optional) [0.0, 360.0-360/num_fibres]

        """
        self.osi = osi
        self.mat = mat
        self.num_fiber = int(num_fiber)
        self.area_fiber = float(area_fiber)
        self.center = center
        self.radius = float(radius)
        self.ang = ang
        self._parameters = [self.op_type, self.mat.tag, self.num_fiber, self.area_fiber, *self.center, self.radius]
        if self.ang is not None:
            self._parameters += self.ang
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)
#
# class PySimple1(UniaxialMaterialBase):
#     op_type = "PySimple1"
#
#     def __init__(self, osi, soil_type: int, p_ult, y50, cd, c):
#         """
#         PySimple1 uniaxial material object
#
#         Parameters
#         ----------
#         osi : opensees_pack.opensees_instance.OpenSeesInstance object
#             An instance of opensees
#         soil_type : int {1, 2}
#             Backbone type for soil
#         p_ult : float
#             Ultimate capacity of the p-y material. Note that “p” or “pult” are distributed loads
#              [force per length of pile] in common design equations, but are both loads for this
#              uniaxialMaterial [i.e., distributed load times the tributary length of the pile].
#         y50 : float
#             Displacement at which 50% of pult is mobilized in monotonic loading.
#         cd : float
#             Variable that sets the drag resistance within a fully-mobilized gap as Cd*pult.
#         c : float
#             The viscous damping term (dashpot) on the far-field (elastic) component
#             of the displacement rate (velocity). (optional Default = 0.0). Nonzero c
#             values are used to represent radiation damping effects
#         """
#         self.soil_type = int(soil_type)
#         self.p_ult = p_ult
#         self.y50 = y50
#         self.cd = cd
#         self.c = c
#         osi.n_mat += 1
#         self._tag = osi.n_mat
#         self._parameters = [self.op_type, self._tag, self.soil_type, self.p_ult, self.y50, self.cd, self.c]
#         self.to_process(osi)
#
#
# class PyLiq1(UniaxialMaterialBase):
#     op_type = "PyLiq1"
#
#     def __init__(self, osi, soil_type, p_ult, y50, cd, c, p_res, ele1=None, ele2=None, time_series=None):
#         """
#
#         Parameters
#         ----------
#         osi : opensees_pack.opensees_instance.OpenSeesInstance object
#             An instance of opensees
#         soil_type
#        p_ult : float
#             Ultimate capacity of the p-y material. Note that “p” or “pult” are distributed loads
#              [force per length of pile] in common design equations, but are both loads for this
#              uniaxialMaterial [i.e., distributed load times the tributary length of the pile].
#         y50 : float
#             Displacement at which 50% of pult is mobilized in monotonic loading.
#         cd : float
#             Variable that sets the drag resistance within a fully-mobilized gap as Cd*pult.
#         c : float
#             The viscous damping term (dashpot) on the far-field (elastic) component
#             of the displacement rate (velocity). (optional Default = 0.0). Nonzero c
#             values are used to represent radiation damping effects
#         p_res : float
#         ele1 : opensees_pack.ElementBase object, optional
#             the eleTag (element numbers) for the two solid element from which PyLiq1
#             will obtain mean effective stresses and excess pore pressures
#         ele2 : opensees_pack.ElementBase object, optional
#             Same as ele1
#         time_series : iterable object, optional
#             A time series of mean effective stress values
#         """
#         self.soil_type = soil_type
#         self.p_ult = p_ult
#         self.y50 = y50
#         self.cd = cd
#         self.c = c
#         self.p_res = p_res
#         self.ele1 = ele1
#         self.ele2 = ele2
#         self.time_series = time_series
#         osi.n_mat += 1
#         self._tag = osi.n_mat
#         self._parameters = [self.op_type, self._tag, self.p_ult, self.y50, self.cd, self.c, self.p_res]
#         if self.ele1 is None:
#             self._parameters.append("-timeSeries")
#             self._parameters.append(self.time_series)
#         else:
#             self._parameters.append(self.ele1)
#             if self.ele2 is None:
#                 self._parameters.append(self.ele2)
#         self.to_process(osi)


class StressDensity(NDMaterialBase):
    op_type = "stressDensity"

    def __init__(self, osi, den, e_init, big_a, n, nu, a1, b1, a2, b2, a3, b3, fd, mu_0, mu_cyc, sc, big_m, p_atm,
                 ssls=None, hsl=None, ps=None):
        """
        Initial method for StressDensity

        Parameters
        ----------
        m_den: float
            Mass density
        e_not: float
            Initial void ratio
        big_a: float
            Constant for elastic shear modulus
        n: float
            Pressure dependency exponent for elastic shear modulus
        nu: float
            Poisson's ratio
        a1: float
            Peak stress ratio coefficient (:math:`etamax = a1 + b1*is`)
        b1: float
            Peak stress ratio coefficient (:math:`etamax = a1 + b1*is`)
        a2: float
            Max shear modulus coefficient (:math:`gn_max = a2 + b2*is`)
        b2: float
            Max shear modulus coefficient (:math:`gn_max = a2 + b2*is`)
        a3: float
            Min shear modulus coefficient (:math:`gn_min = a3 + b3*is`)
        b3: float
            Min shear modulus coefficient (:math:`gn_min = a3 + b3*is`)
        fd: float
            Degradation constant
        mu_not: float
            Dilatancy coefficient (monotonic loading)
        mu_cyc: float
            Dilatancy coefficient (cyclic loading)
        sc: float
            Dilatancy strain
        big_m: float
            Critical state stress ratio
        patm: float
            Atmospheric pressure (in appropriate units)
        ssls: listf
            Void ratio of quasi steady state (qss-line) at pressures [pmin, 10kpa, 30kpa, 50kpa, 100kpa, 200kpa, 400kpa]
            (default = [0.877, 0.877, 0.873, 0.870, 0.860, 0.850, 0.833])
        hsl: float
            Void ratio of upper reference state (ur-line) for all pressures (default = 0.895)
        p1: float
            Pressure corresponding to ssl1 (default = 1.0 kpa)
        """
        self.osi = osi
        self.den = float(den)
        self.e_init = float(e_init)
        self.big_a = float(big_a)
        self.n = float(n)
        self.nu = float(nu)
        self.a1 = float(a1)
        self.b1 = float(b1)
        self.a2 = float(a2)
        self.b2 = float(b2)
        self.a3 = float(a3)
        self.b3 = float(b3)
        self.fd = float(fd)
        self.mu_0 = float(mu_0)
        self.mu_cyc = float(mu_cyc)
        self.sc = float(sc)
        self.big_m = float(big_m)
        self.p_atm = float(p_atm)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat

        self._parameters = [self.op_type, self._tag, self.den, self.e_init, self.big_a, self.n, self.nu, self.a1, self.b1, self.a2,
                            self.b2, self.a3, self.b3, self.fd, self.mu_0, self.mu_cyc, self.sc, self.big_m,
                            self.p_atm]
        if ssls is not None:
            assert len(ssls) == 7, len(ssls)
            self.ssls = [float(x) for x in ssls]
            if hsl is None:
                self.hsl = 0.895
            else:
                self.hsl = float(hsl)
            self._parameters += [*self.ssls, self.hsl]
        if ps is not None:
            self.ps = [float(x) for x in ps]
            self._parameters += [*self.ps]

        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class TFP(ElementBase):
    """
    The TFP Element Class

    This command is used to construct a Triple Friction Pendulum Bearing element object, which is defined by two nodes.
    The element can have zero length or the appropriate bearing height. The bearing has unidirectional (2D) or coupled
    (3D) friction properties (with post-yield stiffening due to the concave sliding surface) for the shear
    deformations, and force-deformation behaviors defined by UniaxialMaterials in the remaining two (2D)
    or four (3D) directions. To capture the uplift behavior of the bearing, the user-specified
    UniaxialMaterial in the axial direction is modified for no-tension behavior. P-Delta
    moments are entirely transferred to the concave sliding surface (iNode). It is
    important to note that rotations of the concave sliding surface (rotations at
    the iNode) affect the shear behavior of the bearing. If the element has
    non-zero length, the local x-axis is determined from the nodal
    geometry unless the optional x-axis vector is specified in
    which case the nodal geometry is ignored and the user-defined orientation is utilized.


    """
    op_type = 'TFP'

    def __init__(self, osi, ele_nodes, r1, r2, r3, r4, db1, db2, db3, db4, d1, d2, d3, d4, mu1, mu2, mu3, mu4, h1, h2,
                 h3, h4, h0, col_load, big_k=None):
        """
        Initial method for TFP

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of two element nodes
        r1: float
            Radius of inner bottom sliding surface
        r2: float
            Radius of inner top sliding surface
        r3: float
            Radius of outer bottom sliding surface
        r4: float
            Radius of outer top sliding surface
        db1: float
            Diameter of inner bottom sliding surface
        db2: float
            Diameter of inner top sliding surface
        db3: float
            Diameter of outer bottom sliding surface
        db4: float
            Diameter of outer top sliding surface
        d1: float
            Diameter of inner slider
        d2: float
            Diameter of inner slider
        d3: float
            Diameter of outer bottom slider
        d4: float
            Diameter of outer top slider
        mu1: float
            Friction coefficient of inner bottom sliding surface
        mu2: float
            Friction coefficient of inner top sliding surface
        mu3: float
            Friction coefficient of outer bottom sliding surface
        mu4: float
            Friction coefficient of outer top sliding surface
        h1: float
            Height from inner bottom sliding surface to center of bearing
        h2: float
            Height from inner top sliding surface to center of bearing
        h3: float
            Height from outer bottom sliding surface to center of bearing
        h4: float
            Height from inner top sliding surface to center of bearing
        h0: float
            Total height of bearing
        col_load: float
            Initial axial load on bearing (only used for first time step then load come from model)
        big_k: float
            Optional, stiffness of spring in vertical dirn (dof 2 if ndm= 2, dof 3 if ndm = 3) (default=1.0e15)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [1, 0]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
        >>> o3.element.TFP(osi, ele_nodes=ele_nodes,
        >>>                r1=1.0, r2=1.0, r3=1.0, r4=1.0,
        >>>                db1=1.0, db2=1.0, db3=1.0, db4=1.0,
        >>>                d1=1.0, d2=1.0, d3=1.0, d4=1.0,
        >>>                mu1=0.3, mu2=0.4, mu3=0.5, mu4=0.5,
        >>>                h1=1.0, h2=1.0, h3=1.0, h4=1.0,
        >>>                h0=1.0, col_load=1.0, big_k=None)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.r1 = float(r1)
        self.r2 = float(r2)
        self.r3 = float(r3)
        self.r4 = float(r4)
        self.db1 = float(db1)
        self.db2 = float(db2)
        self.db3 = float(db3)
        self.db4 = float(db4)
        self.d1 = float(d1)
        self.d2 = float(d2)
        self.d3 = float(d3)
        self.d4 = float(d4)
        self.mu1 = float(mu1)
        self.mu2 = float(mu2)
        self.mu3 = float(mu3)
        self.mu4 = float(mu4)
        self.h1 = float(h1)
        self.h2 = float(h2)
        self.h3 = float(h3)
        self.h4 = float(h4)
        self.h0 = float(h0)
        self.col_load = float(col_load)
        if big_k is not None:
            self.big_k = float(big_k)
        else:
            self.big_k = None
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.r1, self.r2, self.r3, self.r4, self.db1,
                            self.db2, self.db3, self.db4, self.d1, self.d2, self.d3, self.d4, self.mu1, self.mu2,
                            self.mu3, self.mu4, self.h1, self.h2, self.h3, self.h4, self.h0, self.col_load]
        if getattr(self, 'big_k') is not None:
            self._parameters += [self.big_k]
        self.to_process(osi)


class ElastomericBearingBoucWen2D(ElementBase):
    """
    The ElastomericBearingBoucWen2D Element Class

    This command is used to construct an elastomericBearing element object, which is defined by two nodes. The element
    can have zero length or the appropriate bearing height. The bearing has unidirectional (2D) or coupled (3D)
    plasticity properties for the shear deformations, and force-deformation behaviors defined by
    UniaxialMaterials in the remaining two (2D) or four (3D) directions. By default (sDratio =
    0.5) P-Delta moments are equally distributed to the two end-nodes. To avoid the
    introduction of artificial viscous damping in the isolation system (sometimes
    referred to as "damping leakage in the isolation system"), the bearing
    element does not contribute to the Rayleigh damping by default. If
    the element has non-zero length, the local x-axis is determined
    from the nodal geometry unless the optional x-axis vector is
    specified in which case the nodal geometry is ignored and
    the user-defined orientation is utilized.

    For a two-dimensional problem
    """
    op_type = 'elastomericBearingBoucWen'

    def __init__(self, osi, ele_nodes, k_init, qd, alpha1, alpha2, mu, eta, beta, gamma, p_mat=None, mz_mat=None,
                 orient_vals: list = None, shear_dist: float = None, do_rayleigh=False, mass: float = None):
        """
        Initial method for ElastomericBearingBoucWen2D

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of two element nodes
        k_init: float
            Initial elastic stiffness in local shear direction
        qd: float
            Characteristic strength
        alpha1: float
            Post yield stiffness ratio of linear hardening component
        alpha2: float
            Post yield stiffness ratio of non-linear hardening component
        mu: float
            Exponent of non-linear hardening component
        eta: float
            Yielding exponent (sharpness of hysteresis loop corners) (default = 1.0)
        beta: float
            First hysteretic shape parameter (default = 0.5)
        gamma: float
            Second hysteretic shape parameter (default = 0.5)
        p_mat: obj, optional
            Object associated with previously-defined uniaxial_material in axial direction
        mz_mat: obj, optional
            Object associated with previously-defined uniaxial_material in moment direction around local z-axis
        orient_vals: list, optional
            Vector components in global coordinates defining local x-axis , vector components in global coordinates
            defining local y-axis
        shear_dist: float, optional
            Shear distance from inode as a fraction of the element length (optional, default = 0.5)
        do_rayleigh: bool
            To include rayleigh damping from the bearing (optional, default = no rayleigh damping contribution)
        mass: float, optional
            Element mass (optional, default = 0.0)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [1, 0]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
        >>> orient_vals = [1, 0, 0, 1, 0, 1]
        >>> p_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
        >>> mz_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
        >>> o3.element.ElastomericBearingBoucWen2D(osi, ele_nodes=ele_nodes, k_init=1.0, qd=1.0, alpha1=1.0, alpha2=1.0,
        >>>                                        mu=1.0, eta=1.0, beta=1.0, gamma=1.0, p_mat=p_mat, mz_mat=mz_mat,
        >>>                                        orient_vals=orient_vals, shear_dist=1.0, do_rayleigh=False, mass=1.0)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.k_init = float(k_init)
        self.qd = float(qd)
        self.alpha1 = float(alpha1)
        self.alpha2 = float(alpha2)
        self.mu = float(mu)
        self.eta = float(eta)
        self.beta = float(beta)
        self.gamma = float(gamma)
        self.p_mat = p_mat
        self.mz_mat = mz_mat
        self.orient_vals = orient_vals
        if shear_dist is None:
            self.shear_dist = None
        else:
            self.shear_dist = float(shear_dist)
        self.do_rayleigh = do_rayleigh
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.k_init, self.qd, self.alpha1, self.alpha2,
                            self.mu, self.eta, self.beta, self.gamma]
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat.tag]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat.tag]
        if getattr(self, 'orient_vals') is not None:
            self._parameters += ['-orient', *self.orient_vals]
        if getattr(self, 'shear_dist') is not None:
            self._parameters += ['-shearDist', self.shear_dist]
        if getattr(self, 'do_rayleigh'):
            self._parameters += ['-doRayleigh']
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        try:
            self.to_process(osi)
        except ValueError:
            self._parameters[0] = 'ElastomericBearingBoucWen'
            self.to_process(osi)


class ElastomericBearingBoucWen3D(ElementBase):
    """
    The ElastomericBearingBoucWen3D Element Class

    This command is used to construct an elastomericBearing element object, which is defined by two nodes. The element
    can have zero length or the appropriate bearing height. The bearing has unidirectional (2D) or coupled (3D)
    plasticity properties for the shear deformations, and force-deformation behaviors defined by
    UniaxialMaterials in the remaining two (2D) or four (3D) directions. By default (sDratio =
    0.5) P-Delta moments are equally distributed to the two end-nodes. To avoid the
    introduction of artificial viscous damping in the isolation system (sometimes
    referred to as "damping leakage in the isolation system"), the bearing
    element does not contribute to the Rayleigh damping by default. If
    the element has non-zero length, the local x-axis is determined
    from the nodal geometry unless the optional x-axis vector is
    specified in which case the nodal geometry is ignored and
    the user-defined orientation is utilized.

    For a three-dimensional problem
    """
    op_type = 'elastomericBearingBoucWen'

    def __init__(self, osi, ele_nodes, k_init, qd, alpha1, alpha2, mu, eta, beta, gamma, p_mat=None, t_mat=None,
                 my_mat=None, mz_mat=None, orient_vals: list = None, shear_dist: float = None, do_rayleigh=False,
                 mass: float = None):
        """
        Initial method for ElastomericBearingBoucWen3D

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of two element nodes
        k_init: float
            Initial elastic stiffness in local shear direction
        qd: float
            Characteristic strength
        alpha1: float
            Post yield stiffness ratio of linear hardening component
        alpha2: float
            Post yield stiffness ratio of non-linear hardening component
        mu: float
            Exponent of non-linear hardening component
        eta: float
            Yielding exponent (sharpness of hysteresis loop corners) (default = 1.0)
        beta: float
            First hysteretic shape parameter (default = 0.5)
        gamma: float
            Second hysteretic shape parameter (default = 0.5)
        p_mat: obj, optional
            Object associated with previously-defined uniaxial_material in axial direction
        t_mat: obj, optional
            Object associated with previously-defined uniaxial_material in torsional direction
        my_mat: obj, optional
            Object associated with previously-defined uniaxial_material in moment direction around local y-axis
        mz_mat: obj, optional
            Object associated with previously-defined uniaxial_material in moment direction around local z-axis
        orient_vals: list, optional
            Vector components in global coordinates defining local x-axis , vector components in global coordinates
            defining local y-axis
        shear_dist: float, optional
            Shear distance from inode as a fraction of the element length (optional, default = 0.5)
        do_rayleigh: bool
            To include rayleigh damping from the bearing (optional, default = no rayleigh damping contribution)
        mass: float, optional
            Element mass (optional, default = 0.0)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=3, ndf=6)
        >>> coords = [[0, 0, 0], [0, 1, 0]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
        >>> orient_vals = [1, 0, 0]
        >>> p_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
        >>> mz_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
        >>> t_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
        >>> my_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
        >>> o3.element.ElastomericBearingBoucWen3D(osi, ele_nodes=ele_nodes, k_init=1.0, qd=1.0, alpha1=1.0, alpha2=1.0,
        >>>                                        mu=1.0, eta=1.0, beta=1.0, gamma=1.0, p_mat=p_mat, t_mat=t_mat,
        >>>                                        my_mat=my_mat, mz_mat=mz_mat, orient_vals=orient_vals,
        >>>                                        shear_dist=1.0, do_rayleigh=False, mass=1.0)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.k_init = float(k_init)
        self.qd = float(qd)
        self.alpha1 = float(alpha1)
        self.alpha2 = float(alpha2)
        self.mu = float(mu)
        self.eta = float(eta)
        self.beta = float(beta)
        self.gamma = float(gamma)
        self.p_mat = p_mat
        self.t_mat = t_mat
        self.my_mat = my_mat
        self.mz_mat = mz_mat
        self.orient_vals = orient_vals
        if shear_dist is None:
            self.shear_dist = None
        else:
            self.shear_dist = float(shear_dist)
        self.do_rayleigh = do_rayleigh
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.k_init, self.qd, self.alpha1, self.alpha2,
                            self.mu, self.eta, self.beta, self.gamma]
        if getattr(self, 'p_mat') is not None:
            self._parameters += ['-P', self.p_mat.tag]
        if getattr(self, 't_mat') is not None:
            self._parameters += ['-T', self.t_mat.tag]
        if getattr(self, 'my_mat') is not None:
            self._parameters += ['-My', self.my_mat.tag]
        if getattr(self, 'mz_mat') is not None:
            self._parameters += ['-Mz', self.mz_mat.tag]
        if getattr(self, 'orient_vals') is not None:
            self._parameters += ['-orient', *self.orient_vals]
        if getattr(self, 'shear_dist') is not None:
            self._parameters += ['-shearDist', self.shear_dist]
        if getattr(self, 'do_rayleigh'):
            self._parameters += ['-doRayleigh']
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        try:
            self.to_process(osi)
        except ValueError:
            self._parameters[0] = 'ElastomericBearingBoucWen'
            self.to_process(osi)


class Series(UniaxialMaterialBase):
    """
    The Series UniaxialMaterial Class

    This command is used to construct a series material object made up of an arbitrary number of previously-constructed
    UniaxialMaterial objects.
    """
    op_type = 'Series'

    def __init__(self, osi, mats):
        """
        Initial method for Series

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        mats: list
            Identification objects of materials making up the material model

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> mats = [o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None),
        >>>         o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)]
        >>> o3.uniaxial_material.Series(osi, mats=mats)
        """
        self.osi = osi
        self.mats = mats
        self.mat_tags = [x.tag for x in mats]
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, *self.mat_tags]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class PM4Silt(NDMaterialBase):
    """
    The PM4Silt NDMaterial Class

    Code Developed by: **Long Chen** and `pedro <https://www.ce.washington.edu/facultyfinder/pedro-arduino>`_ at
    U.Washington.This command is used to construct a 2-dimensional PM4Silt material.
    """
    op_type = 'PM4Silt'

    def __init__(self, osi, g_o, h_po, den, s_u=None, su_rat=None, su_factor=None, p_atm=None, nu=0.3, n_g=0.75,
                 h0=None, e_init=0.9, lamb=0.06, phicv=32.0, nb_wet=0.8, nb_dry=0.5, nd=0.3, ado=0.8, ru_max=None,
                 zmax=None, cz=100.0, ce=None, cgd=None, ckaf=4.0, m_m=0.01, cg_consol=2.0):
        r"""
        Initial method for PM4Silt

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance

        g_o: float
            Shear modulus constant
        h_po: float
            Contraction rate parameter
        s_u: float
            Undrained shear strength
        su_rat: float
            Undrained shear strength ratio.
        den: float

        su_factor: float
        p_atm: float
            Atmospheric pressure
        nu: float
            Optional: poisson’s ratio. default value is 0.3.
        n_g: float
            Optional: shear modulus exponent. default value is 0.75.
        h0: float
            Optional: variable that adjusts the ratio of plastic modulus to elastic modulus. default value is 0.5.
        e_init: float
            Optional: initial void ratios. default value is 0.90.
        lamb: float
            Optional: the slope of critical state line in e-ln(p) space. default value is 0.060.
        phicv: float
            Optional: critical state effective friction angle. default value is 32 degrees.
        nb_wet: float
            Optional: bounding surface parameter for loose of critical state conditions :math:`1.0 \geq nb_wet \geq
            0.01`. default value is 0.8. in cyclic loading
        nb_dry: float
            Optional: bounding surface parameter for dense of critical state conditions :math:`nb_dry \geq 0`. default
            value is 0.5.
        nd: float
            Optional: dilatancy surface parameter :math:`nd \geq 0`. default value is 0.3.
        ado: float
            Optional: dilatancy parameter. default value is 0.8. with accumulation of fabric
        ru_max: float
            Optional: maximum pore pressure ratio based on p’.
        zmax: None

        cz: float
            Optional: fabric-dilatancy tensor parameter. default value is 100.0.
        ce: float
            Optional: variable that adjusts the rate of strain accumulation in cyclic loading
        cgd: float

        ckaf: float
            Optional: variable that controls the effect that sustained static shear stresses have on plastic modulus.
            default value is 4.0.
        m_m: float
            Optional: yield surface constant (radius of yield surface in stress ratio space). default value is 0.01.
        cg_consol: float
            Optional: reduction factor of elastic modulus for reconsolidation. :math:`cg_consol \geq 1`. default value
            is 2.0.
        """
        self.osi = osi
        if s_u is not None:
            self.s_u = float(s_u)
        else:
            self.s_u = -1.0
        if su_rat is not None:
            self.su_rat = float(su_rat)
        else:
            self.su_rat = -1.0
        self.g_o = float(g_o)
        self.h_po = float(h_po)
        self.den = float(den)
        if su_factor is not None:
            self.su_factor = float(su_factor)
        else:
            self.su_factor = -1.0
        self.p_atm = p_atm
        self.nu = float(nu)
        self.n_g = float(n_g)
        if h0 is None:
            self.h0 = -1.0
        else:
            self.h0 = float(h0)
        self.e_init = float(e_init)
        self.lamb = float(lamb)
        self.phicv = float(phicv)
        self.nb_wet = float(nb_wet)
        self.nb_dry = float(nb_dry)
        self.nd = float(nd)
        self.ado = float(ado)
        if ru_max is None:
            self.ru_max = -1.0
        else:
            self.ru_max = float(ru_max)
        if zmax is None:
            self.zmax = -1.0
        else:
            self.zmax = float(zmax)
        self.cz = float(cz)
        if ce is None:
            self.ce = -1.0
        else:
            self.ce = float(ce)
        if cgd is None:
            self.cgd = -1.0
        else:
            self.cgd = float(cgd)
        self.ckaf = float(ckaf)
        self.m_m = float(m_m)
        self.cg_consol = float(cg_consol)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.s_u, self.su_rat, self.g_o, self.h_po, self.den,
                            self.su_factor, self.p_atm, self.nu, self.n_g, self.h0, self.e_init, self.lamb, self.phicv,
                            self.nb_wet, self.nb_dry, self.nd, self.ado, self.ru_max, self.zmax, self.cz, self.ce,
                            self.cgd, self.ckaf, self.m_m, self.cg_consol]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def update_to_linear(self):
        from o3seespy import update_material_stage
        update_material_stage(self.osi, self, 0)

    def update_to_nonlinear(self):
        from o3seespy import update_material_stage
        update_material_stage(self.osi, self, 1)

    def set_nu(self, nu, ele=None, eles=None):
        from o3seespy import set_parameter
        if ele is not None:
            set_parameter(self.osi, value=nu, eles=[ele], args=['poissonRatio', 1])
        if eles is not None:
            set_parameter(self.osi, value=nu, eles=eles, args=['poissonRatio', 1])

    def set_material_state(self, value, ele=None, eles=None):
        from o3seespy import set_parameter
        if ele is not None:
            set_parameter(self.osi, value=value, eles=[ele], args=['materialState', 1])
        if eles is not None:
            set_parameter(self.osi, value=value, eles=eles, args=['materialState', 1])

    def set_first_call(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, pstr='FirstCall', pval=self.tag, value=value, ele=ele, eles=eles)


class Hysteretic(UniaxialMaterialBase):
    """
    The Hysteretic UniaxialMaterial Class

    This command is used to construct a uniaxial bilinear hysteretic material object with pinching of force and
    deformation, damage due to ductility and energy, and degraded unloading stiffness based on ductility.
    """
    op_type = 'Hysteretic'

    def __init__(self, osi, p1, p2, p3=None, n1=None, n2=None, n3=None, pinch_x=1.0, pinch_y=1.0, damage1=0.0,
                 damage2=0.0, beta=0.0):
        """
        Initial method for Hysteretic

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        p1: list
            ``p1=[s1p, e1p]``, stress and strain (or force & deformation) at first point of the envelope in the positive
            direction
        p2: list
            ``p2=[s2p, e2p]``, stress and strain (or force & deformation) at second point of the envelope in the
            positive direction
        p3: list (default=True), optional
            ``p3=[s3p, e3p]``, stress and strain (or force & deformation) at third point of the envelope in the positive
            direction
        n1: list
            ``n1=[s1n, e1n]``, stress and strain (or force & deformation) at first point of the envelope in the negative
            direction
        n2: list
            ``n2=[s2n, e2n]``, stress and strain (or force & deformation) at second point of the envelope in the
            negative direction
        n3: list (default=True), optional
            ``n3=[s3n, e3n]``, stress and strain (or force & deformation) at third point of the envelope in the negative
            direction
        pinch_x: float
            Pinching factor for strain (or deformation) during reloading
        pinch_y: float
            Pinching factor for stress (or force) during reloading
        damage1: float
            Damage due to ductility: d1(mu-1)
        damage2: float
            Damage due to energy: d2(eii/eult)
        beta: float, optional
            Power used to determine the degraded unloading stiffness based on ductility, mu-beta (optional, default=0.0)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> p1 = [0.5, 0.5]
        >>> p2 = [1.0, 1.0]
        >>> p3 = [0, 1.5]
        >>> n1 = [-0.5, -0.5]
        >>> n2 = [-1.0, -1.0]
        >>> n3 = [0, -1.5]
        >>> o3.uniaxial_material.Hysteretic(osi, p1=p1, p2=p2, p3=p3, n1=n1, n2=n2, n3=n3, pinch_x=1, pinch_y=0, damage1=0, damage2=0)
        """
        self.osi = osi
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        if n1 is None:
            self.n1 = [-p1[0], -p1[1]]
        else:
            self.n1 = n1
        if n2 is None:
            self.n2 = [-p2[0], -p2[1]]
        else:
            self.n2 = n2
        if n3 is None:
            if p3 is None:
                self.n3 = None
            else:
                self.n3 = [-p3[0], -p3[1]]
        else:
            self.n3 = n3
        self.pinch_x = float(pinch_x)
        self.pinch_y = float(pinch_y)
        self.damage1 = float(damage1)
        self.damage2 = float(damage2)
        self.beta = float(beta)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        if self.p3 is None:
            self._parameters = [self.op_type, self._tag, *self.p1, *self.p2, *self.n1, *self.n2]
        else:
            self._parameters = [self.op_type, self._tag, *self.p1, *self.p2, *self.p3, *self.n1, *self.n2, *self.n3]
        self._parameters += [self.pinch_x, self.pinch_y, self.damage1, self.damage2, self.beta]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class Joint2D(ElementBase):
    """
    The Joint2D Element Class

    This command is used to construct a two-dimensional beam-column-joint element object. The two dimensional
    beam-column joint is idealized as a parallelogram shaped shear panel with adjacent elements connected to its
    mid-points. The midpoints of the parallelogram are referred to as external nodes. These nodes are the only
    analysis components that connect the joint element to the surrounding structure.


    """
    op_type = 'Joint2D'

    def __init__(self, osi, ele_nodes, mat1, mat2, mat3, mat4, mat_c, lrg_dsp, dmg, dmg_vals=None):
        """
        Initial method for Joint2D

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of five element nodes = ``[nd1,nd2,nd3,nd4,ndc]``. ``ndc`` is the central node of beam-column joint.
            (the object ``ndc`` is used to generate the internal node, thus, the node should not exist in the domain or be used by
            any other node)
        mat1: int
            Uniaxial material object for interface rotational spring at node 1. use a zero object to indicate the case
            that a beam-column element is rigidly framed to the joint.
        mat2: int
            Uniaxial material object for interface rotational spring at node 2. use a zero object to indicate the case
            that a beam-column element is rigidly framed to the joint.
        mat3: int
            Uniaxial material object for interface rotational spring at node 3. use a zero object to indicate the case
            that a beam-column element is rigidly framed to the joint.
        mat4: int
            Uniaxial material object for interface rotational spring at node 4. use a zero object to indicate the case
            that a beam-column element is rigidly framed to the joint.
        mat_c: int
            Uniaxial material object for rotational spring of the central node that describes shear panel behavior
        lrg_dsp: obj
            An integer indicating the flag for considering large deformations: * ``0`` - for small deformations and
            constant geometry * ``1`` - for large deformations and time varying geometry * ``2`` - for large deformations
            ,time varying geometry and length correction
        dmg: obj
            Damage model object
        dmg1dmg2dmg3dmg4dmg_c: None, optional


        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [1, 0], [1, 1], [0, 1], [0.5, 0.5]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
        >>> o3.element.Joint2D(osi, ele_nodes=ele_nodes, mat1=1, mat2=1, mat3=1, mat4=1, mat_c=1, lrg_dsp='', dmg='', dmg1dmg2dmg3dmg4dmg_c=1)
        """
        self.osi = osi
        self.ele_node_tags = [x.tag for x in ele_nodes]
        self.ele_nodes = ele_nodes
        self.mat1 = int(mat1)
        self.mat2 = int(mat2)
        self.mat3 = int(mat3)
        self.mat4 = int(mat4)
        self.mat_c = int(mat_c)
        self.lrg_dsp = lrg_dsp
        self.dmg = dmg
        self.dmg_vals = dmg_vals
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_node_tags, self.mat1, self.mat2, self.mat3, self.mat4,
                            self.mat_c, self.lrg_dsp.tag]
        if self.dmg is not None:
            self._parameters += ['-damage', self.dmg.tag]
        if getattr(self, 'dmg_vals') is not None:
            self._parameters += ['-damage', self.dmg_vals]
        self.to_process(osi)