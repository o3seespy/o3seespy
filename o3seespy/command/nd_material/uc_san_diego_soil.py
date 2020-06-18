from o3seespy.command.nd_material.base_material import NDMaterialBase


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
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat

        self._parameters = [self.op_type, self._tag, self.nd, self.rho, self.g_mod_ref, self.bulk_mod_ref,
                            self.cohesion, self.peak_strain, self.phi, self.p_ref, self.d]

        if self.yield_surf is not None:
            self._parameters.append(self.n_surf)  # from docs 'add a minus sign in front of noYieldSurf'
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

    def set_nu(self, nu, ele=None, eles=None, adj_g_mod=False):
        if adj_g_mod:
            g_mod = 3 * self.bulk_mod_ref * (1 - 2 * nu) / (2 * (1 + nu))
            self.set_parameter(self.osi, 'shearModulus', g_mod, ele, eles)
        else:
            bulk_mod = 2 * self.g_mod_ref * (1 + nu) / (3 * (1 - 2 * nu))
            self.set_parameter(self.osi, 'bulkModulus', bulk_mod, ele, eles)



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



class PressureDependMultiYield03(NDMaterialBase):
    """
    The PressureDependMultiYield03 NDMaterial Class
    
    The reference for PressureDependMultiYield03 material: Khosravifar, A., Elgamal, A., Lu, J., and Li, J. [2018]. "A
    3D model for earthquake-induced liquefaction triggering and post-liquefaction response." Soil Dynamics and Earthquake
    Engineering, 110, 43-52)PressureDependMultiYield03 is modified from PressureDependMultiYield02 material to comply
    with the established guidelines on the dependence of liquefaction triggering to the number of loading cycles,
    effective overburden stress (Kσ), and static shear stress (Kα).The explanations of parametersSee `notes
    <http://opensees.berkeley.edu/wiki/index.php/PressureDependMultiYield02_Material>`_
    """
    op_type = 'PressureDependMultiYield03'

    def __init__(self, osi, nd, rho, ref_shear_modul, ref_bulk_modul, friction_ang, peak_shear_stra, ref_press, press_depend_coe, pt_ang, ca, cb, cc, cd, ce, da, db, dc, no_yield_surf=20.0, yield_surf: float=None, liquefac1=1, liquefac2=0., pa=101, s0=1.73):
        """
        Initial method for PressureDependMultiYield03

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        nd: None
            
        rho: None
            
        ref_shear_modul: None
            
        ref_bulk_modul: None
            
        friction_ang: None
            
        peak_shear_stra: None
            
        ref_press: None
            
        press_depend_coe: None
            
        pt_ang: None
            
        ca: None
            
        cb: None
            
        cc: None
            
        cd: None
            
        ce: None
            
        da: None
            
        db: None
            
        dc: None
            
        no_yield_surf: None, optional
            
        yield_surf: None (default=True), optional
            
        liquefac1: None, optional
            
        liquefac2: None, optional
            
        pa: None, optional
            
        s0: None, optional
            
        """
        self.osi = osi
        self.nd = nd
        self.rho = rho
        self.ref_shear_modul = ref_shear_modul
        self.ref_bulk_modul = ref_bulk_modul
        self.friction_ang = friction_ang
        self.peak_shear_stra = peak_shear_stra
        self.ref_press = ref_press
        self.press_depend_coe = press_depend_coe
        self.pt_ang = pt_ang
        self.ca = ca
        self.cb = cb
        self.cc = cc
        self.cd = cd
        self.ce = ce
        self.da = da
        self.db = db
        self.dc = dc
        self.no_yield_surf = no_yield_surf
        self.yield_surf = yield_surf
        self.liquefac1 = liquefac1
        self.liquefac2 = liquefac2
        self.pa = pa
        self.s0 = s0
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.nd, self.rho, self.ref_shear_modul, self.ref_bulk_modul, self.friction_ang, self.peak_shear_stra, self.ref_press, self.press_depend_coe, self.pt_ang, self.ca, self.cb, self.cc, self.cd, self.ce, self.da, self.db, self.dc, self.no_yield_surf]
        special_pms = ['yield_surf', 'liquefac1', 'liquefac2', 'pa', 's0']
        packets = [True, False, False, False, False]
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

    def set_update_material_stage(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'updateMaterialStage', value, ele, eles)

    def set_g_mod(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'shearModulus', value, ele, eles)

    def set_bulk_mod(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'bulkModulus', value, ele, eles)

    def set_friction_angle(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'updateMaterialStage', value, ele, eles)

    def set_cohesion(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'shearModulus', value, ele, eles)
