from o3seespy.command.nd_material.base_material import NDMaterialBase


class ElasticIsotropic(NDMaterialBase):
    """
    The ElasticIsotropic NDMaterial Class
    
    This command is used to construct an ElasticIsotropic material object.
    """
    op_type = 'ElasticIsotropic'

    def __init__(self, osi, e_mod, nu, rho=0.0):
        """
        Initial method for ElasticIsotropic

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        e_mod: float
            Elastic modulus
        nu: float
            Poisson's ratio
        rho: float, optional
            Mass density 

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.nd_material.ElasticIsotropic(osi, e_mod=1.0, nu=1.0, rho=0.0)
        """
        self.osi = osi
        self.e_mod = float(e_mod)
        self.nu = float(nu)
        self.rho = float(rho)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.e_mod, self.nu, self.rho]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def set_e_mod(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'E', value, ele, eles)

    def set_v(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'v', value, ele, eles)

    def set_rho(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'rho', value, ele, eles)


class ElasticOrthotropic(NDMaterialBase):
    """
    The ElasticOrthotropic NDMaterial Class
    
    This command is used to construct an ElasticOrthotropic material object.
    """
    op_type = 'ElasticOrthotropic'

    def __init__(self, osi, ex, ey, ez, nu_xy, nu_yz, nu_zx, gxy, gyz, gzx, rho=0.0):
        """
        Initial method for ElasticOrthotropic

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ex: float
            Elastic modulus in x direction
        ey: float
            Elastic modulus in y direction
        ez: float
            Elastic modulus in z direction
        nu_xy: float
            Poisson's ratios in x and y plane
        nu_yz: float
            Poisson's ratios in y and z plane
        nu_zx: float
            Poisson's ratios in z and x plane
        gxy: float
            Shear modulii in x and y plane
        gyz: float
            Shear modulii in y and z plane
        gzx: float
            Shear modulii in z and x plane
        rho: float, optional
            Mass density 

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.nd_material.ElasticOrthotropic(osi, ex=1.0, ey=1.0, ez=1.0, nu_xy=1.0, nu_yz=1.0, nu_zx=1.0, gxy=1.0, gyz=1.0, gzx=1.0, rho=0.0)
        """
        self.osi = osi
        self.ex = float(ex)
        self.ey = float(ey)
        self.ez = float(ez)
        self.nu_xy = float(nu_xy)
        self.nu_yz = float(nu_yz)
        self.nu_zx = float(nu_zx)
        self.gxy = float(gxy)
        self.gyz = float(gyz)
        self.gzx = float(gzx)
        self.rho = float(rho)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.ex, self.ey, self.ez, self.nu_xy, self.nu_yz, self.nu_zx, self.gxy, self.gyz, self.gzx, self.rho]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def set_e_mod_x(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'Ex', value, ele, eles)

    def set_e_mod_y(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'Ey', value, ele, eles)

    def set_e_mod_z(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'Ez', value, ele, eles)

    def set_vyx(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'vyx', value, ele, eles)

    def set_vzy(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'vzy', value, ele, eles)

    def set_vxz(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'vxz', value, ele, eles)

    def set_g_mod_yx(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'Gyx', value, ele, eles)

    def set_g_mod_zy(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'Gzy', value, ele, eles)

    def set_g_mod_xz(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'Gxz', value, ele, eles)

    def set_rho(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'rho', value, ele, eles)


class J2Plasticity(NDMaterialBase):
    """
    The J2Plasticity NDMaterial Class
    
    This command is used to construct an multi dimensional material object that has a von Mises (J2) yield criterium and
    isotropic hardening.
    """
    op_type = 'J2Plasticity'

    def __init__(self, osi, k_mod, g_mod, sig0, sig_inf, delta, big_h):
        """
        Initial method for J2Plasticity

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        k_mod: float
            Bulk modulus
        g_mod: float
            Shear modulus
        sig0: float
            Initial yield stress
        sig_inf: float
            Final saturation yield stress
        delta: float
            Exponential hardening parameter
        big_h: float
            Linear hardening parameter

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.nd_material.J2Plasticity(osi, k_mod=1.0, g_mod=1.0, sig0=1.0, sig_inf=1.0, delta=1.0, big_h=1.0)
        """
        self.osi = osi
        self.k_mod = float(k_mod)
        self.g_mod = float(g_mod)
        self.sig0 = float(sig0)
        self.sig_inf = float(sig_inf)
        self.delta = float(delta)
        self.big_h = float(big_h)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.k_mod, self.g_mod, self.sig0, self.sig_inf, self.delta, self.big_h]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def set_k(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'K', value, ele, eles)

    def set_mu(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'mu', value, ele, eles)

    def set_rho(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'rho', value, ele, eles)


class DruckerPrager(NDMaterialBase):
    """
    The DruckerPrager NDMaterial Class
    
    This command is used to construct an multi dimensional material object that has a Drucker-Prager yield criterium.
    """
    op_type = 'DruckerPrager'

    def __init__(self, osi, k_mod, g_mod, sigma_y, rho, rho_bar, kinf, ko, delta1, delta2, big_h, theta, density, atm_pressure=101e3):
        r"""
        Initial method for DruckerPrager

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        k_mod: float
            Bulk modulus
        g_mod: float
            Shear modulus
        sigma_y: float
            Yield stress
        rho: float
            Frictional strength parameter
        rho_bar: float
            Controls evolution of plastic volume change, :math:`0\le rhobar \le rho`.
        kinf: float
            Nonlinear isotropic strain hardening parameter, :math:`kinf \ge 0`.
        ko: float
            Nonlinear isotropic strain hardening parameter, :math:`ko \ge 0`.
        delta1: float
            Nonlinear isotropic strain hardening parameter, :math:`delta1\ge 0`.
        delta2: float
            Tension softening parameter, :math:`delta2\ge 0`.
        big_h: float
            Linear hardening parameter, :math:`h \ge 0`.
        theta: float
            Controls relative proportions of isotropic and kinematic hardening, :math:`0 \le theta \le 1`.
        density: float
            Mass density of the material
        atm_pressure: float, optional
            Optional atmospheric pressure for update of elastic bulk and shear moduli

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.nd_material.DruckerPrager(osi, k_mod=1.0, g_mod=1.0, sigma_y=1.0, rho=1.0, rho_bar=1.0, kinf=1.0, ko=1.0, delta1=1.0, delta2=1.0, big_h=1.0, theta=1.0, density=1.0, atm_pressure=101e3)
        """
        self.osi = osi
        self.k_mod = float(k_mod)
        self.g_mod = float(g_mod)
        self.sigma_y = float(sigma_y)
        self.rho = float(rho)
        self.rho_bar = float(rho_bar)
        self.kinf = float(kinf)
        self.ko = float(ko)
        self.delta1 = float(delta1)
        self.delta2 = float(delta2)
        self.big_h = float(big_h)
        self.theta = float(theta)
        self.density = float(density)
        self.atm_pressure = float(atm_pressure)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.k_mod, self.g_mod, self.sigma_y, self.rho, self.rho_bar, self.kinf, self.ko, self.delta1, self.delta2, self.big_h, self.theta, self.density, self.atm_pressure]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def set_material_state(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'materialState', value, ele, eles)

    def set_frictional_strength(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'frictionalStrength', value, ele, eles)

    def set_nonassociative_term(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'nonassociativeTerm', value, ele, eles)

    def set_cohesive_intercept(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'cohesiveIntercept', value, ele, eles)

    def set_g_mod(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'shearModulus', value, ele, eles)

    def set_bulk_mod(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'bulkModulus', value, ele, eles)

    def set_update_material_stage(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'updateMaterialStage', value, ele, eles)


class Damage2p(NDMaterialBase):
    """
    The Damage2p NDMaterial Class
    
    This command is used to construct a three-dimensional material object that has a Drucker-Prager plasticity model
    coupled with a two-parameter damage model.
    """
    op_type = 'Damage2p'

    def __init__(self, osi, fcc, fct: float=None, e_mod: float=None, ni: float=None, gt: float=None, gc: float=None, rho_bar: float=None, big_h: float=None, theta: float=None, tangent: float=None):
        r"""
        Initial method for Damage2p

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        fcc: float
            Concrete compressive strength, negative real value (positive input is changed in sign automatically)
        fct: float, optional
            Optional concrete tensile strength, positive real value (for concrete like materials is less than fcc),
            :math:`0.1*abs(fcc)` = :math:`4750*sqrt(abs(fcc))\text{ }if\text{ }abs(fcc)<2000` because fcc is assumed in mpa
            (see aci 318)
        e_mod: float, optional
            Optional young modulus, :math:`57000*sqrt(abs(fcc))` if :math:`abs(fcc)>2000` because fcc is assumed in psi
            (see aci 318)
        ni: float, optional
            Optional poisson coefficient, 0.15 (from comparison with tests by kupfer hilsdorf rusch 1969)
        gt: float, optional
            Optional tension fracture energy density, positive real value (integral of the stress-strain envelope in
            tension), :math:`1840*fct*fct/e` (from comparison with tests by gopalaratnam and shah 1985)
        gc: float, optional
            Optional compression fracture energy density, positive real value (integral of the stress-strain envelope
            after the peak in compression), :math:6250*fcc*fcc/e` (from comparison with tests by karsan and jirsa 1969)
        rho_bar: float, optional
            Optional parameter of plastic volume change, positive real value :math:`0=rhobar< sqrt(2/3)`, 0.2 (from
            comparison with tests by kupfer hilsdorf rusch 1969)
        big_h: float, optional
            Optional linear hardening parameter for plasticity, positive real value (usually less than e),
            :math:`0.25*e` (from comparison with tests by karsan and jirsa 1969 and gopalaratnam and shah 1985)
        theta: float, optional
            Optional ratio between isotropic and kinematic hardening, positive real value :math:`0=theta=1` (with: 0
            hardening kinematic only and 1 hardening isotropic only, 0.5 (from comparison with tests by karsan and jirsa 1969
            and gopalaratnam and shah 1985)
        tangent: float, optional
            Optional integer to choose the computational stiffness matrix, 0: computational tangent; 1: damaged secant
            stiffness (hint: in case of strong nonlinearities use it with krylov-newton algorithm)

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.nd_material.Damage2p(osi, fcc=1.0, fct=1.0, e_mod=1.0, ni=1.0, gt=1.0, gc=1.0, rho_bar=1.0, big_h=1.0, theta=1.0, tangent=1.0)
        """
        self.osi = osi
        self.fcc = float(fcc)
        if fct is None:
            self.fct = None
        else:
            self.fct = float(fct)
        if e_mod is None:
            self.e_mod = None
        else:
            self.e_mod = float(e_mod)
        if ni is None:
            self.ni = None
        else:
            self.ni = float(ni)
        if gt is None:
            self.gt = None
        else:
            self.gt = float(gt)
        if gc is None:
            self.gc = None
        else:
            self.gc = float(gc)
        if rho_bar is None:
            self.rho_bar = None
        else:
            self.rho_bar = float(rho_bar)
        if big_h is None:
            self.big_h = None
        else:
            self.big_h = float(big_h)
        if theta is None:
            self.theta = None
        else:
            self.theta = float(theta)
        if tangent is None:
            self.tangent = None
        else:
            self.tangent = float(tangent)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.fcc]
        if getattr(self, 'fct') is not None:
            self._parameters += ['-fct', self.fct]
        if getattr(self, 'e_mod') is not None:
            self._parameters += ['-E', self.e_mod]
        if getattr(self, 'ni') is not None:
            self._parameters += ['-ni', self.ni]
        if getattr(self, 'gt') is not None:
            self._parameters += ['-Gt', self.gt]
        if getattr(self, 'gc') is not None:
            self._parameters += ['-Gc', self.gc]
        if getattr(self, 'rho_bar') is not None:
            self._parameters += ['-rho_bar', self.rho_bar]
        if getattr(self, 'big_h') is not None:
            self._parameters += ['-H', self.big_h]
        if getattr(self, 'theta') is not None:
            self._parameters += ['-theta', self.theta]
        if getattr(self, 'tangent') is not None:
            self._parameters += ['-tangent', self.tangent]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class PlaneStress(NDMaterialBase):
    """
    The PlaneStress NDMaterial Class
    
    This command is used to construct a plane-stress material wrapper which converts any three-dimensional material into
    a plane stress material via static condensation.
    """
    op_type = 'PlaneStress'

    def __init__(self, osi, mat3d):
        """
        Initial method for PlaneStress

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        mat3d: obj
            Object of perviously defined 3d ndmaterial material

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=3)
        >>> mat_3d = o3.nd_material.ElasticIsotropic(osi, e_mod=1.0, nu=1.0, rho=0.0)
        >>> o3.nd_material.PlaneStress(osi, mat3d=mat_3d)
        """
        self.osi = osi
        self.mat3d = mat3d
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.mat3d.tag]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def set_tangent(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'Tangent', value, ele, eles)

    def set_stress(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'tangent', value, ele, eles)

    def set_stresses(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'stress', value, ele, eles)

    def set_strain(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'stresses', value, ele, eles)

    def set_strains(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'strain', value, ele, eles)


class PlaneStrain(NDMaterialBase):
    """
    The PlaneStrain NDMaterial Class
    
    This command is used to construct a plane-stress material wrapper which converts any three-dimensional material into
    a plane strain material by imposing plain strain conditions on the three-dimensional material.
    """
    op_type = 'PlaneStrain'

    def __init__(self, osi, mat3d):
        """
        Initial method for PlaneStrain

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        mat3d: obj
            Integer object of previously defined 3d ndmaterial material

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> mat_3d = o3.nd_material.ElasticIsotropic(osi, e_mod=1.0, nu=1.0, rho=0.0)
        >>> o3.nd_material.PlaneStrain(osi, mat3d=mat_3d)
        """
        self.osi = osi
        self.mat3d = mat3d
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.mat3d.tag]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def set_tangent(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'Tangent', value, ele, eles)

    def set_stress(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'tangent', value, ele, eles)

    def set_stresses(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'stress', value, ele, eles)

    def set_strain(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'stresses', value, ele, eles)

    def set_strains(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'strain', value, ele, eles)


class MultiaxialCyclicPlasticity(NDMaterialBase):
    """
    The MultiaxialCyclicPlasticity NDMaterial Class
    
    This command is used to construct an multiaxial Cyclic Plasticity model for clays
    """
    op_type = 'MultiaxialCyclicPlasticity'

    def __init__(self, osi, rho, k_mod, g_mod, su, ho, h, m, beta, k_coeff):
        r"""
        Initial method for MultiaxialCyclicPlasticity

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        rho: float
            Density
        k_mod: float
            Buck modulus
        g_mod: float
            Maximum (small strain) shear modulus
        su: float
            Undrained shear strength, size of bounding surface :math:`r=\sqrt{8/3}*su`
        ho: float
            Linear kinematic hardening modulus of bounding surface
        h: float
            Hardening parameter
        m: float
            Hardening parameter
        beta: float
            Integration parameter, usually beta=0.5
        k_coeff: float
            Coefficient of earth pressure, k0

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.nd_material.MultiaxialCyclicPlasticity(osi, rho=1.0, k_mod=1.0, g_mod=1.0, su=1.0, ho=1.0, h=1.0, m=1.0, beta=1.0, k_coeff=1.0)
        """
        self.osi = osi
        self.rho = float(rho)
        self.k_mod = float(k_mod)
        self.g_mod = float(g_mod)
        self.su = float(su)
        self.ho = float(ho)
        self.h = float(h)
        self.m = float(m)
        self.beta = float(beta)
        self.k_coeff = float(k_coeff)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.rho, self.k_mod, self.g_mod, self.su, self.ho, self.h, self.m, self.beta, self.k_coeff]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class BoundingCamClay(NDMaterialBase):
    """
    The BoundingCamClay NDMaterial Class
    
    This command is used to construct a multi-dimensional bounding surface Cam Clay material object after Borja et al.
    (2001).
    """
    op_type = 'BoundingCamClay'

    def __init__(self, osi, mass_density, big_c, bulk_mod, ocr, mu_o, alpha, lamb, h, m):
        """
        Initial method for BoundingCamClay

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        mass_density: float
            Mass density
        big_c: float
            Ellipsoidal axis ratio (defines shape of ellipsoidal loading/bounding surfaces)
        bulk_mod: float
            Initial bulk modulus
        ocr: float
            Overconsolidation ratio
        mu_o: float
            Initial shear modulus
        alpha: float
            Pressure-dependency parameter for modulii (greater than or equal to zero)
        lamb: float
            Soil compressibility index for virgin loading
        h: float
            Hardening parameter for plastic response inside of bounding surface (if h = 0, no hardening)
        m: float
            Hardening parameter (exponent) for plastic response inside of bounding surface (if m = 0, only linear
            hardening)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.nd_material.BoundingCamClay(osi, mass_density=1.0, big_c=1.0, bulk_mod=1.0, ocr=1.0, mu_o=1.0, alpha=1.0, lamb=1.0, h=1.0, m=1.0)
        """
        self.osi = osi
        self.mass_density = float(mass_density)
        self.big_c = float(big_c)
        self.bulk_mod = float(bulk_mod)
        self.ocr = float(ocr)
        self.mu_o = float(mu_o)
        self.alpha = float(alpha)
        self.lamb = float(lamb)
        self.h = float(h)
        self.m = float(m)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.mass_density, self.big_c, self.bulk_mod, self.ocr, self.mu_o, self.alpha, self.lamb, self.h, self.m]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def set_material_state(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'materialState', value, ele, eles)


class PlateFiber(NDMaterialBase):
    """
    The PlateFiber NDMaterial Class
    
    This command is used to construct a plate-fiber material wrapper which converts any three-dimensional material into
    a plate fiber material (by static condensation) appropriate for shell analysis.
    """
    op_type = 'PlateFiber'

    def __init__(self, osi, three_d):
        """
        Initial method for PlateFiber

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        three_d: obj
            Material object for a previously-defined three-dimensional material

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=3)
        >>> mat = o3.nd_material.ElasticIsotropic(osi, e_mod=1.0, nu=1.0, rho=0.0)
        >>> o3.nd_material.PlateFiber(osi, three_d=mat)
        """
        self.osi = osi
        self.three_d = three_d
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.three_d.tag]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def set_tangent(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'Tangent', value, ele, eles)

    def set_stress(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'tangent', value, ele, eles)

    def set_stresses(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'stress', value, ele, eles)

    def set_strain(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'stresses', value, ele, eles)

    def set_strains(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'strain', value, ele, eles)


class FSAM(NDMaterialBase):
    """
    The FSAM NDMaterial Class
    
    This command is used to construct a nDMaterial FSAM (Fixed-Strut-Angle-Model, Figure 1, Kolozvari et al., 2015),
    which is a plane-stress constitutive model for simulating the behavior of RC panel elements under generalized,
    in-plane, reversed-cyclic loading conditions (Ulugtekin, 2010; Orakcal et al., 2012). In the FSAM
    constitutive model, the strain fields acting on concrete and reinforcing steel components of a
    RC panel are assumed to be equal to each other, implying perfect bond assumption between
    concrete and reinforcing steel bars. While the reinforcing steel bars develop uniaxial
    stresses under strains in their longitudinal direction, the behavior of concrete is
    defined using stress-strain relationships in biaxial directions, the orientation
    of which is governed by the state of cracking in concrete. Although the
    concrete stress-strain relationship used in the FSAM is fundamentally
    uniaxial in nature, it also incorporates biaxial softening effects
    including compression softening and biaxial damage. For transfer
    of shear stresses across the cracks, a friction-based
    elasto-plastic shear aggregate interlock model is
    adopted, together with a linear elastic model
    for representing dowel action on the
    reinforcing steel bars (Kolozvari,
    2013). Note that FSAM
    constitutive model
    is implemented to
    be used with Shear-Flexure Interaction model for RC walls (SFI_MVLEM), but it could be also used elsewhere.
    """
    op_type = 'FSAM'

    def __init__(self, osi, rho, s_x, s_y, conc, rou_x, rou_y, nu, alfadow):
        r"""
        Initial method for FSAM

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        rho: float
            Material density
        s_x: obj
            Object of uniaxial_material simulating horizontal (x) reinforcement
        s_y: obj
            Object of uniaxial_material simulating vertical (y) reinforcement
        conc: obj
            Object of uniaxial_material simulating concrete, shall be used with uniaxial_material concretecm
        rou_x: float
            Reinforcing ratio in horizontal (x) direction (:math:`roux = _{s,x}/a_{gross,x}`)
        rou_y: float
            Reinforcing ratio in vertical (x) direction (:math:`rouy = _{s,y}/a_{gross,y}`)
        nu: float
            Concrete friction coefficient (:math:`0.0 < \nu < 1.5`)
        alfadow: float
            Stiffness coefficient of reinforcement dowel action (:math:`0.0 < alfadow < 0.05`)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=3)
        >>> s_x = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
        >>> s_y = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
        >>> conc = o3.uniaxial_material.ConcreteCM(osi, fpcc=1.0, epcc=1.0, ec=1.0, rc=1.0, xcrn=1.0, ft=1.0, et=1.0, rt=1.0, xcrp=1.0,
        >>>                                 gap_close=0)
        >>> o3.nd_material.FSAM(osi, rho=1.0, s_x=s_x, s_y=s_y, conc=conc, rou_x=1.0, rou_y=1.0, nu=1.0, alfadow=1.0)
        """
        self.osi = osi
        self.rho = float(rho)
        self.s_x = s_x
        self.s_y = s_y
        self.conc = conc
        self.rou_x = float(rou_x)
        self.rou_y = float(rou_y)
        self.nu = float(nu)
        self.alfadow = float(alfadow)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.rho, self.s_x.tag, self.s_y.tag, self.conc.tag, self.rou_x, self.rou_y, self.nu, self.alfadow]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class ManzariDafalias(NDMaterialBase):
    """
    The ManzariDafalias NDMaterial Class
    
    This command is used to construct a multi-dimensional Manzari-Dafalias(2004) material.
    """
    op_type = 'ManzariDafalias'

    def __init__(self, osi, g0, nu, e_init, mc, c, lambda_c, e0, ksi, p_atm, m, h0, ch, nb, a0, nd, z_max, cz, den):
        r"""
        Initial method for ManzariDafalias

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        g0: float
            Shear modulus constant
        nu: float
            Poisson ratio
        e_init: float
            Initial void ratio
        mc: float
            Critical state stress ratio
        c: float
            Ratio of critical state stress ratio in extension and compression
        lambda_c: float
            Critical state line constant
        e0: float
            Critical void ratio at p = 0
        ksi: float
            Critical state line constant
        p_atm: float
            Atmospheric pressure
        m: float
            Yield surface constant (radius of yield surface in stress ratio space)
        h0: float
            Constant parameter
        ch: float
            Constant parameter
        nb: float
            Bounding surface parameter, :math:`nb \ge 0`
        a0: float
            Dilatancy parameter
        nd: float
            Dilatancy surface parameter :math:`nd \ge 0`
        z_max: float
            Fabric-dilatancy tensor parameter
        cz: float
            Fabric-dilatancy tensor parameter
        den: float
            Mass density of the material

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.nd_material.ManzariDafalias(osi, g0=1.0, nu=1.0, e_init=1.0, mc=1.0, c=1.0, lambda_c=1.0, e0=1.0, ksi=1.0, p_atm=1.0, m=1.0, h0=1.0, ch=1.0, nb=1.0, a0=1.0, nd=1.0, z_max=1.0, cz=1.0, den=1.0)
        """
        self.osi = osi
        self.g0 = float(g0)
        self.nu = float(nu)
        self.e_init = float(e_init)
        self.mc = float(mc)
        self.c = float(c)
        self.lambda_c = float(lambda_c)
        self.e0 = float(e0)
        self.ksi = float(ksi)
        self.p_atm = float(p_atm)
        self.m = float(m)
        self.h0 = float(h0)
        self.ch = float(ch)
        self.nb = float(nb)
        self.a0 = float(a0)
        self.nd = float(nd)
        self.z_max = float(z_max)
        self.cz = float(cz)
        self.den = float(den)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.g0, self.nu, self.e_init, self.mc, self.c, self.lambda_c, self.e0, self.ksi, self.p_atm, self.m, self.h0, self.ch, self.nb, self.a0, self.nd, self.z_max, self.cz, self.den]
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

    def set_g_mod(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'ShearModulus', value, ele, eles)

    def set_nu(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'poissonRatio', value, ele, eles)

    def set_void_ratio(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'voidRatio', value, ele, eles)

    def set_stress_correction(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'stressCorrection', value, ele, eles)


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



class AcousticMedium(NDMaterialBase):
    """
    The AcousticMedium NDMaterial Class
    
    This command is used to construct an acoustic medium NDMaterial object.
    """
    op_type = 'AcousticMedium'

    def __init__(self, osi, k_mod, rho):
        """
        Initial method for AcousticMedium

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        k_mod: float
            Bulk module of the acoustic medium
        rho: float
            Mass density of the acoustic medium

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.nd_material.AcousticMedium(osi, k_mod=1.0, rho=1.0)
        """
        self.osi = osi
        self.k_mod = float(k_mod)
        self.rho = float(rho)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.k_mod, self.rho]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def set_kf(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'Kf', value, ele, eles)

    def set_rho(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'rho', value, ele, eles)

    def set_gamma(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'gamma', value, ele, eles)
