from o3seespy.command.nd_material.base_material import NDMaterialBase


class ElasticIsotropic(NDMaterialBase):
    """
    The ElasticIsotropic NDMaterial Class
    
    This command is used to construct an ElasticIsotropic material object.
    """
    op_type = 'ElasticIsotropic'

    def __init__(self, osi, e_mod, v, rho=0.0):
        """
        Initial method for ElasticIsotropic

        Parameters
        ----------
        e_mod: float
            Elastic modulus
        v: float
            Poisson's ratio
        rho: float
            Mass density (optional)
        """
        self.e_mod = float(e_mod)
        self.v = float(v)
        self.rho = float(rho)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.e_mod, self.v, self.rho]
        self.to_process(osi)


class ElasticOrthotropic(NDMaterialBase):
    """
    The ElasticOrthotropic NDMaterial Class
    
    This command is used to construct an ElasticOrthotropic material object.
    """
    op_type = 'ElasticOrthotropic'

    def __init__(self, osi, ex, ey, ez, vxy, vyz, vzx, gxy, gyz, gzx, rho=0.0):
        """
        Initial method for ElasticOrthotropic

        Parameters
        ----------
        ex: float
            Elastic modulus in x direction
        ey: float
            Elastic modulus in y direction
        ez: float
            Elastic modulus in z direction
        vxy: float
            Poisson's ratios in x and y plane
        vyz: float
            Poisson's ratios in y and z plane
        vzx: float
            Poisson's ratios in z and x plane
        gxy: float
            Shear modulii in x and y plane
        gyz: float
            Shear modulii in y and z plane
        gzx: float
            Shear modulii in z and x plane
        rho: float
            Mass density (optional)
        """
        self.ex = float(ex)
        self.ey = float(ey)
        self.ez = float(ez)
        self.vxy = float(vxy)
        self.vyz = float(vyz)
        self.vzx = float(vzx)
        self.gxy = float(gxy)
        self.gyz = float(gyz)
        self.gzx = float(gzx)
        self.rho = float(rho)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.ex, self.ey, self.ez, self.vxy, self.vyz, self.vzx, self.gxy, self.gyz, self.gzx, self.rho]
        self.to_process(osi)


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
        """
        self.k_mod = float(k_mod)
        self.g_mod = float(g_mod)
        self.sig0 = float(sig0)
        self.sig_inf = float(sig_inf)
        self.delta = float(delta)
        self.big_h = float(big_h)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.k_mod, self.g_mod, self.sig0, self.sig_inf, self.delta, self.big_h]
        self.to_process(osi)


class DrukerPrager(NDMaterialBase):
    """
    The DrukerPrager NDMaterial Class
    
    This command is used to construct an multi dimensional material object that has a Drucker-Prager yield criterium.
    """
    op_type = 'DrukerPrager'

    def __init__(self, osi, k_mod, g_mod, sigma_y, rho, rho_bar, kinf, ko, delta1, delta2, big_h, theta, density, atm_pressure=101e3):
        """
        Initial method for DrukerPrager

        Parameters
        ----------
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
        atm_pressure: float
            Optional atmospheric pressure for update of elastic bulk and shear moduli
        """
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
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.k_mod, self.g_mod, self.sigma_y, self.rho, self.rho_bar, self.kinf, self.ko, self.delta1, self.delta2, self.big_h, self.theta, self.density, self.atm_pressure]
        self.to_process(osi)


class Damage2p(NDMaterialBase):
    """
    The Damage2p NDMaterial Class
    
    This command is used to construct a three-dimensional material object that has a Drucker-Prager plasticity model
    coupled with a two-parameter damage model.
    """
    op_type = 'Damage2p'

    def __init__(self, osi, fcc, fct: float=None, e_mod: float=None, ni: float=None, gt: float=None, gc: float=None, rho_bar: float=None, big_h: float=None, theta: float=None, tangent: float=None):
        """
        Initial method for Damage2p

        Parameters
        ----------
        fcc: float
            Concrete compressive strength, negative real value (positive input is changed in sign automatically)
        fct: float
            Optional concrete tensile strength, positive real value (for concrete like materials is less than fcc),
            :math:`0.1*abs(fcc)` = :math:`4750*sqrt(abs(fcc))\text{ }if\text{ }abs(fcc)<2000` because fcc is assumed in mpa
            (see aci 318)
        e_mod: float
            Optional young modulus, :math:`57000*sqrt(abs(fcc))` if :math:`abs(fcc)>2000` because fcc is assumed in psi
            (see aci 318)
        ni: float
            Optional poisson coefficient, 0.15 (from comparison with tests by kupfer hilsdorf rusch 1969)
        gt: float
            Optional tension fracture energy density, positive real value (integral of the stress-strain envelope in
            tension), :math:`1840*fct*fct/e` (from comparison with tests by gopalaratnam and shah 1985)
        gc: float
            Optional compression fracture energy density, positive real value (integral of the stress-strain envelope
            after the peak in compression), :math:6250*fcc*fcc/e` (from comparison with tests by karsan and jirsa 1969)
        rho_bar: float
            Optional parameter of plastic volume change, positive real value :math:`0=rhobar< sqrt(2/3)`, 0.2 (from
            comparison with tests by kupfer hilsdorf rusch 1969)
        big_h: float
            Optional linear hardening parameter for plasticity, positive real value (usually less than e),
            :math:`0.25*e` (from comparison with tests by karsan and jirsa 1969 and gopalaratnam and shah 1985)
        theta: float
            Optional ratio between isotropic and kinematic hardening, positive real value :math:`0=theta=1` (with: 0
            hardening kinematic only and 1 hardening isotropic only, 0.5 (from comparison with tests by karsan and jirsa 1969
            and gopalaratnam and shah 1985)
        tangent: float
            Optional integer to choose the computational stiffness matrix, 0: computational tangent; 1: damaged secant
            stiffness (hint: in case of strong nonlinearities use it with krylov-newton algorithm)
        """
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
        mat3d: obj
            Tag of perviously defined 3d ndmaterial material
        """
        self.mat3d = mat3d
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.mat3d.tag]
        self.to_process(osi)


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
        mat3d: obj
            Integer tag of previously defined 3d ndmaterial material
        """
        self.mat3d = mat3d
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.mat3d.tag]
        self.to_process(osi)


class MultiaxialCyclicPlasticity(NDMaterialBase):
    """
    The MultiaxialCyclicPlasticity NDMaterial Class
    
    This command is used to construct an multiaxial Cyclic Plasticity model for clays
    """
    op_type = 'MultiaxialCyclicPlasticity'

    def __init__(self, osi, rho, k_mod, g_mod, su, ho, h, m, beta, k_coeff):
        """
        Initial method for MultiaxialCyclicPlasticity

        Parameters
        ----------
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
        """
        self.rho = float(rho)
        self.k_mod = float(k_mod)
        self.g_mod = float(g_mod)
        self.su = float(su)
        self.ho = float(ho)
        self.h = float(h)
        self.m = float(m)
        self.beta = float(beta)
        self.k_coeff = float(k_coeff)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.rho, self.k_mod, self.g_mod, self.su, self.ho, self.h, self.m, self.beta, self.k_coeff]
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
        """
        self.mass_density = float(mass_density)
        self.big_c = float(big_c)
        self.bulk_mod = float(bulk_mod)
        self.ocr = float(ocr)
        self.mu_o = float(mu_o)
        self.alpha = float(alpha)
        self.lamb = float(lamb)
        self.h = float(h)
        self.m = float(m)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.mass_density, self.big_c, self.bulk_mod, self.ocr, self.mu_o, self.alpha, self.lamb, self.h, self.m]
        self.to_process(osi)


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
        three_d: obj
            Material tag for a previously-defined three-dimensional material
        """
        self.three_d = three_d
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.three_d.tag]
        self.to_process(osi)


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
    defined using stress–strain relationships in biaxial directions, the orientation
    of which is governed by the state of cracking in concrete. Although the
    concrete stress–strain relationship used in the FSAM is fundamentally
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
        """
        Initial method for FSAM

        Parameters
        ----------
        rho: float
            Material density
        s_x: obj
            Tag of uniaxialmaterial simulating horizontal (x) reinforcement
        s_y: obj
            Tag of uniaxialmaterial simulating vertical (y) reinforcement
        conc: obj
            Tag of uniaxialmaterial simulating concrete, shall be used with uniaxialmaterial concretecm
        rou_x: float
            Reinforcing ratio in horizontal (x) direction (:math:`roux = _{s,x}/a_{gross,x}`)
        rou_y: float
            Reinforcing ratio in vertical (x) direction (:math:`rouy = _{s,y}/a_{gross,y}`)
        nu: float
            Concrete friction coefficient (:math:`0.0 < \nu < 1.5`)
        alfadow: float
            Stiffness coefficient of reinforcement dowel action (:math:`0.0 < alfadow < 0.05`)
        """
        self.rho = float(rho)
        self.s_x = s_x
        self.s_y = s_y
        self.conc = conc
        self.rou_x = float(rou_x)
        self.rou_y = float(rou_y)
        self.nu = float(nu)
        self.alfadow = float(alfadow)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.rho, self.s_x.tag, self.s_y.tag, self.conc.tag, self.rou_x, self.rou_y, self.nu, self.alfadow]
        self.to_process(osi)


class ManzariDafalias(NDMaterialBase):
    """
    The ManzariDafalias NDMaterial Class
    
    This command is used to construct a multi-dimensional Manzari-Dafalias(2004) material.
    """
    op_type = 'ManzariDafalias'

    def __init__(self, osi, g0, nu, e_init, mc, c, lambda_c, e0, ksi, p_atm, m, h0, ch, nb, a0, nd, z_max, cz, den):
        """
        Initial method for ManzariDafalias

        Parameters
        ----------
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
        """
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
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.g0, self.nu, self.e_init, self.mc, self.c, self.lambda_c, self.e0, self.ksi, self.p_atm, self.m, self.h0, self.ch, self.nb, self.a0, self.nd, self.z_max, self.cz, self.den]
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
        k_mod: float
            Bulk module of the acoustic medium
        rho: float
            Mass density of the acoustic medium
        """
        self.k_mod = float(k_mod)
        self.rho = float(rho)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.k_mod, self.rho]
        self.to_process(osi)
