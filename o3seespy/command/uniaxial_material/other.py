from o3seespy.command.uniaxial_material.base_material import UniaxialMaterialBase



class Hardening(UniaxialMaterialBase):
    """
    The Hardening UniaxialMaterial Class
    
    This command is used to construct a uniaxial material object with combined linear kinematic and isotropic hardening.
    The model includes optional visco-plasticity using a Perzyna formulation.
    """
    op_type = 'Hardening'

    def __init__(self, osi, e_mod, sigma_y, h_iso, h_kin, eta=0.0):
        """
        Initial method for Hardening

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        e_mod: float
            Tangent stiffness
        sigma_y: float
            Yield stress or force
        h_iso: float
            Isotropic hardening modulus
        h_kin: float
            Kinematic hardening modulus
        eta: float, optional
            Visco-plastic coefficient (optional, default=0.0)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.Hardening(osi, e_mod=1.0, sigma_y=1.0, h_iso=1.0, h_kin=1.0, eta=0.0)
        """
        self.osi = osi
        self.e_mod = float(e_mod)
        self.sigma_y = float(sigma_y)
        self.h_iso = float(h_iso)
        self.h_kin = float(h_kin)
        self.eta = float(eta)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.e_mod, self.sigma_y, self.h_iso, self.h_kin, self.eta]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def set_fy(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'Fy', value, ele, eles)

    def set_e_mod(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'E', value, ele, eles)

    def set_hkin(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'Hkin', value, ele, eles)

    def set_hiso(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'Hiso', value, ele, eles)


class Cast(UniaxialMaterialBase):
    """
    The Cast UniaxialMaterial Class
    
    This command is used to construct a parallel material object made up of an arbitrary number of
    previously-constructed UniaxialMaterial objects.
    """
    op_type = 'Cast'

    def __init__(self, osi, n, bo, h, fy, e_mod, big_l, b, ro, c_r1, c_r2, a1: float=None, a2=1.0, a3: float=None, a4=1.0):
        """
        Initial method for Cast

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        n: int
            Number of yield fingers of the csf-brace
        bo: float
            Width of an individual yielding finger at its base of the csf-brace
        h: float
            Thickness of an individual yielding finger
        fy: float
            Yield strength of the steel material of the yielding finger
        e_mod: float
            Modulus of elasticity of the steel material of the yielding finger
        big_l: float
            Height of an individual yielding finger
        b: float
            Strain hardening ratio
        ro: float
            Parameter that controls the bauschinger effect. recommended values for $ro=between 10 to 30
        c_r1: float
            Parameter that controls the bauschinger effect. recommended value cr1=0.925
        c_r2: float
            Parameter that controls the bauschinger effect. recommended value cr2=0.150
        a1: float (default=True), optional
            Isotropic hardening parameter, increase of compression yield envelope as proportion of yield strength after
            a plastic deformation of a2*(pp/kp)
        a2: float, optional
            Isotropic hardening parameter (see explanation under a1). (optional default = 1.0)
        a3: float (default=True), optional
            Isotropic hardening parameter, increase of tension yield envelope as proportion of yield strength after a
            plastic deformation of a4*(pp/kp)
        a4: float, optional
            Isotropic hardening parameter (see explanation under a3). (optional default = 1.0)

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.Cast(osi, n=1, bo=1.0, h=1.0, fy=1.0, e_mod=1.0, big_l=1.0, b=1.0, ro=1.0, c_r1=1.0, c_r2=1.0, a1=None, a2=1.0, a3=None, a4=1.0)
        """
        self.osi = osi
        self.n = int(n)
        self.bo = float(bo)
        self.h = float(h)
        self.fy = float(fy)
        self.e_mod = float(e_mod)
        self.big_l = float(big_l)
        self.b = float(b)
        self.ro = float(ro)
        self.c_r1 = float(c_r1)
        self.c_r2 = float(c_r2)
        if a1 is None:
            self.a1 = None
        else:
            self.a1 = float(a1)
        self.a2 = float(a2)
        if a3 is None:
            self.a3 = None
        else:
            self.a3 = float(a3)
        self.a4 = float(a4)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.n, self.bo, self.h, self.fy, self.e_mod, self.big_l, self.b, self.ro, self.c_r1, self.c_r2]
        special_pms = ['a1', 'a2', 'a3', 'a4']
        packets = [False, False, False, False]
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


class ViscousDamper(UniaxialMaterialBase):
    """
    The ViscousDamper UniaxialMaterial Class
    
    This command is used to construct a ViscousDamper material, which represents the Maxwell Model (linear spring and
    nonlinear dashpot in series). The ViscousDamper material simulates the hysteretic response of nonlinear viscous
    dampers. An adaptive iterative algorithm has been implemented and validated to solve numerically the
    constitutive equations within a nonlinear viscous damper with a high-precision accuracy.
    """
    op_type = 'ViscousDamper'

    def __init__(self, osi, k_el, cd, alpha, l_gap=0.0, nm=1, rel_tol=1e-6, abs_tol=1e-10, max_half=15):
        """
        Initial method for ViscousDamper

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        k_el: float
            Elastic stiffness of linear spring to model the axial flexibility of a viscous damper (e.g. combined
            stiffness of the supporting brace and internal damper portion)
        cd: float
            Damping coefficient
        alpha: float
            Velocity exponent
        l_gap: float, optional
            Gap length to simulate the gap length due to the pin tolerance
        nm: int, optional
            Employed adaptive numerical algorithm (default value nm = 1; * ``1`` = dormand-prince54, * ``2`` = 6th order
            adams-bashforth-moulton, * ``3`` = modified rosenbrock triple)
        rel_tol: float, optional
            Tolerance for absolute relative error control of the adaptive iterative algorithm (default value 10^-6)
        abs_tol: float, optional
            Tolerance for absolute error control of adaptive iterative algorithm (default value 10^-10)
        max_half: int, optional
            Maximum number of sub-step iterations within an integration step (default value 15)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.ViscousDamper(osi, k_el=1.0, cd=1.0, alpha=1.0, l_gap=0.0, nm=1, rel_tol=1e-6, abs_tol=1e-10, max_half=15)
        """
        self.osi = osi
        self.k_el = float(k_el)
        self.cd = float(cd)
        self.alpha = float(alpha)
        self.l_gap = float(l_gap)
        self.nm = int(nm)
        self.rel_tol = float(rel_tol)
        self.abs_tol = float(abs_tol)
        self.max_half = int(max_half)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.k_el, self.cd, self.alpha, self.l_gap, self.nm, self.rel_tol, self.abs_tol, self.max_half]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def set_e_mod(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'E', value, ele, eles)

    def set_eta(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'eta', value, ele, eles)


class BilinearOilDamper(UniaxialMaterialBase):
    """
    The BilinearOilDamper UniaxialMaterial Class
    
    This command is used to construct a BilinearOilDamper material, which simulates the hysteretic response of bilinear
    oil dampers with relief valve. Two adaptive iterative algorithms have been implemented and validated to solve
    numerically the constitutive equations within a bilinear oil damper with a high-precision accuracy.
    """
    op_type = 'BilinearOilDamper'

    def __init__(self, osi, k_el, cd, fr=1.0, p=1.0, l_gap=0.0, nm=1, rel_tol=1e-6, abs_tol=1e-10, max_half=15):
        """
        Initial method for BilinearOilDamper

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        k_el: float
            Elastic stiffness of linear spring to model the axial flexibility of a viscous damper (e.g. combined
            stiffness of the supporting brace and internal damper portion)
        cd: float
            Damping coefficient
        fr: float, optional
            Damper relief load (default=1.0, damper property)
        p: float, optional
            Post-relief viscous damping coefficient ratio (default=1.0, linear oil damper)
        l_gap: float, optional
            Gap length to simulate the gap length due to the pin tolerance
        nm: int, optional
            Employed adaptive numerical algorithm (default value nm = 1; * ``1`` = dormand-prince54, * ``2`` = 6th order
            adams-bashforth-moulton, * ``3`` = modified rosenbrock triple)
        rel_tol: float, optional
            Tolerance for absolute relative error control of the adaptive iterative algorithm (default value 10^-6)
        abs_tol: float, optional
            Tolerance for absolute error control of adaptive iterative algorithm (default value 10^-10)
        max_half: int, optional
            Maximum number of sub-step iterations within an integration step (default value 15)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.BilinearOilDamper(osi, k_el=1.0, cd=1.0, fr=1.0, p=1.0, l_gap=0.0, nm=1, rel_tol=1e-6, abs_tol=1e-10, max_half=15)
        """
        self.osi = osi
        self.k_el = float(k_el)
        self.cd = float(cd)
        self.fr = float(fr)
        self.p = float(p)
        self.l_gap = float(l_gap)
        self.nm = int(nm)
        self.rel_tol = float(rel_tol)
        self.abs_tol = float(abs_tol)
        self.max_half = int(max_half)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.k_el, self.cd, self.fr, self.p, self.l_gap, self.nm, self.rel_tol, self.abs_tol, self.max_half]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def set_k(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'K', value, ele, eles)

    def set_c(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'C', value, ele, eles)

    def set_fr(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'Fr', value, ele, eles)

    def set_p(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'p', value, ele, eles)

    def set_l_gap(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'LGap', value, ele, eles)


class Bilin(UniaxialMaterialBase):
    """
    The Bilin UniaxialMaterial Class
    
    This command is used to construct a bilin material. The bilin material simulates the modified Ibarra-Krawinkler
    deterioration model with bilinear hysteretic response. Note that the hysteretic response of this material has been
    calibrated with respect to more than 350 experimental data of steel beam-to-column connections and multivariate
    regression formulas are provided to estimate the deterioration parameters of the model for different
    connection types. These relationships were developed by Lignos and Krawinkler (2009, 2011) and have
    been adopted by PEER/ATC (2010). The input parameters for this component model can be computed
    interactively from this `link <http://dimitrios-lignos.research.mcgill.ca/databases/>`_.
    **Use the module Component Model.**
    """
    op_type = 'Bilin'

    def __init__(self, osi, k0, as__plus, as__neg, my__plus, my__neg, lamda_s, lamda_c, lamda_a, lamda_k, c_s, c_c, c_a, c_k, theta_p__plus, theta_p__neg, theta_pc__plus, theta_pc__neg, res__pos, res__neg, theta_u__plus, theta_u__neg, d__plus, d__neg, n_factor=0.0):
        """
        Initial method for Bilin

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        k0: float
            Elastic stiffness
        as__plus: float
            Strain hardening ratio for positive loading direction
        as__neg: float
            Strain hardening ratio for negative loading direction
        my__plus: float
            Effective yield strength for positive loading direction
        my__neg: float
            Effective yield strength for negative loading direction (negative value)
        lamda_s: float
            Cyclic deterioration parameter for strength deterioration [e_t=lamda_s*m_y; set lamda_s = 0 to disable this
            mode of deterioration]
        lamda_c: float
            Cyclic deterioration parameter for post-capping strength deterioration [e_t=lamda_c*m_y; set lamda_c = 0 to
            disable this mode of deterioration]
        lamda_a: float
            Cyclic deterioration parameter for acceleration reloading stiffness deterioration (is not a deterioration
            mode for a component with bilinear hysteretic response) [input value is required, but not used; set lamda_a = 0].
        lamda_k: float
            Cyclic deterioration parameter for unloading stiffness deterioration [e_t=lamda_k*m_y; set lamda_k = 0 to
            disable this mode of deterioration]
        c_s: float
            Rate of strength deterioration. the default value is 1.0.
        c_c: float
            Rate of post-capping strength deterioration. the default value is 1.0.
        c_a: float
            Rate of accelerated reloading deterioration. the default value is 1.0.
        c_k: float
            Rate of unloading stiffness deterioration. the default value is 1.0.
        theta_p__plus: float
            Pre-capping rotation for positive loading direction (often noted as plastic rotation capacity)
        theta_p__neg: float
            Pre-capping rotation for negative loading direction (often noted as plastic rotation capacity) (positive
            value)
        theta_pc__plus: float
            Post-capping rotation for positive loading direction
        theta_pc__neg: float
            Post-capping rotation for negative loading direction (positive value)
        res__pos: float
            Residual strength ratio for positive loading direction
        res__neg: float
            Residual strength ratio for negative loading direction (positive value)
        theta_u__plus: float
            Ultimate rotation capacity for positive loading direction
        theta_u__neg: float
            Ultimate rotation capacity for negative loading direction (positive value)
        d__plus: float
            Rate of cyclic deterioration in the positive loading direction (this parameter is used to create assymetric
            hysteretic behavior for the case of a composite beam). for symmetric hysteretic response use 1.0.
        d__neg: float
            Rate of cyclic deterioration in the negative loading direction (this parameter is used to create assymetric
            hysteretic behavior for the case of a composite beam). for symmetric hysteretic response use 1.0.
        n_factor: float, optional
            Elastic stiffness amplification factor, mainly for use with concentrated plastic hinge elements (optional,
            default = 0).

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.Bilin(osi, k0=1.0, as__plus=1.0, as__neg=1.0, my__plus=1.0, my__neg=1.0, lamda_s=1.0, lamda_c=1.0, lamda_a=1.0, lamda_k=1.0, c_s=1.0, c_c=1.0, c_a=1.0, c_k=1.0, theta_p__plus=1.0, theta_p__neg=1.0, theta_pc__plus=1.0, theta_pc__neg=1.0, res__pos=1.0, res__neg=1.0, theta_u__plus=1.0, theta_u__neg=1.0, d__plus=1.0, d__neg=1.0, n_factor=0.0)
        """
        self.osi = osi
        self.k0 = float(k0)
        self.as__plus = float(as__plus)
        self.as__neg = float(as__neg)
        self.my__plus = float(my__plus)
        self.my__neg = float(my__neg)
        self.lamda_s = float(lamda_s)
        self.lamda_c = float(lamda_c)
        self.lamda_a = float(lamda_a)
        self.lamda_k = float(lamda_k)
        self.c_s = float(c_s)
        self.c_c = float(c_c)
        self.c_a = float(c_a)
        self.c_k = float(c_k)
        self.theta_p__plus = float(theta_p__plus)
        self.theta_p__neg = float(theta_p__neg)
        self.theta_pc__plus = float(theta_pc__plus)
        self.theta_pc__neg = float(theta_pc__neg)
        self.res__pos = float(res__pos)
        self.res__neg = float(res__neg)
        self.theta_u__plus = float(theta_u__plus)
        self.theta_u__neg = float(theta_u__neg)
        self.d__plus = float(d__plus)
        self.d__neg = float(d__neg)
        self.n_factor = float(n_factor)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.k0, self.as__plus, self.as__neg, self.my__plus, self.my__neg, self.lamda_s, self.lamda_c, self.lamda_a, self.lamda_k, self.c_s, self.c_c, self.c_a, self.c_k, self.theta_p__plus, self.theta_p__neg, self.theta_pc__plus, self.theta_pc__neg, self.res__pos, self.res__neg, self.theta_u__plus, self.theta_u__neg, self.d__plus, self.d__neg, self.n_factor]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class ModIMKPeakOriented(UniaxialMaterialBase):
    """
    The ModIMKPeakOriented UniaxialMaterial Class
    
    This command is used to construct a ModIMKPeakOriented material. This material simulates the modified
    Ibarra-Medina-Krawinkler deterioration model with peak-oriented hysteretic response. Note that the
    hysteretic response of this material has been calibrated with respect to 200 experimental data of
    RC beams in order to estimate the deterioration parameters of the model. This information was
    developed by Lignos and Krawinkler (2012). NOTE: before you use this material make sure that
    you have downloaded the latest OpenSees version. A youtube video presents a summary of this
    model including the way to be used within openSees `youtube link <http://youtu.be/YHBHQ-xuybE>`_.
    """
    op_type = 'ModIMKPeakOriented'

    def __init__(self, osi, k0, as__plus, as__neg, my__plus, my__neg, lamda_s, lamda_c, lamda_a, lamda_k, c_s, c_c, c_a, c_k, theta_p__plus, theta_p__neg, theta_pc__plus, theta_pc__neg, res__pos, res__neg, theta_u__plus, theta_u__neg, d__plus, d__neg):
        """
        Initial method for ModIMKPeakOriented

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        k0: float
            Elastic stiffness
        as__plus: float
            Strain hardening ratio for positive loading direction
        as__neg: float
            Strain hardening ratio for negative loading direction
        my__plus: float
            Effective yield strength for positive loading direction
        my__neg: float
            Effective yield strength for negative loading direction (negative value)
        lamda_s: float
            Cyclic deterioration parameter for strength deterioration [e_t=lamda_s*m_y, see lignos and krawinkler
            (2011); set lamda_s = 0 to disable this mode of deterioration]
        lamda_c: float
            Cyclic deterioration parameter for post-capping strength deterioration [e_t=lamda_c*m_y, see lignos and
            krawinkler (2011); set lamda_c = 0 to disable this mode of deterioration]
        lamda_a: float
            Cyclic deterioration parameter for accelerated reloading stiffness deterioration [e_t=lamda_a*m_y, see
            lignos and krawinkler (2011); set lamda_a = 0 to disable this mode of deterioration]
        lamda_k: float
            Cyclic deterioration parameter for unloading stiffness deterioration [e_t=lamda_k*m_y, see lignos and
            krawinkler (2011); set lamda_k = 0 to disable this mode of deterioration]
        c_s: float
            Rate of strength deterioration. the default value is 1.0.
        c_c: float
            Rate of post-capping strength deterioration. the default value is 1.0.
        c_a: float
            Rate of accelerated reloading deterioration. the default value is 1.0.
        c_k: float
            Rate of unloading stiffness deterioration. the default value is 1.0.
        theta_p__plus: float
            Pre-capping rotation for positive loading direction (often noted as plastic rotation capacity)
        theta_p__neg: float
            Pre-capping rotation for negative loading direction (often noted as plastic rotation capacity) (must be
            defined as a positive value)
        theta_pc__plus: float
            Post-capping rotation for positive loading direction
        theta_pc__neg: float
            Post-capping rotation for negative loading direction (must be defined as a positive value)
        res__pos: float
            Residual strength ratio for positive loading direction
        res__neg: float
            Residual strength ratio for negative loading direction (must be defined as a positive value)
        theta_u__plus: float
            Ultimate rotation capacity for positive loading direction
        theta_u__neg: float
            Ultimate rotation capacity for negative loading direction (must be defined as a positive value)
        d__plus: float
            Rate of cyclic deterioration in the positive loading direction (this parameter is used to create assymetric
            hysteretic behavior for the case of a composite beam). for symmetric hysteretic response use 1.0.
        d__neg: float
            Rate  of cyclic deterioration in the negative loading direction (this parameter is used to create assymetric
            hysteretic behavior for the case of a composite beam). for symmetric hysteretic response use 1.0.

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.ModIMKPeakOriented(osi, k0=1.0, as__plus=1.0, as__neg=1.0, my__plus=1.0, my__neg=1.0, lamda_s=1.0, lamda_c=1.0, lamda_a=1.0, lamda_k=1.0, c_s=1.0, c_c=1.0, c_a=1.0, c_k=1.0, theta_p__plus=1.0, theta_p__neg=1.0, theta_pc__plus=1.0, theta_pc__neg=1.0, res__pos=1.0, res__neg=1.0, theta_u__plus=1.0, theta_u__neg=1.0, d__plus=1.0, d__neg=1.0)
        """
        self.osi = osi
        self.k0 = float(k0)
        self.as__plus = float(as__plus)
        self.as__neg = float(as__neg)
        self.my__plus = float(my__plus)
        self.my__neg = float(my__neg)
        self.lamda_s = float(lamda_s)
        self.lamda_c = float(lamda_c)
        self.lamda_a = float(lamda_a)
        self.lamda_k = float(lamda_k)
        self.c_s = float(c_s)
        self.c_c = float(c_c)
        self.c_a = float(c_a)
        self.c_k = float(c_k)
        self.theta_p__plus = float(theta_p__plus)
        self.theta_p__neg = float(theta_p__neg)
        self.theta_pc__plus = float(theta_pc__plus)
        self.theta_pc__neg = float(theta_pc__neg)
        self.res__pos = float(res__pos)
        self.res__neg = float(res__neg)
        self.theta_u__plus = float(theta_u__plus)
        self.theta_u__neg = float(theta_u__neg)
        self.d__plus = float(d__plus)
        self.d__neg = float(d__neg)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.k0, self.as__plus, self.as__neg, self.my__plus, self.my__neg, self.lamda_s, self.lamda_c, self.lamda_a, self.lamda_k, self.c_s, self.c_c, self.c_a, self.c_k, self.theta_p__plus, self.theta_p__neg, self.theta_pc__plus, self.theta_pc__neg, self.res__pos, self.res__neg, self.theta_u__plus, self.theta_u__neg, self.d__plus, self.d__neg]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class ModIMKPinching(UniaxialMaterialBase):
    """
    The ModIMKPinching UniaxialMaterial Class
    
    This command is used to construct a ModIMKPinching material. This material simulates the modified
    Ibarra-Medina-Krawinkler deterioration model with pinching hysteretic response. NOTE: **before you
    use this material make sure that you have downloaded the latest OpenSees version**. A youtube
    video presents a summary of this model including the way to be used within openSees `youtube
    link <http://youtu.be/YHBHQ-xuybE>`_.
    """
    op_type = 'ModIMKPinching'

    def __init__(self, osi, k0, as__plus, as__neg, my__plus, my__neg, fpr_pos, fpr_neg, a_pinch, lamda_s, lamda_c, lamda_a, lamda_k, c_s, c_c, c_a, c_k, theta_p__plus, theta_p__neg, theta_pc__plus, theta_pc__neg, res__pos, res__neg, theta_u__plus, theta_u__neg, d__plus, d__neg):
        """
        Initial method for ModIMKPinching

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        k0: float
            Elastic stiffness
        as__plus: float
            Strain hardening ratio for positive loading direction
        as__neg: float
            Strain hardening ratio for negative loading direction
        my__plus: float
            Effective yield strength for positive loading direction
        my__neg: float
            Effective yield strength for negative loading direction (must be defined as a negative value)
        fpr_pos: float
            Ratio of the force at which reloading begins to force corresponding to the maximum historic deformation
            demand (positive loading direction)
        fpr_neg: float
            Ratio of the force at which reloading begins to force corresponding to the absolute maximum historic
            deformation demand (negative loading direction)
        a_pinch: float
            Ratio of reloading stiffness
        lamda_s: float
            Cyclic deterioration parameter for strength deterioration [e_t=lamda_s*m_y, see lignos and krawinkler
            (2011); set lamda_s = 0 to disable this mode of deterioration]
        lamda_c: float
            Cyclic deterioration parameter for post-capping strength deterioration [e_t=lamda_c*m_y, see lignos and
            krawinkler (2011); set lamda_c = 0 to disable this mode of deterioration]
        lamda_a: float
            Cyclic deterioration parameter for accelerated reloading stiffness deterioration [e_t=lamda_a*m_y, see
            lignos and krawinkler (2011); set lamda_a = 0 to disable this mode of deterioration]
        lamda_k: float
            Cyclic deterioration parameter for unloading stiffness deterioration [e_t=lamda_k*m_y, see lignos and
            krawinkler (2011); set lamda_k = 0 to disable this mode of deterioration]
        c_s: float
            Rate of strength deterioration. the default value is 1.0.
        c_c: float
            Rate of post-capping strength deterioration. the default value is 1.0.
        c_a: float
            Rate of accelerated reloading deterioration. the default value is 1.0.
        c_k: float
            Rate of unloading stiffness deterioration. the default value is 1.0.
        theta_p__plus: float
            Pre-capping rotation for positive loading direction (often noted as plastic rotation capacity)
        theta_p__neg: float
            Pre-capping rotation for negative loading direction (often noted as plastic rotation capacity) (must be
            defined as a positive value)
        theta_pc__plus: float
            Post-capping rotation for positive loading direction
        theta_pc__neg: float
            Post-capping rotation for negative loading direction (must be defined as a positive value)
        res__pos: float
            Residual strength ratio for positive loading direction
        res__neg: float
            Residual strength ratio for negative loading direction (must be defined as a positive value)
        theta_u__plus: float
            Ultimate rotation capacity for positive loading direction
        theta_u__neg: float
            Ultimate rotation capacity for negative loading direction (must be defined as a positive value)
        d__plus: float
            Rate of cyclic deterioration in the positive loading direction (this parameter is used to create assymetric
            hysteretic behavior for the case of a composite beam). for symmetric hysteretic response use 1.0.
        d__neg: float
            Rate of cyclic deterioration in the negative loading direction (this parameter is used to create assymetric
            hysteretic behavior for the case of a composite beam). for symmetric hysteretic response use 1.0.

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.ModIMKPinching(osi, k0=1.0, as__plus=1.0, as__neg=1.0, my__plus=1.0, my__neg=1.0, fpr_pos=1.0, fpr_neg=1.0, a_pinch=1.0, lamda_s=1.0, lamda_c=1.0, lamda_a=1.0, lamda_k=1.0, c_s=1.0, c_c=1.0, c_a=1.0, c_k=1.0, theta_p__plus=1.0, theta_p__neg=1.0, theta_pc__plus=1.0, theta_pc__neg=1.0, res__pos=1.0, res__neg=1.0, theta_u__plus=1.0, theta_u__neg=1.0, d__plus=1.0, d__neg=1.0)
        """
        self.osi = osi
        self.k0 = float(k0)
        self.as__plus = float(as__plus)
        self.as__neg = float(as__neg)
        self.my__plus = float(my__plus)
        self.my__neg = float(my__neg)
        self.fpr_pos = float(fpr_pos)
        self.fpr_neg = float(fpr_neg)
        self.a_pinch = float(a_pinch)
        self.lamda_s = float(lamda_s)
        self.lamda_c = float(lamda_c)
        self.lamda_a = float(lamda_a)
        self.lamda_k = float(lamda_k)
        self.c_s = float(c_s)
        self.c_c = float(c_c)
        self.c_a = float(c_a)
        self.c_k = float(c_k)
        self.theta_p__plus = float(theta_p__plus)
        self.theta_p__neg = float(theta_p__neg)
        self.theta_pc__plus = float(theta_pc__plus)
        self.theta_pc__neg = float(theta_pc__neg)
        self.res__pos = float(res__pos)
        self.res__neg = float(res__neg)
        self.theta_u__plus = float(theta_u__plus)
        self.theta_u__neg = float(theta_u__neg)
        self.d__plus = float(d__plus)
        self.d__neg = float(d__neg)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.k0, self.as__plus, self.as__neg, self.my__plus, self.my__neg, self.fpr_pos, self.fpr_neg, self.a_pinch, self.lamda_s, self.lamda_c, self.lamda_a, self.lamda_k, self.c_s, self.c_c, self.c_a, self.c_k, self.theta_p__plus, self.theta_p__neg, self.theta_pc__plus, self.theta_pc__neg, self.res__pos, self.res__neg, self.theta_u__plus, self.theta_u__neg, self.d__plus, self.d__neg]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class SAWS(UniaxialMaterialBase):
    """
    The SAWS UniaxialMaterial Class
    
    This file contains the class definition for SAWSMaterial. SAWSMaterial provides the implementation of a
    one-dimensional hysteretic model develeped as part of the CUREe Caltech wood frame project.
    """
    op_type = 'SAWS'

    def __init__(self, osi, f0, fi, du, s0, r1, r2, r3, r4, alpha, beta):
        """
        Initial method for SAWS

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        f0: float
            Intercept strength of the shear wall spring element for the asymtotic line to the envelope curve f0 > fi > 0
        fi: float
            Intercept strength of the spring element for the pinching branch of the hysteretic curve. (fi > 0).
        du: float
            Spring element displacement at ultimate load. (du > 0).
        s0: float
            Initial stiffness of the shear wall spring element (s0 > 0).
        r1: float
            Stiffness ratio of the asymptotic line to the spring element envelope curve. the slope of this line is r1
            s0. (0 < r1 < 1.0).
        r2: float
            Stiffness ratio of the descending branch of the spring element envelope curve. the slope of this line is r2
            s0. ( r2 < 0).
        r3: float
            Stiffness ratio of the unloading branch off the spring element envelope curve. the slope of this line is r3
            s0. ( r3 1).
        r4: float
            Stiffness ratio of the pinching branch for the spring element. the slope of this line is r4 s0. ( r4 > 0).
        alpha: float
            Stiffness degradation parameter for the shear wall spring element. (alpha > 0).
        beta: float
            Stiffness degradation parameter for the spring element. (beta > 0).

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.SAWS(osi, f0=1.0, fi=1.0, du=1.0, s0=1.0, r1=1.0, r2=1.0, r3=1.0, r4=1.0, alpha=1.0, beta=1.0)
        """
        self.osi = osi
        self.f0 = float(f0)
        self.fi = float(fi)
        self.du = float(du)
        self.s0 = float(s0)
        self.r1 = float(r1)
        self.r2 = float(r2)
        self.r3 = float(r3)
        self.r4 = float(r4)
        self.alpha = float(alpha)
        self.beta = float(beta)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.f0, self.fi, self.du, self.s0, self.r1, self.r2, self.r3, self.r4, self.alpha, self.beta]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class BarSlip(UniaxialMaterialBase):
    """
    The BarSlip UniaxialMaterial Class
    
    This command is used to construct a uniaxial material that simulates the bar force versus slip response of a
    reinforcing bar anchored in a beam-column joint. The model exhibits degradation under cyclic loading. Cyclic
    degradation of strength and stiffness occurs in three ways: unloading stiffness degradation, reloading
    stiffness degradation, strength degradation.
    """
    op_type = 'BarSlip'

    def __init__(self, osi, fc, fy, es, fu, eh, db, ld, nb, depth, height, bs_flag, otype, anc_lratio=1.0, damage='Damage', unit='psi'):
        """
        Initial method for BarSlip

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        fc: float
            Positive floating point value defining the compressive strength of the concrete in which the reinforcing bar
            is anchored
        fy: float
            Positive floating point value defining the yield strength of the reinforcing steel
        es: float
            Floating point value defining the modulus of elasticity of the reinforcing steel
        fu: float
            Positive floating point value defining the ultimate strength of the reinforcing steel
        eh: float
            Floating point value defining the hardening modulus of the reinforcing steel
        db: float
            Point value defining the diameter of reinforcing steel
        ld: float
            Floating point value defining the development length of the reinforcing steel
        nb: int
            An integer defining the number of anchored bars
        depth: float
            Floating point value defining the dimension of the member (beam or column) perpendicular to the dimension of
            the plane of the paper
        height: float
            Floating point value defining the height of the flexural member, perpendicular to direction in which the
            reinforcing steel is placed, but in the plane of the paper
        bs_flag: str
            String indicating relative bond strength for the anchored reinforcing bar (options: ``'strong'`` or
            ``'weak'``)
        otype: str
            String indicating where the reinforcing bar is placed. (options: ``'beamtop'``, ``'beambot'`` or
            ``'column'``)
        anc_lratio: float, optional
            Floating point value defining the ratio of anchorage length used for the reinforcing bar to the dimension of
            the joint in the direction of the reinforcing bar (optional, default: 1.0)
        damage: str, optional
            String indicating type of damage:whether there is full damage in the material or no damage (optional,
            options: ``'damage'``, ``'nodamage'`` ; default: ``'damage'``)
        unit: str, optional
            String indicating the type of unit system used (optional, options: ``'psi'``, ``'mpa'``, ``'pa'``,
            ``'psf'``, ``'ksi'``, ``'ksf'``) (default: ``'psi'`` / ``'mpa'``)

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.BarSlip(osi, fc=1.0, fy=1.0, es=1.0, fu=1.0, eh=1.0, db=1.0, ld=1.0, nb=1, depth=1.0, height=1.0, anc_lratio=1.0, bs_flag="Strong", otype="beamtop", damage='Damage', unit='psi')
        """
        self.osi = osi
        self.fc = float(fc)
        self.fy = float(fy)
        self.es = float(es)
        self.fu = float(fu)
        self.eh = float(eh)
        self.db = float(db)
        self.ld = float(ld)
        self.nb = int(nb)
        self.depth = float(depth)
        self.height = float(height)
        self.anc_lratio = float(anc_lratio)
        self.bs_flag = bs_flag
        self.otype = otype
        self.damage = damage
        self.unit = unit
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.fc, self.fy, self.es, self.fu, self.eh, self.db, self.ld, self.nb, self.depth, self.height, self.anc_lratio, self.bs_flag, self.otype, self.damage, self.unit]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class BondSP01(UniaxialMaterialBase):
    """
    The BondSP01 UniaxialMaterial Class
    
    This command is used to construct a uniaxial material object for capturing strain penetration effects at the
    column-to-footing, column-to-bridge bent caps, and wall-to-footing intersections. In these cases, the bond slip
    associated with strain penetration typically occurs along a portion of the anchorage length. This model can
    also be applied to the beam end regions, where the strain penetration may include slippage of the bar
    along the entire anchorage length, but the model parameters should be chosen appropriately.This
    model is for fully anchored steel reinforcement bars that experience bond slip along a portion
    of the anchorage length due to strain penetration effects, which are usually the case for
    column and wall longitudinal bars anchored into footings or bridge joints
    """
    op_type = 'Bond_SP01'

    def __init__(self, osi, fy, sy, fu, su, b, big_r):
        """
        Initial method for BondSP01

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        fy: float
            Yield strength of the reinforcement steel
        sy: float
            Rebar slip at member interface under yield stress. (see notes below)
        fu: float
            Ultimate strength of the reinforcement steel
        su: float
            Rebar slip at the loaded end at the bar fracture strength
        b: float
            Initial hardening ratio in the monotonic slip vs. bar stress response (0.3~0.5)
        big_r: float
            Pinching factor for the cyclic slip vs. bar response (0.5~1.0)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.BondSP01(osi, fy=1.0, sy=1.0, fu=1.0, su=1.0, b=1.0, big_r=1.0)
        """
        self.osi = osi
        self.fy = float(fy)
        self.sy = float(sy)
        self.fu = float(fu)
        self.su = float(su)
        self.b = float(b)
        self.big_r = float(big_r)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.fy, self.sy, self.fu, self.su, self.b, self.big_r]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class Fatigue(UniaxialMaterialBase):
    """
    The Fatigue UniaxialMaterial Class
    
    The fatigue material uses a modified rainflow cycle counting algorithm to accumulate damage in a material using
    Miner's Rule. Element stress/strain relationships become zero when fatigue life is exhausted.
    """
    op_type = 'Fatigue'

    def __init__(self, osi, other, e0: float=None, m: float=None, min: float=None, max: float=None):
        """
        Initial method for Fatigue

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        other: obj
            Unique material object integer object for the material that is being wrapped
        e0: float, optional
            Value of strain at which one cycle will cause failure (default 0.191)
        m: float, optional
            Slope of coffin-manson curve in log-log space (default -0.458)
        min: float, optional
            Global minimum value for strain or deformation (default -1e16)
        max: float, optional
            Global maximum value for strain or deformation (default 1e16)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> other = o3.uniaxial_material.Hardening(osi, e_mod=1.0, sigma_y=1.0, h_iso=1.0, h_kin=1.0, eta=0.0)
        >>> o3.uniaxial_material.Fatigue(osi, other=other, e0=0.191, m=-0.458, min=-1e16, max=1e16)
        """
        self.osi = osi
        self.other = other
        if e0 is None:
            self.e0 = None
        else:
            self.e0 = float(e0)
        if m is None:
            self.m = None
        else:
            self.m = float(m)
        if min is None:
            self.min = None
        else:
            self.min = float(min)
        if max is None:
            self.max = None
        else:
            self.max = float(max)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.other.tag]
        if getattr(self, 'e0') is not None:
            self._parameters += ['-E0', self.e0]
        if getattr(self, 'm') is not None:
            self._parameters += ['-m', self.m]
        if getattr(self, 'min') is not None:
            self._parameters += ['-min', self.min]
        if getattr(self, 'max') is not None:
            self._parameters += ['-max', self.max]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class ImpactMaterial(UniaxialMaterialBase):
    """
    The ImpactMaterial UniaxialMaterial Class
    
    This command is used to construct an impact material object
    """
    op_type = 'ImpactMaterial'

    def __init__(self, osi, k1, k2, sigy, gap):
        """
        Initial method for ImpactMaterial

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        k1: float
            Initial stiffness
        k2: float
            Secondary stiffness
        sigy: float
            Yield displacement
        gap: float
            Initial gap

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.ImpactMaterial(osi, k1=1.0, k2=1.0, sigy=1.0, gap=1.0)
        """
        self.osi = osi
        self.k1 = float(k1)
        self.k2 = float(k2)
        self.sigy = float(sigy)
        self.gap = float(gap)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.k1, self.k2, self.sigy, self.gap]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class HyperbolicGapMaterial(UniaxialMaterialBase):
    """
    The HyperbolicGapMaterial UniaxialMaterial Class
    
    This command is used to construct a hyperbolic gap material object.
    """
    op_type = 'HyperbolicGapMaterial'

    def __init__(self, osi, kmax, kur, rf, fult, gap):
        """
        Initial method for HyperbolicGapMaterial

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        kmax: float
            Initial stiffness
        kur: float
            Unloading/reloading stiffness
        rf: float
            Failure ratio
        fult: float
            Ultimate (maximum) passive resistance
        gap: float
            Initial gap

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.HyperbolicGapMaterial(osi, kmax=1.0, kur=1.0, rf=1.0, fult=1.0, gap=1.0)
        """
        self.osi = osi
        self.kmax = float(kmax)
        self.kur = float(kur)
        self.rf = float(rf)
        self.fult = float(fult)
        self.gap = float(gap)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.kmax, self.kur, self.rf, self.fult, self.gap]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class LimitState(UniaxialMaterialBase):
    """
    The LimitState UniaxialMaterial Class
    
    This command is used to construct a uniaxial hysteretic material object with pinching of force and deformation,
    damage due to ductility and energy, and degraded unloading stiffness based on ductility. Failure of the material
    is defined by the associated Limit Curve.
    """
    op_type = 'LimitState'

    def __init__(self, osi, s1p, e1p, s2p, e2p, s3p, e3p, s1n, e1n, s2n, e2n, s3n, e3n, pinch_x, pinch_y, damage1, damage2, beta, curve, curve_type):
        """
        Initial method for LimitState

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        s1p: float
            Stress and strain (or force & deformation) at first point of the envelope in the positive direction
        e1p: float
            Stress and strain (or force & deformation) at first point of the envelope in the positive direction
        s2p: float
            Stress and strain (or force & deformation) at second point of the envelope in the positive direction
        e2p: float
            Stress and strain (or force & deformation) at second point of the envelope in the positive direction
        s3p: float
            Stress and strain (or force & deformation) at third point of the envelope in the positive direction
        e3p: float
            Stress and strain (or force & deformation) at third point of the envelope in the positive direction
        s1n: float
            Stress and strain (or force & deformation) at first point of the envelope in the negative direction
        e1n: float
            Stress and strain (or force & deformation) at first point of the envelope in the negative direction
        s2n: float
            Stress and strain (or force & deformation) at second point of the envelope in the negative direction
        e2n: float
            Stress and strain (or force & deformation) at second point of the envelope in the negative direction
        s3n: float
            Stress and strain (or force & deformation) at third point of the envelope in the negative direction
        e3n: float
            Stress and strain (or force & deformation) at third point of the envelope in the negative direction
        pinch_x: float
            Pinching factor for strain (or deformation) during reloading
        pinch_y: float
            Pinching factor for stress (or force) during reloading
        damage1: float
            Damage due to ductility: d1(m-1)
        damage2: float
            Damage due to energy: d2(ei/eult)
        beta: float
            Power used to determine the degraded unloading stiffness based on ductility, m-b (optional, default=0.0)
        curve: obj
            An integer object for the limit curve defining the limit surface
        curve_type: int
            An integer defining the type of limitcurve (0 = no curve, 1 = axial curve, all other curves can be any other
            integer)

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> curve = 1
        >>> o3.uniaxial_material.LimitState(osi, s1p=1.0, e1p=1.0, s2p=1.0, e2p=1.0, s3p=1.0, e3p=1.0, s1n=1.0, e1n=1.0, s2n=1.0, e2n=1.0, s3n=1.0, e3n=1.0, pinch_x=1.0, pinch_y=1.0, damage1=1.0, damage2=1.0, beta=1.0, curve=curve, curve_type=1)
        """
        self.osi = osi
        self.s1p = float(s1p)
        self.e1p = float(e1p)
        self.s2p = float(s2p)
        self.e2p = float(e2p)
        self.s3p = float(s3p)
        self.e3p = float(e3p)
        self.s1n = float(s1n)
        self.e1n = float(e1n)
        self.s2n = float(s2n)
        self.e2n = float(e2n)
        self.s3n = float(s3n)
        self.e3n = float(e3n)
        self.pinch_x = float(pinch_x)
        self.pinch_y = float(pinch_y)
        self.damage1 = float(damage1)
        self.damage2 = float(damage2)
        self.beta = float(beta)
        self.curve = curve
        self.curve_type = int(curve_type)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.s1p, self.e1p, self.s2p, self.e2p, self.s3p, self.e3p, self.s1n, self.e1n, self.s2n, self.e2n, self.s3n, self.e3n, self.pinch_x, self.pinch_y, self.damage1, self.damage2, self.beta, self.curve.tag, self.curve_type]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def set_state_flag(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'stateFlag', value, ele, eles)


class MinMax(UniaxialMaterialBase):
    """
    The MinMax UniaxialMaterial Class
    
    This command is used to construct a MinMax material object. This stress-strain behaviour for this material is
    provided by another material. If however the strain ever falls below or above certain threshold values, the
    other material is assumed to have failed. From that point on, values of 0.0 are returned for the tangent and stress.
    """
    op_type = 'MinMax'

    def __init__(self, osi, other, min_strain: float=None, max_strain: float=None):
        """
        Initial method for MinMax

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        other: obj
            Object of the other material
        min_strain: float, optional
            Minimum value of strain. optional default = -1.0e16.
        max_strain: float, optional
            Max value of strain. optional default = 1.0e16.

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> other = o3.uniaxial_material.Hardening(osi, e_mod=1.0, sigma_y=1.0, h_iso=1.0, h_kin=1.0, eta=0.0)
        >>> o3.uniaxial_material.MinMax(osi, other=other, min_strain=1e-16, max_strain=1e16)
        """
        self.osi = osi
        self.other = other
        if min_strain is None:
            self.min_strain = None
        else:
            self.min_strain = float(min_strain)
        if max_strain is None:
            self.max_strain = None
        else:
            self.max_strain = float(max_strain)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.other.tag]
        if getattr(self, 'min_strain') is not None:
            self._parameters += ['-min', self.min_strain]
        if getattr(self, 'max_strain') is not None:
            self._parameters += ['-max', self.max_strain]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class ElasticBilin(UniaxialMaterialBase):
    """
    The ElasticBilin UniaxialMaterial Class
    
    This command is used to construct an elastic bilinear uniaxial material object. Unlike all other bilinear materials,
    the unloading curve follows the loading curve exactly.
    """
    op_type = 'ElasticBilin'

    def __init__(self, osi, ep1, ep2, eps_p2, en1: float=None, en2: float=None, eps_n2: float=None):
        """
        Initial method for ElasticBilin

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ep1: float
            Tangent in tension for stains: 0 <= strains <=    ``epsp2``
        ep2: float
            Tangent when material in tension with strains >    ``epsp2``
        eps_p2: float
            Strain at which material changes tangent in tension.
        en1: float (default=True), optional
            Optional, default =    ``ep1``. tangent in compression for stains: 0 < strains <=    ``epsn2``
        en2: float (default=True), optional
            Optional, default =    ``ep2``. tangent in compression with strains <    ``epsn2``
        eps_n2: float (default=True), optional
            Optional, default = ``-epsp2``. strain at which material changes tangent in compression.

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.ElasticBilin(osi, ep1=1.0, ep2=1.0, eps_p2=1.0, en1=None, en2=None, eps_n2=None)
        """
        self.osi = osi
        self.ep1 = float(ep1)
        self.ep2 = float(ep2)
        self.eps_p2 = float(eps_p2)
        if en1 is None:
            self.en1 = None
        else:
            self.en1 = float(en1)
        if en2 is None:
            self.en2 = None
        else:
            self.en2 = float(en2)
        if eps_n2 is None:
            self.eps_n2 = None
        else:
            self.eps_n2 = float(eps_n2)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.ep1, self.ep2, self.eps_p2]
        special_pms = ['en1', 'en2', 'eps_n2']
        packets = [False, False, False]
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


class ElasticMultiLinear(UniaxialMaterialBase):
    """
    The ElasticMultiLinear UniaxialMaterial Class
    
    This command is used to construct a multi-linear elastic uniaxial material object. The nonlinear stress-strain
    relationship is given by a multi-linear curve that is define by a set of points. The behavior is nonlinear but it
    is elastic. This means that the material loads and unloads along the same curve, and no energy is dissipated.
    The slope given by the last two specified points on the positive strain axis is extrapolated to infinite
    positive strain. Similarly, the slope given by the last two specified points on the negative strain
    axis is extrapolated to infinite negative strain. The number of provided strain points needs to be
    equal to the number of provided stress points.
    """
    op_type = 'ElasticMultiLinear'

    def __init__(self, osi, eta=0.0, strain: list=None, stress: list=None):
        """
        Initial method for ElasticMultiLinear

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        eta: float, optional
            Damping tangent (optional, default=0.0)
        strain: list, optional
            List of strain points along stress-strain curve
        stress: list, optional
            List of stress points along stress-strain curve

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> strain = [1.0, 1.0]
        >>> stress = [1.0, 1.0]
        >>> o3.uniaxial_material.ElasticMultiLinear(osi, eta=0.0, strain=strain, stress=stress)
        """
        self.osi = osi
        self.eta = float(eta)
        self.strain = strain
        self.stress = stress
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.eta]
        if getattr(self, 'strain') is not None:
            self._parameters += ['-strain', *self.strain]
        if getattr(self, 'stress') is not None:
            self._parameters += ['-stress', *self.stress]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class MultiLinear(UniaxialMaterialBase):
    """
    The MultiLinear UniaxialMaterial Class
    
    This command is used to construct a uniaxial multilinear material object.
    """
    op_type = 'MultiLinear'

    def __init__(self, osi, pts):
        """
        Initial method for MultiLinear

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        pts: list
            A list of strain and stress points ``pts = [strain1, stress1, strain2, stress2, ..., ]``

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> pts = [0.01, 1.0, 0.02, 2.0]
        >>> o3.uniaxial_material.MultiLinear(osi, pts=pts)
        """
        self.osi = osi
        self.pts = pts
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, *self.pts]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class InitStrainMaterial(UniaxialMaterialBase):
    """
    The InitStrainMaterial UniaxialMaterial Class
    
    This command is used to construct an Initial Strain material object. The stress-strain behaviour for this material
    is defined by another material. Initial Strain Material enables definition of initial strains for the material under
    consideration. The stress that corresponds to the initial strain will be calculated from the other material.
    """
    op_type = 'InitStrainMaterial'

    def __init__(self, osi, other, init_strain):
        """
        Initial method for InitStrainMaterial

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        other: obj
            Object of the other material
        init_strain: float
            Initial strain

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> other = o3.uniaxial_material.ModIMKPinching(osi, k0=1.0, as__plus=1.0, as__neg=1.0, my__plus=1.0, my__neg=1.0, fpr_pos=1.0, fpr_neg=1.0, a_pinch=1.0, lamda_s=1.0, lamda_c=1.0, lamda_a=1.0, lamda_k=1.0, c_s=1.0, c_c=1.0, c_a=1.0, c_k=1.0, theta_p__plus=1.0, theta_p__neg=1.0, theta_pc__plus=1.0, theta_pc__neg=1.0, res__pos=1.0, res__neg=1.0, theta_u__plus=1.0, theta_u__neg=1.0, d__plus=1.0, d__neg=1.0)
        >>> o3.uniaxial_material.InitStrainMaterial(osi, other=other, init_strain=1.0)
        """
        self.osi = osi
        self.other = other
        self.init_strain = float(init_strain)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.other.tag, self.init_strain]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class InitStressMaterial(UniaxialMaterialBase):
    """
    The InitStressMaterial UniaxialMaterial Class
    
    This command is used to construct an Initial Stress material object. The stress-strain behaviour for this material
    is defined by another material. Initial Stress Material enables definition of initial stress for the material under
    consideration. The strain that corresponds to the initial stress will be calculated from the other material.
    """
    op_type = 'InitStressMaterial'

    def __init__(self, osi, other, init_stress):
        """
        Initial method for InitStressMaterial

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        other: obj
            Object of the other material
        init_stress: float
            Initial stress

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> other = o3.uniaxial_material.ModIMKPinching(osi, k0=1.0, as__plus=1.0, as__neg=1.0, my__plus=1.0, my__neg=1.0,
        >>>                                             fpr_pos=1.0, fpr_neg=1.0, a_pinch=1.0, lamda_s=1.0, lamda_c=1.0,
        >>>                                             lamda_a=1.0, lamda_k=1.0, c_s=1.0, c_c=1.0, c_a=1.0, c_k=1.0,
        >>>                                             theta_p__plus=1.0, theta_p__neg=1.0, theta_pc__plus=1.0,
        >>>                                             theta_pc__neg=1.0, res__pos=1.0, res__neg=1.0, theta_u__plus=1.0,
        >>>                                             theta_u__neg=1.0, d__plus=1.0, d__neg=1.0)
        >>> o3.uniaxial_material.InitStressMaterial(osi, other=other, init_stress=1.0)
        """
        self.osi = osi
        self.other = other
        self.init_stress = float(init_stress)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.other.tag, self.init_stress]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def set_f0(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'F0', value, ele, eles)


class PathIndependent(UniaxialMaterialBase):
    """
    The PathIndependent UniaxialMaterial Class
    
    This command is to create a PathIndependent material
    """
    op_type = 'PathIndependent'

    def __init__(self, osi, other):
        """
        Initial method for PathIndependent

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        other: obj
            A pre-defined material

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> other = o3.uniaxial_material.Hardening(osi, e_mod=1.0, sigma_y=1.0, h_iso=1.0, h_kin=1.0, eta=0.0)
        >>> o3.uniaxial_material.PathIndependent(osi, other=other)
        """
        self.osi = osi
        self.other = other
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.other.tag]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class ECC01(UniaxialMaterialBase):
    """
    The ECC01 UniaxialMaterial Class
    
    This command is used to construct a uniaxial Engineered Cementitious Composites (ECC)material object based on the
    ECC material model of Han, et al. (see references). Reloading in tension and compression is linear.
    """
    op_type = 'ECC01'

    def __init__(self, osi, sigt0, epst0, sigt1, epst1, epst2, sigc0, epsc0, epsc1, alpha_t1, alpha_t2, alpha_c, alpha_cu, beta_t, beta_c):
        """
        Initial method for ECC01

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        sigt0: float
            Tensile cracking stress
        epst0: float
            Strain at tensile cracking stress
        sigt1: float
            Peak tensile stress
        epst1: float
            Strain at peak tensile stress
        epst2: float
            Ultimate tensile strain
        sigc0: float
            Compressive strength (see notes)
        epsc0: float
            Strain at compressive strength (see notes)
        epsc1: float
            Ultimate compressive strain (see notes)
        alpha_t1: float
            Exponent of the unloading curve in tensile strain hardening region
        alpha_t2: float
            Exponent of the unloading curve in tensile softening region
        alpha_c: float
            Exponent of the unloading curve in the compressive softening
        alpha_cu: float
            Exponent of the compressive softening curve (use 1 for linear softening)
        beta_t: float
            Parameter to determine permanent strain in tension
        beta_c: float
            Parameter to determine permanent strain in compression

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.ECC01(osi, sigt0=1.0, epst0=1.0, sigt1=1.0, epst1=1.0, epst2=1.0, sigc0=1.0, epsc0=1.0, epsc1=1.0, alpha_t1=1.0, alpha_t2=1.0, alpha_c=1.0, alpha_cu=1.0, beta_t=1.0, beta_c=1.0)
        """
        self.osi = osi
        self.sigt0 = float(sigt0)
        self.epst0 = float(epst0)
        self.sigt1 = float(sigt1)
        self.epst1 = float(epst1)
        self.epst2 = float(epst2)
        self.sigc0 = float(sigc0)
        self.epsc0 = float(epsc0)
        self.epsc1 = float(epsc1)
        self.alpha_t1 = float(alpha_t1)
        self.alpha_t2 = float(alpha_t2)
        self.alpha_c = float(alpha_c)
        self.alpha_cu = float(alpha_cu)
        self.beta_t = float(beta_t)
        self.beta_c = float(beta_c)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.sigt0, self.epst0, self.sigt1, self.epst1, self.epst2, self.sigc0, self.epsc0, self.epsc1, self.alpha_t1, self.alpha_t2, self.alpha_c, self.alpha_cu, self.beta_t, self.beta_c]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class SelfCentering(UniaxialMaterialBase):
    """
    The SelfCentering UniaxialMaterial Class
    
    This command is used to construct a uniaxial self-centering (flag-shaped) material object with optional
    non-recoverable slip behaviour and an optional stiffness increase at high strains (bearing behaviour).
    """
    op_type = 'SelfCentering'

    def __init__(self, osi, k1, k2, sig_act, beta, eps_slip=0, eps_bear=0, r_bear: float=None):
        """
        Initial method for SelfCentering

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        k1: float
            Initial stiffness
        k2: float
            Post-activation stiffness (0<   ``k2``<   ``k1``)
        sig_act: float
            Forward activation stress/force
        beta: float
            Ratio of forward to reverse activation stress/force
        eps_slip: float, optional
            Slip strain/deformation (if    ``epsslip`` = 0, there will be no slippage)
        eps_bear: float, optional
            Bearing strain/deformation (if    ``epsbear`` = 0, there will be no bearing)
        r_bear: float (default=True), optional
            Ratio of bearing stiffness to initial stiffness    ``k1``

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.SelfCentering(osi, k1=1.0, k2=1.0, sig_act=1.0, beta=1.0, eps_slip=0, eps_bear=0, r_bear=None)
        """
        self.osi = osi
        self.k1 = float(k1)
        self.k2 = float(k2)
        self.sig_act = float(sig_act)
        self.beta = float(beta)
        self.eps_slip = float(eps_slip)
        self.eps_bear = float(eps_bear)
        if r_bear is None:
            self.r_bear = None
        else:
            self.r_bear = float(r_bear)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.k1, self.k2, self.sig_act, self.beta, self.eps_slip, self.eps_bear]
        special_pms = ['r_bear']
        packets = [False]
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


class Viscous(UniaxialMaterialBase):
    """
    The Viscous UniaxialMaterial Class
    
    This command is used to construct a uniaxial viscous material object. stress =C(strain-rate)^alpha
    """
    op_type = 'Viscous'

    def __init__(self, osi, big_c, alpha):
        """
        Initial method for Viscous

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        big_c: float
            Damping coeficient
        alpha: float
            Power factor (=1 means linear damping)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.Viscous(osi, big_c=1.0, alpha=1.0)
        """
        self.osi = osi
        self.big_c = float(big_c)
        self.alpha = float(alpha)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.big_c, self.alpha]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def set_c(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'C', value, ele, eles)

    def set_alpha(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'Alpha', value, ele, eles)

    def set_min_vel(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'minVel', value, ele, eles)


class BoucWen(UniaxialMaterialBase):
    """
    The BoucWen UniaxialMaterial Class
    
    This command is used to construct a uniaxial Bouc-Wen smooth hysteretic material object. This material model is an
    extension of the original Bouc-Wen model that includes stiffness and strength degradation (Baber and Noori (1985)).
    """
    op_type = 'BoucWen'

    def __init__(self, osi, alpha, ko, n, gamma, beta, ao, delta_a, delta_nu, delta_eta):
        """
        Initial method for BoucWen

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        alpha: float
            Ratio of post-yield stiffness to the initial elastic stiffenss (0< alpha <1)
        ko: float
            Initial elastic stiffness
        n: float
            Parameter that controls transition from linear to nonlinear range (as n increases the transition becomes
            sharper; n is usually grater or equal to 1)
        gamma: float
            Parameters that control shape of hysteresis loop; depending on the values of gamma and beta softening,
            hardening or quasi-linearity can be simulated (look at the notes)
        beta: float
            Parameters that control shape of hysteresis loop; depending on the values of gamma and beta softening,
            hardening or quasi-linearity can be simulated (look at the notes)
        ao: float
            Parameters that control tangent stiffness
        delta_a: float
            Parameters that control tangent stiffness
        delta_nu: float
            Parameters that control material degradation
        delta_eta: float
            Parameters that control material degradation

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.BoucWen(osi, alpha=1.0, ko=1.0, n=1.0, gamma=1.0, beta=1.0, ao=1.0, delta_a=1.0, delta_nu=1.0, delta_eta=1.0)
        """
        self.osi = osi
        self.alpha = float(alpha)
        self.ko = float(ko)
        self.n = float(n)
        self.gamma = float(gamma)
        self.beta = float(beta)
        self.ao = float(ao)
        self.delta_a = float(delta_a)
        self.delta_nu = float(delta_nu)
        self.delta_eta = float(delta_eta)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.alpha, self.ko, self.n, self.gamma, self.beta, self.ao, self.delta_a, self.delta_nu, self.delta_eta]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def set_alpha(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'alpha', value, ele, eles)

    def set_ko(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'ko', value, ele, eles)

    def set_n(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'n', value, ele, eles)

    def set_gamma(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'gamma', value, ele, eles)

    def set_beta(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'beta', value, ele, eles)

    def set_ao(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'Ao', value, ele, eles)

    def set_delta_a(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'deltaA', value, ele, eles)

    def set_delta_nu(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'deltaNu', value, ele, eles)

    def set_delta_eta(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'deltaEta', value, ele, eles)


class BWBN(UniaxialMaterialBase):
    """
    The BWBN UniaxialMaterial Class
    
    This command is used to construct a uniaxial Bouc-Wen pinching hysteretic material object. This material model is an
    extension of the original Bouc-Wen model that includes pinching (Baber and Noori (1986) and Foliente (1995)).
    """
    op_type = 'BWBN'

    def __init__(self, osi, alpha, ko, n, gamma, beta, ao, q, zetas, p, shi, delta_shi, lamb, tol, max_iter):
        """
        Initial method for BWBN

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        alpha: float
            Ratio of post-yield stiffness to the initial elastic stiffenss (0< alpha <1)
        ko: float
            Initial elastic stiffness
        n: float
            Parameter that controls transition from linear to nonlinear range (as n increases the transition becomes
            sharper; n is usually grater or equal to 1)
        gamma: float
            Parameters that control shape of hysteresis loop; depending on the values of gamma and beta softening,
            hardening or quasi-linearity can be simulated (look at the boucwen material)
        beta: float
            Parameters that control shape of hysteresis loop; depending on the values of gamma and beta softening,
            hardening or quasi-linearity can be simulated (look at the boucwen material)
        ao: float
            Parameter that controls tangent stiffness
        q: float
            Parameters that control pinching
        zetas: float
            Parameters that control pinching
        p: float
            Parameters that control pinching
        shi: float
            Parameters that control pinching
        delta_shi: float
            Parameters that control pinching
        lamb: float
            Parameters that control pinching
        tol: float
            Tolerance
        max_iter: float
            Maximum iterations

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.BWBN(osi, alpha=1.0, ko=1.0, n=1.0, gamma=1.0, beta=1.0, ao=1.0, q=1.0, zetas=1.0, p=1.0, shi=1.0, delta_shi=1.0, lamb=1.0, tol=1.0, max_iter=1.0)
        """
        self.osi = osi
        self.alpha = float(alpha)
        self.ko = float(ko)
        self.n = float(n)
        self.gamma = float(gamma)
        self.beta = float(beta)
        self.ao = float(ao)
        self.q = float(q)
        self.zetas = float(zetas)
        self.p = float(p)
        self.shi = float(shi)
        self.delta_shi = float(delta_shi)
        self.lamb = float(lamb)
        self.tol = float(tol)
        self.max_iter = float(max_iter)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.alpha, self.ko, self.n, self.gamma, self.beta, self.ao, self.q, self.zetas, self.p, self.shi, self.delta_shi, self.lamb, self.tol, self.max_iter]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class AxialSp(UniaxialMaterialBase):
    """
    The AxialSp UniaxialMaterial Class
    
    This command is used to construct a uniaxial AxialSp material object. This material model produces axial
    stress-strain curve of elastomeric bearings.
    """
    op_type = 'AxialSp'

    def __init__(self, osi, sce, fty, fcy, bte, bty, bcy, fcr):
        """
        Initial method for AxialSp

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        sce: float
            Compressive modulus
        fty: float
            Yield stress under tension (   ``fty``) and compression (   ``fcy``) (see note 1)
        fcy: float
            Yield stress under tension (   ``fty``) and compression (   ``fcy``) (see note 1)
        bte: float
            Reduction rate for tensile elastic range (   ``bte``), tensile yielding (   ``bty``) and compressive
            yielding (   ``bcy``) (see note 1)
        bty: float
            Reduction rate for tensile elastic range (   ``bte``), tensile yielding (   ``bty``) and compressive
            yielding (   ``bcy``) (see note 1)
        bcy: float
            Reduction rate for tensile elastic range (   ``bte``), tensile yielding (   ``bty``) and compressive
            yielding (   ``bcy``) (see note 1)
        fcr: float
            Target point stress (see note 1)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.AxialSp(osi, sce=1.0, fty=1.0, fcy=1.0, bte=1.0, bty=1.0, bcy=1.0, fcr=1.0)
        """
        self.osi = osi
        self.sce = float(sce)
        self.fty = float(fty)
        self.fcy = float(fcy)
        self.bte = float(bte)
        self.bty = float(bty)
        self.bcy = float(bcy)
        self.fcr = float(fcr)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.sce, self.fty, self.fcy, self.bte, self.bty, self.bcy, self.fcr]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class AxialSpHD(UniaxialMaterialBase):
    """
    The AxialSpHD UniaxialMaterial Class
    
    This command is used to construct a uniaxial AxialSpHD material object. This material model produces axial
    stress-strain curve of elastomeric bearings including hardening behavior.
    """
    op_type = 'AxialSpHD'

    def __init__(self, osi, sce, fty, fcy, bte, bty, bth, bcy, fcr, ath):
        """
        Initial method for AxialSpHD

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        sce: float
            Compressive modulus
        fty: float
            Yield stress under tension (``fty``) and compression (``fcy``) (see note 1)
        fcy: float
            Yield stress under tension (``fty``) and compression (``fcy``) (see note 1)
        bte: float
            Reduction rate for tensile elastic range (``bte``), tensile yielding (``bty``), tensile hardening (  
            ``bth``) and compressive yielding (``bcy``) (see note 1)
        bty: float
            Reduction rate for tensile elastic range (``bte``), tensile yielding (``bty``), tensile hardening (  
            ``bth``) and compressive yielding (``bcy``) (see note 1)
        bth: float
            Reduction rate for tensile elastic range (``bte``), tensile yielding (``bty``), tensile hardening (  
            ``bth``) and compressive yielding (``bcy``) (see note 1)
        bcy: float
            Reduction rate for tensile elastic range (``bte``), tensile yielding (``bty``), tensile hardening (  
            ``bth``) and compressive yielding (``bcy``) (see note 1)
        fcr: float
            Target point stress (see note 1)
        ath: float
            Hardening strain ratio to yield strain

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.AxialSpHD(osi, sce=1.0, fty=1.0, fcy=1.0, bte=1.0, bty=1.0, bth=1.0, bcy=1.0, fcr=1.0, ath=1.0)
        """
        self.osi = osi
        self.sce = float(sce)
        self.fty = float(fty)
        self.fcy = float(fcy)
        self.bte = float(bte)
        self.bty = float(bty)
        self.bth = float(bth)
        self.bcy = float(bcy)
        self.fcr = float(fcr)
        self.ath = float(ath)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.sce, self.fty, self.fcy, self.bte, self.bty, self.bth, self.bcy, self.fcr, self.ath]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class CFSWSWP(UniaxialMaterialBase):
    """
    The CFSWSWP UniaxialMaterial Class
    
    This command is used to construct a uniaxialMaterial model that simulates the hysteresis response (Shear
    strength-Lateral displacement) of a wood-sheathed cold-formed steel shear wall panel (CFS-SWP). The
    hysteresis model has smooth curves and takes into account the strength and stiffness degradation,
    as well as pinching effect.This uniaxialMaterial gives results in Newton and Meter units, for
    strength and displacement, respectively.
    """
    op_type = 'CFSWSWP'

    def __init__(self, osi, height, width, fut, tf, ife, ifi, ts, np, ds, vs, sc, nc, otype, opening_area, opening_length):
        """
        Initial method for CFSWSWP

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        height: float
            Swp's height (mm)
        width: float
            Swp's width (mm)
        fut: float
            Tensile strength of framing members (mpa)
        tf: float
            Framing thickness (mm)
        ife: float
            Moment of inertia of the double end-stud (mm4)
        ifi: float
            Moment of inertia of the intermediate stud (mm4)
        ts: float
            Sheathing thickness (mm)
        np: float
            Sheathing number (one or two sides sheathed)
        ds: float
            Screws diameter (mm)
        vs: float
            Screws shear strength (n)
        sc: float
            Screw spacing on the swp perimeter (mm)
        nc: float
            Total number of screws located on the swp perimeter
        otype: int
            Integer identifier used to define wood sheathing type (dfp=1, osb=2, csp=3)
        opening_area: float
            Total area of openings (mm2)
        opening_length: float
            Cumulative length of openings (mm)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.CFSWSWP(osi, height=1.0, width=1.0, fut=1.0, tf=1.0, ife=1.0, ifi=1.0, ts=1.0, np=1.0, ds=1.0, vs=1.0, sc=1.0, nc=1.0, otype=1, opening_area=1.0, opening_length=1.0)
        """
        self.osi = osi
        self.height = float(height)
        self.width = float(width)
        self.fut = float(fut)
        self.tf = float(tf)
        self.ife = float(ife)
        self.ifi = float(ifi)
        self.ts = float(ts)
        self.np = float(np)
        self.ds = float(ds)
        self.vs = float(vs)
        self.sc = float(sc)
        self.nc = float(nc)
        self.otype = int(otype)
        self.opening_area = float(opening_area)
        self.opening_length = float(opening_length)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.height, self.width, self.fut, self.tf, self.ife, self.ifi, self.ts, self.np, self.ds, self.vs, self.sc, self.nc, self.otype, self.opening_area, self.opening_length]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class CFSSSWP(UniaxialMaterialBase):
    """
    The CFSSSWP UniaxialMaterial Class
    
    This command is used to construct a uniaxialMaterial model that simulates the hysteresis response (Shear
    strength-lateral Displacement) of a Steel-Sheathed Cold-Formed Steel Shear Wall Panel (CFS-SWP). The
    hysteresis model has smooth curves and takes into account the strength and stiffness degradation,
    as well as pinching effect.This uniaxialMaterial gives results in Newton and Meter units, for
    strength and displacement, respectively.
    """
    op_type = 'CFSSSWP'

    def __init__(self, osi, height, width, fuf, fyf, tf, af, fus, fys, ts, np, ds, vs, sc, dt, opening_area, opening_length):
        """
        Initial method for CFSSSWP

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        height: float
            Swp's height (mm)
        width: float
            Swp's width (mm)
        fuf: float
            Tensile strength of framing members (mpa)
        fyf: float
            Yield strength of framing members (mpa)
        tf: float
            Framing thickness (mm)
        af: float
            Framing cross section area (mm2)
        fus: float
            Tensile strength of steel sheet sheathing (mpa)
        fys: float
            Yield strength of steel sheet sheathing (mpa)
        ts: float
            Sheathing thickness (mm)
        np: float
            Sheathing number (one or two sides sheathed)
        ds: float
            Screws diameter (mm)
        vs: float
            Screws shear strength (n)
        sc: float
            Screw spacing on the swp perimeter (mm)
        dt: float
            Anchor bolt's diameter (mm)
        opening_area: float
            Total area of openings (mm2)
        opening_length: float
            Cumulative length of openings (mm)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.CFSSSWP(osi, height=1.0, width=1.0, fuf=1.0, fyf=1.0, tf=1.0, af=1.0, fus=1.0, fys=1.0, ts=1.0, np=1.0, ds=1.0, vs=1.0, sc=1.0, dt=1.0, opening_area=1.0, opening_length=1.0)
        """
        self.osi = osi
        self.height = float(height)
        self.width = float(width)
        self.fuf = float(fuf)
        self.fyf = float(fyf)
        self.tf = float(tf)
        self.af = float(af)
        self.fus = float(fus)
        self.fys = float(fys)
        self.ts = float(ts)
        self.np = float(np)
        self.ds = float(ds)
        self.vs = float(vs)
        self.sc = float(sc)
        self.dt = float(dt)
        self.opening_area = float(opening_area)
        self.opening_length = float(opening_length)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.height, self.width, self.fuf, self.fyf, self.tf, self.af, self.fus, self.fys, self.ts, self.np, self.ds, self.vs, self.sc, self.dt, self.opening_area, self.opening_length]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)
