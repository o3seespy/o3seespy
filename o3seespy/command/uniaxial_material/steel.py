from o3seespy.command.uniaxial_material.base_material import UniaxialMaterialBase



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



class Steel02(UniaxialMaterialBase):
    """
    The Steel02 UniaxialMaterial Class
    
    This command is used to construct a uniaxial Giuffre-Menegotto-Pinto steel material object with isotropic strain
    hardening.
    """
    op_type = 'Steel02'

    def __init__(self, osi, fy, e0, b, params, a1: float=None, a2=1.0, a3: float=None, a4=1.0, sig_init=0.0):
        r"""
        Initial method for Steel02

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        fy: float
            Yield strength
        e0: float
            Initial elastic tangent
        b: float
            Strain-hardening ratio (ratio between post-yield tangent and initial elastic tangent)
        params: list
            Parameters to control the transition from elastic to plastic branches. ``params=[r0,cr1,cr2]``. recommended
            values: r0=between 10 and 20, cr1=0.925, cr2=0.15
        a1: float (default=True), optional
            Isotropic hardening parameter, increase of compression yield envelope as proportion of yield strength after
            a plastic strain of :math:`a_2*(f_y/e_0)` 
        a2: float, optional
            Isotropic hardening parameter
        a3: float (default=True), optional
            Isotropic hardening parameter, increase of tension yield envelope as proportion of yield strength after a
            plastic strain of :math:`a_4*(f_y/e_0)`. 
        a4: float, optional
            Isotropic hardening parameter (see explanation
        sig_init: float, optional
            Initial stress value (optional, default: 0.0) the strain is calculated from ``epsp=siginit/e`` :: if
            (siginit!= 0.0) { double epsinit = siginit/e; eps = trialstrain+epsinit; } else { eps = trialstrain; }

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.Steel02(osi, fy=1.0, e0=1.0, b=1.0, params=[15, 0.925, 0.15])
        """
        self.osi = osi
        self.fy = float(fy)
        self.e0 = float(e0)
        self.b = float(b)
        self.params = params
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
        self.sig_init = float(sig_init)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.fy, self.e0, self.b, *self.params]
        special_pms = ['a1', 'a2', 'a3', 'a4', 'sig_init']
        packets = [False, False, False, False, False]
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


class Hysteretic(UniaxialMaterialBase):
    """
    The Hysteretic UniaxialMaterial Class
    
    This command is used to construct a uniaxial bilinear hysteretic material object with pinching of force and
    deformation, damage due to ductility and energy, and degraded unloading stiffness based on ductility.
    """
    op_type = 'Hysteretic'

    def __init__(self, osi, p1, p2, p3, n1, n2, n3, pinch_x, pinch_y, damage1, damage2, beta=0.0):
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
        self.n1 = n1
        self.n2 = n2
        self.n3 = n3
        self.pinch_x = float(pinch_x)
        self.pinch_y = float(pinch_y)
        self.damage1 = float(damage1)
        self.damage2 = float(damage2)
        self.beta = float(beta)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, *self.p1, *self.p2]
        special_pms = ['p3', 'n1', 'n2', 'n3', 'pinch_x', 'pinch_y', 'damage1', 'damage2', 'beta']
        packets = [True, True, True, True, False, False, False, False, False]
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


class ReinforcingSteelGABuck(UniaxialMaterialBase):
    """
    The ReinforcingSteelGABuck UniaxialMaterial Class
    
    This command is used to construct a ReinforcingSteel uniaxial material object. This object is intended to be used in
    a reinforced concrete fiber section as the steel reinforcing material.
    """
    op_type = 'ReinforcingSteel'

    def __init__(self, osi, fy, fu, es, esh, eps_sh, eps_ult, lsr, beta, r, gamma):
        """
        Initial method for ReinforcingSteelGABuck

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        fy: float
            Yield stress in tension
        fu: float
            Ultimate stress in tension
        es: float
            Initial elastic tangent
        esh: float
            Tangent at initial strain hardening
        eps_sh: float
            Strain corresponding to initial strain hardening
        eps_ult: float
            Strain at peak stress
        lsr: float
            Slenderness ratio
        beta: float
            Amplification factor for the buckled stress strain curve.
        r: float
            Buckling reduction factor r can be a real number between [0.0 and 1.0] r=1.0 full reduction (no buckling)
            r=0.0 no reduction 0.0<r<1.0 linear interpolation between buckled and unbuckled curves
        gamma: float
            Buckling constant

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.ReinforcingSteelGABuck(osi, fy=1.0, fu=1.0, es=1.0, esh=1.0, eps_sh=1.0, eps_ult=1.0, lsr=1.0, beta=1.0, r=1.0, gamma=1.0)
        """
        self.osi = osi
        self.fy = float(fy)
        self.fu = float(fu)
        self.es = float(es)
        self.esh = float(esh)
        self.eps_sh = float(eps_sh)
        self.eps_ult = float(eps_ult)
        self.lsr = float(lsr)
        self.beta = float(beta)
        self.r = float(r)
        self.gamma = float(gamma)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.fy, self.fu, self.es, self.esh, self.eps_sh, self.eps_ult, '-GABuck', self.lsr, self.beta, self.r, self.gamma]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

class ReinforcingSteelDMBuck(UniaxialMaterialBase):
    """
    The ReinforcingSteelDMBuck UniaxialMaterial Class
    
    This command is used to construct a ReinforcingSteel uniaxial material object. This object is intended to be used in
    a reinforced concrete fiber section as the steel reinforcing material.
    """
    op_type = 'ReinforcingSteel'

    def __init__(self, osi, fy, fu, es, esh, eps_sh, eps_ult, lsr_2, alpha=1.0):
        """
        Initial method for ReinforcingSteelDMBuck

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        fy: float
            Yield stress in tension
        fu: float
            Ultimate stress in tension
        es: float
            Initial elastic tangent
        esh: float
            Tangent at initial strain hardening
        eps_sh: float
            Strain corresponding to initial strain hardening
        eps_ult: float
            Strain at peak stress
        lsr_2: None
            
        alpha: float, optional
            Coffin-manson constant a

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.ReinforcingSteelDMBuck(osi, fy=1.0, fu=1.0, es=1.0, esh=1.0, eps_sh=1.0, eps_ult=1.0, lsr_2=1, alpha=1.0)
        """
        self.osi = osi
        self.fy = float(fy)
        self.fu = float(fu)
        self.es = float(es)
        self.esh = float(esh)
        self.eps_sh = float(eps_sh)
        self.eps_ult = float(eps_ult)
        self.lsr_2 = lsr_2
        self.alpha = float(alpha)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.fy, self.fu, self.es, self.esh, self.eps_sh, self.eps_ult, '-DMBuck', self.lsr_2, self.alpha]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

class ReinforcingSteelCMFatigue(UniaxialMaterialBase):
    """
    The ReinforcingSteelCMFatigue UniaxialMaterial Class
    
    This command is used to construct a ReinforcingSteel uniaxial material object. This object is intended to be used in
    a reinforced concrete fiber section as the steel reinforcing material.
    """
    op_type = 'ReinforcingSteel'

    def __init__(self, osi, fy, fu, es, esh, eps_sh, eps_ult, cf, alpha_2, cd):
        """
        Initial method for ReinforcingSteelCMFatigue

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        fy: float
            Yield stress in tension
        fu: float
            Ultimate stress in tension
        es: float
            Initial elastic tangent
        esh: float
            Tangent at initial strain hardening
        eps_sh: float
            Strain corresponding to initial strain hardening
        eps_ult: float
            Strain at peak stress
        cf: float
            Coffin-manson constant c
        alpha_2: None
            
        cd: float
            Cyclic strength reduction constant

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.ReinforcingSteelCMFatigue(osi, fy=1.0, fu=1.0, es=1.0, esh=1.0, eps_sh=1.0, eps_ult=1.0, cf=1.0, alpha_2=1, cd=1.0)
        """
        self.osi = osi
        self.fy = float(fy)
        self.fu = float(fu)
        self.es = float(es)
        self.esh = float(esh)
        self.eps_sh = float(eps_sh)
        self.eps_ult = float(eps_ult)
        self.cf = float(cf)
        self.alpha_2 = alpha_2
        self.cd = float(cd)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.fy, self.fu, self.es, self.esh, self.eps_sh, self.eps_ult, '-CMFatigue', self.cf, self.alpha_2, self.cd]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

class ReinforcingSteelIsoHard(UniaxialMaterialBase):
    """
    The ReinforcingSteelIsoHard UniaxialMaterial Class
    
    This command is used to construct a ReinforcingSteel uniaxial material object. This object is intended to be used in
    a reinforced concrete fiber section as the steel reinforcing material.
    """
    op_type = 'ReinforcingSteel'

    def __init__(self, osi, fy, fu, es, esh, eps_sh, eps_ult, a1=4.3, limit=1.0):
        """
        Initial method for ReinforcingSteelIsoHard

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        fy: float
            Yield stress in tension
        fu: float
            Ultimate stress in tension
        es: float
            Initial elastic tangent
        esh: float
            Tangent at initial strain hardening
        eps_sh: float
            Strain corresponding to initial strain hardening
        eps_ult: float
            Strain at peak stress
        a1: float, optional
            Hardening constant (default = 4.3)
        limit: float, optional
            Limit for the reduction of the yield plateau. % of original plateau length to remain (0.01 < limit < 1.0 )
            limit =1.0, then no reduction takes place (default =0.01)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.ReinforcingSteelIsoHard(osi, fy=1.0, fu=1.0, es=1.0, esh=1.0, eps_sh=1.0, eps_ult=1.0, a1=4.3, limit=1.0)
        """
        self.osi = osi
        self.fy = float(fy)
        self.fu = float(fu)
        self.es = float(es)
        self.esh = float(esh)
        self.eps_sh = float(eps_sh)
        self.eps_ult = float(eps_ult)
        self.a1 = float(a1)
        self.limit = float(limit)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.fy, self.fu, self.es, self.esh, self.eps_sh, self.eps_ult, '-IsoHard', self.a1, self.limit]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

class ReinforcingSteelMPCurveParams(UniaxialMaterialBase):
    """
    The ReinforcingSteelMPCurveParams UniaxialMaterial Class
    
    This command is used to construct a ReinforcingSteel uniaxial material object. This object is intended to be used in
    a reinforced concrete fiber section as the steel reinforcing material.
    """
    op_type = 'ReinforcingSteel'

    def __init__(self, osi, fy, fu, es, esh, eps_sh, eps_ult, r1=0.333, r2=18.0, r3=4.0):
        """
        Initial method for ReinforcingSteelMPCurveParams

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        fy: float
            Yield stress in tension
        fu: float
            Ultimate stress in tension
        es: float
            Initial elastic tangent
        esh: float
            Tangent at initial strain hardening
        eps_sh: float
            Strain corresponding to initial strain hardening
        eps_ult: float
            Strain at peak stress
        r1: float, optional
            (default = 0.333)
        r2: float, optional
            (default = 18)
        r3: float, optional
            (default = 4)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.ReinforcingSteelMPCurveParams(osi, fy=1.0, fu=1.0, es=1.0, esh=1.0, eps_sh=1.0, eps_ult=1.0, r1=0.333, r2=18.0, r3=4.0)
        """
        self.osi = osi
        self.fy = float(fy)
        self.fu = float(fu)
        self.es = float(es)
        self.esh = float(esh)
        self.eps_sh = float(eps_sh)
        self.eps_ult = float(eps_ult)
        self.r1 = float(r1)
        self.r2 = float(r2)
        self.r3 = float(r3)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.fy, self.fu, self.es, self.esh, self.eps_sh, self.eps_ult, '-MPCurveParams', self.r1, self.r2, self.r3]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class DoddRestrepo(UniaxialMaterialBase):
    """
    The DoddRestrepo UniaxialMaterial Class
    
    This command is used to construct a Dodd-Restrepo steel material
    """
    op_type = 'Dodd_Restrepo'

    def __init__(self, osi, fy, fsu, esh, esu, youngs, eshi, fshi, omega_fac=1.0):
        """
        Initial method for DoddRestrepo

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        fy: float
            Yield strength
        fsu: float
            Ultimate tensile strength (uts)
        esh: float
            Tensile strain at initiation of strain hardening
        esu: float
            Tensile strain at the uts
        youngs: float
            Modulus of elasticity
        eshi: float
            Tensile strain for a point on strain hardening curve, recommended range of values for eshi: [ (esu +
            5*esh)/6, (esu + 3*esh)/4]
        fshi: float
            Tensile stress at point on strain hardening curve corresponding to eshi
        omega_fac: float, optional
            Roundedness factor for bauschinger curve in cycle reversals from the strain hardening curve. range: [0.75,
            1.15]. largest value tends to near a bilinear bauschinger curve. default = 1.0.

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.DoddRestrepo(osi, fy=1.0, fsu=1.0, esh=1.0, esu=1.0, youngs=1.0, eshi=1.0, fshi=1.0, omega_fac=1.0)
        """
        self.osi = osi
        self.fy = float(fy)
        self.fsu = float(fsu)
        self.esh = float(esh)
        self.esu = float(esu)
        self.youngs = float(youngs)
        self.eshi = float(eshi)
        self.fshi = float(fshi)
        self.omega_fac = float(omega_fac)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.fy, self.fsu, self.esh, self.esu, self.youngs, self.eshi, self.fshi, self.omega_fac]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class RambergOsgoodSteel(UniaxialMaterialBase):
    """
    The RambergOsgoodSteel UniaxialMaterial Class
    
    This command is used to construct a Ramberg-Osgood steel material object.
    """
    op_type = 'RambergOsgoodSteel'

    def __init__(self, osi, fy, e0, a, n):
        """
        Initial method for RambergOsgoodSteel

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        fy: float
            Yield strength
        e0: float
            Initial elastic tangent
        a: float
            "yield offset" and the commonly used value for a is 0.002
        n: float
            Parameters to control the transition from elastic to plastic branches. and controls the hardening of the
            material by increasing the "n" hardening ratio will be decreased. commonly used values for n are ~5 or greater.

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.RambergOsgoodSteel(osi, fy=1.0, e0=1.0, a=1.0, n=1.0)
        """
        self.osi = osi
        self.fy = float(fy)
        self.e0 = float(e0)
        self.a = float(a)
        self.n = float(n)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.fy, self.e0, self.a, self.n]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class SteelMPF(UniaxialMaterialBase):
    """
    The SteelMPF UniaxialMaterial Class
    
    This command is used to construct a uniaxialMaterial SteelMPF (Kolozvari et al., 2015), which represents the
    well-known uniaxial constitutive nonlinear hysteretic material model for steel proposed by Menegotto and Pinto
    (1973), and extended by Filippou et al. (1983) to include isotropic strain hardening effects.
    """
    op_type = 'SteelMPF'

    def __init__(self, osi, fyp, fyn, e0, bp, bn, params, a1=0.0, a2=1.0, a3=0.0, a4=1.0):
        """
        Initial method for SteelMPF

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        fyp: float
            Yield strength in tension (positive loading direction)
        fyn: float
            Yield strength in compression (negative loading direction)
        e0: float
            Initial tangent modulus
        bp: float
            Strain hardening ratio in tension (positive loading direction)
        bn: float
            Strain hardening ratio in compression (negative loading direction)
        params: list
            Parameters to control the transition from elastic to plastic branches. ``params=[r0,cr1,cr2]``. recommended
            values: ``r0=20``, ``cr1=0.925``, ``cr2=0.15`` or ``cr2=0.0015``
        a1: float, optional
            Isotropic hardening in compression parameter (optional, default = 0.0). shifts compression yield envelope by
            a proportion of compressive yield strength after a maximum plastic tensile strain of a2(fyp/e0)
        a2: float, optional
            Isotropic hardening in compression parameter (optional, default = 1.0).
        a3: float, optional
            Isotropic hardening in tension parameter (optional, default = 0.0). shifts tension yield envelope by a
            proportion of tensile yield strength after a maximum plastic compressive strain of a3(fyn/e0).
        a4: float, optional
            Isotropic hardening in tension parameter (optional, default = 1.0). see explanation of a3.

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.SteelMPF(osi, fyp=1.0, fyn=1.0, e0=1.0, bp=1.0, bn=1.0, params=[1.0, 1.0, 1.0], a1=0.0, a2=1.0, a3=0.0, a4=1.0)
        """
        self.osi = osi
        self.fyp = float(fyp)
        self.fyn = float(fyn)
        self.e0 = float(e0)
        self.bp = float(bp)
        self.bn = float(bn)
        self.params = params
        self.a1 = float(a1)
        self.a2 = float(a2)
        self.a3 = float(a3)
        self.a4 = float(a4)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.fyp, self.fyn, self.e0, self.bp, self.bn, *self.params, self.a1, self.a2, self.a3, self.a4]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class Steel01Thermal(UniaxialMaterialBase):
    """
    The Steel01Thermal UniaxialMaterial Class
    
    
    """
    op_type = 'Steel01Thermal'

    def __init__(self, osi, fy, e0, b, a1, a2, a3, a4):
        r"""
        Initial method for Steel01Thermal

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        fy: float
            Yield strength
        e0: float
            Initial elastic tangent
        b: float
            Strain-hardening ratio (ratio between post-yield tangent and initial elastic tangent)
        a1: float
            Isotropic hardening parameter, increase of compression yield envelope as proportion of yield strength after
            a plastic strain of :math:`a_2*(f_y/e_0)` 
        a2: float
            Isotropic hardening parameter
        a3: float
            Isotropic hardening parameter, increase of tension yield envelope as proportion of yield strength after a
            plastic strain of :math:`a_4*(f_y/e_0)`. 
        a4: float
            Isotropic hardening parameter (see explanation

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.Steel01Thermal(osi, fy=1.0, e0=1.0, b=1.0, a1=1.0, a2=1.0, a3=1.0, a4=1.0)
        """
        self.osi = osi
        self.fy = float(fy)
        self.e0 = float(e0)
        self.b = float(b)
        self.a1 = float(a1)
        self.a2 = float(a2)
        self.a3 = float(a3)
        self.a4 = float(a4)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.fy, self.e0, self.b, self.a1, self.a2, self.a3, self.a4]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def set_fy(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'fy', value, ele, eles)

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
