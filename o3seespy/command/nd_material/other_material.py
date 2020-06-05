from o3seespy.command.nd_material.base_material import NDMaterialBase
from o3seespy.command.common import update_material_stage


class PM4Sand(NDMaterialBase):
    op_type = "PM4Sand"

    def __init__(self, osi, d_r, g_o, h_po, den, p_atm, h_o=-1.0, e_max=0.8, e_min=0.5, n_b=0.5, n_d=0.1, a_do=None,
                 z_max=None, c_z=250.0, c_e=None, phi_cv=33.0, nu=0.3, g_degr=2.0, c_dr=None, c_kaf=None, q_bolt=10.0,
                 r_bolt=1.5, m_par=0.01, f_sed=None, p_sed=None):
        """
        This command is used to construct a 2-dimensional PM4Sand material.

        Parameters
        ==========
        d_r: float
            Relative density, in fraction
        g_o: float
            Shear modulus constant
        h_po: float
            Contraction rate parameter
        den: float
            Mass density of the material
        p_atm: float, optional
            Atmospheric pressure
        h_o: float, optional
            Variable that adjusts the ratio of plastic modulus to elastic modulus
        e_max: float, optional
            Maximum and minimum void ratios
        e_min: float, optional
            Maximum and minimum void ratios
        n_b: float, optional
            Bounding surface parameter, :math:`nb \ge 0`
        n_d: float, optional
            Dilatancy surface parameter :math:`nd \ge 0`
        a_do: float, optional
            Dilatancy parameter, will be computed at the time of initialization if input value is negative
        z_max: float, optional
            Fabric-dilatancy tensor parameter
        c_z: float, optional
            Fabric-dilatancy tensor parameter
        c_e: float, optional
            Variable that adjusts the rate of strain accumulation in cyclic loading
        phi_cv: float, optional
            Critical state effective friction angle
        nu: float, optional
            Poisson's ratio
        c_gd: float, optional
            Variable that adjusts degradation of elastic modulus with accumulation of fabric
        c_dr: float,   optional
            Variable that controls the rotated dilatancy surface
        c_kaf: float,  optional
            Variable that controls the effect that sustained
                                          con shear stresses have on plastic modulus
        q_bolt: float,     optional
            Critical state line parameter
        r_bolt: float,     optional
            Critical state line parameter
        m_par: float,     optional
            Yield surface constant (radius of yield surface in stress ratio space)
        f_sed: float, optional,
            Variable that controls the minimum value the reduction factor of the elastic moduli
            can get during reconsolidation
        p_sed: float, optional
            Mean effective stress up to which reconsolidation strains are enhanced
        """
        self.osi = osi
        self.d_r = float(d_r)
        self.g_o = float(g_o)
        self.h_po = float(h_po)
        self.den = float(den)
        self.p_atm = float(p_atm)
        self.h_o = float(h_o)
        self.e_max = float(e_max)
        self.e_min = float(e_min)
        self.n_b = float(n_b)
        self.n_d = float(n_d)
        if a_do is None:
            self.a_do = -1.0
        else:
            self.a_do = float(a_do)
        if z_max is None:
            self.z_max = -1.0
        else:
            self.z_max = float(z_max)
        self.c_z = float(c_z)
        if c_e is None:
            self.c_e = -1.0
        else:
            self.c_e = float(c_e)
        self.phi_cv = float(phi_cv)
        self.nu = float(nu)
        self.g_degr = float(g_degr)
        if c_dr is None:
            self.c_dr = -1.0
        else:
            self.c_dr = float(c_dr)
        if c_kaf is None:
            self.c_kaf = -1.0
        else:
            self.c_kaf = float(c_kaf)
        self.q_bolt = float(q_bolt)
        self.r_bolt = float(r_bolt)
        self.m_par = float(m_par)
        if f_sed is None:
            self.f_sed = -1.0
        else:
            self.f_sed = float(f_sed)
        if p_sed is None:
            self.p_sed = -1.0
        else:
            self.p_sed = float(p_sed)

        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat

        self._parameters = [self.op_type, self._tag, self.d_r, self.g_o, self.h_po, self.den, self.p_atm, self.h_o,
                            self.e_max, self.e_min, self.n_b, self.n_d, self.a_do, self.z_max, self.c_z, self.c_e, self.phi_cv,
                            self.nu, self.g_degr, self.c_dr, self.c_kaf, self.q_bolt, self.r_bolt, self.m_par, self.f_sed,
                            self.p_sed]

        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def pre_dynamic(self, osi):  # deprecated
        update_material_stage(osi, self, stage=1)
        # opw.set_parameter(osi, value=0, eles=[ele], args=['FirstCall', 1])

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





