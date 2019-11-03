from o3seespy.command.nd_material.base_material import NDMaterialBase
import numpy as np


class PressureIndependMultiYield(NDMaterialBase):
    op_type = "PressureIndependMultiYield"

    def __init__(self, osi,  nd, rho, ref_shear_modul, ref_bulk_modul, cohesi, peak_shear_stra, friction_ang=0.,
                 ref_press=100., press_depend_coe=0., no_yield_surf=20, strains=None, ratios=None):
        """
        PressureIndependMultiYield material
        """
        self.nd = nd
        self.rho = float(rho)
        self.ref_shear_modul = float(ref_shear_modul)
        self.ref_bulk_modul = float(ref_bulk_modul)
        self.cohesi = float(cohesi)
        self.peak_shear_stra = float(peak_shear_stra)
        self.friction_ang = float(friction_ang)
        self.ref_press = float(ref_press)
        self.press_depend_coe = float(press_depend_coe)
        assert no_yield_surf < 40
        self.no_yield_surf = int(no_yield_surf)
        if strains is not None:
            assert len(strains) == len(ratios)
            yield_surf = []
            for i in range(len(strains)):
                yield_surf.append(ratios[i])
                yield_surf.append(strains[i])
            self.yield_surf = yield_surf
        else:
            self.yield_surf = None

        osi.n_mats += 1
        self._tag = osi.n_mats

        self._parameters = [self.op_type, self._tag, self.nd, self.rho, self.ref_shear_modul, self.ref_bulk_modul,
                            self.cohesi, self.peak_shear_stra, self.friction_ang, self.ref_press, self.press_depend_coe]

        if self.yield_surf is not None:
            self._parameters.append(self.no_yield_surf)  # from docs 'add a minus sign in front of noYieldSurf'
            self._parameters += list(self.yield_surf)

        else:
            # self._keyword_args['noYieldSurf'] = self.no_yield_surf
            self._parameters.append(self.no_yield_surf)
        self.to_process(osi)


class PM4Sand(NDMaterialBase):
    op_type = "PM4Sand"

    def __init__(self, osi, d_r, g0, hpo, den, p_atm, h0=-1.0, emax=0.8, emin=0.5, nb=0.5, nd=0.1, ado=None,
                 z_max=None, cz=250.0, ce=None, phic=33.0, nu=0.3, cgd=2.0, cdr=None, ckaf=None, big_q=10.0,
                 big_r=1.5, m=0.01, fsed_min=None, p_sedo=None):
        """
        This command is used to construct a 2-dimensional PM4Sand material.

        ================================   ===========================================================================
        ``matTag`` |int|                   integer tag identifying material
        ``d_r`` |float|                     Relative density, in fraction
        ``g0`` |float|                     Shear modulus constant
        ``hpo`` |float|                    Contraction rate parameter
        ``Den`` |float|                    Mass density of the material
        ``P_atm`` |float|                  Optional, Atmospheric pressure
        ``h0`` |float|                     Optional, Variable that adjusts the ratio of plastic modulus
                                          to elastic modulus
        ``emax`` |float|                   Optional, Maximum and minimum void ratios
        ``emin`` |float|                   Optional, Maximum and minimum void ratios
        ``nb`` |float|                     Optional, Bounding surface parameter, :math:`nb \ge 0`
        ``nd`` |float|                     Optional, Dilatancy surface parameter :math:`nd \ge 0`
        ``ado`` |float|                    Optional, Dilatancy parameter, will be computed at the time
                                          of initialization if input value is negative
        ``z_max`` |float|                  Optional, Fabric-dilatancy tensor parameter
        ``cz`` |float|                     Optional, Fabric-dilatancy tensor parameter
        ``ce`` |float|                     Optional, Variable that adjusts the rate of strain accumulation
                                          in cyclic loading
        ``phic`` |float|                   Optional, Critical state effective friction angle
        ``nu`` |float|                     Optional, Poisson's ratio
        ``cgd`` |float|                    Optional, Variable that adjusts degradation of elastic modulus
                                          with accumulation of fabric
        ``cdr`` |float|                    Optional, Variable that controls the rotated dilatancy surface
        ``ckaf`` |float|                   Optional, Variable that controls the effect that sustained
                                          con shear stresses have on plastic modulus
        ``big_q`` |float|                      Optional, Critical state line parameter
        ``big_r`` |float|                      Optional, Critical state line parameter
        ``m`` |float|                      Optional, Yield surface constant (radius of yield surface
                                          in stress ratio space)
        ``fsed_min`` |float|               Optional, Variable that controls the minimum value the
                                          reduction factor of the elastic moduli can get during reconsolidation
        ``p_sedo`` |float|                 Optional, Mean effective stress up to which reconsolidation
                                          strains are enhanced
        ================================   =======================================================================
        """
        self.d_r = float(d_r)
        self.g0 = float(g0)
        self.hpo = float(hpo)
        self.den = float(den)
        self.p_atm = float(p_atm)
        self.h0 = float(h0)
        self.emax = float(emax)
        self.emin = float(emin)
        self.nb = float(nb)
        self.nd = float(nd)
        if ado is None:
            self.ado = -1.0
        else:
            self.ado = float(ado)
        if z_max is None:
            self.z_max = -1.0
        else:
            self.z_max = float(z_max)
        self.cz = float(cz)
        if ce is None:
            self.ce = -1.0
        else:
            self.ce = float(ce)
        self.p_atm = float(p_atm)
        self.phic = float(phic)
        self.nu = float(nu)
        self.cgd = float(cgd)
        if cdr is None:
            self.cdr = -1.0
        else:
            self.cdr = float(cdr)
        if ckaf is None:
            self.ckaf = -1.0
        else:
            self.ckaf = float(ckaf)
        self.big_q = float(big_q)
        self.big_r = float(big_r)
        self.m = float(m)
        if fsed_min is None:
            self.fsed_min = -1.0
        else:
            self.fsed_min = float(fsed_min)
        if p_sedo is None:
            self.p_sedo = -1.0
        else:
            self.p_sedo = float(p_sedo)

        osi.n_mats += 1
        self._tag = osi.n_mats

        self._parameters = [self.op_type, self._tag, self.d_r, self.g0, self.hpo, self.den, self.p_atm, self.h0,
                            self.emax, self.emin, self.nb, self.nd, self.ado, self.z_max, self.cz, self.ce, self.phic,
                            self.nu, self.cgd, self.cdr, self.ckaf, self.big_q, self.big_r, self.m, self.fsed_min,
                            self.p_sedo]

        self.to_process(osi)


class StressDensityModel(NDMaterialBase):
    op_type = "StressDensityModel"

    def __init__(self, osi, den, e_init, big_a, n, nu, a1, b1, a2, b2, a3, b3, fd, mu_not, mu_cyc, sc, big_m, p_atm,
                 ssls=None, hsl=None, ps=None):
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
        self.mu_not = float(mu_not)
        self.mu_cyc = float(mu_cyc)
        self.sc = float(sc)
        self.big_m = float(big_m)
        self.p_atm = float(p_atm)
        if ssls is None:
            self.ssls = np.array([1.0, 10.0, 30.0, 50.0, 100.0, 200.0, 400.0, 400.0, 400.0, 400.0])
        else:
            assert len(ssls) == 10
            self.ssls = np.array(ssls, dtype=np.float)
        if hsl is None:
            self.hsl = 0.895
        else:
            self.hsl = float(hsl)
        if ps is None:
            self.ps = np.array([0.877, 0.877, 0.873, 0.870, 0.860, 0.850, 0.833, 0.833, 0.833, 0.833])
        else:
            self.ps = np.array(ps, dtype=np.float)

        osi.n_mats += 1
        self._tag = osi.n_mats

        self._parameters = [self.osi, self.den, self.e_init, self.big_a, self.n, self.nu, self.a1, self.b1, self.a2,
                            self.b2, self.a3, self.b3, self.fd, self.mu_not, self.mu_cyc, self.sc, self.big_m,
                            self.p_atm, *self.ssls, self.hsl, *self.ps]

        self.to_process(osi)

