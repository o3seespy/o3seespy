from o3seespy.command.nd_material.base_material import NDMaterialBase


class ElasticIsotropic(NDMaterialBase):

    def __init__(self, osi, big_e, v, rho=0.0):
        self.big_e = float(big_e)
        self.v = float(v)
        self.rho = float(rho)
        osi.n_mats += 1
        self._tag = osi.mats
        self._parameters = [self.op_type, self._tag, self.big_e, self.v, self.rho]
        self.to_process(osi)


class ElasticOrthotropic(NDMaterialBase):

    def __init__(self, osi, ex, ey, ez, vxy, vyz, vzx, gxy, gyz, gzx, rho=0.0):
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
        osi.n_mats += 1
        self._tag = osi.mats
        self._parameters = [self.op_type, self._tag, self.ex, self.ey, self.ez, self.vxy, self.vyz, self.vzx, self.gxy, self.gyz, self.gzx, self.rho]
        self.to_process(osi)


class J2Plasticity(NDMaterialBase):

    def __init__(self, osi, big_k, big_g, sig0, sig_inf, delta, big_h):
        self.big_k = float(big_k)
        self.big_g = float(big_g)
        self.sig0 = float(sig0)
        self.sig_inf = float(sig_inf)
        self.delta = float(delta)
        self.big_h = float(big_h)
        osi.n_mats += 1
        self._tag = osi.mats
        self._parameters = [self.op_type, self._tag, self.big_k, self.big_g, self.sig0, self.sig_inf, self.delta, self.big_h]
        self.to_process(osi)


class DrukerPrager(NDMaterialBase):

    def __init__(self, osi, big_k, big_g, sigma_y, rho, rho_bar, kinf, ko, delta1, delta2, big_h, theta, density, atm_pressure=101e3):
        self.big_k = float(big_k)
        self.big_g = float(big_g)
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
        osi.n_mats += 1
        self._tag = osi.mats
        self._parameters = [self.op_type, self._tag, self.big_k, self.big_g, self.sigma_y, self.rho, self.rho_bar, self.kinf, self.ko, self.delta1, self.delta2, self.big_h, self.theta, self.density, self.atm_pressure]
        self.to_process(osi)


class Damage2p(NDMaterialBase):

    def __init__(self, osi, fcc, fct=None, big_e=None, ni=None, gt=None, gc=None, rho_bar=None, big_h=None, theta=None, tangent=None):
        self.fcc = float(fcc)
        self.fct = float(fct)
        self.big_e = float(big_e)
        self.ni = float(ni)
        self.gt = float(gt)
        self.gc = float(gc)
        self.rho_bar = float(rho_bar)
        self.big_h = float(big_h)
        self.theta = float(theta)
        self.tangent = float(tangent)
        osi.n_mats += 1
        self._tag = osi.mats
        self._parameters = [self.op_type, self._tag, self.fcc]
        if getattr(self, 'fct') is not None:
            self._parameters += ['-fct', self.fct]
        if getattr(self, 'big_e') is not None:
            self._parameters += ['-E', self.big_e]
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

    def __init__(self, osi, three_dtag):
        self.three_dtag = int(three_dtag)
        osi.n_mats += 1
        self._tag = osi.mats
        self._parameters = [self.op_type, self._tag, self.three_dtag]
        self.to_process(osi)


class PlaneStrain(NDMaterialBase):

    def __init__(self, osi, three_dtag):
        self.three_dtag = int(three_dtag)
        osi.n_mats += 1
        self._tag = osi.mats
        self._parameters = [self.op_type, self._tag, self.three_dtag]
        self.to_process(osi)


class MultiaxialCyclicPlasticity(NDMaterialBase):

    def __init__(self, osi, rho, big_k, big_g, su, ho, h, m, beta, k_coeff):
        self.rho = float(rho)
        self.big_k = float(big_k)
        self.big_g = float(big_g)
        self.su = float(su)
        self.ho = float(ho)
        self.h = float(h)
        self.m = float(m)
        self.beta = float(beta)
        self.k_coeff = float(k_coeff)
        osi.n_mats += 1
        self._tag = osi.mats
        self._parameters = [self.op_type, self._tag, self.rho, self.big_k, self.big_g, self.su, self.ho, self.h, self.m, self.beta, self.k_coeff]
        self.to_process(osi)


class BoundingCamClay(NDMaterialBase):

    def __init__(self, osi, mass_density, big_c, bulk_mod, ocr, mu_o, alpha, lamb, h, m):
        self.mass_density = float(mass_density)
        self.big_c = float(big_c)
        self.bulk_mod = float(bulk_mod)
        self.ocr = float(ocr)
        self.mu_o = float(mu_o)
        self.alpha = float(alpha)
        self.lamb = float(lamb)
        self.h = float(h)
        self.m = float(m)
        osi.n_mats += 1
        self._tag = osi.mats
        self._parameters = [self.op_type, self._tag, self.mass_density, self.big_c, self.bulk_mod, self.ocr, self.mu_o, self.alpha, self.lamb, self.h, self.m]
        self.to_process(osi)


class PlateFiber(NDMaterialBase):

    def __init__(self, osi, three_d):
        self.three_d = three_d.tag
        osi.n_mats += 1
        self._tag = osi.mats
        self._parameters = [self.op_type, self._tag, self.three_d.tag]
        self.to_process(osi)


class FSAM(NDMaterialBase):

    def __init__(self, osi, rho, s_x, s_y, conc, rou_x, rou_y, nu, alfadow):
        self.rho = float(rho)
        self.s_x = float(s_x)
        self.s_y = float(s_y)
        self.conc = float(conc)
        self.rou_x = float(rou_x)
        self.rou_y = float(rou_y)
        self.nu = float(nu)
        self.alfadow = float(alfadow)
        osi.n_mats += 1
        self._tag = osi.mats
        self._parameters = [self.op_type, self._tag, self.rho, self.s_x, self.s_y, self.conc, self.rou_x, self.rou_y, self.nu, self.alfadow]
        self.to_process(osi)


class ManzariDafalias(NDMaterialBase):

    def __init__(self, osi, g0, nu, e_init, mc, c, lambda_c, e0, ksi, p_atm, m, h0, ch, nb, a0, nd, z_max, cz, den):
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
        osi.n_mats += 1
        self._tag = osi.mats
        self._parameters = [self.op_type, self._tag, self.g0, self.nu, self.e_init, self.mc, self.c, self.lambda_c, self.e0, self.ksi, self.p_atm, self.m, self.h0, self.ch, self.nb, self.a0, self.nd, self.z_max, self.cz, self.den]
        self.to_process(osi)


class PM4Sand(NDMaterialBase):

    def __init__(self, osi, d_r, g_o, h_po, den, p_atm, h_o, e_max, e_min, n_b, n_d, a_do, z_max, c_z, c_e, phi_cv, nu, c_degr, c_dr, c_kaf, big_q, big_r, m, f_sed, p_sedo):
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
        self.a_do = float(a_do)
        self.z_max = float(z_max)
        self.c_z = float(c_z)
        self.c_e = float(c_e)
        self.phi_cv = float(phi_cv)
        self.nu = float(nu)
        self.c_degr = float(c_degr)
        self.c_dr = float(c_dr)
        self.c_kaf = float(c_kaf)
        self.big_q = float(big_q)
        self.big_r = float(big_r)
        self.m = float(m)
        self.f_sed = float(f_sed)
        self.p_sedo = float(p_sedo)
        osi.n_mats += 1
        self._tag = osi.mats
        self._parameters = [self.op_type, self._tag, self.d_r, self.g_o, self.h_po, self.den, self.p_atm, self.h_o, self.e_max, self.e_min, self.n_b, self.n_d, self.a_do, self.z_max, self.c_z, self.c_e, self.phi_cv, self.nu, self.c_degr, self.c_dr, self.c_kaf, self.big_q, self.big_r, self.m, self.f_sed, self.p_sedo]
        self.to_process(osi)


class StressDensityModel(NDMaterialBase):

    def __init__(self, osi, m_den, e_not, big_a, n, nu, a1, b1, a2, b2, a3, b3, fd, mu_not, mu_cyc, sc, big_m, patm, ssl1, ssl2, ssl3, ssl4, ssl5, ssl6, ssl7, ssl8, ssl9, ssl10, hsl, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10):
        self.m_den = float(m_den)
        self.e_not = float(e_not)
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
        self.patm = float(patm)
        self.ssl1 = float(ssl1)
        self.ssl2 = float(ssl2)
        self.ssl3 = float(ssl3)
        self.ssl4 = float(ssl4)
        self.ssl5 = float(ssl5)
        self.ssl6 = float(ssl6)
        self.ssl7 = float(ssl7)
        self.ssl8 = float(ssl8)
        self.ssl9 = float(ssl9)
        self.ssl10 = float(ssl10)
        self.hsl = float(hsl)
        self.p1 = float(p1)
        self.p2 = float(p2)
        self.p3 = float(p3)
        self.p4 = float(p4)
        self.p5 = float(p5)
        self.p6 = float(p6)
        self.p7 = float(p7)
        self.p8 = float(p8)
        self.p9 = float(p9)
        self.p10 = float(p10)
        osi.n_mats += 1
        self._tag = osi.mats
        self._parameters = [self.op_type, self._tag, self.m_den, self.e_not, self.big_a, self.n, self.nu, self.a1, self.b1, self.a2, self.b2, self.a3, self.b3, self.fd, self.mu_not, self.mu_cyc, self.sc, self.big_m, self.patm, self.ssl1, self.ssl2, self.ssl3, self.ssl4, self.ssl5, self.ssl6, self.ssl7, self.ssl8, self.ssl9, self.ssl10, self.hsl, self.p1, self.p2, self.p3, self.p4, self.p5, self.p6, self.p7, self.p8, self.p9, self.p10]
        self.to_process(osi)


class AcousticMedium(NDMaterialBase):

    def __init__(self, osi, big_k, rho):
        self.big_k = float(big_k)
        self.rho = float(rho)
        osi.n_mats += 1
        self._tag = osi.mats
        self._parameters = [self.op_type, self._tag, self.big_k, self.rho]
        self.to_process(osi)
