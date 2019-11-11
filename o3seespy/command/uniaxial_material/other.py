from o3seespy.command.uniaxial_material.base_material import UniaxialMaterialBase



class Hardening(UniaxialMaterialBase):
    op_type = 'Hardening'

    def __init__(self, osi, big_e, sigma_y, h_iso, h_kin, eta=0.0):
        self.big_e = float(big_e)
        self.sigma_y = float(sigma_y)
        self.h_iso = float(h_iso)
        self.h_kin = float(h_kin)
        self.eta = float(eta)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.big_e, self.sigma_y, self.h_iso, self.h_kin, self.eta]
        self.to_process(osi)


class Cast(UniaxialMaterialBase):
    op_type = 'Cast'

    def __init__(self, osi, n, bo, h, fy, big_e, big_l, b, ro, c_r1, c_r2, a1=None, a2=1.0, a3=None, a4=1.0):
        self.n = int(n)
        self.bo = float(bo)
        self.h = float(h)
        self.fy = float(fy)
        self.big_e = float(big_e)
        self.big_l = float(big_l)
        self.b = float(b)
        self.ro = float(ro)
        self.c_r1 = float(c_r1)
        self.c_r2 = float(c_r2)
        self.a1 = float(a1)
        self.a2 = float(a2)
        self.a3 = float(a3)
        self.a4 = float(a4)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.n, self.bo, self.h, self.fy, self.big_e, self.big_l, self.b, self.ro, self.c_r1, self.c_r2]
        special_pms = ['a1', 'a2', 'a3', 'a4']
        packets = [False, False, False, False]
        for i, pm in enumerate(special_pms):
            if getattr(self, pm) is not None:
                if packets[i]:
                    self._parameters += [*getattr(self, pm)]
                else:
                    self._parameters += [getattr(self, pm)]
        self.to_process(osi)


class ViscousDamper(UniaxialMaterialBase):
    op_type = 'ViscousDamper'

    def __init__(self, osi, big_k, cd, alpha, l_gap=0.0, nm=1, rel_tol=1e-6, abs_tol=1e-10, max_half=15):
        self.big_k = float(big_k)
        self.cd = float(cd)
        self.alpha = float(alpha)
        self.l_gap = float(l_gap)
        self.nm = int(nm)
        self.rel_tol = float(rel_tol)
        self.abs_tol = float(abs_tol)
        self.max_half = int(max_half)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.big_k, self.cd, self.alpha, self.l_gap, self.nm, self.rel_tol, self.abs_tol, self.max_half]
        self.to_process(osi)


class BilinearOilDamper(UniaxialMaterialBase):
    op_type = 'BilinearOilDamper'

    def __init__(self, osi, big_k, cd, fr=1.0, p=1.0, l_gap=0.0, nm=1, rel_tol=1e-6, abs_tol=1e-10, max_half=15):
        self.big_k = float(big_k)
        self.cd = float(cd)
        self.fr = float(fr)
        self.p = float(p)
        self.l_gap = float(l_gap)
        self.nm = int(nm)
        self.rel_tol = float(rel_tol)
        self.abs_tol = float(abs_tol)
        self.max_half = int(max_half)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.big_k, self.cd, self.fr, self.p, self.l_gap, self.nm, self.rel_tol, self.abs_tol, self.max_half]
        self.to_process(osi)


class Bilin(UniaxialMaterialBase):
    op_type = 'Bilin'

    def __init__(self, osi, k0, as__plus, as__neg, my__plus, my__neg, lamda_s, lamda_c, lamda_a, lamda_k, c_s, c_c, c_a, c_k, theta_p__plus, theta_p__neg, theta_pc__plus, theta_pc__neg, res__pos, res__neg, theta_u__plus, theta_u__neg, d__plus, d__neg, n_factor=0.0):
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
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.k0, self.as__plus, self.as__neg, self.my__plus, self.my__neg, self.lamda_s, self.lamda_c, self.lamda_a, self.lamda_k, self.c_s, self.c_c, self.c_a, self.c_k, self.theta_p__plus, self.theta_p__neg, self.theta_pc__plus, self.theta_pc__neg, self.res__pos, self.res__neg, self.theta_u__plus, self.theta_u__neg, self.d__plus, self.d__neg, self.n_factor]
        self.to_process(osi)


class ModIMKPeakOriented(UniaxialMaterialBase):
    op_type = 'ModIMKPeakOriented'

    def __init__(self, osi, k0, as__plus, as__neg, my__plus, my__neg, lamda_s, lamda_c, lamda_a, lamda_k, c_s, c_c, c_a, c_k, theta_p__plus, theta_p__neg, theta_pc__plus, theta_pc__neg, res__pos, res__neg, theta_u__plus, theta_u__neg, d__plus, d__neg):
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
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.k0, self.as__plus, self.as__neg, self.my__plus, self.my__neg, self.lamda_s, self.lamda_c, self.lamda_a, self.lamda_k, self.c_s, self.c_c, self.c_a, self.c_k, self.theta_p__plus, self.theta_p__neg, self.theta_pc__plus, self.theta_pc__neg, self.res__pos, self.res__neg, self.theta_u__plus, self.theta_u__neg, self.d__plus, self.d__neg]
        self.to_process(osi)


class ModIMKPinching(UniaxialMaterialBase):
    op_type = 'ModIMKPinching'

    def __init__(self, osi, k0, as__plus, as__neg, my__plus, my__neg, fpr_pos, fpr_neg, a_pinch, lamda_s, lamda_c, lamda_a, lamda_k, c_s, c_c, c_a, c_k, theta_p__plus, theta_p__neg, theta_pc__plus, theta_pc__neg, res__pos, res__neg, theta_u__plus, theta_u__neg, d__plus, d__neg):
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
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.k0, self.as__plus, self.as__neg, self.my__plus, self.my__neg, self.fpr_pos, self.fpr_neg, self.a_pinch, self.lamda_s, self.lamda_c, self.lamda_a, self.lamda_k, self.c_s, self.c_c, self.c_a, self.c_k, self.theta_p__plus, self.theta_p__neg, self.theta_pc__plus, self.theta_pc__neg, self.res__pos, self.res__neg, self.theta_u__plus, self.theta_u__neg, self.d__plus, self.d__neg]
        self.to_process(osi)


class SAWS(UniaxialMaterialBase):
    op_type = 'SAWS'

    def __init__(self, osi, f0, fi, du, s0, r1, r2, r3, r4, alpha, beta):
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
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.f0, self.fi, self.du, self.s0, self.r1, self.r2, self.r3, self.r4, self.alpha, self.beta]
        self.to_process(osi)


class BarSlip(UniaxialMaterialBase):
    op_type = 'BarSlip'

    def __init__(self, osi, fc, fy, es, fu, eh, db, ld, nb, depth, height, bs_flag, otype, anc_lratio=1.0, damage='Damage', unit='psi'):
        self.fc = float(fc)
        self.fy = float(fy)
        self.es = float(es)
        self.fu = float(fu)
        self.eh = float(eh)
        self.db = float(db)
        self.ld = float(ld)
        self.nb = float(nb)
        self.depth = float(depth)
        self.height = float(height)
        self.anc_lratio = float(anc_lratio)
        self.bs_flag = bs_flag
        self.otype = otype
        self.damage = damage
        self.unit = unit
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.fc, self.fy, self.es, self.fu, self.eh, self.db, self.ld, self.nb, self.depth, self.height, self.anc_lratio, self.bs_flag, self.otype, self.damage, self.unit]
        self.to_process(osi)


class BondSP01(UniaxialMaterialBase):
    op_type = 'Bond_SP01'

    def __init__(self, osi, fy, sy, fu, su, b, big_r):
        self.fy = float(fy)
        self.sy = float(sy)
        self.fu = float(fu)
        self.su = float(su)
        self.b = float(b)
        self.big_r = float(big_r)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.fy, self.sy, self.fu, self.su, self.b, self.big_r]
        self.to_process(osi)


class Fatigue(UniaxialMaterialBase):
    op_type = 'Fatigue'

    def __init__(self, osi, other, e0=None, m=None, min=None, max=None):
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
        self.to_process(osi)


class ImpactMaterial(UniaxialMaterialBase):
    op_type = 'ImpactMaterial'

    def __init__(self, osi, k1, k2, sigy, gap):
        self.k1 = float(k1)
        self.k2 = float(k2)
        self.sigy = float(sigy)
        self.gap = float(gap)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.k1, self.k2, self.sigy, self.gap]
        self.to_process(osi)


class HyperbolicGapMaterial(UniaxialMaterialBase):
    op_type = 'HyperbolicGapMaterial'

    def __init__(self, osi, kmax, kur, rf, fult, gap):
        self.kmax = float(kmax)
        self.kur = float(kur)
        self.rf = float(rf)
        self.fult = float(fult)
        self.gap = float(gap)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.kmax, self.kur, self.rf, self.fult, self.gap]
        self.to_process(osi)


class LimitState(UniaxialMaterialBase):
    op_type = 'LimitState'

    def __init__(self, osi, s1p, e1p, s2p, e2p, s3p, e3p, s1n, e1n, s2n, e2n, s3n, e3n, pinch_x, pinch_y, damage1, damage2, beta, curve, curve_type):
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
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.s1p, self.e1p, self.s2p, self.e2p, self.s3p, self.e3p, self.s1n, self.e1n, self.s2n, self.e2n, self.s3n, self.e3n, self.pinch_x, self.pinch_y, self.damage1, self.damage2, self.beta, self.curve.tag, self.curve_type]
        self.to_process(osi)


class MinMax(UniaxialMaterialBase):
    op_type = 'MinMax'

    def __init__(self, osi, other, min_strain=None, max_strain=None):
        self.other = other
        if min_strain is None:
            self.min_strain = None
        else:
            self.min_strain = float(min_strain)
        if max_strain is None:
            self.max_strain = None
        else:
            self.max_strain = float(max_strain)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.other.tag]
        if getattr(self, 'min_strain') is not None:
            self._parameters += ['-min', self.min_strain]
        if getattr(self, 'max_strain') is not None:
            self._parameters += ['-max', self.max_strain]
        self.to_process(osi)


class ElasticBilin(UniaxialMaterialBase):
    op_type = 'ElasticBilin'

    def __init__(self, osi, ep1, ep2, eps_p2, en1=None, en2=None, eps_n2=None):
        self.ep1 = float(ep1)
        self.ep2 = float(ep2)
        self.eps_p2 = float(eps_p2)
        self.en1 = float(en1)
        self.en2 = float(en2)
        self.eps_n2 = float(eps_n2)
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
        self.to_process(osi)


class ElasticMultiLinear(UniaxialMaterialBase):
    op_type = 'ElasticMultiLinear'

    def __init__(self, osi, strain=None, stress=None, eta=0.0):
        self.eta = float(eta)
        self.strain = strain
        self.stress = stress
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.eta]
        if getattr(self, 'strain') is not None:
            self._parameters += ['-strain', *self.strain]
        if getattr(self, 'stress') is not None:
            self._parameters += ['-stress', *self.stress]
        self.to_process(osi)


class MultiLinear(UniaxialMaterialBase):
    op_type = 'MultiLinear'

    def __init__(self, osi, pts):
        self.pts = pts
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, *self.pts]
        self.to_process(osi)


class InitStrainMaterial(UniaxialMaterialBase):
    op_type = 'InitStrainMaterial'

    def __init__(self, osi, other, init_strain):
        self.other = other
        self.init_strain = float(init_strain)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.other.tag, self.init_strain]
        self.to_process(osi)


class InitStressMaterial(UniaxialMaterialBase):
    op_type = 'InitStressMaterial'

    def __init__(self, osi, other, init_stress):
        self.other = other
        self.init_stress = float(init_stress)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.other.tag, self.init_stress]
        self.to_process(osi)


class PathIndependent(UniaxialMaterialBase):
    op_type = 'PathIndependent'

    def __init__(self, osi, tag):
        self.tag = int(tag)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.tag]
        self.to_process(osi)


class ECC01(UniaxialMaterialBase):
    op_type = 'ECC01'

    def __init__(self, osi, sigt0, epst0, sigt1, epst1, epst2, sigc0, epsc0, epsc1, alpha_t1, alpha_t2, alpha_c, alpha_cu, beta_t, beta_c):
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
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.sigt0, self.epst0, self.sigt1, self.epst1, self.epst2, self.sigc0, self.epsc0, self.epsc1, self.alpha_t1, self.alpha_t2, self.alpha_c, self.alpha_cu, self.beta_t, self.beta_c]
        self.to_process(osi)


class SelfCentering(UniaxialMaterialBase):
    op_type = 'SelfCentering'

    def __init__(self, osi, k1, k2, sig_act, beta, eps_slip=0, eps_bear=0, r_bear=None):
        self.k1 = float(k1)
        self.k2 = float(k2)
        self.sig_act = float(sig_act)
        self.beta = float(beta)
        self.eps_slip = float(eps_slip)
        self.eps_bear = float(eps_bear)
        self.r_bear = float(r_bear)
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
        self.to_process(osi)


class Viscous(UniaxialMaterialBase):
    op_type = 'Viscous'

    def __init__(self, osi, big_c, alpha):
        self.big_c = float(big_c)
        self.alpha = float(alpha)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.big_c, self.alpha]
        self.to_process(osi)


class BoucWen(UniaxialMaterialBase):
    op_type = 'BoucWen'

    def __init__(self, osi, alpha, ko, n, gamma, beta, ao, delta_a, delta_nu, delta_eta):
        self.alpha = float(alpha)
        self.ko = float(ko)
        self.n = float(n)
        self.gamma = float(gamma)
        self.beta = float(beta)
        self.ao = float(ao)
        self.delta_a = float(delta_a)
        self.delta_nu = float(delta_nu)
        self.delta_eta = float(delta_eta)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.alpha, self.ko, self.n, self.gamma, self.beta, self.ao, self.delta_a, self.delta_nu, self.delta_eta]
        self.to_process(osi)


class BWBN(UniaxialMaterialBase):
    op_type = 'BWBN'

    def __init__(self, osi, alpha, ko, n, gamma, beta, ao, q, zetas, p, shi, delta_shi, lamb, tol, max_iter):
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
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.alpha, self.ko, self.n, self.gamma, self.beta, self.ao, self.q, self.zetas, self.p, self.shi, self.delta_shi, self.lamb, self.tol, self.max_iter]
        self.to_process(osi)


class AxialSp(UniaxialMaterialBase):
    op_type = 'AxialSp'

    def __init__(self, osi, sce, fty, fcy, bte, bty, bcy, fcr):
        self.sce = float(sce)
        self.fty = float(fty)
        self.fcy = float(fcy)
        self.bte = float(bte)
        self.bty = float(bty)
        self.bcy = float(bcy)
        self.fcr = float(fcr)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.sce, self.fty, self.fcy, self.bte, self.bty, self.bcy, self.fcr]
        self.to_process(osi)


class AxialSpHD(UniaxialMaterialBase):
    op_type = 'AxialSpHD'

    def __init__(self, osi, sce, fty, fcy, bte, bty, bth, bcy, fcr, ath):
        self.sce = float(sce)
        self.fty = float(fty)
        self.fcy = float(fcy)
        self.bte = float(bte)
        self.bty = float(bty)
        self.bth = float(bth)
        self.bcy = float(bcy)
        self.fcr = float(fcr)
        self.ath = float(ath)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.sce, self.fty, self.fcy, self.bte, self.bty, self.bth, self.bcy, self.fcr, self.ath]
        self.to_process(osi)


class CFSWSWP(UniaxialMaterialBase):
    op_type = 'CFSWSWP'

    def __init__(self, osi, height, width, fut, tf, ife, ifi, ts, np, ds, vs, sc, nc, otype, opening_area, opening_length):
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
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.height, self.width, self.fut, self.tf, self.ife, self.ifi, self.ts, self.np, self.ds, self.vs, self.sc, self.nc, self.otype, self.opening_area, self.opening_length]
        self.to_process(osi)


class CFSSSWP(UniaxialMaterialBase):
    op_type = 'CFSSSWP'

    def __init__(self, osi, height, width, fuf, fyf, tf, af, fus, fys, ts, np, ds, vs, sc, dt, opening_area, opening_length):
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
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.height, self.width, self.fuf, self.fyf, self.tf, self.af, self.fus, self.fys, self.ts, self.np, self.ds, self.vs, self.sc, self.dt, self.opening_area, self.opening_length]
        self.to_process(osi)
