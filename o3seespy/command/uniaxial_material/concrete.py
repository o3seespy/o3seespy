from o3seespy.command.uniaxial_material.base_material import UniaxialMaterialBase



class Concrete01(UniaxialMaterialBase):

    def __init__(self, osi, fpc, epsc0, fpcu, eps_u):
        self.fpc = float(fpc)
        self.epsc0 = float(epsc0)
        self.fpcu = float(fpcu)
        self.eps_u = float(eps_u)
        osi.n_mats += 1
        self._tag = osi.mats
        self._parameters = [self.op_type, self._tag, self.fpc, self.epsc0, self.fpcu, self.eps_u]
        self.to_process(osi)


class Concrete02(UniaxialMaterialBase):

    def __init__(self, osi, fpc, epsc0, fpcu, eps_u, lamb, ft, ets):
        self.fpc = float(fpc)
        self.epsc0 = float(epsc0)
        self.fpcu = float(fpcu)
        self.eps_u = float(eps_u)
        self.lamb = float(lamb)
        self.ft = float(ft)
        self.ets = float(ets)
        osi.n_mats += 1
        self._tag = osi.mats
        self._parameters = [self.op_type, self._tag, self.fpc, self.epsc0, self.fpcu, self.eps_u, self.lamb, self.ft, self.ets]
        self.to_process(osi)


class Concrete04(UniaxialMaterialBase):

    def __init__(self, osi, fc, epsc, epscu, ec, fct, et, beta):
        self.fc = float(fc)
        self.epsc = float(epsc)
        self.epscu = float(epscu)
        self.ec = float(ec)
        self.fct = float(fct)
        self.et = float(et)
        self.beta = float(beta)
        osi.n_mats += 1
        self._tag = osi.mats
        self._parameters = [self.op_type, self._tag, self.fc, self.epsc, self.epscu, self.ec, self.fct, self.et, self.beta]
        self.to_process(osi)


class Concrete06(UniaxialMaterialBase):

    def __init__(self, osi, fc, e0, n, k, alpha1, fcr, ecr, b, alpha2):
        self.fc = float(fc)
        self.e0 = float(e0)
        self.n = float(n)
        self.k = float(k)
        self.alpha1 = float(alpha1)
        self.fcr = float(fcr)
        self.ecr = float(ecr)
        self.b = float(b)
        self.alpha2 = float(alpha2)
        osi.n_mats += 1
        self._tag = osi.mats
        self._parameters = [self.op_type, self._tag, self.fc, self.e0, self.n, self.k, self.alpha1, self.fcr, self.ecr, self.b, self.alpha2]
        self.to_process(osi)


class Concrete07(UniaxialMaterialBase):

    def __init__(self, osi, fc, epsc, ec, ft, et, xp, xn, r):
        self.fc = float(fc)
        self.epsc = float(epsc)
        self.ec = float(ec)
        self.ft = float(ft)
        self.et = float(et)
        self.xp = float(xp)
        self.xn = float(xn)
        self.r = float(r)
        osi.n_mats += 1
        self._tag = osi.mats
        self._parameters = [self.op_type, self._tag, self.fc, self.epsc, self.ec, self.ft, self.et, self.xp, self.xn, self.r]
        self.to_process(osi)


class Concrete01WithSITC(UniaxialMaterialBase):

    def __init__(self, osi, fpc, epsc0, fpcu, eps_u, end_strain_sitc=0.01):
        self.fpc = float(fpc)
        self.epsc0 = float(epsc0)
        self.fpcu = float(fpcu)
        self.eps_u = float(eps_u)
        self.end_strain_sitc = float(end_strain_sitc)
        osi.n_mats += 1
        self._tag = osi.mats
        self._parameters = [self.op_type, self._tag, self.fpc, self.epsc0, self.fpcu, self.eps_u, self.end_strain_sitc]
        self.to_process(osi)


class ConfinedConcrete01(UniaxialMaterialBase):

    def __init__(self, osi, sec_type, fpc, ec, epscu_type, epscu_val, nu, l1, l2, l3, phis, big_s, fyh, es0, ha_ratio, mu, phi_lon, internal_args=None, wrap_args=None, gravel=False, silica=False, tol=None, max_num_iter=None, epscu_limit=None, st_ratio=None):
        self.sec_type = sec_type
        self.fpc = float(fpc)
        self.ec = float(ec)
        self.epscu_type = epscu_type
        self.epscu_val = float(epscu_val)
        self.nu = nu
        self.l1 = float(l1)
        self.l2 = float(l2)
        self.l3 = float(l3)
        self.phis = float(phis)
        self.big_s = float(big_s)
        self.fyh = float(fyh)
        self.es0 = float(es0)
        self.ha_ratio = float(ha_ratio)
        self.mu = float(mu)
        self.phi_lon = float(phi_lon)
        self.internal_args = internal_args
        self.wrap_args = wrap_args
        self.gravel = gravel
        self.silica = silica
        self.tol = float(tol)
        self.max_num_iter = int(max_num_iter)
        self.epscu_limit = float(epscu_limit)
        self.st_ratio = st_ratio
        osi.n_mats += 1
        self._tag = osi.mats
        self._parameters = [self.op_type, self._tag, self.sec_type, self.fpc, self.ec, self.epscu_type, self.epscu_val, self.nu, self.l1, self.l2, self.l3, self.phis, self.big_s, self.fyh, self.es0, self.ha_ratio, self.mu, self.phi_lon]
        if getattr(self, 'internal_args') is not None:
            self._parameters += ['-internal', *self.internal_args]
        if getattr(self, 'wrap_args') is not None:
            self._parameters += ['-wrap', *self.wrap_args]
        if getattr(self, 'gravel') is not None:
            self._parameters += ['-gravel']
        if getattr(self, 'silica') is not None:
            self._parameters += ['-silica']
        if getattr(self, 'tol') is not None:
            self._parameters += ['-tol', self.tol]
        if getattr(self, 'max_num_iter') is not None:
            self._parameters += ['-maxNumIter', self.max_num_iter]
        if getattr(self, 'epscu_limit') is not None:
            self._parameters += ['-epscuLimit', self.epscu_limit]
        if getattr(self, 'st_ratio') is not None:
            self._parameters += ['-stRatio', self.st_ratio]
        self.to_process(osi)


class ConcreteD(UniaxialMaterialBase):

    def __init__(self, osi, fc, epsc, ft, epst, ec, alphac, alphat, cesp=0.25, etap=1.15):
        self.fc = float(fc)
        self.epsc = float(epsc)
        self.ft = float(ft)
        self.epst = float(epst)
        self.ec = float(ec)
        self.alphac = float(alphac)
        self.alphat = float(alphat)
        self.cesp = float(cesp)
        self.etap = float(etap)
        osi.n_mats += 1
        self._tag = osi.mats
        self._parameters = [self.op_type, self._tag, self.fc, self.epsc, self.ft, self.epst, self.ec, self.alphac, self.alphat, self.cesp, self.etap]
        self.to_process(osi)


class FRPConfinedConcrete(UniaxialMaterialBase):

    def __init__(self, osi, fpc1, fpc2, epsc0, big_d, c, ej, sj, tj, eju, big_s, fyl, fyh, dlong, dtrans, es, vo, k, use_buck):
        self.fpc1 = float(fpc1)
        self.fpc2 = float(fpc2)
        self.epsc0 = float(epsc0)
        self.big_d = float(big_d)
        self.c = float(c)
        self.ej = float(ej)
        self.sj = float(sj)
        self.tj = float(tj)
        self.eju = float(eju)
        self.big_s = float(big_s)
        self.fyl = float(fyl)
        self.fyh = float(fyh)
        self.dlong = float(dlong)
        self.dtrans = float(dtrans)
        self.es = float(es)
        self.vo = float(vo)
        self.k = float(k)
        self.use_buck = float(use_buck)
        osi.n_mats += 1
        self._tag = osi.mats
        self._parameters = [self.op_type, self._tag, self.fpc1, self.fpc2, self.epsc0, self.big_d, self.c, self.ej, self.sj, self.tj, self.eju, self.big_s, self.fyl, self.fyh, self.dlong, self.dtrans, self.es, self.vo, self.k, self.use_buck]
        self.to_process(osi)


class ConcreteCM(UniaxialMaterialBase):

    def __init__(self, osi, fpcc, epcc, ec, rc, xcrn, ft, et, rt, xcrp, gap_close=0):
        self.fpcc = float(fpcc)
        self.epcc = float(epcc)
        self.ec = float(ec)
        self.rc = float(rc)
        self.xcrn = float(xcrn)
        self.ft = float(ft)
        self.et = float(et)
        self.rt = float(rt)
        self.xcrp = float(xcrp)
        self.gap_close = float(gap_close)
        osi.n_mats += 1
        self._tag = osi.mats
        self._parameters = [self.op_type, self._tag, self.fpcc, self.epcc, self.ec, self.rc, self.xcrn, self.ft, self.et, self.rt, self.xcrp]
        if getattr(self, 'gap_close') is not None:
            self._parameters += ['-GapClose', self.gap_close]
        self.to_process(osi)
