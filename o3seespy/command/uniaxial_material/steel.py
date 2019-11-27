from o3seespy.command.uniaxial_material.base_material import UniaxialMaterialBase



class Steel01(UniaxialMaterialBase):
    op_type = 'Steel01'

    def __init__(self, osi, fy, e0, b, a1, a2, a3, a4):
        self.fy = float(fy)
        self.e0 = float(e0)
        self.b = float(b)
        self.a1 = float(a1)
        self.a2 = float(a2)
        self.a3 = float(a3)
        self.a4 = float(a4)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.fy, self.e0, self.b, self.a1, self.a2, self.a3, self.a4]
        self.to_process(osi)


class Steel02(UniaxialMaterialBase):
    op_type = 'Steel02'

    def __init__(self, osi, fy, e0, b, params, a1=None, a2=1.0, a3=None, a4=1.0, sig_init=0.0):
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
        self.to_process(osi)


class Hysteretic(UniaxialMaterialBase):
    op_type = 'Hysteretic'

    def __init__(self, osi, p1, p2, p3, n1, n2, n3, pinch_x, pinch_y, damage1, damage2, beta):
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
        self.to_process(osi)


class ReinforcingSteelGABuck(UniaxialMaterialBase):
    op_type = 'ReinforcingSteel'

    def __init__(self, osi, fy, fu, es, esh, eps_sh, eps_ult, lsr, beta, r, gamma):
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
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.fy, self.fu, self.es, self.esh, self.eps_sh, self.eps_ult, '-GABuck', self.lsr, self.beta, self.r, self.gamma]
        self.to_process(osi)

class ReinforcingSteelDMBuck(UniaxialMaterialBase):
    op_type = 'ReinforcingSteel'

    def __init__(self, osi, fy, fu, es, esh, eps_sh, eps_ult, lsr_2, alpha=1.0):
        self.fy = float(fy)
        self.fu = float(fu)
        self.es = float(es)
        self.esh = float(esh)
        self.eps_sh = float(eps_sh)
        self.eps_ult = float(eps_ult)
        self.lsr_2 = lsr_2
        self.alpha = float(alpha)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.fy, self.fu, self.es, self.esh, self.eps_sh, self.eps_ult, '-DMBuck', self.lsr_2, self.alpha]
        self.to_process(osi)

class ReinforcingSteelCMFatigue(UniaxialMaterialBase):
    op_type = 'ReinforcingSteel'

    def __init__(self, osi, fy, fu, es, esh, eps_sh, eps_ult, cf, alpha_2, cd):
        self.fy = float(fy)
        self.fu = float(fu)
        self.es = float(es)
        self.esh = float(esh)
        self.eps_sh = float(eps_sh)
        self.eps_ult = float(eps_ult)
        self.cf = float(cf)
        self.alpha_2 = alpha_2
        self.cd = float(cd)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.fy, self.fu, self.es, self.esh, self.eps_sh, self.eps_ult, '-CMFatigue', self.cf, self.alpha_2, self.cd]
        self.to_process(osi)

class ReinforcingSteelIsoHard(UniaxialMaterialBase):
    op_type = 'ReinforcingSteel'

    def __init__(self, osi, fy, fu, es, esh, eps_sh, eps_ult, a1=4.3, limit=1.0):
        self.fy = float(fy)
        self.fu = float(fu)
        self.es = float(es)
        self.esh = float(esh)
        self.eps_sh = float(eps_sh)
        self.eps_ult = float(eps_ult)
        self.a1 = float(a1)
        self.limit = float(limit)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.fy, self.fu, self.es, self.esh, self.eps_sh, self.eps_ult, '-IsoHard', self.a1, self.limit]
        self.to_process(osi)

class ReinforcingSteelMPCurveParams(UniaxialMaterialBase):
    op_type = 'ReinforcingSteel'

    def __init__(self, osi, fy, fu, es, esh, eps_sh, eps_ult, r1=0.333, r2=18.0, r3=4.0):
        self.fy = float(fy)
        self.fu = float(fu)
        self.es = float(es)
        self.esh = float(esh)
        self.eps_sh = float(eps_sh)
        self.eps_ult = float(eps_ult)
        self.r1 = float(r1)
        self.r2 = float(r2)
        self.r3 = float(r3)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.fy, self.fu, self.es, self.esh, self.eps_sh, self.eps_ult, '-MPCurveParams', self.r1, self.r2, self.r3]
        self.to_process(osi)


class DoddRestrepo(UniaxialMaterialBase):
    op_type = 'Dodd_Restrepo'

    def __init__(self, osi, fy, fsu, esh, esu, youngs, eshi, fshi, omega_fac=1.0):
        self.fy = float(fy)
        self.fsu = float(fsu)
        self.esh = float(esh)
        self.esu = float(esu)
        self.youngs = float(youngs)
        self.eshi = float(eshi)
        self.fshi = float(fshi)
        self.omega_fac = float(omega_fac)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.fy, self.fsu, self.esh, self.esu, self.youngs, self.eshi, self.fshi, self.omega_fac]
        self.to_process(osi)


class RambergOsgoodSteel(UniaxialMaterialBase):
    op_type = 'RambergOsgoodSteel'

    def __init__(self, osi, fy, e0, a, n):
        self.fy = float(fy)
        self.e0 = float(e0)
        self.a = float(a)
        self.n = float(n)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.fy, self.e0, self.a, self.n]
        self.to_process(osi)


class SteelMPF(UniaxialMaterialBase):
    op_type = 'SteelMPF'

    def __init__(self, osi, fyp, fyn, e0, bp, bn, r0, c_r1, c_r2, a1=0.0, a2=1.0, a3=0.0, a4=1.0):
        self.fyp = float(fyp)
        self.fyn = float(fyn)
        self.e0 = float(e0)
        self.bp = float(bp)
        self.bn = float(bn)
        self.r0 = float(r0)
        self.c_r1 = float(c_r1)
        self.c_r2 = float(c_r2)
        self.a1 = float(a1)
        self.a2 = float(a2)
        self.a3 = float(a3)
        self.a4 = float(a4)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.fyp, self.fyn, self.e0, self.bp, self.bn, self.r0, self.c_r1, self.c_r2, self.a1, self.a2, self.a3, self.a4]
        self.to_process(osi)


class Steel01Thermal(UniaxialMaterialBase):
    op_type = 'Steel01Thermal'

    def __init__(self, osi, fy, e0, b, a1, a2, a3, a4):
        self.fy = float(fy)
        self.e0 = float(e0)
        self.b = float(b)
        self.a1 = float(a1)
        self.a2 = float(a2)
        self.a3 = float(a3)
        self.a4 = float(a4)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.fy, self.e0, self.b, self.a1, self.a2, self.a3, self.a4]
        self.to_process(osi)
