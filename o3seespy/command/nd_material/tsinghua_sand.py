from o3seespy.command.nd_material.base_material import NDMaterialBase


class CycLiqCP(NDMaterialBase):
    op_type = 'CycLiqCP'

    def __init__(self, osi, g0, kappa, h, mfc, dre1, mdc, dre2, rdr, alpha, dir, ein, rho):
        self.g0 = float(g0)
        self.kappa = float(kappa)
        self.h = float(h)
        self.mfc = float(mfc)
        self.dre1 = float(dre1)
        self.mdc = float(mdc)
        self.dre2 = float(dre2)
        self.rdr = float(rdr)
        self.alpha = float(alpha)
        self.dir = float(dir)
        self.ein = float(ein)
        self.rho = float(rho)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.g0, self.kappa, self.h, self.mfc, self.dre1, self.mdc, self.dre2, self.rdr, self.alpha, self.dir, self.ein, self.rho]
        self.to_process(osi)


class CycLiqCPSP(NDMaterialBase):
    op_type = 'CycLiqCPSP'

    def __init__(self, osi, g0, kappa, h, big_m, dre1, dre2, rdr, alpha, dir, lambdac, ksi, e0, np, nd, ein, rho):
        self.g0 = float(g0)
        self.kappa = float(kappa)
        self.h = float(h)
        self.big_m = float(big_m)
        self.dre1 = float(dre1)
        self.dre2 = float(dre2)
        self.rdr = float(rdr)
        self.alpha = float(alpha)
        self.dir = float(dir)
        self.lambdac = float(lambdac)
        self.ksi = float(ksi)
        self.e0 = float(e0)
        self.np = float(np)
        self.nd = float(nd)
        self.ein = float(ein)
        self.rho = float(rho)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.g0, self.kappa, self.h, self.big_m, self.dre1, self.dre2, self.rdr, self.alpha, self.dir, self.lambdac, self.ksi, self.e0, self.np, self.nd, self.ein, self.rho]
        self.to_process(osi)
