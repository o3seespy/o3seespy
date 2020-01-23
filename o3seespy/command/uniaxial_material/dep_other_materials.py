# from o3seespy.command.uniaxial_material.base_material import UniaxialMaterialBase
#
#
# class Viscous(UniaxialMaterialBase):
#     op_type = "Viscous"
#
#     def __init__(self, osi, c, alpha):
#         self.c = c
#         self.alpha = float(alpha)
#         osi.n_mat += 1
#         self._tag = osi.n_mat
#         self._parameters = [self.op_type, self._tag, self.c, self.alpha]
#         self.to_process(osi)
#
#
#
# class MultiLinear(UniaxialMaterialBase):
#     op_type = "MultiLinear"
#
#     def __init__(self, osi, points):
#         self.points = points
#         osi.n_mat += 1
#         self._tag = osi.n_mat
#         self._parameters = [self.op_type, self._tag, self.points]
#         self.to_process(osi)
#
#
# class BoucWen(UniaxialMaterialBase):
#     op_type = "BoucWen"
#
#     def __init__(self, osi, alpha, ko, n, gamma, beta, ao, delta_a, delta_nu, delta_eta):
#         self.alpha = float(alpha)
#         self.ko = float(ko)
#         self.n = float(n)
#         self.gamma = float(gamma)
#         self.beta = float(beta)
#         self.ao = float(ao)
#         self.delta_a = float(delta_a)
#         self.delta_nu = float(delta_nu)
#         self.delta_eta = float(delta_eta)
#         osi.n_mat += 1
#         self._tag = osi.n_mat
#         self._parameters = [self.op_type, self._tag, self.alpha, self.ko, self.n, self.gamma, self.beta, self.ao, self.delta_a, self.delta_nu, self.delta_eta]
#         self.to_process(osi)
#
#
# class BondSP01(UniaxialMaterialBase):
#     op_type = "Bond_SP01"
#
#     def __init__(self, osi, fy, sy, fu, su, b, big_r):
#         self.fy = float(fy)
#         self.sy = float(sy)
#         self.fu = float(fu)
#         self.su = float(su)
#         self.b = float(b)
#         self.big_r = float(big_r)
#         osi.n_mat += 1
#         self._tag = osi.n_mat
#         self._parameters = [self.op_type, self._tag, self.fy, self.sy, self.fu, self.su, self.b, self.big_r]
#         self.to_process(osi)
#
#
# class BilinearOilDamper(UniaxialMaterialBase):
#     op_type = "BilinearOilDamper"
#
#     def __init__(self, osi, k, cd, fr=1.0, p=1.0, l_gap=0.0, nm=1, rel_tol=1e-6, abs_tol=1e-10, max_half=15):
#         self.k = float(k)
#         self.cd = float(cd)
#         self.fr = float(fr)
#         self.p = float(p)
#         self.l_gap = float(l_gap)
#         self.nm = int(nm)
#         self.rel_tol = float(rel_tol)
#         self.abs_tol = float(abs_tol)
#         self.max_half = int(max_half)
#         osi.n_mat += 1
#         self._tag = osi.n_mat
#         self._parameters = [self.op_type, self._tag, self.k, self.cd, self.fr, self.p, self.l_gap, self.nm, self.rel_tol, self.abs_tol, self.max_half]
#         self.to_process(osi)
#
#
# class InitStrainMaterial(UniaxialMaterialBase):
#     op_type = "InitStrainMaterial"
#
#     def __init__(self, osi, other, init_strain):
#         self.other = other
#         self.init_strain = float(init_strain)
#         osi.n_mat += 1
#         self._tag = osi.n_mat
#         self._parameters = [self.op_type, self._tag, self.other.tag, self.init_strain]
#         self.to_process(osi)
