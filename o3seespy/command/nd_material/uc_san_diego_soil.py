from o3seespy.command.nd_material.base_material import NDMaterialBase


class PressureIndependMultiYield(NDMaterialBase):

    def __init__(self, osi, nd, rho, ref_shear_modul, ref_bulk_modul, cohesi, peak_shear_stra, friction_ang=0., ref_press=100., press_depend_coe=0., no_yield_surf=20, yield_surf):
        self.nd = float(nd)
        self.rho = float(rho)
        self.ref_shear_modul = float(ref_shear_modul)
        self.ref_bulk_modul = float(ref_bulk_modul)
        self.cohesi = float(cohesi)
        self.peak_shear_stra = float(peak_shear_stra)
        self.friction_ang = float(friction_ang)
        self.ref_press = float(ref_press)
        self.press_depend_coe = float(press_depend_coe)
        self.no_yield_surf = float(no_yield_surf)
        self.yield_surf = yield_surf
        osi.n_mats += 1
        self._tag = osi.mats
        self._parameters = [self.op_type, self._tag, self.nd, self.rho, self.ref_shear_modul, self.ref_bulk_modul, self.cohesi, self.peak_shear_stra, self.friction_ang, self.ref_press, self.press_depend_coe, self.no_yield_surf, *self.yield_surf]
        self.to_process(osi)
