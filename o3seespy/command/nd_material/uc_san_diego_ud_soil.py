from o3seespy.command.nd_material.base_material import NDMaterialBase


class FluidSolidPorousMaterial(NDMaterialBase):

    def __init__(self, osi, nd, soil_mat, combined_bulk_modul, pa=101.0):
        self.nd = float(nd)
        self.soil_mat = soil_mat.tag
        self.combined_bulk_modul = float(combined_bulk_modul)
        self.pa = float(pa)
        osi.n_mats += 1
        self._tag = osi.mats
        self._parameters = [self.op_type, self._tag, self.nd, self.soil_mat.tag, self.combined_bulk_modul, self.pa]
        self.to_process(osi)
