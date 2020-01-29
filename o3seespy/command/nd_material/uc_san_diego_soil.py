from o3seespy.command.nd_material.base_material import NDMaterialBase

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

        osi.n_mat += 1
        self._tag = osi.n_mat

        self._parameters = [self.op_type, self._tag, self.nd, self.rho, self.ref_shear_modul, self.ref_bulk_modul,
                            self.cohesi, self.peak_shear_stra, self.friction_ang, self.ref_press, self.press_depend_coe]

        if self.yield_surf is not None:
            self._parameters.append(self.no_yield_surf)  # from docs 'add a minus sign in front of noYieldSurf'
            self._parameters += list(self.yield_surf)

        else:
            # self._keyword_args['noYieldSurf'] = self.no_yield_surf
            self._parameters.append(self.no_yield_surf)
        self.to_process(osi)

