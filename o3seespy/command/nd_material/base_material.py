from o3seespy.base_model import OpenSeesObject


class NDMaterialBase(OpenSeesObject):
    op_base_type = "nDMaterial"
    op_type = None
    built = 1

    def set_parameter(self, osi, pstr, value, ele, eles, pval=1):
        from o3seespy import set_parameter
        if ele is not None:
            set_parameter(osi, value=value, eles=[ele], args=[pstr, pval])
        if eles is not None:
            set_parameter(osi, value=value, eles=eles, args=[pstr, pval])

    def build(self, osi):
        self.osi = osi
        osi.n_mat += 1
        self._tag = osi.n_mat
        ind = self.parameters.index(None)
        self.parameters[ind] = self._tag
        self.to_process(osi)
        self.built = 1
