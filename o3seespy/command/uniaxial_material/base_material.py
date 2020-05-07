from o3seespy.base_model import OpenSeesObject


class UniaxialMaterialBase(OpenSeesObject):
    op_base_type = "uniaxialMaterial"
    op_type = None

    def set_parameter(self, osi, pstr, value, ele, eles):
        from o3seespy import set_parameter
        if ele is not None:
            set_parameter(osi, value=value, eles=[ele], args=[pstr, 1])
        if eles is not None:
            set_parameter(osi, value=value, eles=eles, args=[pstr, 1])

    def build(self, osi):
        self.osi = osi
        osi.n_mat += 1
        self._tag = osi.n_mat
        self.to_process(osi)
        self.built = 1
