from o3seespy.base_model import OpenSeesObject


class NDMaterialBase(OpenSeesObject):
    op_base_type = "nDMaterial"
    op_type = None

    def update_parameter(self, osi, pstr, value, ele, eles):
        from o3seespy import set_parameter
        if ele is not None:
            set_parameter(osi, value=value, eles=[ele], args=[pstr, 1])
        if eles is not None:
            set_parameter(osi, value=value, eles=eles, args=[pstr, 1])

