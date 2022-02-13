from o3seespy.base_model import OpenSeesObject


class NDMaterialBase(OpenSeesObject):
    op_base_type = "nDMaterial"
    op_type = None
    built = 1

    def set_parameter(self, osi, pstr, value, ele, eles, ele_tag_range=None, pval=1):
        from o3seespy import set_parameter
        if pval is None:
            args = [pstr]
        else:
            args = [pstr, pval]
        if ele is not None:
            set_parameter(osi, value=value, eles=[ele], args=args)
        elif eles is not None:
            set_parameter(osi, value=value, eles=eles, args=args)
        elif ele_tag_range is not None:
            set_parameter(osi, value=value, ele_tag_range=ele_tag_range, args=args)
        else:
            set_parameter(osi, value=value, args=args)

    def build(self, osi):
        self.osi = osi
        osi.n_mat += 1
        self._tag = osi.n_mat
        # ind = self.parameters.index(None)
        self.parameters[1] = self._tag
        self.to_process(osi)
        self.built = 1
