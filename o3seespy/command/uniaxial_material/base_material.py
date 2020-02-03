from o3seespy.base_model import OpenSeesObject


class UniaxialMaterialBase(OpenSeesObject):
    op_base_type = "uniaxialMaterial"
    op_type = None

