from o3seespy.command.element.base_element import ElementBase


class ZeroLength(ElementBase):
    """
        The ZeroLength Element Class

        This command is used to construct a zeroLength element object, which is defined by two nodes at the same location.
        The nodes are connected by multiple UniaxialMaterial objects to represent the force-deformation relationship for the
        element.
        """
    op_type = 'zeroLength'

    def __init__(self, osi, node_i, node_j, mat=None, mat_x=None, mat_y=None, mat_z=None,
                 mat_rot_x=None, mat_rot_y=None, mat_rot_z=None, r_flag=0):
        self.node_i = node_i
        self.node_j = node_j
        self.mat = mat
        self.mats = [mat_x, mat_y, mat_z, mat_rot_x, mat_rot_y, mat_rot_z]
        self.r_flag = r_flag
        osi.n_ele += 1
        self._tag = osi.n_ele

        self._parameters = [self.op_type, self._tag, *[self.node_i.tag, self.node_j.tag]]
        mats = []
        dirs = []
        for i, mat in enumerate(self.mats):
            if mat is not None:
                mats.append(mat.tag)
                dirs.append(i + 1)
        if len(mats) == 1:
            mats = mats[0]
            dirs = dirs[0]
        elif len(mats) == 0:
            if self.mat is None:
                raise ValueError("Must set either mat or a mat direction")
            mats = self.mat.tag
            dirs = 1
        self._parameters.append("-mat")
        self._parameters.append(mats)
        self._parameters.append("-dir")
        self._parameters.append(dirs)
        self._parameters.append("-doRayleigh")
        self._parameters.append(self.r_flag)
        self.to_process(osi)
