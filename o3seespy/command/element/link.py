from o3seespy.command.element.base_element import ElementBase


class TwoNodeLink(ElementBase):
    """
    The TwoNodeLink Element Class
    
    This command is used to construct a twoNodeLink element object, which is defined by two nodes. The element can have
    zero or non-zero length. This element can have 1 to 6 degrees of freedom, where only the transverse and rotational
    degrees of freedom are coupled as long as the element has non-zero length. In addition, if the element length is
    larger than zero, the user can optionally specify how the P-Delta moments around the local x- and y-axis are
    distributed among a moment at node i, a moment at node j, and a shear couple. The sum of these three ratios
    is always equal to 1. In addition the shear center can be specified as a fraction of the element length
    from the iNode. The element does not contribute to the Rayleigh damping by default. If the element has
    non-zero length, the local x-axis is determined from the nodal geometry unless the optional x-axis
    vector is specified in which case the nodal geometry is ignored and the user-defined orientation
    is utilized. It is important to recognize that if this element has zero length, it does not
    consider the geometry as given by the nodal coordinates, but utilizes the user-defined
    orientation vectors to determine the directions of the springs.

    
    """
    op_type = 'twoNodeLink'

    def __init__(self, osi, ele_nodes, mats: list=None, dir=None, p_delta_vals: list=None, shear_dist=None, do_rayleigh=False, orient: list=None, mass: float=None):
        """
        Initial method for TwoNodeLink

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        ele_nodes: list
            A list of two element nodes
        mats: list, optional
            A list of objects associated with previously-defined uniaxial_material objects
        dir: None, optional
            
        p_delta_vals: list, optional
            P-delta moment contribution ratios, size of ratio vector is 2 for 2d-case and 4 for 3d-case (entries:
            ``[my_inode, my_jnode, mz_inode, mz_jnode]``) ``my_inode`` + ``my_jnode`` <= 1.0, ``mz_inode`` + ``mz_jnode`` <=
            1.0. remaining p-delta moments are resisted by shear couples. 
        shear_dist: None, optional
            
        do_rayleigh: bool
            To include rayleigh damping from the element (optional, default = no rayleigh damping contribution)
        orient: list, optional
            
        mass: float, optional
            Element mass (optional, default = 0.0)

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [0, 1]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(2)]
        >>> mats = [o3.uniaxial_material.Elastic(osi, 1.0),
        >>>         o3.uniaxial_material.Elastic(osi, 1.0)]
        >>> p_delta_vals = [1.0, 1.0]
        >>> o3.element.TwoNodeLink(osi, ele_nodes=ele_nodes, mats=mats, dir=[1, 1], p_delta_vals=p_delta_vals)
        """
        self.osi = osi
        self.ele_nodes = [x.tag for x in ele_nodes]
        if mats is None:
            self.mats = None
        else:
            self.mats = [x.tag for x in mats]
        self.dir = dir
        self.p_delta_vals = p_delta_vals
        self.shear_dist = shear_dist
        self.do_rayleigh = do_rayleigh
        self.orient = orient
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes]
        if getattr(self, 'mats') is not None:
            self._parameters += ['-mat', *self.mats]
        if getattr(self, 'dir') is not None:
            self._parameters += ['-dir', *self.dir]
        if getattr(self, 'p_delta_vals') is not None:
            self._parameters += ['-pDelta', *self.p_delta_vals]
        if getattr(self, 'shear_dist') is not None:
            self._parameters += ['-shearDist', *self.shear_dist]
        if getattr(self, 'do_rayleigh'):
            self._parameters += ['-doRayleigh']
        if getattr(self, 'orient') is not None:
            self._parameters += ['--orient', *self.orient]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        self.to_process(osi)
