from o3seespy.base_model import OpenSeesObject


class Node(OpenSeesObject):
    op_base_type = "node"
    op_type = "node"

    def __init__(self, osi, x: float, y=None, z=None, vel=None, acc=None, mass: list=None,
                 x_mass=None, y_mass=None, z_mass=None, x_rot_mass=None, y_rot_mass=None, z_rot_mass=None, tag=None):
        """
        An OpenSees node

        Parameters
        ----------
        osi : o3seespy.opensees_instance.OpenSeesInstance object
            An instance of OpenSees
        x : float
            x-coordinate
        y : float, optional
            y-coordinate
        z : float, optional
            z-coordinate
        vel : iterable object, optional
            nodal velocities (x, y, z)
        acc : iterable object, optional
        mass: iterable object, option
            nodal masses
        """

        self.x = float(x)
        if y is not None:
            self.y = float(y)
        if z is not None:
            self.z = float(z)
        self.vel = vel
        self.acc = acc
        if tag is None:
            osi.n_node += 1
            self._tag = osi.n_node
        else:
            self._tag = int(tag)
        if osi.ndm == 1:
            self._parameters = [self._tag, self.x]
            poss_mass = [x_mass]
        elif osi.ndm == 2:
            self._parameters = [self._tag, self.x, self.y]
            poss_mass = [x_mass, y_mass, z_rot_mass]
        elif osi.ndm == 3:
            self._parameters = [self._tag, self.x, self.y, self.z]
            poss_mass = [x_mass, y_mass, z_mass, x_rot_mass, y_rot_mass, z_rot_mass]
        else:
            raise NotImplementedError(f"Currently only supports 1-3D analyses, ndm={osi.ndm}")
        if mass is None:
            mass = []
            if poss_mass[0] is not None:
                none_found = False
                for mval in poss_mass:
                    if mval is None:
                        none_found = True
                    else:
                        if not none_found:
                            mass.append(float(mval))
                        else:
                            mstr = ','.join([str(x) for x in poss_mass])
                            raise ValueError(f'Cannot set mass, since None in mass=[{mstr}]')
        else:
            mass = [float(x) for x in mass]
        if len(mass):
            self.mass = mass
            self._parameters += ["-mass", *mass]
        if self.vel is not None:
            self._parameters += ["-vel", self.vel]
        if self.acc is not None:
            self._parameters += ["-accel", self.acc]
        self.to_process(osi)


def build_regular_node_mesh(osi, xs, ys, zs=None, active=None, tags=None):
    """
    Creates an array of nodes that are in a regular mesh

    The mesh has len(xs) nodes in the x-direction and len(ys) in the y-direction.
    If zs is not None then has len(zs) in the z-direction.

    Parameters
    ----------
    osi
    xs
    ys
    zs
    active
    tags: array_like
        array of node tags

    Returns
    -------
    np.array
        axis-0 = x-direction
        axis-1 = y-direction
        axis-2 = z  # not included if len(zs)=1 or zs=None
    """
    # axis-0 = x  # unless x or y are singular
    # axis-1 = y
    # axis-2 = z  # not included if len(zs)=1 or
    tag = None
    from numpy import array
    if not hasattr(zs, '__len__'):
        zs = [zs]
    sn = []
    for xx in range(len(xs)):
        sn.append([])
        for yy in range(len(ys)):

            if len(zs) == 1:
                if tags is not None:
                    tag = tags[xx][yy]
                if active is None or active[xx][yy]:
                    if osi.ndm == 2:
                        pms = [osi, xs[xx], ys[yy]]
                    else:
                        pms = [osi, xs[xx], ys[yy], zs[0]]
                    sn[xx].append(Node(*pms, tag=tag))
                else:
                    sn[xx].append(None)
            else:
                sn[xx].append([])
                for zz in range(len(zs)):
                    if tags is not None:
                        tag = tags[xx][yy][zz]
                    # Establish left and right nodes
                    if active is None or active[xx][yy][zz]:
                        sn[xx][yy].append(Node(osi, xs[xx], ys[yy], zs[zz], tag=tag))
                    else:
                        sn[xx][yy].append(None)
    # if len(zs) == 1:
    #     return sn[0]
    return array(sn)


def build_varied_y_node_mesh(osi, xs, ys, zs=None, active=None):
    """
    Creates an array of nodes that in vertical lines, but vary in height

    The mesh has len(xs)=ln(ys) nodes in the x-direction and len(ys[0]) in the y-direction.
    If zs is not None then has len(zs) in the z-direction.

    Parameters
    ----------
    osi
    xs
    ys
    zs
    active

    Returns
    -------
    np.array
        axis-0 = x-direction
        axis-1 = y-direction
        axis-2 = z  # not included if len(zs)=1 or zs=None
    """
    # axis-0 = x  # unless x or y are singular
    # axis-1 = y
    # axis-2 = z  # not included if len(zs)=1 or
    from numpy import array
    if not hasattr(zs, '__len__'):
        zs = [zs]
    sn = []
    for xx in range(len(xs)):
        sn.append([])
        for yy in range(len(ys[xx])):

            if len(zs) == 1:
                if active is None or active[xx][yy]:
                    if osi.ndm == 2:
                        pms = [osi, xs[xx], ys[xx][yy]]
                    else:
                        pms = [osi, xs[xx], ys[xx][yy], zs[0]]
                    sn[xx].append(Node(*pms))
                else:
                    sn[xx].append(None)
            else:
                sn[xx].append([])
                for zz in range(len(zs)):
                    # Establish left and right nodes
                    if active is None or active[xx][yy][zz]:
                        sn[xx][yy].append(Node(osi, xs[xx], ys[xx][yy], zs[zz]))
                    else:
                        sn[xx][yy].append(None)
    # if len(zs) == 1:
    #     return sn[0]
    return array(sn)


# UNUSED?
def duplicate_node(osi, node):
    """
    Copy a node to initialise in another processor in parallel mode

    Note: It has the same node number
    """
    if osi.ndm == 1:
        _parameters = [node.tag, *[node.x]]
        pms = ['mass']
    elif osi.ndm == 2:
        _parameters = [node.tag, *[node.x, node.y]]
        pms = ['mass']
    elif osi.ndm == 3:
        _parameters = [node.tag, *[node.x, node.y, node.z]]
        pms = ['mass']
    else:
        raise NotImplementedError("Currently only supports 1-3D analyses")
    masses = []
    none_found = 0
    if hasattr(node, 'mass') and node.mass is not None:
        _parameters += ["-mass", *node.mass]
    if node.vel is not None:
        _parameters += ["-vel", node.vel]
    if node.acc is not None:
        _parameters += ["-accel", node.acc]
    osi.to_process('node', _parameters)
    

def repeat_node(osi, node):
    """
    Copy a node to initialise in another processor in parallel mode

    Note: It has the same node number
    """
    
    if osi.ndm == 1:
        return Node(osi, node.x)
    elif osi.ndm == 2:
        return Node(osi, node.x, node.y)
    elif osi.ndm == 3:
        return Node(osi, node.x, node.y, node.z)
    else:
        raise NotImplementedError("Currently only supports 1-3D analyses")


# def build_node_if_within_segment(osi, coords, segment):
#     """
#
#     Parameters
#     ----------
#     coords
#     segment: array_like
#         if osi.ndm = 1, segment is [[min_x], [max_x]]
#         if osi.num = 2, segment is [[
#
#     Returns
#     -------
#
#     """
#     pass


def hash_coords(coords, nsf=8):
    return '_'.join(['{n:.{nsf}}'.format(n=float(x), nsf=nsf) for x in coords])


def build_node_tag_hash_dict_from_mesh(ndm, xs, ys, zs=None, active=None, init_tag=0, nsf=8):
    """
    Creates an array of nodes that in vertical lines, but vary in height

    The mesh has len(xs) nodes in the x-direction and len(ys[0]) in the y-direction.
    If zs is not None then has len(zs) in the z-direction.

    Parameters
    ----------
    osi
    xs
    ys
    zs
    active

    Returns
    -------
    np.array
        axis-0 = x-direction
        axis-1 = y-direction
        axis-2 = z  # not included if len(zs)=1 or zs=None
        """
    # axis-0 = x  # unless x or y are singular
    # axis-1 = y
    # axis-2 = z  # not included if len(zs)=1 or
    import numpy as np
    ys = np.array(ys)
    nx = xs.shape[0]
    ny = ys.shape[min(len(ys.shape), 1)]
    if not hasattr(zs, '__len__'):
        zs = [zs]
    if len(active.shape) == 3:
        azi = np.arange(len(zs))
    else:
        azi = [None] * len(xs)
    if len(xs.shape) >= 2:
        yind4x = np.arange(ny)
    else:
        yind4x = [None] * ny
    if len(ys.shape) >= 2:
        xind4y = np.arange(nx)
    else:
        xind4y = [None] * nx

    nz = zs.shape[-1]
    sd = {}
    tag = init_tag
    for xx in range(nx):
        for yy in range(ny):
            for zz in range(nz):
                if active is None or active[xx, yy, azi[zz]]:
                    if ndm == 2:
                        pms = [xs[xx, yind4x[yy]], ys[xind4y[xx], yy]]
                    else:
                        pms = [xs[xx, yind4x[yy]], ys[xind4y[xx], yy], zs[zz]]

                    fstr = hash_coords(pms, nsf)
                    if fstr not in sd:
                        sd[fstr] = []
                    sd[fstr].append((tag, *pms))
                    tag += 1

    return sd, tag
