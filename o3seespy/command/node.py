from o3seespy.base_model import OpenSeesObject


class Node(OpenSeesObject):
    op_base_type = "node"
    op_type = "node"
    # x_con = None
    # y_con = None
    # z_con = None
    # x_rot_con = None
    # y_rot_con = None
    # z_rot_con = None
    # x_mass = 0.0
    # y_mass = 0.0
    # z_rot_mass = 0.0

    def __init__(self, osi, x: float, y=None, z=None, vel=None, acc=None,
                 x_mass=None, y_mass=None, z_mass=None, x_rot_mass=None, y_rot_mass=None, z_rot_mass=None):
        """
        An OpenSEES node

        Parameters
        ----------
        osi : opensees_pack.opensees_instance.OpenSeesInstance object
            An instance of opensees
        x : float
            x-coordinate
        y : float, optional
            y-coordinate
        z : float, optional
            z-coordinate
        vel : iterable object, optional
            nodal velocities (x, y, z)
        acc : iterable object, optional
        """

        self.x = float(x)
        if y is not None:
            self.y = float(y)
        if z is not None:
            self.z = float(z)
        self.x_mass = x_mass
        self.y_mass = y_mass
        self.z_mass = z_mass
        self.x_rot_mass = x_rot_mass
        self.y_rot_mass = y_rot_mass
        self.z_rot_mass = z_rot_mass
        self.vel = vel
        self.acc = acc
        osi.n_node += 1
        self._tag = osi.n_node
        if osi.ndm == 1:
            self._parameters = [self._tag, *[self.x]]
            pms = ['x_mass']
        elif osi.ndm == 2:
            self._parameters = [self._tag, *[self.x, self.y]]
            pms = ['x_mass', 'y_mass', 'z_rot_mass']
        elif osi.ndm == 3:
            self._parameters = [self._tag, *[self.x, self.y, self.z]]
            pms = ['x_mass', 'y_mass', 'z_mass', 'x_rot_mass', 'y_rot_mass', 'z_rot_mass']
        else:
            raise NotImplementedError("Currently only supports 1-3D analyses")
        masses = []
        none_found = 0
        for pm in pms:
            val = getattr(self, pm)
            if val is None:
                none_found = pm
            else:
                setattr(self, pm, float(val))
                if not none_found:
                    masses.append(float(val))
                else:
                    raise ValueError(f'Cannot set {pm} since {none_found} is None')
        if len(masses):
            self._parameters += ["-mass", *[self.x_mass]]
        if self.vel is not None:
            self._parameters += ["-vel", self.vel]
        if self.acc is not None:
            self._parameters += ["-accel", self.acc]
        self.to_process(osi)


def build_regular_node_mesh(osi, xs, ys, zs=None, active=None):
    # axis-2 = x  # unless x or y are singlar
    # axis-1 = y
    # axis-1 = z
    from numpy import array
    if not hasattr(zs, '__len__'):
        zs = [zs]
    sn = []
    for xx in range(len(xs)):
        sn.append([])
        for yy in range(len(ys)):

            if len(zs) == 1:
                if active is None or active[xx][yy]:
                    sn[xx].append(Node(osi, xs[xx], ys[yy], zs[0]))
                else:
                    sn[xx].append(None)
            else:
                sn[xx].append([])
                for zz in range(len(zs)):
                    # Establish left and right nodes
                    if active is None or active[xx][yy][zz]:
                        sn[xx][yy].append(Node(osi, xs[xx], ys[yy], zs[zz]))
                    else:
                        sn[xx][yy].append(None)
    # if len(zs) == 1:
    #     return sn[0]
    return array(sn)

