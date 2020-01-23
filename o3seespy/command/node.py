from o3seespy.base_model import OpenseesObject
from o3seespy.command import common
import openseespy.opensees as opy


class Node(OpenseesObject):
    op_base_type = "node"
    op_type = "node"
    x_con = None
    y_con = None
    z_con = None
    x_rot_con = None
    y_rot_con = None
    z_rot_con = None
    x_mass = 0.0
    y_mass = 0.0
    z_rot_mass = 0.0

    def __init__(self, osi, x: float, y=None, z=None, vel=None, acc=None, x_mass=0.0, y_mass=0.0, z_rot_mass=0.0):
        """
        An Opensees node

        Parameters
        ----------
        osi : opensees_pack.opensees_instance.OpenseesInstance object
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
        self.vel = vel
        self.acc = acc
        self.x_mass = x_mass
        self.y_mass = y_mass
        self.z_rot_mass = z_rot_mass
        osi.n_node += 1
        self._tag = osi.n_node
        if osi.dimensions == 2:
            self._parameters = [self._tag, *[self.x, self.y]]
            self._parameters += ["-mass", *[self.x_mass, self.y_mass, self.z_rot_mass]]
        else:
            raise NotImplementedError("Currently only supports 2D analyses")
        if self.vel is not None:
            self._parameters += ["-vel", self.vel]
        if self.acc is not None:
            self._parameters += ["-accel", self.acc]
        self.to_process(osi)
