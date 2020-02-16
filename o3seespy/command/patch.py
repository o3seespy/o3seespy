from o3seespy.base_model import OpenSeesObject


class PatchBase(OpenSeesObject):
    op_base_type = "patch"


class Quad(PatchBase):
    """
    The Quad Patch Class
    
    The patch command is used to generate a number of fibers over a cross-sectional area. Currently there are three
    types of cross-section that fibers can be generated: quadrilateral, rectangular and circular.
    """
    op_type = 'quad'

    def __init__(self, osi, mat, num_subdiv_ij, num_subdiv_jk, crds_i, crds_j, crds_k, crds_l):
        """
        Initial method for Quad

        Parameters
        ----------
        mat: obj
            Material tag associated with this fiber (uniaxialmaterial tag for a fibersection and ndmaterial tag for use
            in an ndfibersection).
        num_subdiv_ij: int
            Number of subdivisions (fibers) in the ij direction.
        num_subdiv_jk: int
            Number of subdivisions (fibers) in the jk direction.
        crds_i: list
            Y & z-coordinates of vertex i (local coordinate system)
        crds_j: list
            Y & z-coordinates of vertex j (local coordinate system)
        crds_k: list
            Y & z-coordinates of vertex k (local coordinate system)
        crds_l: list
            Y & z-coordinates of vertex l (local coordinate system)
        """
        self.mat = mat
        self.num_subdiv_ij = int(num_subdiv_ij)
        self.num_subdiv_jk = int(num_subdiv_jk)
        self.crds_i = crds_i
        self.crds_j = crds_j
        self.crds_k = crds_k
        self.crds_l = crds_l
        self._parameters = [self.op_type, self.mat.tag, self.num_subdiv_ij, self.num_subdiv_jk, *self.crds_i, *self.crds_j, *self.crds_k, *self.crds_l]
        self.to_process(osi)

class Rect(PatchBase):
    """
    The Rect Patch Class
    
    This is the command to generate a rectangular patch. The geometry of the patch is defined by coordinates of
    vertices: I and J. The first vertex, I, is the bottom-left point and the second vertex, J, is the top-right
    point, having as a reference the local y-z plane.
    """
    op_type = 'rect'

    def __init__(self, osi, mat, num_subdiv_y, num_subdiv_z, crds_i, crds_j):
        """
        Initial method for Rect

        Parameters
        ----------
        mat: obj
            Material tag associated with this fiber (uniaxialmaterial tag for a fibersection and ndmaterial tag for use
            in an ndfibersection).
        num_subdiv_y: int
            Number of subdivisions (fibers) in local y direction.
        num_subdiv_z: int
            Number of subdivisions (fibers) in local z direction.
        crds_i: list
            Y & z-coordinates of vertex i (local coordinate system)
        crds_j: list
            Y & z-coordinates of vertex j (local coordinate system)
        """
        self.mat = mat
        self.num_subdiv_y = int(num_subdiv_y)
        self.num_subdiv_z = int(num_subdiv_z)
        self.crds_i = crds_i
        self.crds_j = crds_j
        self._parameters = [self.op_type, self.mat.tag, self.num_subdiv_y, self.num_subdiv_z, *self.crds_i, *self.crds_j]
        self.to_process(osi)

class Circ(PatchBase):
    """
    The Circ Patch Class
    
    This is the command to generate a circular shaped patch
    """
    op_type = 'circ'

    def __init__(self, osi, mat, num_subdiv_circ, num_subdiv_rad, center, rad, ang):
        """
        Initial method for Circ

        Parameters
        ----------
        mat: obj
            Material tag associated with this fiber (uniaxialmaterial tag for a fibersection and ndmaterial tag for use
            in an ndfibersection).
        num_subdiv_circ: int
            Number of subdivisions (fibers) in the circumferential direction (number of wedges)
        num_subdiv_rad: int
            Number of subdivisions (fibers) in the radial direction (number of rings)
        center: list
            Y & z-coordinates of the center of the circle
        rad: list
            Internal & external radius
        ang: list
            Starting & ending-coordinates angles (degrees)
        """
        self.mat = mat
        self.num_subdiv_circ = int(num_subdiv_circ)
        self.num_subdiv_rad = int(num_subdiv_rad)
        self.center = center
        self.rad = rad
        self.ang = ang
        self._parameters = [self.op_type, self.mat.tag, self.num_subdiv_circ, self.num_subdiv_rad, *self.center, *self.rad, *self.ang]
        self.to_process(osi)