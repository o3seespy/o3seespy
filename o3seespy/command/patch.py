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
        osi: o3seespy.OpenSeesInstance
        mat: obj
            Material object associated with this fiber (uniaxial_material object for a fibersection and ndmaterial
            object for use in an ndfibersection).
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

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> crds_i = [1.0, 1.0]
        >>> crds_j = [1.0, 1.0]
        >>> crds_k = [1.0, 1.0]
        >>> crds_l = [1.0, 1.0]
        >>> conc_conf = o3.uniaxial_material.Concrete01(osi, fpc=-5.0, epsc0=-0.005, fpcu=-3.5, eps_u=-0.02)
        >>> o3.patch.Quad(osi, conc_conf, num_subdiv_ij=1, num_subdiv_jk=1, crds_i=crds_i, crds_j=crds_j, crds_k=crds_k, crds_l=crds_l)
        """
        self.osi = osi
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
        osi: o3seespy.OpenSeesInstance
        mat: obj
            Material object associated with this fiber (uniaxial_material object for a fibersection and ndmaterial
            object for use in an ndfibersection).
        num_subdiv_y: int
            Number of subdivisions (fibers) in local y direction.
        num_subdiv_z: int
            Number of subdivisions (fibers) in local z direction.
        crds_i: list
            Y & z-coordinates of vertex i (local coordinate system)
        crds_j: list
            Y & z-coordinates of vertex j (local coordinate system)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> crds_i = [1.0, 1.0]
        >>> crds_j = [1.0, 1.0]
        >>> conc_conf = o3.uniaxial_material.Concrete01(osi, fpc=-5.0, epsc0=-0.005, fpcu=-3.5, eps_u=-0.02)
        >>> o3.patch.Rect(osi, conc_conf, num_subdiv_y=1, num_subdiv_z=1, crds_i=crds_i, crds_j=crds_j)
        """
        self.osi = osi
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
        osi: o3seespy.OpenSeesInstance
        mat: obj
            Material object associated with this fiber (uniaxial_material object for a fibersection and ndmaterial
            object for use in an ndfibersection).
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

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> center = [1.0, 1.0]
        >>> rad = [1.0, 1.0]
        >>> ang = [1.0, 1.0]
        >>> conc_conf = o3.uniaxial_material.Concrete01(osi, fpc=-5.0, epsc0=-0.005, fpcu=-3.5, eps_u=-0.02)
        >>> o3.patch.Circ(osi, conc_conf, num_subdiv_circ=1, num_subdiv_rad=1, center=center, rad=rad, ang=ang)
        """
        self.osi = osi
        self.mat = mat
        self.num_subdiv_circ = int(num_subdiv_circ)
        self.num_subdiv_rad = int(num_subdiv_rad)
        self.center = center
        self.rad = rad
        self.ang = ang
        self._parameters = [self.op_type, self.mat.tag, self.num_subdiv_circ, self.num_subdiv_rad, *self.center, *self.rad, *self.ang]
        self.to_process(osi)
