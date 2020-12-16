from o3seespy.base_model import OpenSeesObject


class MeshBase(OpenSeesObject):
    op_base_type = "mesh"


class Line(MeshBase):
    """
    The Line Mesh Class
    
    Create a line mesh object. 
    """
    op_type = 'line'

    def __init__(self, osi, numnodes, ndtags, id, ndf, meshsize, ele_type='', ele_args: float=None):
        """
        Initial method for Line

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        numnodes: int
            Number of nodes for defining consective lines.
        ndtags: list
            The node objects
        id: int
            Mesh id. meshes with same id are considered as same structure of fluid identity. * ``id`` = 0 : not in fsi *
            ``id`` > 0 : structure * ``id`` < 0 : fluid
        ndf: int
            Ndf for nodes to be created.
        meshsize: float
            Mesh size.
        ele_type: str, optional
            The type of the element,  * :ref:`elasticbeamcolumn` * :ref:`forcebeamcolumn-element` *
            :ref:`dispbeamcolumn-element` if no type is given, only nodes are created
        ele_args: list (default=True), optional
            A list of element arguments. the arguments are same as in the element commands, but without element object,
            and node objects.  for example, ``eleargs = ['elasticbeamcolumn', a, e, iz, transfobject]``

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> coords = [[0, 0], [0, 0]]
        >>> ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(2)]
        >>> o3.mesh.Line(osi, numnodes=1, ndtags=ele_nodes, id=1, ndf=1, meshsize=1.0, ele_type='', ele_args=None)
        """
        self.osi = osi
        self.numnodes = int(numnodes)
        self.ndtags = ndtags
        self.id = int(id)
        self.ndf = int(ndf)
        self.meshsize = float(meshsize)
        self.ele_type = ele_type
        self.ele_args = ele_args
        osi.n_mesh += 1
        self._tag = osi.n_mesh
        self._parameters = [self.op_type, self._tag, self.numnodes, *self.ndtags, self.id, self.ndf, self.meshsize, self.ele_type]
        special_pms = ['ele_args']
        packets = [True]
        for i, pm in enumerate(special_pms):
            if getattr(self, pm) is not None:
                if packets[i]:
                    self._parameters += [*getattr(self, pm)]
                else:
                    self._parameters += [getattr(self, pm)]
            else:
                break
        self.to_process(osi)


class Tri(MeshBase):
    """
    The Tri Mesh Class
    
    Create a triangular mesh object.
    """
    op_type = 'tri'

    def __init__(self, osi, numlines, ltags, id, ndf, meshsize, ele_type='', ele_args: float=None):
        """
        Initial method for Tri

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        numlines: int
            Number of lines (:ref:`linemesh`) for defining a polygon.
        ltags: list
            The :ref:`linemesh` objects
        id: int
            Mesh id. meshes with same id are considered as same structure of fluid identity. * ``id`` = 0 : not in fsi *
            ``id`` > 0 : structure * ``id`` < 0 : fluid
        ndf: int
            Ndf for nodes to be created.
        meshsize: float
            Mesh size.
        ele_type: str, optional
            The element type,  * :doc:`pfemelementbubble` * :doc:`pfemelementcompressible` * :doc:`tri31` *
            :doc:`elasticbeamcolumn` * :doc:`forcebeamcolumn` * :doc:`dispbeamcolumn` if no type is given, only nodes
            are created. if beam elements are given, beams are created instead of triangular elements.
        ele_args: list (default=True), optional
            A list of element arguments. the arguments are same as in the element commands, but without element object,
            and node objects.  for example, ``eleargs = ['pfemelementbubble', rho, mu, b1, b2, thickness, kappa]``

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> ltags = [1, 1]
        >>> o3.mesh.Tri(osi, numlines=1, ltags=ltags, id=1, ndf=1, meshsize=1.0, ele_type='', ele_args=None)
        """
        self.osi = osi
        self.numlines = int(numlines)
        self.ltags = ltags
        self.id = int(id)
        self.ndf = int(ndf)
        self.meshsize = float(meshsize)
        self.ele_type = ele_type
        self.ele_args = ele_args
        osi.n_mesh += 1
        self._tag = osi.n_mesh
        self._parameters = [self.op_type, self._tag, self.numlines, *self.ltags, self.id, self.ndf, self.meshsize, self.ele_type]
        special_pms = ['ele_args']
        packets = [True]
        for i, pm in enumerate(special_pms):
            if getattr(self, pm) is not None:
                if packets[i]:
                    self._parameters += [*getattr(self, pm)]
                else:
                    self._parameters += [getattr(self, pm)]
            else:
                break
        self.to_process(osi)


class Quad(MeshBase):
    """
    The Quad Mesh Class
    
    Create a quad mesh object. The number of lines must be 4. These lines are continuousto form a loop.
    """
    op_type = 'quad'

    def __init__(self, osi, numlines, ltags, id, ndf, meshsize, ele_type='', ele_args: float=None):
        """
        Initial method for Quad

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        numlines: int
            Number of lines (:ref:`linemesh`) for defining a polygon.
        ltags: list
            The :ref:`linemesh` objects
        id: int
            Mesh id. meshes with same id are considered as same structure of fluid identity. * ``id`` = 0 : not in fsi *
            ``id`` > 0 : structure * ``id`` < 0 : fluid
        ndf: int
            Ndf for nodes to be created.
        meshsize: float
            Mesh size.
        ele_type: str, optional
            The element type,  * :doc:`pfemelementbubble` * :doc:`pfemelementcompressible` * :doc:`tri31` *
            :doc:`elasticbeamcolumn` * :doc:`forcebeamcolumn` * :doc:`dispbeamcolumn` * :doc:`shellmitc4` if no type
            is given, only nodes are created. if beam elements are given, beams are created instead of quad
            elements. if triangular elements are given, they are created by dividing one quad to two triangles.
        ele_args: list (default=True), optional
            A list of element arguments. the arguments are same as in the element commands, but without element object,
            and node objects.  for example, ``eleargs = ['pfemelementbubble', rho, mu, b1, b2, thickness, kappa]``

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> ltags = [1, 1]
        >>> o3.mesh.Quad(osi, numlines=1, ltags=ltags, id=1, ndf=1, meshsize=1.0, ele_type='', ele_args=None)
        """
        self.osi = osi
        self.numlines = int(numlines)
        self.ltags = ltags
        self.id = int(id)
        self.ndf = int(ndf)
        self.meshsize = float(meshsize)
        self.ele_type = ele_type
        self.ele_args = ele_args
        osi.n_mesh += 1
        self._tag = osi.n_mesh
        self._parameters = [self.op_type, self._tag, self.numlines, *self.ltags, self.id, self.ndf, self.meshsize, self.ele_type]
        special_pms = ['ele_args']
        packets = [True]
        for i, pm in enumerate(special_pms):
            if getattr(self, pm) is not None:
                if packets[i]:
                    self._parameters += [*getattr(self, pm)]
                else:
                    self._parameters += [getattr(self, pm)]
            else:
                break
        self.to_process(osi)


class Tet(MeshBase):
    """
    The Tet Mesh Class
    
    Create a 3D tetrahedron mesh object.
    """
    op_type = 'tet'

    def __init__(self, osi, nummesh, mtags, id, ndf, meshsize, ele_type='', ele_args: float=None):
        """
        Initial method for Tet

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        nummesh: int
            Number of 2d mesh for defining a 3d body.
        mtags: list
            The mesh objects
        id: int
            Mesh id. meshes with same id are considered as same structure of fluid identity. * ``id`` = 0 : not in fsi *
            ``id`` > 0 : structure * ``id`` < 0 : fluid
        ndf: int
            Ndf for nodes to be created.
        meshsize: float
            Mesh size.
        ele_type: str, optional
            The element type,  * :doc:`fournodetetrahedron` if no type is given, only nodes are created.
        ele_args: list (default=True), optional
            A list of element arguments. the arguments are same as in the element commands, but without element object,
            and node objects. 

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> mtags = [1, 1]
        >>> o3.mesh.Tet(osi, nummesh=1, mtags=mtags, id=1, ndf=1, meshsize=1.0, ele_type='', ele_args=None)
        """
        self.osi = osi
        self.nummesh = int(nummesh)
        self.mtags = mtags
        self.id = int(id)
        self.ndf = int(ndf)
        self.meshsize = float(meshsize)
        self.ele_type = ele_type
        self.ele_args = ele_args
        osi.n_mesh += 1
        self._tag = osi.n_mesh
        self._parameters = [self.op_type, self._tag, self.nummesh, *self.mtags, self.id, self.ndf, self.meshsize, self.ele_type]
        special_pms = ['ele_args']
        packets = [True]
        for i, pm in enumerate(special_pms):
            if getattr(self, pm) is not None:
                if packets[i]:
                    self._parameters += [*getattr(self, pm)]
                else:
                    self._parameters += [getattr(self, pm)]
            else:
                break
        self.to_process(osi)


class Partvel(MeshBase):
    """
    The Partvel Mesh Class
    
    Create or return a group of particles which will be used for background mesh.
    """
    op_type = 'part'

    def __init__(self, osi, otype, p_args, vel0, ele_type='', ele_args: float=None):
        """
        Initial method for Partvel

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        otype: str
            Type of the mesh
        p_args: list
            Coordinates of points defining the mesh region and z directions * ``'quad'`` : [x1, y1, x2, y2, x3, y3, x4,
            y4, nx, ny] coordinates of four corners in counter-clock wise order. * ``'cube'`` : [x1, y1, z1, x2, y2, z2, x3, y3,
            z3, x4, y4, z4, x5, y5, z5, x6, y6, z6, x7, y7, z7, x8, y8, z8, nx, ny, nz] coordinates of four corners at bottom
            and at top in counter-clock wise order * ``'tri'`` : [x1, y1, x2, y2, x3, y3, nx, ny] coordinates of three
            corners in counter-clock wise order * ``'line'`` : [x1, y1, x2, y2, nx] coordinates of two ends in
            counter-clock wise order * ``'pointlist'`` : [x1, y1, <z1>, vx1, vy1, <vz1>, ax1, ay1, <az1>, p1,
             x2, y2, <z2>, vx2, vy2, <vz2>, ax2, ay2, <az2>, p2, ..] input particles' data in a list, in the
            order of coordinates of last time step, current coordinates, velocity, acceleration, and
            pressure. * ``'pointlist'`` without list return a list of current particles' data in
            this mesh [object1, x1, y1, <z1>, vx1, vy1, <vz1>, ax1, ay1, <az1>, p1, object2,
            x2, y2, <z2>, vx2, vy2, <vz2>, ax2, ay2, <az2>, p2, ..] the format is similar
            to the input list, but with an additional object for each particle.
        vel0: list
            A list of initial velocities. 
        ele_type: str, optional
            The element type,  * :doc:`pfemelementbubble` * :doc:`pfemelementcompressible` * :doc:`tri31` if no type is
            given, only nodes are created
        ele_args: list (default=True), optional
            A list of element arguments. (optional, see :doc:`linemesh` and :doc:`trimesh`)

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> p_args = [1.0, 1.0]
        >>> vel0 = [1.0, 1.0]
        >>> o3.mesh.Partvel(osi, otype="string", p_args=p_args, ele_type='', ele_args=None, vel0=vel0)
        """
        self.osi = osi
        self.otype = otype
        self.p_args = p_args
        self.ele_type = ele_type
        self.ele_args = ele_args
        self.vel0 = vel0
        osi.n_mesh += 1
        self._tag = osi.n_mesh
        self._parameters = [self.op_type, self._tag, self.otype, *self.p_args, self.ele_type]
        special_pms = ['ele_args', 'vel0', 'p0']
        packets = [True, True, False]
        for i, pm in enumerate(special_pms):
            if getattr(self, pm) is not None:
                if packets[i]:
                    self._parameters += [*getattr(self, pm)]
                else:
                    self._parameters += [getattr(self, pm)]
            else:
                break
        self.to_process(osi)

class Partpressure(MeshBase):
    """
    The Partpressure Mesh Class
    
    Create or return a group of particles which will be used for background mesh.
    """
    op_type = 'part'

    def __init__(self, osi, otype, p_args, p0, ele_type='', ele_args: float=None):
        """
        Initial method for Partpressure

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        otype: str
            Type of the mesh
        p_args: list
            Coordinates of points defining the mesh region and z directions * ``'quad'`` : [x1, y1, x2, y2, x3, y3, x4,
            y4, nx, ny] coordinates of four corners in counter-clock wise order. * ``'cube'`` : [x1, y1, z1, x2, y2, z2, x3, y3,
            z3, x4, y4, z4, x5, y5, z5, x6, y6, z6, x7, y7, z7, x8, y8, z8, nx, ny, nz] coordinates of four corners at bottom
            and at top in counter-clock wise order * ``'tri'`` : [x1, y1, x2, y2, x3, y3, nx, ny] coordinates of three
            corners in counter-clock wise order * ``'line'`` : [x1, y1, x2, y2, nx] coordinates of two ends in
            counter-clock wise order * ``'pointlist'`` : [x1, y1, <z1>, vx1, vy1, <vz1>, ax1, ay1, <az1>, p1,
             x2, y2, <z2>, vx2, vy2, <vz2>, ax2, ay2, <az2>, p2, ..] input particles' data in a list, in the
            order of coordinates of last time step, current coordinates, velocity, acceleration, and
            pressure. * ``'pointlist'`` without list return a list of current particles' data in
            this mesh [object1, x1, y1, <z1>, vx1, vy1, <vz1>, ax1, ay1, <az1>, p1, object2,
            x2, y2, <z2>, vx2, vy2, <vz2>, ax2, ay2, <az2>, p2, ..] the format is similar
            to the input list, but with an additional object for each particle.
        p0: float
            Initial pressure. 
        ele_type: str, optional
            The element type,  * :doc:`pfemelementbubble` * :doc:`pfemelementcompressible` * :doc:`tri31` if no type is
            given, only nodes are created
        ele_args: list (default=True), optional
            A list of element arguments. (optional, see :doc:`linemesh` and :doc:`trimesh`)

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> p_args = [1.0, 1.0]
        >>> o3.mesh.Partpressure(osi, otype="string", p_args=p_args, ele_type='', ele_args=None, p0=1.0)
        """
        self.osi = osi
        self.otype = otype
        self.p_args = p_args
        self.ele_type = ele_type
        self.ele_args = ele_args
        self.p0 = float(p0)
        osi.n_mesh += 1
        self._tag = osi.n_mesh
        self._parameters = [self.op_type, self._tag, self.otype, *self.p_args, self.ele_type]
        special_pms = ['ele_args', 'vel0', 'p0']
        packets = [True, True, False]
        for i, pm in enumerate(special_pms):
            if getattr(self, pm) is not None:
                if packets[i]:
                    self._parameters += [*getattr(self, pm)]
                else:
                    self._parameters += [getattr(self, pm)]
            else:
                break
        self.to_process(osi)


class Bgwave(MeshBase):
    """
    The Bgwave Mesh Class
    
    Create a background mesh. 
    """
    op_type = 'bg'

    def __init__(self, osi, lower, upper, wavefilename, numl, locations, tol: float=None, meshtol: float=None, numsub=None):
        """
        Initial method for Bgwave

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        lower: list
            A list of coordinates of the lower point of the background region.
        upper: list
            A list of coordinates of the uuper point of the background region.
        wavefilename: str
            A filename to record wave heights and velocities 
        numl: int
            Number of locations to record wave 
        locations: list
            Coordinates of the locations 
        tol: float, optional
            Tolerance for intri check. (optional, default 1e-10)
        meshtol: float, optional
            Tolerance for cell boundary check. (optional, default 0.1)
        numsub: None, optional
            

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> lower = [1.0, 1.0]
        >>> upper = [1.0, 1.0]
        >>> locations = [1.0, 1.0]
        >>> o3.mesh.Bgwave(osi, lower=lower, upper=upper, tol=1.0, meshtol=1.0, wavefilename="string", numl=1, locations=locations, numsub=1)
        """
        self.osi = osi
        self.lower = lower
        self.upper = upper
        if tol is None:
            self.tol = None
        else:
            self.tol = float(tol)
        if meshtol is None:
            self.meshtol = None
        else:
            self.meshtol = float(meshtol)
        self.wavefilename = wavefilename
        self.numl = int(numl)
        self.locations = locations
        self.numsub = numsub
        osi.n_mesh += 1
        self._tag = osi.n_mesh
        self._parameters = [self.op_type, self._tag, *self.lower, *self.upper, '-wave', self.wavefilename, self.numl, *self.locations]
        if getattr(self, 'tol') is not None:
            self._parameters += ['-tol', self.tol]
        if getattr(self, 'meshtol') is not None:
            self._parameters += ['-meshtol', self.meshtol]
        if getattr(self, 'numsub') is not None:
            self._parameters += ['-numsub', self.numsub]
        self.to_process(osi)

class Bgstructure(MeshBase):
    """
    The Bgstructure Mesh Class
    
    Create a background mesh. 
    """
    op_type = 'bg'

    def __init__(self, osi, lower, upper, id, numnodes, snodes, tol: float=None, meshtol: float=None):
        """
        Initial method for Bgstructure

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        lower: list
            A list of coordinates of the lower point of the background region.
        upper: list
            A list of coordinates of the uuper point of the background region.
        id: int
            Structural id > 0, same meaning as :doc:`trimesh` 
        numnodes: None
            
        snodes: None
            
        tol: float, optional
            Tolerance for intri check. (optional, default 1e-10)
        meshtol: float, optional
            Tolerance for cell boundary check. (optional, default 0.1)

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> lower = [1.0, 1.0]
        >>> upper = [1.0, 1.0]
        >>> o3.mesh.Bgstructure(osi, lower=lower, upper=upper, tol=1.0, meshtol=1.0, id=1, numnodes=1, snodes=1)
        """
        self.osi = osi
        self.lower = lower
        self.upper = upper
        if tol is None:
            self.tol = None
        else:
            self.tol = float(tol)
        if meshtol is None:
            self.meshtol = None
        else:
            self.meshtol = float(meshtol)
        self.id = int(id)
        self.numnodes = numnodes
        self.snodes = snodes
        osi.n_mesh += 1
        self._tag = osi.n_mesh
        self._parameters = [self.op_type, self._tag, *self.lower, *self.upper, '-structure', self.id, self.numnodes, *self.snodes]
        if getattr(self, 'tol') is not None:
            self._parameters += ['-tol', self.tol]
        if getattr(self, 'meshtol') is not None:
            self._parameters += ['-meshtol', self.meshtol]
        self.to_process(osi)

class BglargeSize(MeshBase):
    """
    The BglargeSize Mesh Class
    
    Create a background mesh. 
    """
    op_type = 'bg'

    def __init__(self, osi, lower, upper, level, llower, lupper, tol: float=None, meshtol: float=None):
        """
        Initial method for BglargeSize

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        lower: list
            A list of coordinates of the lower point of the background region.
        upper: list
            A list of coordinates of the uuper point of the background region.
        level: int
            Some regions can have larger mesh size with larger ``level``. ``level = 1`` means same as basic mesh size.
        llower: list
            A list of coordinates of the lower point of the region with larger mesh size 
        lupper: list
            A list of coordinates of the upper point of the region with larger mesh size
        tol: float, optional
            Tolerance for intri check. (optional, default 1e-10)
        meshtol: float, optional
            Tolerance for cell boundary check. (optional, default 0.1)

        Examples
        --------
        >>> import o3seespy as o3
        >>> # Example is currently not working
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> lower = [1.0, 1.0]
        >>> upper = [1.0, 1.0]
        >>> llower = [1.0, 1.0]
        >>> lupper = [1.0, 1.0]
        >>> o3.mesh.BglargeSize(osi, lower=lower, upper=upper, tol=1.0, meshtol=1.0, level=1, llower=llower, lupper=lupper)
        """
        self.osi = osi
        self.lower = lower
        self.upper = upper
        if tol is None:
            self.tol = None
        else:
            self.tol = float(tol)
        if meshtol is None:
            self.meshtol = None
        else:
            self.meshtol = float(meshtol)
        self.level = int(level)
        self.llower = llower
        self.lupper = lupper
        osi.n_mesh += 1
        self._tag = osi.n_mesh
        self._parameters = [self.op_type, self._tag, *self.lower, *self.upper, '-largeSize', self.level, *self.llower, *self.lupper]
        if getattr(self, 'tol') is not None:
            self._parameters += ['-tol', self.tol]
        if getattr(self, 'meshtol') is not None:
            self._parameters += ['-meshtol', self.meshtol]
        self.to_process(osi)
