from o3seespy.base_model import OpenSeesObject


class GeomTransfBase(OpenSeesObject):
    op_base_type = "geomTransf"


class Linear2D(GeomTransfBase):
    """
    The Linear2D GeomTransf Class
    
    
    """
    op_type = 'Linear'

    def __init__(self, osi, d_i: list=None, d_j: list=None):
        """
        Initial method for Linear2D

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        d_i: list, optional
            Joint offset values -- offsets specified with respect to the global coordinate system for element-end node i
            (the number of arguments depends on the dimensions of the current model).
        d_j: list, optional
            Joint offset values -- offsets specified with respect to the global coordinate system for element-end node j
            (the number of arguments depends on the dimensions of the current model).

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> d_i = [1.0, 1.0]
        >>> d_j = [1.0, 1.0]
        >>> o3.geom_transf.Linear2D(osi, d_i=d_i, d_j=d_j)
        """
        self.osi = osi
        self.d_i = d_i
        self.d_j = d_j
        osi.n_transformation += 1
        self._tag = osi.n_transformation
        self._parameters = [self.op_type, self._tag]
        if getattr(self, 'd_i') is not None:
            self._parameters += ['-jntOffset', *self.d_i]
        if getattr(self, 'd_j') is not None:
            if getattr(self, 'd_i') is None:
                raise ValueError('Cannot set: d_j and not: d_i')
            self._parameters += [*self.d_j]
        self.to_process(osi)


class Linear3D(GeomTransfBase):
    """
    The Linear3D GeomTransf Class
    
    This command is used to construct a linear coordinate transformation (LinearCrdTransf) object, which performs a
    linear geometric transformation of beam stiffness and resisting force from the basic system to the
    global-coordinate system.
    """
    op_type = 'Linear'

    def __init__(self, osi, vecxz, d_i: list=None, d_j: list=None):
        """
        Initial method for Linear3D

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        vecxz: list
            X, y, and z components of vecxz, the vector used to define the local x-z plane of the local-coordinate
            system. the local y-axis is defined by taking the cross product of the vecxz vector and the x-axis. these
            components are specified in the global-coordinate system x,y,z and define a vector that is in a plane
            parallel to the x-z plane of the local-coordinate system. these items need to be specified for the
            three-dimensional problem.
        d_i: list, optional
            Joint offset values -- offsets specified with respect to the global coordinate system for element-end node i
            (the number of arguments depends on the dimensions of the current model).
        d_j: list, optional
            Joint offset values -- offsets specified with respect to the global coordinate system for element-end node j
            (the number of arguments depends on the dimensions of the current model).

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> vecxz = [1.0, 1.0]
        >>> d_i = [1.0, 1.0]
        >>> d_j = [1.0, 1.0]
        >>> o3.geom_transf.Linear3D(osi, vecxz=vecxz, d_i=d_i, d_j=d_j)
        """
        self.osi = osi
        self.vecxz = vecxz
        self.d_i = d_i
        self.d_j = d_j
        osi.n_transformation += 1
        self._tag = osi.n_transformation
        self._parameters = [self.op_type, self._tag, *self.vecxz]
        if getattr(self, 'd_i') is not None:
            self._parameters += ['-jntOffset', *self.d_i]
        if getattr(self, 'd_j') is not None:
            if getattr(self, 'd_i') is None:
                raise ValueError('Cannot set: d_j and not: d_i')
            self._parameters += [*self.d_j]
        self.to_process(osi)


class PDelta2D(GeomTransfBase):
    """
    The PDelta2D GeomTransf Class
    
    
    """
    op_type = 'PDelta'

    def __init__(self, osi, d_i: list=None, d_j: list=None):
        """
        Initial method for PDelta2D

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        d_i: list, optional
            Joint offset values -- offsets specified with respect to the global coordinate system for element-end node i
            (the number of arguments depends on the dimensions of the current model).
        d_j: list, optional
            Joint offset values -- offsets specified with respect to the global coordinate system for element-end node j
            (the number of arguments depends on the dimensions of the current model).

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> d_i = [1.0, 1.0]
        >>> d_j = [1.0, 1.0]
        >>> o3.geom_transf.PDelta2D(osi, d_i=d_i, d_j=d_j)
        """
        self.osi = osi
        self.d_i = d_i
        self.d_j = d_j
        osi.n_transformation += 1
        self._tag = osi.n_transformation
        self._parameters = [self.op_type, self._tag]
        if getattr(self, 'd_i') is not None:
            self._parameters += ['-jntOffset', *self.d_i]
        if getattr(self, 'd_j') is not None:
            if getattr(self, 'd_i') is None:
                raise ValueError('Cannot set: d_j and not: d_i')
            self._parameters += [*self.d_j]
        self.to_process(osi)


class PDelta3D(GeomTransfBase):
    """
    The PDelta3D GeomTransf Class
    
    This command is used to construct the P-Delta Coordinate Transformation (PDeltaCrdTransf) object, which performs a
    linear geometric transformation of beam stiffness and resisting force from the basic system to the global coordinate
    system, considering second-order P-Delta effects.
    """
    op_type = 'PDelta'

    def __init__(self, osi, vecxz, d_i: list=None, d_j: list=None):
        """
        Initial method for PDelta3D

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        vecxz: list
            X, y, and z components of vecxz, the vector used to define the local x-z plane of the local-coordinate
            system. the local y-axis is defined by taking the cross product of the vecxz vector and the x-axis. these
            components are specified in the global-coordinate system x,y,z and define a vector that is in a plane
            parallel to the x-z plane of the local-coordinate system. these items need to be specified for the
            three-dimensional problem.
        d_i: list, optional
            Joint offset values -- offsets specified with respect to the global coordinate system for element-end node i
            (the number of arguments depends on the dimensions of the current model).
        d_j: list, optional
            Joint offset values -- offsets specified with respect to the global coordinate system for element-end node j
            (the number of arguments depends on the dimensions of the current model).

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> vecxz = [1.0, 1.0]
        >>> d_i = [1.0, 1.0]
        >>> d_j = [1.0, 1.0]
        >>> o3.geom_transf.PDelta3D(osi, vecxz=vecxz, d_i=d_i, d_j=d_j)
        """
        self.osi = osi
        self.vecxz = vecxz
        self.d_i = d_i
        self.d_j = d_j
        osi.n_transformation += 1
        self._tag = osi.n_transformation
        self._parameters = [self.op_type, self._tag, *self.vecxz]
        if getattr(self, 'd_i') is not None:
            self._parameters += ['-jntOffset', *self.d_i]
        if getattr(self, 'd_j') is not None:
            if getattr(self, 'd_i') is None:
                raise ValueError('Cannot set: d_j and not: d_i')
            self._parameters += [*self.d_j]
        self.to_process(osi)


class Corotational2D(GeomTransfBase):
    """
    The Corotational2D GeomTransf Class
    
    
    """
    op_type = 'Corotational'

    def __init__(self, osi, d_i: list=None, d_j: list=None):
        """
        Initial method for Corotational2D

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        d_i: list, optional
            Joint offset values -- offsets specified with respect to the global coordinate system for element-end node i
            (the number of arguments depends on the dimensions of the current model).
        d_j: list, optional
            Joint offset values -- offsets specified with respect to the global coordinate system for element-end node j
            (the number of arguments depends on the dimensions of the current model).

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> d_i = [1.0, 1.0]
        >>> d_j = [1.0, 1.0]
        >>> o3.geom_transf.Corotational2D(osi, d_i=d_i, d_j=d_j)
        """
        self.osi = osi
        self.d_i = d_i
        self.d_j = d_j
        osi.n_transformation += 1
        self._tag = osi.n_transformation
        self._parameters = [self.op_type, self._tag]
        if getattr(self, 'd_i') is not None:
            self._parameters += ['-jntOffset', *self.d_i]
        if getattr(self, 'd_j') is not None:
            if getattr(self, 'd_i') is None:
                raise ValueError('Cannot set: d_j and not: d_i')
            self._parameters += [*self.d_j]
        self.to_process(osi)


class Corotational3D(GeomTransfBase):
    """
    The Corotational3D GeomTransf Class
    
    This command is used to construct the Corotational Coordinate Transformation (CorotCrdTransf) object. Corotational
    transformation can be used in large displacement-small strain problems.
    """
    op_type = 'Corotational'

    def __init__(self, osi, vecxz):
        """
        Initial method for Corotational3D

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        vecxz: list
            X, y, and z components of vecxz, the vector used to define the local x-z plane of the local-coordinate
            system. the local y-axis is defined by taking the cross product of the vecxz vector and the x-axis. these
            components are specified in the global-coordinate system x,y,z and define a vector that is in a plane
            parallel to the x-z plane of the local-coordinate system. these items need to be specified for the
            three-dimensional problem.

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> vecxz = [1.0, 1.0]
        >>> o3.geom_transf.Corotational3D(osi, vecxz=vecxz)
        """
        self.osi = osi
        self.vecxz = vecxz
        osi.n_transformation += 1
        self._tag = osi.n_transformation
        self._parameters = [self.op_type, self._tag, *self.vecxz]
        self.to_process(osi)
