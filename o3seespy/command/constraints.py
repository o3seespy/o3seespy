from o3seespy.base_model import OpenSeesObject


class ConstraintsBase(OpenSeesObject):
    op_base_type = "constraints"


class Plain(ConstraintsBase):
    """
    The Plain Constraints Class
    
    This command is used to construct a Plain constraint handler. A plain constraint handler can only enforce
    homogeneous single point constraints (fix command) and multi-point constraints constructed where the
    constraint matrix is equal to the identity (equalDOF command). The following is the command to
    construct a plain constraint handler:.. note::As mentioned, this constraint handler can only
    enforce homogeneous single point constraints (fix command) and multi-pont constraints where
    the constraint matrix is equal to the identity (equalDOF command).
    """
    op_type = 'Plain'

    def __init__(self, osi):
        """
        Initial method for Plain

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.constraints.Plain(osi)
        """
        self.osi = osi
        self._parameters = [self.op_type]
        self.to_process(osi)


class Lagrange(ConstraintsBase):
    """
    The Lagrange Constraints Class
    
    This command is used to construct a LagrangeMultiplier constraint handler, which enforces the constraints by
    introducing Lagrange multiplies to the system of equation. The following is the command to construct a plain
    constraint handler:
    """
    op_type = 'Lagrange'

    def __init__(self, osi, alpha_s=1.0, alpha_m=1.0):
        r"""
        Initial method for Lagrange

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        alpha_s: float, optional
            :math:`\alpha_s` factor on single points.
        alpha_m: float, optional
            :math:`\alpha_m` factor on multi-points.

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.constraints.Lagrange(osi, alpha_m=1.0)
        """
        self.osi = osi
        self.alpha_s = float(alpha_s)
        self.alpha_m = float(alpha_m)
        self._parameters = [self.op_type, self.alpha_s, self.alpha_m]
        self.to_process(osi)


class Penalty(ConstraintsBase):
    """
    The Penalty Constraints Class
    
    This command is used to construct a Penalty constraint handler, which enforces the constraints using the penalty
    method. The following is the command to construct a penalty constraint handler:
    """
    op_type = 'Penalty'

    def __init__(self, osi, alpha_s=1.0, alpha_m=1.0):
        r"""
        Initial method for Penalty

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        alpha_s: float, optional
            :math:`\alpha_s` factor on single points.
        alpha_m: float, optional
            :math:`\alpha_m` factor on multi-points.

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.constraints.Penalty(osi, alpha_m=1.0)
        """
        self.osi = osi
        self.alpha_s = float(alpha_s)
        self.alpha_m = float(alpha_m)
        self._parameters = [self.op_type, self.alpha_s, self.alpha_m]
        self.to_process(osi)


class Transformation(ConstraintsBase):
    """
    The Transformation Constraints Class
    
    This command is used to construct a transformation constraint handler, which enforces the constraints using the
    transformation method. The following is the command to construct a transformation constraint handler.. note::* The
    single-point constraints when using the transformation method are done directly. The matrix equation is not
    manipulated to enforce them, rather the trial displacements are set directly at the nodes at the start of
    each analysis step.* Great care must be taken when multiple constraints are being enforced as the
    transformation method does not follow constraints:#. If a node is fixed, constrain it with the
    fix command and not equalDOF or other type of constraint.#. If multiple nodes are
    constrained, make sure that the retained node is not constrained in any other
    constraint.And remember if a node is constrained to multiple nodes in your model it probably means you have messed up.
    """
    op_type = 'Transformation'

    def __init__(self, osi):
        """
        Initial method for Transformation

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.constraints.Transformation(osi)
        """
        self.osi = osi
        self._parameters = [self.op_type]
        self.to_process(osi)
