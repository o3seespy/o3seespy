from o3seespy.base_model import OpenSeesObject


class NumbererBase(OpenSeesObject):
    op_base_type = "numberer"


class Plain(NumbererBase):
    op_type = "Plain"

    def __init__(self, osi):
        """
        This command is used to construct a Plain degree-of-freedom numbering object to provide the mapping between
        the degrees-of-freedom at the nodes and the equation numbers.
        A Plain numberer just takes whatever order the domain gives it nodes and numbers them, this ordering is both
        dependent on node numbering and size of the model.

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        """
        self.osi = osi
        self._parameters = [self.op_type]
        self.to_process(osi)


class RCM(NumbererBase):
    op_type = "RCM"

    def __init__(self, osi):
        """
        This command is used to construct an RCM degree-of-freedom numbering object to provide the mapping
        between the degrees-of-freedom at the nodes and the equation numbers. An RCM numberer uses the reverse
        Cuthill-McKee scheme to order the matrix equations.

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        """
        self.osi = osi
        self._parameters = [self.op_type]
        self.to_process(osi)


class AMD(NumbererBase):
    op_type = "AMD"

    def __init__(self, osi):
        """
        This command is used to construct an AMD degree-of-freedom numbering object to provide the mapping between
        the degrees-of-freedom at the nodes and the equation numbers. An AMD numberer uses the approximate
        minimum degree scheme to order the matrix equations.

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        """
        self.osi = osi
        self._parameters = [self.op_type]
        self.to_process(osi)


class ParallelPlain(NumbererBase):
    op_type = "ParallelPlain"

    def __init__(self, osi):
        """
        This command is used to construct a parallel version of Plain degree-of-freedom numbering object to
        provide the mapping between the degrees-of-freedom at the nodes and the equation numbers.
        A Plain numberer just takes whatever order the domain gives it nodes and numbers them,
        this ordering is both dependent on node numbering and size of the model.

        Use this command only for parallel model.

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        """
        self.osi = osi
        self._parameters = [self.op_type]
        self.to_process(osi)


class ParallelRCM(NumbererBase):
    op_type = "ParallelRCM"

    def __init__(self, osi):
        """
        This command is used to construct a parallel version of RCM degree-of-freedom numbering object to
        provide the mapping between the degrees-of-freedom at the nodes and the equation numbers.
        A Plain numberer just takes whatever order the domain gives it nodes and numbers them,
        this ordering is both dependent on node numbering and size of the model.

        Use this command only for parallel model.

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        """
        self.osi = osi
        self._parameters = [self.op_type]
        self.to_process(osi)
