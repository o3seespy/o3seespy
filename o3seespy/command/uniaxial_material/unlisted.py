from o3seespy.command.uniaxial_material.base_material import UniaxialMaterialBase


class PySimple2(UniaxialMaterialBase):
    """
    The PySimple2 UniaxialMaterial Class

    This command is used to construct a PySimple2 uniaxial material object.
    """
    op_type = 'PySimple2'

    def __init__(self, osi, soil_type, pult, y50, cd, c=0.0):
        """
        Initial method for PySimple2

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        soil_type: int
            Soiltype = 1 backbone of p-y curve approximates matlock (1970) soft clay relation. soiltype = 2 backbone of
            p-y curve approximates api (1993) sand relation.
        pult: float
            Ultimate capacity of the p-y material. note that "p" or "pult" are distributed loads [force per length of
            pile] in common design equations, but are both loads for this uniaxial_material [i.e., distributed load times the
            tributary length of the pile].
        y50: float
            Displacement at which 50% of pult is mobilized in monotonic loading.
        cd: float
            Variable that sets the drag resistance within a fully-mobilized gap as cd*pult.
        c: float, optional
            The viscous damping term (dashpot) on the far-field (elastic) component of the displacement rate (velocity).
            (optional default = 0.0). nonzero c values are used to represent radiation damping effects

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.PySimple2(osi, soil_type=1, pult=1.0, y50=1.0, cd=1.0, c=0.0)
        """
        self.osi = osi
        self.soil_type = int(soil_type)
        self.pult = float(pult)
        self.y50 = float(y50)
        self.cd = float(cd)
        self.c = float(c)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.soil_type, self.pult, self.y50, self.cd, self.c]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class QzSimple2(UniaxialMaterialBase):
    """
    The QzSimple2 UniaxialMaterial Class

    This command is used to construct a QzSimple2 uniaxial material object.
    """
    op_type = 'QzSimple2'

    def __init__(self, osi, qz_type, qult, z50, suction=0.0, c=0.0):
        """
        Initial method for QzSimple2

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        qz_type: int
            Qztype = 1 backbone of q-z curve approximates reese and o'neill's (1987) relation for drilled shafts in
            clay. qztype = 2 backbone of q-z curve approximates vijayvergiya's (1977) relation for piles in sand.
        qult: float
            Ultimate capacity of the q-z material. see note 1.
        z50: float
            Displacement at which 50% of qult is mobilized in monotonic loading. see note 2.
        suction: float, optional
            Uplift resistance is equal to suction*qult. default = 0.0. the value of suction must be 0.0 to 0.1.*
        c: float, optional
            The viscous damping term (dashpot) on the far-field (elastic) component of the displacement rate (velocity).
            default = 0.0. nonzero c values are used to represent radiation damping effects.*

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.QzSimple2(osi, qz_type=1, qult=1.0, z50=1.0, suction=0.0, c=0.0)
        """
        self.osi = osi
        self.qz_type = int(qz_type)
        self.qult = float(qult)
        self.z50 = float(z50)
        self.suction = float(suction)
        self.c = float(c)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.qz_type, self.qult, self.z50, self.suction, self.c]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class TzSimple2(UniaxialMaterialBase):
    """
    The TzSimple2 UniaxialMaterial Class

    This command is used to construct a TzSimple2 uniaxial material object.
    """
    op_type = 'TzSimple2'

    def __init__(self, osi, soil_type, tult, z50, c=0.0):
        """
        Initial method for TzSimple2

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        soil_type: int
            Soiltype = 1 backbone of t-z curve approximates reese and o'neill (1987). soiltype = 2 backbone of t-z curve
            approximates mosher (1984) relation.
        tult: float
            Ultimate capacity of the t-z material. see note 1.
        z50: float
            Displacement at which 50% of tult is mobilized in monotonic loading.
        c: float, optional
            The viscous damping term (dashpot) on the far-field (elastic) component of the displacement rate (velocity).
            (optional default = 0.0). see note 2.

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.TzSimple2(osi, soil_type=1, tult=1.0, z50=1.0, c=0.0)
        """
        self.osi = osi
        self.soil_type = int(soil_type)
        self.tult = float(tult)
        self.z50 = float(z50)
        self.c = float(c)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.soil_type, self.tult, self.z50, self.c]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)
