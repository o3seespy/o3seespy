from o3seespy.command.uniaxial_material.base_material import UniaxialMaterialBase


class Elastic(UniaxialMaterialBase):
    type = "Elastic"

    def __init__(self, osi, e_mod, eta=0.0, e_mod_comp=None):
        """
        Elastic uniaxial material

        Parameters
        ----------
        osi : opensees_pack.opensees_instance.OpenseesInstance object
            An instance of opensees
        e_mod : float
            Tangent
        eta : float, optional
            Damping tangent
        e_mod_comp: float, optional
            Tangent in compression
        """
        self.e_mod = e_mod
        self.eta = eta
        if e_mod_comp is None:
            e_mod_comp = e_mod
        self.e_mod_comp = e_mod_comp
        osi.n_mats += 1
        self._tag = osi.n_mats
        self._parameters = [self.type, self._tag, self.e_mod, self.eta, self.e_mod_comp]
        self.to_process(osi)


class ElasticPP(UniaxialMaterialBase):
    type = "ElasticPP"

    def __init__(self, osi, e_mod, epsy_p, epsy_n=None, eps0=0.0):
        """
        Elastic perfectly plastic behaviour

        Parameters
        ----------
        osi : opensees_pack.opensees_instance.OpenseesInstance object
            An instance of opensees
        e_mod : float
            Tangent
        epsy_p : float
            Strain or deformation at which material reaches plastic state in tension
        epsy_n : float, optional
            Strain or deformation at which material reaches plastic state in compression.
        eps0 : float, optional
            Initial strain (default is 0.0)
        """
        self.e_mod = e_mod
        self.epsy_p = epsy_p
        if epsy_n is None:
            epsy_n = epsy_p
        self.epsy_n = epsy_n
        self.eps0 = eps0
        osi.n_mats += 1
        self._tag = osi.n_mats
        self._parameters = [self.type, self._tag, self.e_mod, self.epsy_p, self.epsy_n, self.eps0]
        self.to_process(osi)


class PySimple1(UniaxialMaterialBase):
    type = "PySimple1"

    def __init__(self, osi, soil_type: int, p_ult, y50, cd, c):
        """
        PySimple1 uniaxial material object

        Parameters
        ----------
        osi : opensees_pack.opensees_instance.OpenseesInstance object
            An instance of opensees
        soil_type : int {1, 2}
            Backbone type for soil
        p_ult : float
            Ultimate capacity of the p-y material. Note that “p” or “pult” are distributed loads
             [force per length of pile] in common design equations, but are both loads for this
             uniaxialMaterial [i.e., distributed load times the tributary length of the pile].
        y50 : float
            Displacement at which 50% of pult is mobilized in monotonic loading.
        cd : float
            Variable that sets the drag resistance within a fully-mobilized gap as Cd*pult.
        c : float
            The viscous damping term (dashpot) on the far-field (elastic) component
            of the displacement rate (velocity). (optional Default = 0.0). Nonzero c
            values are used to represent radiation damping effects
        """
        self.soil_type = int(soil_type)
        self.p_ult = p_ult
        self.y50 = y50
        self.cd = cd
        self.c = c
        osi.n_mats += 1
        self._tag = osi.n_mats
        self._parameters = [self.type, self._tag, self.p_ult, self.y50, self.cd, self.c]
        self.to_process(osi)


class PyLiq1(UniaxialMaterialBase):
    op_type = "PyLiq1"

    def __init__(self, osi, soil_type, p_ult, y50, cd, c, p_res, ele1=None, ele2=None, time_series=None):
        """

        Parameters
        ----------
        osi : opensees_pack.opensees_instance.OpenseesInstance object
            An instance of opensees
        soil_type
       p_ult : float
            Ultimate capacity of the p-y material. Note that “p” or “pult” are distributed loads
             [force per length of pile] in common design equations, but are both loads for this
             uniaxialMaterial [i.e., distributed load times the tributary length of the pile].
        y50 : float
            Displacement at which 50% of pult is mobilized in monotonic loading.
        cd : float
            Variable that sets the drag resistance within a fully-mobilized gap as Cd*pult.
        c : float
            The viscous damping term (dashpot) on the far-field (elastic) component
            of the displacement rate (velocity). (optional Default = 0.0). Nonzero c
            values are used to represent radiation damping effects
        p_res : float
        ele1 : opensees_pack.ElementBase object, optional
            the eleTag (element numbers) for the two solid element from which PyLiq1
            will obtain mean effective stresses and excess pore pressures
        ele2 : opensees_pack.ElementBase object, optional
            Same as ele1
        time_series : iterable object, optional
            A time series of mean effective stress values
        """
        self.soil_type = soil_type
        self.p_ult = p_ult
        self.y50 = y50
        self.cd = cd
        self.c = c
        self.p_res = p_res
        self.ele1 = ele1
        self.ele2 = ele2
        self.time_series = time_series
        osi.n_mats += 1
        self._tag = osi.n_mats
        self._parameters = [self.op_type, self._tag, self.p_ult, self.y50, self.cd, self.c, self.p_res]
        if self.ele1 is None:
            self._parameters.append("-timeSeries")
            self._parameters.append(self.time_series)
        else:
            self._parameters.append(self.ele1)
            if self.ele2 is None:
                self._parameters.append(self.ele2)
        self.to_process(osi)
