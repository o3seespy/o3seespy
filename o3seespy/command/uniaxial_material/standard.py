from o3seespy.command.uniaxial_material.base_material import UniaxialMaterialBase



class Elastic(UniaxialMaterialBase):
    """
    The Elastic UniaxialMaterial Class
    
    This command is used to construct an elastic uniaxial material object.
    """
    op_type = 'Elastic'

    def __init__(self, osi, e_mod, eta=0.0, eneg: float=None):
        """
        Initial method for Elastic

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        e_mod: float
            Tangent
        eta: float, optional
            Damping tangent (optional, default=0.0)
        eneg: float (default=True), optional
            Tangent in compression (optional, default=e)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
        """
        self.osi = osi
        self.e_mod = float(e_mod)
        self.eta = float(eta)
        if eneg is None:
            self.eneg = None
        else:
            self.eneg = float(eneg)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.e_mod, self.eta]
        special_pms = ['eneg']
        packets = [False]
        for i, pm in enumerate(special_pms):
            if getattr(self, pm) is not None:
                if packets[i]:
                    self._parameters += [*getattr(self, pm)]
                else:
                    self._parameters += [getattr(self, pm)]
            else:
                break
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def set_e_mod(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'E', value, ele, eles)

    def set_epos(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'Epos', value, ele, eles)

    def set_eneg(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'Eneg', value, ele, eles)

    def set_eta(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'eta', value, ele, eles)


class ElasticPP(UniaxialMaterialBase):
    """
    The ElasticPP UniaxialMaterial Class
    
    This command is used to construct an elastic perfectly-plastic uniaxial material object.
    """
    op_type = 'ElasticPP'

    def __init__(self, osi, e_mod, epsy_p, epsy_n: float=None, eps0=0.0):
        """
        Initial method for ElasticPP

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        e_mod: float
            Tangent
        epsy_p: float
            Strain or deformation at which material reaches plastic state in tension
        epsy_n: float (default=True), optional
            Strain or deformation at which material reaches plastic state in compression. (optional, default is tension
            value)
        eps0: float, optional
            Initial strain (optional, default: zero)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.ElasticPP(osi, e_mod=1.0, epsy_p=1.0, epsy_n=None, eps0=0.0)
        """
        self.osi = osi
        self.e_mod = float(e_mod)
        self.epsy_p = float(epsy_p)
        if epsy_n is None:
            self.epsy_n = None
        else:
            self.epsy_n = float(epsy_n)
        self.eps0 = float(eps0)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.e_mod, self.epsy_p]
        special_pms = ['epsy_n', 'eps0']
        packets = [False, False]
        for i, pm in enumerate(special_pms):
            if getattr(self, pm) is not None:
                if packets[i]:
                    self._parameters += [*getattr(self, pm)]
                else:
                    self._parameters += [getattr(self, pm)]
            else:
                break
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def set_fy(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'Fy', value, ele, eles)

    def set_e_mod(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'E', value, ele, eles)

    def set_ep(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'ep', value, ele, eles)


class ElasticPPGap(UniaxialMaterialBase):
    """
    The ElasticPPGap UniaxialMaterial Class
    
    This command is used to construct an elastic perfectly-plastic gap uniaxial material object.
    """
    op_type = 'ElasticPPGap'

    def __init__(self, osi, e_mod, fy, gap, eta=0.0, damage='noDamage'):
        """
        Initial method for ElasticPPGap

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        e_mod: float
            Tangent
        fy: float
            Stress or force at which material reaches plastic state
        gap: float
            Initial gap (strain or deformation)
        eta: float, optional
            Hardening ratio (=eh/e), which can be negative
        damage: str, optional
            An optional string to specify whether to accumulate damage or not in the material. with the default
            re-center on load reversal. is provided this recentering will not occur and gap will grow.

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.ElasticPPGap(osi, e_mod=1.0, fy=1.0, gap=1.0, eta=0.0, damage='noDamage')
        """
        self.osi = osi
        self.e_mod = float(e_mod)
        self.fy = float(fy)
        self.gap = float(gap)
        self.eta = float(eta)
        self.damage = damage
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.e_mod, self.fy, self.gap, self.eta, self.damage]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class ENT(UniaxialMaterialBase):
    """
    The ENT UniaxialMaterial Class
    
    This command is used to construct a uniaxial elastic-no tension material object.
    """
    op_type = 'ENT'

    def __init__(self, osi, e_mod):
        """
        Initial method for ENT

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        e_mod: float
            Tangent

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.uniaxial_material.ENT(osi, e_mod=1.0)
        """
        self.osi = osi
        self.e_mod = float(e_mod)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.e_mod]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def set_e_mod(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'E', value, ele, eles)


class Parallel(UniaxialMaterialBase):
    """
    The Parallel UniaxialMaterial Class
    
    This command is used to construct a parallel material object made up of an arbitrary number of
    previously-constructed UniaxialMaterial objects.
    """
    op_type = 'Parallel'

    def __init__(self, osi, mats, factor_args: list=None):
        """
        Initial method for Parallel

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        mats: list
            Identification objects of materials making up the material model
        factor_args: list, optional
            Factors to create a linear combination of the specified materials. factors can be negative to subtract one
            material from an other. (optional, default = 1.0)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> mats = [o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None),
        >>>         o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)]
        >>> factor_args = [1.0, 1.0]
        >>> o3.uniaxial_material.Parallel(osi, mats=mats, factor_args=factor_args)
        """
        self.osi = osi
        self.mats = [x.tag for x in mats]
        self.factor_args = factor_args
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, *self.mats]
        if getattr(self, 'factor_args') is not None:
            self._parameters += ['-factors', *self.factor_args]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class Series(UniaxialMaterialBase):
    """
    The Series UniaxialMaterial Class
    
    This command is used to construct a series material object made up of an arbitrary number of previously-constructed
    UniaxialMaterial objects.
    """
    op_type = 'Series'

    def __init__(self, osi, mats):
        """
        Initial method for Series

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        mats: list
            Identification objects of materials making up the material model

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> mats = [o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None),
        >>>         o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)]
        >>> o3.uniaxial_material.Series(osi, mats=mats)
        """
        self.osi = osi
        self.mats = [x.tag for x in mats]
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, *self.mats]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)
