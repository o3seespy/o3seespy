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
        e_mod: float
            Tangent
        eta: float
            Damping tangent (optional, default=0.0)
        eneg: float (default=True)
            Tangent in compression (optional, default=e)
        """
        self.e_mod = float(e_mod)
        self.eta = float(eta)
        if eneg is None:
            self.eneg = None
        else:
            self.eneg = float(eneg)
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
        self.to_process(osi)


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
        e_mod: float
            Tangent
        epsy_p: float
            Strain or deformation at which material reaches plastic state in tension
        epsy_n: float (default=True)
            Strain or deformation at which material reaches plastic state in compression. (optional, default is tension
            value)
        eps0: float
            Initial strain (optional, default: zero)
        """
        self.e_mod = float(e_mod)
        self.epsy_p = float(epsy_p)
        if epsy_n is None:
            self.epsy_n = None
        else:
            self.epsy_n = float(epsy_n)
        self.eps0 = float(eps0)
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
        self.to_process(osi)


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
        e_mod: float
            Tangent
        fy: float
            Stress or force at which material reaches plastic state
        gap: float
            Initial gap (strain or deformation)
        eta: float
            Hardening ratio (=eh/e), which can be negative
        damage: str
            An optional string to specify whether to accumulate damage or not in the material. with the default
            re-center on load reversal. is provided this recentering will not occur and gap will grow.
        """
        self.e_mod = float(e_mod)
        self.fy = float(fy)
        self.gap = float(gap)
        self.eta = float(eta)
        self.damage = damage
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.e_mod, self.fy, self.gap, self.eta, self.damage]
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
        e_mod: float
            Tangent
        """
        self.e_mod = float(e_mod)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.e_mod]
        self.to_process(osi)


class Parallel(UniaxialMaterialBase):
    """
    The Parallel UniaxialMaterial Class
    
    This command is used to construct a parallel material object made up of an arbitrary number of
    previously-constructed UniaxialMaterial objects.
    """
    op_type = 'Parallel'

    def __init__(self, osi, tags, factor_args=None):
        """
        Initial method for Parallel

        Parameters
        ----------
        tags: listi
            Identification tags of materials making up the material model
        factor_args: listf
            Factors to create a linear combination of the specified materials. factors can be negative to subtract one
            material from an other. (optional, default = 1.0)
        """
        self.tags = tags
        self.factor_args = factor_args
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, *self.tags]
        if getattr(self, 'factor_args') is not None:
            self._parameters += ['-factor', *self.factor_args]
        self.to_process(osi)


class Series(UniaxialMaterialBase):
    """
    The Series UniaxialMaterial Class
    
    This command is used to construct a series material object made up of an arbitrary number of previously-constructed
    UniaxialMaterial objects.
    """
    op_type = 'Series'

    def __init__(self, osi, tags):
        """
        Initial method for Series

        Parameters
        ----------
        tags: listi
            Identification tags of materials making up the material model
        """
        self.tags = tags
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, *self.tags]
        self.to_process(osi)
