from o3seespy.command.nd_material.base_material import NDMaterialBase


class PlaneStressUserMaterial(NDMaterialBase):
    """
    The PlaneStressUserMaterial NDMaterial Class
    
    This command is used to create the multi-dimensional concrete material model that is based on the damage mechanism
    and smeared crack model.
    """
    op_type = 'PlaneStressUserMaterial'

    def __init__(self, osi, nstatevs, nprops, fc, ft, fcu, epsc0, epscu, epstu, stc):
        """
        Initial method for PlaneStressUserMaterial

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        nstatevs: None
            
        nprops: int
            Number of material properties (usually 7)
        fc: float
            Concrete compressive strength at 28 days (positive)
        ft: float
            Concrete tensile strength (positive)
        fcu: float
            Concrete crushing strength (negative)
        epsc0: float
            Concrete strain at maximum strength (negative)
        epscu: float
            Concrete strain at crushing strength (negative)
        epstu: float
            Ultimate tensile strain (positive)
        stc: float
            Shear retention factor

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> o3.nd_material.PlaneStressUserMaterial(osi, nstatevs=1, nprops=1, fc=1.0, ft=1.0, fcu=1.0, epsc0=1.0, epscu=1.0, epstu=1.0, stc=1.0)
        """
        self.osi = osi
        self.nstatevs = nstatevs
        self.nprops = int(nprops)
        self.fc = float(fc)
        self.ft = float(ft)
        self.fcu = float(fcu)
        self.epsc0 = float(epsc0)
        self.epscu = float(epscu)
        self.epstu = float(epstu)
        self.stc = float(stc)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.nstatevs, self.nprops, self.fc, self.ft, self.fcu, self.epsc0, self.epscu, self.epstu, self.stc]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class PlateFromPlaneStress(NDMaterialBase):
    """
    The PlateFromPlaneStress NDMaterial Class
    
    This command is used to create the multi-dimensional concrete material model that is based on the damage mechanism
    and smeared crack model.
    """
    op_type = 'PlateFromPlaneStress'

    def __init__(self, osi, pre_def_mat, outof_plane_modulus):
        """
        Initial method for PlateFromPlaneStress

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        pre_def_mat: obj
            Integer object identifying planestress material
        outof_plane_modulus: float
            Shear modulus for out of plane stresses

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> mat = o3.nd_material.ElasticIsotropic(osi, e_mod=1.0, nu=1.0, rho=0.0)
        >>> o3.nd_material.PlateFromPlaneStress(osi, pre_def_mat=mat, outof_plane_modulus=1.0)
        """
        self.osi = osi
        self.pre_def_mat = pre_def_mat
        self.outof_plane_modulus = float(outof_plane_modulus)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.pre_def_mat.tag, self.outof_plane_modulus]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class PlateRebar(NDMaterialBase):
    """
    The PlateRebar NDMaterial Class
    
    This command is used to create the multi-dimensional reinforcement material.
    """
    op_type = 'PlateRebar'

    def __init__(self, osi, pre_def_mat, sita):
        """
        Initial method for PlateRebar

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        pre_def_mat: obj
            Integer object identifying uniaxial material
        sita: float
            Define the angle of reinforcement layer, 90 (longitudinal), 0 (tranverse)

        Examples
        --------
        >>> import o3seespy as o3
        >>> osi = o3.OpenSeesInstance(ndm=2)
        >>> mat = o3.uniaxial_material.Elastic(osi, 1.0)
        >>> o3.nd_material.PlateRebar(osi, pre_def_mat=mat, sita=1.0)
        """
        self.osi = osi
        self.pre_def_mat = pre_def_mat
        self.sita = float(sita)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.pre_def_mat.tag, self.sita]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)


class PlasticDamageConcretePlaneStress(NDMaterialBase):
    """
    The PlasticDamageConcretePlaneStress NDMaterial Class
    
    No documentation is available yet. If you have the manual, please let me know.
    """
    op_type = 'PlasticDamageConcretePlaneStress'

    def __init__(self, osi, e_mod, nu, ft, fc, beta, ap, an, bn):
        """
        Initial method for PlasticDamageConcretePlaneStress

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        e_mod: None
            
        nu: None
            
        ft: None
            
        fc: None
            
        beta: None
            
        ap: None
            
        an: None
            
        bn: None
            
        """
        self.osi = osi
        self.e_mod = e_mod
        self.nu = nu
        self.ft = ft
        self.fc = fc
        self.beta = beta
        self.ap = ap
        self.an = an
        self.bn = bn
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.e_mod, self.nu, self.ft, self.fc, self.beta, self.ap, self.an, self.bn]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)
