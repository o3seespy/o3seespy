from o3seespy.command.nd_material.base_material import NDMaterialBase


class PlaneStressUserMaterial(NDMaterialBase):
    """
    The PlaneStressUserMaterial NDMaterial Class
    
    This command is used to create the multi-dimensional concrete material model that is based on the damage mechanism
    and smeared crack model.
    """
    op_type = 'PlaneStressUserMaterial'

    def __init__(self, osi, fc, ft, fcu, epsc0, epscu, epstu, stc):
        """
        Initial method for PlaneStressUserMaterial

        Parameters
        ----------
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
        """
        self.fc = float(fc)
        self.ft = float(ft)
        self.fcu = float(fcu)
        self.epsc0 = float(epsc0)
        self.epscu = float(epscu)
        self.epstu = float(epstu)
        self.stc = float(stc)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.fc, self.ft, self.fcu, self.epsc0, self.epscu, self.epstu, self.stc]
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
        pre_def_mat: obj
            Integer tag identifying planestress material
        outof_plane_modulus: float
            Shear modulus for out of plane stresses
        """
        self.pre_def_mat = pre_def_mat
        self.outof_plane_modulus = float(outof_plane_modulus)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.pre_def_mat.tag, self.outof_plane_modulus]
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
        pre_def_mat: obj
            Integer tag identifying uniaxial material
        sita: float
            Define the angle of reinforcement layer, 90 (longitudinal), 0 (tranverse)
        """
        self.pre_def_mat = pre_def_mat
        self.sita = float(sita)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.pre_def_mat.tag, self.sita]
        self.to_process(osi)
