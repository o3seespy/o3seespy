from o3seespy.command.nd_material.base_material import NDMaterialBase


class FluidSolidPorous(NDMaterialBase):
    """
    The FluidSolidPorous NDMaterial Class
    
    FluidSolidPorous material couples the responses of two phases: fluid and solid. The fluid phase response is only
    volumetric and linear elastic. The solid phase can be any NDMaterial. This material is developed to simulate the
    response of saturated porous media under fully undrained condition.
    """
    op_type = 'FluidSolidPorous'

    def __init__(self, osi, nd, soil_mat, combined_bulk_modul, pa=101.0):
        r"""
        Initial method for FluidSolidPorous

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        nd: float
            Number of dimensions, 2 for plane-strain, and 3 for 3d analysis.
        soil_mat: obj
            The material number for the solid phase material (previously defined).
        combined_bulk_modul: float
            Combined undrained bulk modulus :math:`b_c` relating changes in pore pressure and volumetric strain, may be
            approximated by: :math:`b_c \approx b_f /n` where :math:`b_f` is the bulk modulus of fluid phase (2.2x106 kpa (or
            3.191x105 psi) for water), and :math:`n` the initial porosity.
        pa: float, optional
            Optional atmospheric pressure for normalization (typically 101 kpa in si units, or 14.65 psi in english
            units)
        """
        self.osi = osi
        self.nd = float(nd)
        self.soil_mat = soil_mat
        self.combined_bulk_modul = float(combined_bulk_modul)
        self.pa = float(pa)
        if osi is not None:
            osi.n_mat += 1
            self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.nd, self.soil_mat.tag, self.combined_bulk_modul, self.pa]
        if osi is None:
            self.built = 0
        if osi is not None:
            self.to_process(osi)

    def set_update_material_stage(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'updateMaterialStage', value, ele, eles)

    def set_combined_bulk_modulus(self, value, ele=None, eles=None):
        self.set_parameter(self.osi, 'combinedBulkModulus', value, ele, eles)
