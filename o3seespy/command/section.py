from o3seespy.base_model import OpenSeesObject


class SectionBase(OpenSeesObject):
    op_base_type = "section"


class Elastic2D(SectionBase):
    """
    The Elastic2D Section Class
    
    
    """
    op_type = 'Elastic'

    def __init__(self, osi, e_mod, area, iz, g_mod: float=None, alpha_y: float=None):
        """
        Initial method for Elastic2D

        Parameters
        ----------
        e_mod: float
            Young's modulus
        area: float
            Cross-sectional area of section
        iz: float
            Second moment of area about the local z-axis
        g_mod: float (default=True)
            Shear modulus (optional for 2d analysis, required for 3d analysis)
        alpha_y: float (default=True)
            Shear shape factor along the local y-axis (optional)
        """
        self.e_mod = float(e_mod)
        self.area = float(area)
        self.iz = float(iz)
        if g_mod is None:
            self.g_mod = None
        else:
            self.g_mod = float(g_mod)
        if alpha_y is None:
            self.alpha_y = None
        else:
            self.alpha_y = float(alpha_y)
        osi.n_sect += 1
        self._tag = osi.n_sect
        self._parameters = [self.op_type, self._tag, self.e_mod, self.area, self.iz]
        special_pms = ['g_mod', 'alpha_y']
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


class Elastic3D(SectionBase):
    """
    The Elastic3D Section Class
    
    This command allows the user to construct an ElasticSection. The inclusion of shear deformations is optional. The
    dofs for 2D elastic section are ``[P, Mz]``,for 3D are ``[P,Mz,My,T]``.
    """
    op_type = 'Elastic'

    def __init__(self, osi, e_mod, area, iz, iy, g_mod, jxx, alpha_y: float=None, alpha_z: float=None):
        """
        Initial method for Elastic3D

        Parameters
        ----------
        e_mod: float
            Young's modulus
        area: float
            Cross-sectional area of section
        iz: float
            Second moment of area about the local z-axis
        iy: float
            Second moment of area about the local y-axis (required for 3d analysis)
        g_mod: float
            Shear modulus (optional for 2d analysis, required for 3d analysis)
        jxx: float
            Torsional moment of inertia of section (required for 3d analysis)
        alpha_y: float (default=True)
            Shear shape factor along the local y-axis (optional)
        alpha_z: float (default=True)
            Shear shape factor along the local z-axis (optional)
        """
        self.e_mod = float(e_mod)
        self.area = float(area)
        self.iz = float(iz)
        self.iy = float(iy)
        self.g_mod = float(g_mod)
        self.jxx = float(jxx)
        if alpha_y is None:
            self.alpha_y = None
        else:
            self.alpha_y = float(alpha_y)
        if alpha_z is None:
            self.alpha_z = None
        else:
            self.alpha_z = float(alpha_z)
        osi.n_sect += 1
        self._tag = osi.n_sect
        self._parameters = [self.op_type, self._tag, self.e_mod, self.area, self.iz, self.iy, self.g_mod, self.jxx]
        special_pms = ['alpha_y', 'alpha_z']
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

class Fiber(SectionBase):
    """
    The Fiber Section Class

    This command allows the user to construct a FiberSection object. Each FiberSection object is composed of Fibers,
    with each fiber containing a UniaxialMaterial, an area and a location (y,z). The dofs for 2D section are ``[P,
    Mz]``,for 3D are ``[P,Mz,My,T]``.
    """
    op_type = 'Fiber'

    def __init__(self, osi, gj: float = None, torsion_mat=None):
        """
        Initial method for Fiber

        Parameters
        ----------
        gj: float
            Linear-elastic torsional stiffness assigned to the section
        torsion_mat: obj
            Uniaxialmaterial tag assigned to the section for torsional response (can be nonlinear)
        """
        if gj is None:
            self.gj = None
        else:
            self.gj = float(gj)
        self.torsion_mat = torsion_mat
        osi.n_sect += 1
        self._tag = osi.n_sect
        self._parameters = [self.op_type, self._tag]
        if getattr(self, 'gj') is not None:
            self._parameters += ['-GJ', self.gj]
        if getattr(self, 'torsion_mat') is not None:
            self._parameters += ['-torsion', self.torsion_mat.tag]
        self.to_process(osi)

class Fiber(SectionBase):
    """
    The Fiber Section Class

    This command allows the user to construct a FiberSection object. Each FiberSection object is composed of Fibers,
    with each fiber containing a UniaxialMaterial, an area and a location (y,z). The dofs for 2D section are ``[P,
    Mz]``,for 3D are ``[P,Mz,My,T]``.
    """
    op_type = 'Fiber'

    def __init__(self, osi, gj: float = None, torsion_mat=None):
        """
        Initial method for Fiber

        Parameters
        ----------
        gj: float
            Linear-elastic torsional stiffness assigned to the section
        torsion_mat: obj
            Uniaxialmaterial tag assigned to the section for torsional response (can be nonlinear)
        """
        if gj is None:
            self.gj = None
        else:
            self.gj = float(gj)
        self.torsion_mat = torsion_mat
        osi.n_sect += 1
        self._tag = osi.n_sect
        self._parameters = [self.op_type, self._tag]
        if getattr(self, 'gj') is not None:
            self._parameters += ['-GJ', self.gj]
        if getattr(self, 'torsion_mat') is not None:
            self._parameters += ['-torsion', self.torsion_mat.tag]
        self.to_process(osi)

class Fiber(SectionBase):
    """
    The Fiber Section Class

    This command allows the user to construct a FiberSection object. Each FiberSection object is composed of Fibers,
    with each fiber containing a UniaxialMaterial, an area and a location (y,z). The dofs for 2D section are ``[P,
    Mz]``,for 3D are ``[P,Mz,My,T]``.
    """
    op_type = 'Fiber'

    def __init__(self, osi, gj: float = None, torsion_mat=None):
        """
        Initial method for Fiber

        Parameters
        ----------
        gj: float
            Linear-elastic torsional stiffness assigned to the section
        torsion_mat: obj
            Uniaxialmaterial tag assigned to the section for torsional response (can be nonlinear)
        """
        if gj is None:
            self.gj = None
        else:
            self.gj = float(gj)
        self.torsion_mat = torsion_mat
        osi.n_sect += 1
        self._tag = osi.n_sect
        self._parameters = [self.op_type, self._tag]
        if getattr(self, 'gj') is not None:
            self._parameters += ['-GJ', self.gj]
        if getattr(self, 'torsion_mat') is not None:
            self._parameters += ['-torsion', self.torsion_mat.tag]
        self.to_process(osi)



class FiberThermal(SectionBase):
    """
    The FiberThermal Section Class
    
    This command create a FiberSectionThermal object.The dofs for 2D section are ``[P, Mz]``,for 3D are ``[P,Mz,My]``...
    note::#. The commands below should be called after the section command to generate all the fibers in the section.#. The
    patch and layer commands can be used to generate multiple fibers in a single command... toctree:::maxdepth: 2:caption:
    Commands to generate all fibersfiberpatchlayer
    """
    op_type = 'FiberThermal'

    def __init__(self, osi, gj=None):
        """
        Initial method for FiberThermal

        Parameters
        ----------
        gj: None
            
        """
        self.gj = gj
        osi.n_sect += 1
        self._tag = osi.n_sect
        self._parameters = [self.op_type, self._tag]
        if getattr(self, 'gj') is not None:
            self._parameters += ['-GJ', self.gj]
        self.to_process(osi)
