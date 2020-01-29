from o3seespy.base_model import OpenseesObject
from o3seespy.command.section import SectionBase


class LayeredShell(SectionBase):
    """
    The LayeredShell Section Class

    This command will create the section of the multi-layer shell element, including the multi-dimensional concrete,
    reinforcement material and the corresponding thickness.
    """
    op_type = 'LayeredShell'

    def __init__(self, osi, mats):
        """
        Initial method for LayeredShell

        Parameters
        ----------
        mats: list
            A list of material objs and thickness, ``[[mat1,thk1], ..., [mat2,thk2]]``
        """
        self.mats = []
        for i, mat in enumerate(mats):
            # self.mats.append([mats[i][0].tag, mats[i][1]])
            self.mats.append(mats[i][0].tag)
            self.mats.append(mats[i][1])

        self.n_layers = int(len(self.mats))
        osi.n_sect += 1
        self._tag = osi.n_sect
        self._parameters = [self.op_type, self._tag, self.n_layers, *self.mats]
        self.to_process(osi)


class Aggregator(SectionBase):
    """
    The Aggregator Section Class

    This command is used to construct a SectionAggregator object which aggregates groups previously-defined
    UniaxialMaterial objects into a single section force-deformation model. Each UniaxialMaterial object
    represents the section force-deformation response for a particular section degree-of-freedom (dof).
    There is no interaction between responses in different dof directions. The aggregation can include
    one previously defined section.
    """
    op_type = 'Aggregator'

    def __init__(self, osi, mats, section=None):
        """
        Initial method for Aggregator

        Parameters
        ----------
        mats: list
            List of mat objs and dofs of previously-defined uniaxialmaterial objects, ``mats =
            [[mattag1,dof1],[mattag2,dof2],...]`` the force-deformation quantity to be modeled by this
            section object. one of the following section dof may be used: * ``'p'`` axial
            force-deformation * ``'mz'`` moment-curvature about section local z-axis *
            ``'vy'`` shear force-deformation along section local y-axis * ``'my'``
            moment-curvature about section local y-axis * ``'vz'`` shear
            force-deformation along section local z-axis * ``'t'`` torsion force-deformation
        section: obj
            Tag of previously-defined section object to which the uniaxialmaterial objects are aggregated as additional
            force-deformation relationships (optional)
        """
        self.mats = []
        for i, mat in enumerate(mats):
            self.mats.append(mats[i][0].tag)
            self.mats.append(mats[i][1])
        self.section = section
        osi.n_sect += 1
        self._tag = osi.n_sect
        self._parameters = [self.op_type, self._tag, *self.mats]
        if getattr(self, 'section') is not None:
            self._parameters += ['-section', self.section.tag]
        self.to_process(osi)


from o3seespy.command.nd_material.base_material import NDMaterialBase
from o3seespy.command.uniaxial_material.base_material import UniaxialMaterialBase


class PressureIndependMultiYield(NDMaterialBase):
    op_type = "PressureIndependMultiYield"

    def __init__(self, osi,  nd, rho, ref_shear_modul, ref_bulk_modul, cohesi, peak_shear_stra, friction_ang=0.,
                 ref_press=100., press_depend_coe=0., no_yield_surf=20, strains=None, ratios=None):
        """
        PressureIndependMultiYield material
        """
        self.nd = nd
        self.rho = float(rho)
        self.ref_shear_modul = float(ref_shear_modul)
        self.ref_bulk_modul = float(ref_bulk_modul)
        self.cohesi = float(cohesi)
        self.peak_shear_stra = float(peak_shear_stra)
        self.friction_ang = float(friction_ang)
        self.ref_press = float(ref_press)
        self.press_depend_coe = float(press_depend_coe)
        assert no_yield_surf < 40
        self.no_yield_surf = int(no_yield_surf)
        if strains is not None:
            assert len(strains) == len(ratios)
            yield_surf = []
            for i in range(len(strains)):
                yield_surf.append(ratios[i])
                yield_surf.append(strains[i])
            self.yield_surf = yield_surf
        else:
            self.yield_surf = None

        osi.n_mat += 1
        self._tag = osi.n_mat

        self._parameters = [self.op_type, self._tag, self.nd, self.rho, self.ref_shear_modul, self.ref_bulk_modul,
                            self.cohesi, self.peak_shear_stra, self.friction_ang, self.ref_press, self.press_depend_coe]

        if self.yield_surf is not None:
            self._parameters.append(self.no_yield_surf)  # from docs 'add a minus sign in front of noYieldSurf'
            self._parameters += list(self.yield_surf)

        else:
            # self._keyword_args['noYieldSurf'] = self.no_yield_surf
            self._parameters.append(self.no_yield_surf)
        self.to_process(osi)


class Steel01(UniaxialMaterialBase):
    """
    The Steel01 UniaxialMaterial Class

    This command is used to construct a uniaxial bilinear steel material object with kinematic hardening and optional
    isotropic hardening described by a non-linear evolution equation (REF: Fedeas).
    """
    op_type = "Steel01"

    def __init__(self, osi, fy: float, e0: float, b: float, a1=None, a2=None, a3=None, a4=None):
        """
        Initial method for Steel01

        Parameters
        ----------
        fy: float
            Yield strength
        e0: float
            Initial elastic tangent
        b: float
            Strain-hardening ratio (ratio between post-yield tangent and initial elastic tangent)
        a1: float
            Isotropic hardening parameter, increase of compression yield envelope as proportion of yield strength after
            a plastic strain of :math:`a_2*(f_y/e_0)` (optional)
        a2: float
            Isotropic hardening parameter
        a3: float
            Isotropic hardening parameter, increase of tension yield envelope as proportion of yield strength after a
            plastic strain of :math:`a_4*(f_y/e_0)`. (optional)
        a4: float
            Isotropic hardening parameter (see explanation
        """
        self.fy = float(fy)
        self.e0 = float(e0)
        self.b = float(b)
        self.a_values = [a1, a2, a3, a4]
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.a4 = a4
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self.tag, self.fy, self.e0, self.b]
        for a in self.a_values:
            if a is None:
                break
            self._parameters.append(a)
        self.to_process(osi)


#
# class PySimple1(UniaxialMaterialBase):
#     op_type = "PySimple1"
#
#     def __init__(self, osi, soil_type: int, p_ult, y50, cd, c):
#         """
#         PySimple1 uniaxial material object
#
#         Parameters
#         ----------
#         osi : opensees_pack.opensees_instance.OpenseesInstance object
#             An instance of opensees
#         soil_type : int {1, 2}
#             Backbone type for soil
#         p_ult : float
#             Ultimate capacity of the p-y material. Note that “p” or “pult” are distributed loads
#              [force per length of pile] in common design equations, but are both loads for this
#              uniaxialMaterial [i.e., distributed load times the tributary length of the pile].
#         y50 : float
#             Displacement at which 50% of pult is mobilized in monotonic loading.
#         cd : float
#             Variable that sets the drag resistance within a fully-mobilized gap as Cd*pult.
#         c : float
#             The viscous damping term (dashpot) on the far-field (elastic) component
#             of the displacement rate (velocity). (optional Default = 0.0). Nonzero c
#             values are used to represent radiation damping effects
#         """
#         self.soil_type = int(soil_type)
#         self.p_ult = p_ult
#         self.y50 = y50
#         self.cd = cd
#         self.c = c
#         osi.n_mat += 1
#         self._tag = osi.n_mat
#         self._parameters = [self.op_type, self._tag, self.soil_type, self.p_ult, self.y50, self.cd, self.c]
#         self.to_process(osi)
#
#
# class PyLiq1(UniaxialMaterialBase):
#     op_type = "PyLiq1"
#
#     def __init__(self, osi, soil_type, p_ult, y50, cd, c, p_res, ele1=None, ele2=None, time_series=None):
#         """
#
#         Parameters
#         ----------
#         osi : opensees_pack.opensees_instance.OpenseesInstance object
#             An instance of opensees
#         soil_type
#        p_ult : float
#             Ultimate capacity of the p-y material. Note that “p” or “pult” are distributed loads
#              [force per length of pile] in common design equations, but are both loads for this
#              uniaxialMaterial [i.e., distributed load times the tributary length of the pile].
#         y50 : float
#             Displacement at which 50% of pult is mobilized in monotonic loading.
#         cd : float
#             Variable that sets the drag resistance within a fully-mobilized gap as Cd*pult.
#         c : float
#             The viscous damping term (dashpot) on the far-field (elastic) component
#             of the displacement rate (velocity). (optional Default = 0.0). Nonzero c
#             values are used to represent radiation damping effects
#         p_res : float
#         ele1 : opensees_pack.ElementBase object, optional
#             the eleTag (element numbers) for the two solid element from which PyLiq1
#             will obtain mean effective stresses and excess pore pressures
#         ele2 : opensees_pack.ElementBase object, optional
#             Same as ele1
#         time_series : iterable object, optional
#             A time series of mean effective stress values
#         """
#         self.soil_type = soil_type
#         self.p_ult = p_ult
#         self.y50 = y50
#         self.cd = cd
#         self.c = c
#         self.p_res = p_res
#         self.ele1 = ele1
#         self.ele2 = ele2
#         self.time_series = time_series
#         osi.n_mat += 1
#         self._tag = osi.n_mat
#         self._parameters = [self.op_type, self._tag, self.p_ult, self.y50, self.cd, self.c, self.p_res]
#         if self.ele1 is None:
#             self._parameters.append("-timeSeries")
#             self._parameters.append(self.time_series)
#         else:
#             self._parameters.append(self.ele1)
#             if self.ele2 is None:
#                 self._parameters.append(self.ele2)
#         self.to_process(osi)
