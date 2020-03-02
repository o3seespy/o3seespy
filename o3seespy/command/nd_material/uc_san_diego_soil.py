from o3seespy.command.nd_material.base_material import NDMaterialBase


class PressureIndependMultiYield(NDMaterialBase):
    op_type = "PressureIndependMultiYield"

    def __init__(self, osi,  nd, rho, g_mod_ref, bulk_mod_ref, cohesion, peak_strain, phi=0.,
                 p_ref=100., m=0., n_surf=20, strains=None, ratios=None):
        """
        PressureIndependMultiYield material

        Parameters
        ==========
        matTag: int
            integer tag identifying material
        nd: float
            Number of dimensions, 2 for plane-strain, and 3 for 3D analysis.
       rho: float
            Saturated soil mass density.
       g_mod_ref, float
            (:math:`G_0`) Reference low-strain shear modulus,
            specified at a reference mean effective confining pressure (`p_ref`).
       bulk_mod_ref: float
            (:math:`B_r`) Reference bulk modulus,
            specified at a reference mean effective confining pressure (`p_ref`)).
       cohesion: float
            (:math:`c`) Apparent cohesion at zero effective confinement.
       peak_strain: float
            (:math:`\\gamma_{max}`) An octahedral shear strain at
            which the maximum shear strength is reached,
            specified at a reference mean effective confining
            pressure refPress of p'r (see below).
       phi: float
            (:math:`phi`) Friction angle at peak shear
            strength in degrees, optional (default is 0.0).
       p_ref: float
            (:math:`p'_ref`) Reference mean effective confining pressure at which
                          :math:`G_r`, :math:`B_r`, and :math:`\gamma_{max}`
                          are defined, optional (default is 100. kPa).
       m: float
            (:math:`d`) A positive constant defining variations
                        of :math:`G` and :math:`B` as a function of
                          instantaneous effective
                          confinement :math:`p'` (default is 0.0)

                          :math:`G=G_r(\frac{p'}{p'_ref})^d`

                          :math:`B=B_r(\frac{p'}{p'_ref})^d`

                          If :math:`\phi=0`, :math:`d` is reset to 0.0.

       n_surf: float, optional
            Number of yield surfaces, optional (must be less
            than 40, default is 20). The surfaces are generated
            based on the hyperbolic relation.
       strains: list
            Instead of automatic surfaces generation,
            you can define yield surfaces directly based on
            desired shear modulus reduction curve.
            provide a list of shear strains and corresponding shear modulus reduction ratios (`ratios`)
        ratios: list
            Instead of automatic surfaces generation,
            you can define yield surfaces directly based on
            desired shear modulus reduction curve.
            provide a list of shear strains and corresponding shear modulus reduction ratios (`ratios`)
        """
        self.nd = nd
        self.rho = float(rho)
        self.g_mod_ref = float(g_mod_ref)
        self.bulk_mod_ref = float(bulk_mod_ref)
        self.cohesion = float(cohesion)
        self.peak_strain = float(peak_strain)
        self.phi = float(phi)
        self.p_ref = float(p_ref)
        self.m = float(m)
        assert n_surf < 40
        self.n_surf = int(n_surf)
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

        self._parameters = [self.op_type, self._tag, self.nd, self.rho, self.g_mod_ref, self.bulk_mod_ref,
                            self.cohesion, self.peak_strain, self.phi, self.p_ref, self.m]

        if self.yield_surf is not None:
            self._parameters.append(self.n_surf)  # from docs 'add a minus sign in front of noYieldSurf'
            self._parameters += list(self.yield_surf)

        else:
            # self._keyword_args['noYieldSurf'] = self.no_yield_surf
            self._parameters.append(self.n_surf)
        self.to_process(osi)

