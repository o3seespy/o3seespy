import o3seespy as o3


def test_pressure_depend_multi_yield():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.nd_material.PressureDependMultiYield(osi, nd=2, rho=1.7, g_mod_ref=6.0e4, bulk_mod_ref=16.0e4, phi=31., peak_strain=0.1, p_ref=101,
                 d=0.5, pt_ang=31., con_rate=0.21, dil_rates=[0., 0.], liquefac=[10, 0.02, 1], n_surf=20.0,
                 e_init=0.85, cs_params=None, c=0.3)


def test_pressure_depend_multi_yield02():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.nd_material.PressureDependMultiYield02(osi, nd=2, rho=1.7, g_mod_ref=6.0e4, bulk_mod_ref=16.0e4, phi=31.,
                                            peak_strain=0.1, p_ref=101,
                                            d=0.5, pt_ang=31., con_rates=[0.08, 5, 0.18], dil_rates=[0., 3., 0.],
                                            liquefac=[1, 0.], n_surf=20.0,
                                            e_init=0.85, cs_params=None, c=0.3)