import openseespy.opensees as opy
import o3seespy as o3


def test_uniaxial_steel01():
    osi = o3.OpenseesInstance(ndm=2)

    # Define material
    bilinear_mat = o3.uniaxial_material.Steel01(osi, fy=1, e0=1, b=1)


# def test_bond_sp01():
#     osi = o3.OpenseesInstance(ndm=2)
#     o3.uniaxial_material.BondSP01(osi, fy=1, sy=1.0, fu=1.0, su=1.0, b=1.0, big_r=1.0)
#
#
# def test_bouc_wen():
#     osi = o3.OpenseesInstance(ndm=2)
#     o3.uniaxial_material.BoucWen(osi, alpha=1.0, ko=1.0, n=1.0, gamma=1.0, beta=1.0, ao=1.0, delta_a=1.0, delta_nu=1.0, delta_eta=1.0)
#
#
# def test_bilinear_oil_damper():
#     osi = o3.OpenseesInstance(ndm=2)
#     o3.uniaxial_material.BilinearOilDamper(osi, big_k=1.0, cd=1.0, fr=1.0, p=1.0, l_gap=0.0, nm=1, rel_tol=1e-6,
#                                             abs_tol=1e-10, max_half=15)
#
#
# def test_init_strain_material():
#     osi = o3.OpenseesInstance(ndm=2)
#     bilinear_mat = o3.uniaxial_material.Steel01(osi, fy=1, e0=1, b=1)
#     o3.uniaxial_material.InitStrainMaterial(osi, other=bilinear_mat, init_strain=1.0)
