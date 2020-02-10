import o3seespy as o3  # for testing only


def test_quad():
    osi = o3.OpenSeesInstance(ndm=2)
    crds_i = [1.0, 1.0]
    crds_j = [1.0, 1.0]
    crds_k = [1.0, 1.0]
    crds_l = [1.0, 1.0]
    conc_conf = o3.uniaxial_material.Concrete01(osi, fpc=-5.0, epsc0=-0.005, fpcu=-3.5, eps_u=-0.02)
    o3.patch.Quad(osi, conc_conf, num_subdiv_ij=1, num_subdiv_jk=1, crds_i=crds_i, crds_j=crds_j, crds_k=crds_k, crds_l=crds_l)


def test_rect():
    osi = o3.OpenSeesInstance(ndm=2)
    crds_i = [1.0, 1.0]
    crds_j = [1.0, 1.0]
    conc_conf = o3.uniaxial_material.Concrete01(osi, fpc=-5.0, epsc0=-0.005, fpcu=-3.5, eps_u=-0.02)
    o3.patch.Rect(osi, conc_conf, num_subdiv_y=1, num_subdiv_z=1, crds_i=crds_i, crds_j=crds_j)


def test_circ():
    osi = o3.OpenSeesInstance(ndm=2)
    center = [1.0, 1.0]
    rad = [1.0, 1.0]
    ang = [1.0, 1.0]
    conc_conf = o3.uniaxial_material.Concrete01(osi, fpc=-5.0, epsc0=-0.005, fpcu=-3.5, eps_u=-0.02)
    o3.patch.Circ(osi, conc_conf, num_subdiv_circ=1, num_subdiv_rad=1, center=center, rad=rad, ang=ang)

