import o3seespy as o3  # for testing only
import pytest


def test_elastomeric_bearing_plasticity2d():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [0, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    p_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    mz_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    o3.element.ElastomericBearingPlasticity2D(osi, ele_nodes=ele_nodes, k_init=1.0, qd=1.0, alpha1=1.0, alpha2=1.0,
                                              mu=1.0, p_mat=p_mat, mz_mat=mz_mat)


def test_elastomeric_bearing_plasticity3d():
    osi = o3.OpenSeesInstance(ndm=3, ndf=6)
    coords = [[0, 0, 0], [0, 1, 0]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    p_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    mz_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    t_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    my_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    orient_vals = [1, 0, 0]
    o3.element.ElastomericBearingPlasticity3D(osi, ele_nodes=ele_nodes, k_init=1.0, qd=1.0, alpha1=1.0, alpha2=1.0,
                                              mu=1.0, p_mat=p_mat, t_mat=t_mat, my_mat=my_mat, mz_mat=mz_mat,
                                              orient=orient_vals)


def test_elastomeric_bearing_bouc_wen2d():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    orient_vals = [1, 0, 0, 1, 0, 1]
    p_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    mz_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    o3.element.ElastomericBearingBoucWen2D(osi, ele_nodes=ele_nodes, k_init=1.0, qd=1.0, alpha1=1.0, alpha2=1.0,
                                           mu=1.0, eta=1.0, beta=1.0, gamma=1.0, p_mat=p_mat, mz_mat=mz_mat,
                                           orient_vals=orient_vals, shear_dist=1.0, do_rayleigh=False, mass=1.0)


def test_elastomeric_bearing_bouc_wen3d():
    osi = o3.OpenSeesInstance(ndm=3, ndf=6)
    coords = [[0, 0, 0], [0, 1, 0]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    orient_vals = [1, 0, 0]
    p_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    mz_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    t_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    my_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    o3.element.ElastomericBearingBoucWen3D(osi, ele_nodes=ele_nodes, k_init=1.0, qd=1.0, alpha1=1.0, alpha2=1.0,
                                           mu=1.0, eta=1.0, beta=1.0, gamma=1.0, p_mat=p_mat, t_mat=t_mat,
                                           my_mat=my_mat, mz_mat=mz_mat, orient_vals=orient_vals,
                                           shear_dist=1.0, do_rayleigh=False, mass=1.0)


def test_flat_slider_bearing2d():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    p_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    mz_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    frn1 = o3.friction_model.Coulomb(osi, mu=1.0)
    o3.element.FlatSliderBearing2D(osi, ele_nodes=ele_nodes, frn_mdl=frn1, k_init=1.0, p_mat=p_mat, mz_mat=mz_mat,
                                   do_rayleigh=False, max_iter=1, tol=1.0, orient=None, mass=1.0, shear_dist=1.0)


def test_flat_slider_bearing3d():
    osi = o3.OpenSeesInstance(ndm=3, ndf=6)
    coords = [[0, 0, 0], [0, 1, 0]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    p_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    mz_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    t_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    my_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    frn1 = o3.friction_model.Coulomb(osi, mu=1.0)
    orient_vals = [1, 0, 0]
    o3.element.FlatSliderBearing3D(osi, ele_nodes=ele_nodes, frn_mdl=frn1, k_init=1.0, p_mat=p_mat, t_mat=t_mat,
                                   my_mat=my_mat, mz_mat=mz_mat, do_rayleigh=False, max_iter=None, tol=None,
                                   mass=1.0, shear_dist=1.0, orient=orient_vals)


def test_single_fp_bearing2d():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    p_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    mz_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    frn1 = o3.friction_model.Coulomb(osi, mu=1.0)
    o3.element.SingleFPBearing2D(osi, ele_nodes=ele_nodes, frn_mdl=frn1, reff=1.0, k_init=1.0, p_mat=p_mat,
                                 mz_mat=mz_mat, do_rayleigh=False, max_iter=1, tol=1.0, orient=None,
                                 mass=1.0, shear_dist=1.0)


def test_single_fp_bearing3d():
    osi = o3.OpenSeesInstance(ndm=3, ndf=6)
    coords = [[0, 0, 0], [0, 1, 0]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    p_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    mz_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    t_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    my_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    frn1 = o3.friction_model.Coulomb(osi, mu=1.0)
    orient_vals = [1, 0, 0]
    o3.element.SingleFPBearing3D(osi, ele_nodes=ele_nodes, frn_mdl=frn1, reff=1.0, k_init=1.0, p_mat=p_mat, t_mat=t_mat,
                                 my_mat=my_mat, mz_mat=mz_mat, do_rayleigh=False, max_iter=None, tol=None,
                                 orient=orient_vals, mass=1.0, shear_dist=1.0)


def test_tfp():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    o3.element.TFP(osi, ele_nodes=ele_nodes,
                   r1=1.0, r2=1.0, r3=1.0, r4=1.0,
                   db1=1.0, db2=1.0, db3=1.0, db4=1.0,
                   d1=1.0, d2=1.0, d3=1.0, d4=1.0,
                   mu1=0.3, mu2=0.4, mu3=0.5, mu4=0.5,
                   h1=1.0, h2=1.0, h3=1.0, h4=1.0,
                   h0=1.0, col_load=1.0, big_k=None)


@pytest.mark.skip()
def test_triple_friction_pendulum():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    frn1 = o3.friction_model.Coulomb(osi, mu=1.0)
    frn2 = o3.friction_model.Coulomb(osi, mu=1.0)
    frn3 = o3.friction_model.Coulomb(osi, mu=1.0)
    vert_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    rot_z_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    rot_x_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    rot_y_mat = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    o3.element.TripleFrictionPendulum(osi, ele_nodes=ele_nodes, frn1=frn1, frn2=frn2, frn3=frn3, vert_mat=vert_mat,
                                      rot_z_mat=rot_z_mat, rot_x_mat=rot_x_mat, rot_y_mat=rot_y_mat, l1=1.0, l2=1.0,
                                      l3=1.0, d1=1.0, d2=1.0, d3=1.0, big_w=1.0, uy=1.0, kvt=1.0, min_fv=None, tol=1.0)


def test_multiple_shear_spring():
    osi = o3.OpenSeesInstance(ndm=3, ndf=6)
    coords = [[0, 0, 0], [1, 0, 0]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    mat = o3.uniaxial_material.Elastic(osi, 1.0)
    o3.element.MultipleShearSpring(osi, ele_nodes=ele_nodes, n_spring=1, mat=mat, lim=1.0, mass=1.0, orient=None)


@pytest.mark.skip()
def test_kikuchi_bearingadjust_pd_output():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    mat_mss = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    mat_mns = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    o3.element.KikuchiBearingadjustPDOutput(osi, ele_nodes=ele_nodes, shape=1.0, size=1.0, total_rubber=1.0, total_height=1.0, n_mss=1, mat_mss=mat_mss, lim_disp=1.0, n_mns=1, mat_mns=mat_mns, lamb=1.0, no_pd_input="string", no_tilt="string", ci=1.0, cj=1.0, orient=[0.0, 0.0], mass=1.0)


@pytest.mark.skip()
def test_kikuchi_bearingdo_balance():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    mat_mss = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    mat_mns = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    o3.element.KikuchiBearingdoBalance(osi, ele_nodes=ele_nodes, shape=1.0, size=1.0, total_rubber=1.0, total_height=1.0, n_mss=1, mat_mss=mat_mss, lim_disp=1.0, n_mns=1, mat_mns=mat_mns, lamb=1.0, no_pd_input="string", no_tilt="string", lim_fo=1.0, lim_fi=1.0, n_iter=1.0, orient=[0.0, 0.0], mass=1.0)


@pytest.mark.skip()
def test_yamamoto_biaxial_hd_rco_rs():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    o3.element.YamamotoBiaxialHDRcoRS(osi, ele_nodes=ele_nodes, tp=1, d_do=1.0, d_di=1.0, hr=1.0, cr=1.0, cs=1.0, orient=[0.0, 0.0], mass=1.0)


@pytest.mark.skip()
def test_elastomeric_x():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    o3.element.ElastomericX(osi, ele_nodes=ele_nodes, fy=1.0, alpha=1.0, gr=1.0, kbulk=1.0, d1=1.0, d2=1.0, ts=1.0,
                            tr=1.0, n=1, x1=1.0, x2=1.0, x3=1.0, y1=1.0, y2=1.0, y3=1.0, kc=1.0, phi_m=1.0, ac=1.0,
                            s_dratio=1.0, m=1.0, cd=1.0, tc=1.0, tag1=1.0, tag2=1.0, tag3=1.0, tag4=1.0)


def test_rj_watson_eqs_bearing2d():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [0, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    p_mat = o3.uniaxial_material.Elastic(osi, 1, 1)
    vy_mat = o3.uniaxial_material.Elastic(osi, 1, 1)
    mz_mat = o3.uniaxial_material.Elastic(osi, 1, 1)
    frn_mdl = o3.friction_model.Coulomb(osi, mu=1.0)
    o3.element.RJWatsonEqsBearing2D(osi, ele_nodes=ele_nodes, frn_mdl=frn_mdl, k_init=1.0, p_mat=p_mat, vy_mat=vy_mat,
                                    mz_mat=mz_mat, do_rayleigh=False, max_iter=1, tol=1.0, orient=None, mass=1.0,
                                    shear_dist=1.0)


def test_rj_watson_eqs_bearing3d():
    osi = o3.OpenSeesInstance(ndm=3, ndf=6)
    coords = [[0, 0, 0], [0, 1, 0]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(2)]
    p_mat = o3.uniaxial_material.Elastic(osi, 1, 1)
    vy_mat = o3.uniaxial_material.Elastic(osi, 1, 1)
    vz_mat = o3.uniaxial_material.Elastic(osi, 1, 1)
    t_mat = o3.uniaxial_material.Elastic(osi, 1, 1)
    my_mat = o3.uniaxial_material.Elastic(osi, 1, 1)
    mz_mat = o3.uniaxial_material.Elastic(osi, 1, 1)
    orient_vals = [1, 0, 0]
    frn_mdl = o3.friction_model.Coulomb(osi, mu=1.0)
    o3.element.RJWatsonEqsBearing3D(osi, ele_nodes=ele_nodes, frn_mdl=frn_mdl, k_init=1.0, p_mat=p_mat,
                                    vy_mat=vy_mat, vz_mat=vz_mat, t_mat=t_mat, my_mat=my_mat, mz_mat=mz_mat,
                                    do_rayleigh=False, max_iter=1, tol=1.0, orient=orient_vals, mass=1.0, shear_dist=1.0)


@pytest.mark.skip()
def test_lead_rubber_x():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    o3.element.LeadRubberX(osi, ele_nodes=ele_nodes, fy=1.0, alpha=1.0, gr=1.0, kbulk=1.0, d1=1.0, d2=1.0, ts=1.0,
                           tr=1.0, n=1, x1=1.0, x2=1.0, x3=1.0, y1=1.0, y2=1.0, y3=1.0, kc=1.0, phi_m=1.0, ac=1.0,
                           s_dratio=1.0, m=1.0, cd=1.0, tc=1.0, q_l=1.0, c_l=1.0, k_s=1.0, a_s=1.0,
                           tag1=1, tag2=1, tag3=1, tag4=1, tag5=1)


@pytest.mark.skip()
def test_hdr():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    o3.element.HDR(osi, ele_nodes=ele_nodes, gr=1.0, kbulk=1.0, d1=1.0, d2=1.0, ts=1.0, tr=1.0, n=1, a1=1.0, a2=1.0, a3=1.0, b1=1.0, b2=1.0, b3=1.0, c1=1.0, c2=1.0, c3=1.0, c4=1.0, x1=1.0, x2=1.0, x3=1.0, y1=1.0, y2=1.0, y3=1.0, kc=1.0, phi_m=1.0, ac=1.0, s_dratio=1.0, m=1.0, tc=1.0)


@pytest.mark.skip()
def test_fp_bearing_ptv():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]
    o3.element.FPBearingPTV(osi, ele_nodes=ele_nodes, mu_ref=1.0, is_pressure_dependent=1, p_ref=1.0, is_temperature_dependent=1, diffusivity=1.0, conductivity=1.0, is_velocity_dependent=1, rate_parameter=1.0, reffective_fp=1.0, radius__contact=1.0, k_initial=1.0, the_material_a=1, the_material_b=1, the_material_c=1, the_material_d=1, x1=1.0, x2=1.0, x3=1.0, y1=1.0, y2=1.0, y3=1.0, shear_dist=1.0, do_rayleigh=1, mass=1.0, max_iter=1, tol=1.0, unit=1)


if __name__ == '__main__':
    test_elastomeric_bearing_bouc_wen2d()
