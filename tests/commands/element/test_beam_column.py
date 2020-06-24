import o3seespy as o3  # for testing only
import pytest


def test_elastic_beam_column2d():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
    transf = o3.geom_transf.Linear2D(osi, [])
    o3.element.ElasticBeamColumn2D(osi, ele_nodes=ele_nodes, area=1.0, e_mod=1.0, iz=1.0, transf=transf, mass=1.0)

def test_elastic_beam_column3d():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
    transf = o3.geom_transf.Linear2D(osi, [])
    o3.element.ElasticBeamColumn3D(osi, ele_nodes=ele_nodes, area=1.0, e_mod=1.0, g_mod=1.0, jxx=1.0, iy=1.0, iz=1.0, transf=transf, mass=1.0)


def test_mod_elastic_beam2dmass():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
    transf = o3.geom_transf.Linear2D(osi, [])
    o3.element.ModElasticBeam2D(osi, ele_nodes=ele_nodes, area=1.0, e_mod=1.0, iz=1.0, k11=1.0, k33=1.0, k44=1.0, transf=transf, mass=1.0)


def test_elastic_timoshenko_beam2d():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
    transf = o3.geom_transf.Linear2D(osi, [])
    o3.element.ElasticTimoshenkoBeam2D(osi, ele_nodes=ele_nodes, e_mod=1.0, g_mod=1.0, area=1.0, iz=1.0, avy=1.0, transf=transf, mass=1.0)


def test_elastic_timoshenko_beam3d():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
    transf = o3.geom_transf.Linear2D(osi, [])
    o3.element.ElasticTimoshenkoBeam3D(osi, ele_nodes=ele_nodes, e_mod=1.0, g_mod=1.0, area=1.0, iz=1.0, jxx=1.0, iy=1.0, iz_2=1, avy=1.0, avz=1.0, transf=transf, mass=1.0, c_mass=True)


def test_disp_beam_column():
    osi = o3.OpenSeesInstance(ndm=2)
    i_node = o3.node.Node(osi, 0.0, 0.0)
    j_node = o3.node.Node(osi, 0.0, 1.0)
    transf = o3.geom_transf.Linear2D(osi, [])
    sec = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
    integration = o3.beam_integration.Lobatto(osi, sec, 5)
    o3.element.DispBeamColumn(osi, ele_nodes=[i_node, j_node], transf=transf, integration=integration, c_mass=1, mass=0.0)


def test_force_beam_column():
    osi = o3.OpenSeesInstance(ndm=2)
    i_node = o3.node.Node(osi, 0.0, 0.0)
    j_node = o3.node.Node(osi, 0.0, 1.0)
    transf = o3.geom_transf.Linear2D(osi, [])
    sec = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
    integration = o3.beam_integration.Lobatto(osi, sec, 5)
    o3.element.ForceBeamColumn(osi, ele_nodes=[i_node, j_node], transf=transf, integration=integration, max_iter=10, tol=1e-12, mass=0.0)


@pytest.mark.skip()  # can't find int_type that works with this!
def test_nonlinear_beam_column():
    osi = o3.OpenSeesInstance(ndm=2)
    i_node = o3.node.Node(osi, 0.0, 0.0)
    j_node = o3.node.Node(osi, 0.0, 1.0)
    sec = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
    transf = o3.geom_transf.Linear2D(osi, [])
    o3.element.NonlinearBeamColumn(osi, ele_nodes=[i_node, j_node], num_intgr_pts=1, sec=sec, transf=transf,
                                   max_iter=10, tol=1e-12, mass=0.0, int_type="radau")


@pytest.mark.skip()
def test_disp_beam_column_int():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(2)]
    sec = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
    transf = o3.geom_transf.Linear2D(osi, [])
    o3.element.DispBeamColumnInt(osi, ele_nodes=ele_nodes, num_intgr_pts=4, sec=sec, transf=transf, c_rot=1.0, mass=1.0)


@pytest.mark.skip()
def test_mvlem():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
    mat_conc = [o3.uniaxial_material.Concrete01(osi, fpc=1.0, epsc0=1.0, fpcu=1.0, eps_u=1.0),
                o3.uniaxial_material.Concrete01(osi, fpc=1.0, epsc0=1.0, fpcu=1.0, eps_u=1.0)]
    mat_steel = [o3.uniaxial_material.Steel02(osi, fy=1.0, e0=1.0, b=1.0, params=[15, 0.925, 0.15])]
    mat_shear = o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)
    o3.element.MVLEM(osi, dens=1.0, ele_nodes=ele_nodes, m=1, c=1.0, thick=[1.0, 1.0], widths=[1, 1], rho=[1., 1.],
                     mat_concretes=mat_conc, mat_steels=mat_steel, mat_shear=mat_shear)


@pytest.mark.skip()
def test_sfimvlem():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]
    mats = [o3.uniaxial_material.Elastic(osi, 1.0, 1.0), o3.uniaxial_material.Elastic(osi, 1.0, 1.0)]
    mat_tags = [x.tag for x in mats]  # TODO: should pass in mats not mat tags
    o3.element.SFIMVLEM(osi, ele_nodes=ele_nodes, m=1, c=1.0, thick=[1.0, 1.0], widths=[1.0, 1.0], mat_tags=mat_tags)

