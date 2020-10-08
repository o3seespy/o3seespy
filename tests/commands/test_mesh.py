import o3seespy as o3  # for testing only
import pytest


@pytest.mark.skip()
def test_line():
    osi = o3.OpenSeesInstance(ndm=2)
    coords = [[0, 0], [0, 0]]
    ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(2)]
    o3.mesh.Line(osi, numnodes=1, ndtags=ele_nodes, id=1, ndf=1, meshsize=1.0, ele_type='', ele_args=None)


@pytest.mark.skip()
def test_tri():
    osi = o3.OpenSeesInstance(ndm=2)
    ltags = [1, 1]
    o3.mesh.Tri(osi, numlines=1, ltags=ltags, id=1, ndf=1, meshsize=1.0, ele_type='', ele_args=None)


@pytest.mark.skip()
def test_quad():
    osi = o3.OpenSeesInstance(ndm=2)
    ltags = [1, 1]
    o3.mesh.Quad(osi, numlines=1, ltags=ltags, id=1, ndf=1, meshsize=1.0, ele_type='', ele_args=None)


@pytest.mark.skip()
def test_tet():
    osi = o3.OpenSeesInstance(ndm=2)
    mtags = [1, 1]
    o3.mesh.Tet(osi, nummesh=1, mtags=mtags, id=1, ndf=1, meshsize=1.0, ele_type='', ele_args=None)


@pytest.mark.skip()
def test_partvel():
    osi = o3.OpenSeesInstance(ndm=2)
    p_args = [1.0, 1.0]
    vel0 = [1.0, 1.0]
    o3.mesh.Partvel(osi, otype="string", p_args=p_args, ele_type='', ele_args=None, vel0=vel0)


@pytest.mark.skip()
def test_partpressure():
    osi = o3.OpenSeesInstance(ndm=2)
    p_args = [1.0, 1.0]
    o3.mesh.Partpressure(osi, otype="string", p_args=p_args, ele_type='', ele_args=None, p0=1.0)


@pytest.mark.skip()
def test_bgwave():
    osi = o3.OpenSeesInstance(ndm=2)
    lower = [1.0, 1.0]
    upper = [1.0, 1.0]
    locations = [1.0, 1.0]
    o3.mesh.Bgwave(osi, lower=lower, upper=upper, tol=1.0, meshtol=1.0, wavefilename="string", numl=1, locations=locations, numsub=1)


@pytest.mark.skip()
def test_bgstructure():
    osi = o3.OpenSeesInstance(ndm=2)
    lower = [1.0, 1.0]
    upper = [1.0, 1.0]
    o3.mesh.Bgstructure(osi, lower=lower, upper=upper, tol=1.0, meshtol=1.0, id=1, numnodes=1, snodes=1)


@pytest.mark.skip()
def test_bglarge_size():
    osi = o3.OpenSeesInstance(ndm=2)
    lower = [1.0, 1.0]
    upper = [1.0, 1.0]
    llower = [1.0, 1.0]
    lupper = [1.0, 1.0]
    o3.mesh.BglargeSize(osi, lower=lower, upper=upper, tol=1.0, meshtol=1.0, level=1, llower=llower, lupper=lupper)

