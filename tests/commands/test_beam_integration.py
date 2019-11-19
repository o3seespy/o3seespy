import o3seespy as o3  # for testing only
import pytest


@pytest.mark.skip()
def test_user_hinge():
    osi = o3.OpenseesInstance(dimensions=2)
    sec_e = o3.section.Elastic(osi, 10.0, 1.0, 1.0)
    o3.beam_integration.UserHinge(osi, sec_e=sec_e, np_l=1, secs_l=1, locs_l=1, wts_l=1, np_r=1, secs_r=1, locs_r=1, wts_r=1)


def test_hinge_midpoint():
    osi = o3.OpenseesInstance(dimensions=2)
    sec_i = o3.section.Elastic(osi, 10.0, 1.0, 1.0)
    sec_j = o3.section.Elastic(osi, 10.0, 1.0, 1.0)
    sec_e = o3.section.Elastic(osi, 10.0, 1.0, 1.0)
    o3.beam_integration.HingeMidpoint(osi, sec_i=sec_i, lp_i=1.0, sec_j=sec_j, lp_j=1.0, sec_e=sec_e)


def test_hinge_radau():
    osi = o3.OpenseesInstance(dimensions=2)
    sec_i = o3.section.Elastic(osi, 10.0, 1.0, 1.0)
    sec_j = o3.section.Elastic(osi, 10.0, 1.0, 1.0)
    sec_e = o3.section.Elastic(osi, 10.0, 1.0, 1.0)
    o3.beam_integration.HingeRadau(osi, sec_i=sec_i, lp_i=1, sec_j=sec_j, lp_j=1, sec_e=sec_e)


def test_hinge_radau_two():
    osi = o3.OpenseesInstance(dimensions=2)
    sec_i = o3.section.Elastic(osi, 10.0, 1.0, 1.0)
    sec_j = o3.section.Elastic(osi, 10.0, 1.0, 1.0)
    sec_e = o3.section.Elastic(osi, 10.0, 1.0, 1.0)
    o3.beam_integration.HingeRadauTwo(osi, sec_i=sec_i, lp_i=1, sec_j=sec_j, lp_j=1, sec_e=sec_e)


def test_beamhinge_endpoint():
    osi = o3.OpenseesInstance(dimensions=2)
    sec_j = o3.section.Elastic(osi, 10.0, 1.0, 1.0)
    sec_e = o3.section.Elastic(osi, 10.0, 1.0, 1.0)
    o3.beam_integration.BeamhingeEndpoint(osi, lp_i=1.0, sec_j=sec_j, lp_j=1.0, sec_e=sec_e)

