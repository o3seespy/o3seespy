import o3seespy as o3  # for testing only
import pytest


def test_lobatto():
    osi = o3.OpenSeesInstance(ndm=2)
    sec = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
    o3.beam_integration.Lobatto(osi, sec=sec, big_n=1)


def test_legendre():
    osi = o3.OpenSeesInstance(ndm=2)
    sec = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
    o3.beam_integration.Legendre(osi, sec=sec, big_n=1)


def test_newton_cotes():
    osi = o3.OpenSeesInstance(ndm=2)
    sec = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
    o3.beam_integration.NewtonCotes(osi, sec=sec, big_n=1)


def test_radau():
    osi = o3.OpenSeesInstance(ndm=2)
    sec = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
    o3.beam_integration.Radau(osi, sec=sec, big_n=1)


def test_trapezoidal():
    osi = o3.OpenSeesInstance(ndm=2)
    sec = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
    o3.beam_integration.Trapezoidal(osi, sec=sec, big_n=1)


def test_composite_simpson():
    osi = o3.OpenSeesInstance(ndm=2)
    sec = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
    o3.beam_integration.CompositeSimpson(osi, sec=sec, big_n=1)



def test_user_defined():
    osi = o3.OpenSeesInstance(ndm=2)
    secs = [o3.section.Elastic2D(osi, e_mod=1.0, area=1.0, iz=1.0),
             o3.section.Elastic2D(osi, e_mod=1.0, area=1.0, iz=1.0)]
    o3.beam_integration.UserDefined(osi, big_n=2, secs=secs, locs=[0.2, 0.9], wts=[0.5, 0.5])


def test_fixed_location():
    osi = o3.OpenSeesInstance(ndm=2)
    secs = [o3.section.Elastic2D(osi, e_mod=1.0, area=1.0, iz=1.0),
            o3.section.Elastic2D(osi, e_mod=1.0, area=1.0, iz=1.0)]
    o3.beam_integration.FixedLocation(osi, big_n=2, secs=secs, locs=[0.2, 0.9])



def test_low_order():
    osi = o3.OpenSeesInstance(ndm=2)
    secs = [o3.section.Elastic2D(osi, e_mod=1.0, area=1.0, iz=1.0),
            o3.section.Elastic2D(osi, e_mod=1.0, area=1.0, iz=1.0)]
    o3.beam_integration.LowOrder(osi, big_n=2, secs=secs, locs=[0.2, 0.9], wts=[0.5, 0.5])


def test_mid_distance():
    osi = o3.OpenSeesInstance(ndm=2)
    secs = [o3.section.Elastic2D(osi, e_mod=1.0, area=1.0, iz=1.0),
            o3.section.Elastic2D(osi, e_mod=1.0, area=1.0, iz=1.0)]
    o3.beam_integration.MidDistance(osi, big_n=2, secs=secs, locs=[0.2, 0.9])


def test_user_hinge():
    osi = o3.OpenSeesInstance(ndm=2)
    sec_e = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
    secs_l = [o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)]
    secs_r = [o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)]
    o3.beam_integration.UserHinge(osi, sec_e=sec_e, np_l=1, secs_ls=secs_l, locs_l=[1], wts_l=[1], np_r=1,
                                  secs_rs=secs_r, locs_r=[1], wts_r=[1])


def test_hinge_midpoint():
    osi = o3.OpenSeesInstance(ndm=2)
    sec_i = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
    sec_j = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
    sec_e = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
    o3.beam_integration.HingeMidpoint(osi, sec_i=sec_i, lp_i=1.0, sec_j=sec_j, lp_j=1.0, sec_e=sec_e)


def test_hinge_radau():
    osi = o3.OpenSeesInstance(ndm=2)
    sec_i = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
    sec_j = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
    sec_e = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
    o3.beam_integration.HingeRadau(osi, sec_i=sec_i, lp_i=1.0, sec_j=sec_j, lp_j=1.0, sec_e=sec_e)


def test_hinge_radau_two():
    osi = o3.OpenSeesInstance(ndm=2)
    sec_i = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
    sec_j = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
    sec_e = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
    o3.beam_integration.HingeRadauTwo(osi, sec_i=sec_i, lp_i=1.0, sec_j=sec_j, lp_j=1.0, sec_e=sec_e)


@pytest.mark.skip()
def test_beamhinge_endpoint():
    osi = o3.OpenSeesInstance(ndm=2)
    sec_j = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
    sec_e = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)
    o3.beam_integration.BeamhingeEndpoint(osi, lp_i=1.0, sec_j=sec_j, lp_j=1.0, sec_e=sec_e)

