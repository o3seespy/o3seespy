import o3seespy as o3  # for testing only
import pytest

def test_central_difference():
    osi = o3.OpenseesInstance(dimensions=2)
    o3.integrator.CentralDifference(osi)


def test_newmark():
    osi = o3.OpenseesInstance(dimensions=2)
    o3.integrator.Newmark(osi, gamma=1.0, beta=1.0, form=1)


def test_hht():
    osi = o3.OpenseesInstance(dimensions=2)
    o3.integrator.HHT(osi, alpha=1.0, gamma=1.0, beta=1.0)


def test_generalized_alpha():
    osi = o3.OpenseesInstance(dimensions=2)
    o3.integrator.GeneralizedAlpha(osi, alpha_m=1.0, alpha_f=1.0, gamma=1.0, beta=1.0)


def test_trbdf2():
    osi = o3.OpenseesInstance(dimensions=2)
    o3.integrator.TRBDF2(osi)


def test_explicit_difference():
    osi = o3.OpenseesInstance(dimensions=2)
    o3.integrator.ExplicitDifference(osi)

