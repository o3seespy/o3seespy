import o3seespy as o3  # for testing only
import pytest


def test_elastic():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.Elastic(osi, e_mod=1.0, eta=0.0, eneg=None)


def test_elastic_pp():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.ElasticPP(osi, e_mod=1.0, epsy_p=1.0, epsy_n=None, eps0=0.0)


def test_elastic_pp_gap():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.ElasticPPGap(osi, e_mod=1.0, fy=1.0, gap=1.0, eta=0.0, damage='noDamage')


def test_ent():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.ENT(osi, e_mod=1.0)


@pytest.mark.skip()
def test_parallel():
    osi = o3.OpenSeesInstance(ndm=2)
    tags = [1, 1]
    factor_args = [1.0, 1.0]
    o3.uniaxial_material.Parallel(osi, tags=tags, factor_args=factor_args)


@pytest.mark.skip()
def test_series():
    osi = o3.OpenSeesInstance(ndm=2)
    tags = [1, 1]
    o3.uniaxial_material.Series(osi, tags=tags)

