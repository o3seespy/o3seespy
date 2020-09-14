import o3seespy as o3  # for testing only
import pytest


def test_band_gen():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.system.BandGen(osi)


def test_band_spd():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.system.BandSPD(osi)


def test_profile_spd():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.system.ProfileSPD(osi)


def test_super_lu():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.system.SuperLU(osi)


def test_umf_pack():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.system.UmfPack(osi)


def test_full_general():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.system.FullGeneral(osi)


@pytest.mark.skip()
def test_sparse_sym():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.system.SparseSYM(osi)


#
# def test_mumps():
#     osi = o3.OpenSeesInstance(ndm=2)
#     o3.system.Mumps(osi, icntl14=20.0, icntl7=7)

