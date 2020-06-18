import o3seespy as o3  # for testing only
import pytest


def test_plane_stress_user_material():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.nd_material.PlaneStressUserMaterial(osi, nstatevs=1, nprops=1, fc=1.0, ft=1.0, fcu=1.0, epsc0=1.0, epscu=1.0, epstu=1.0, stc=1.0)


def test_plate_from_plane_stress():
    osi = o3.OpenSeesInstance(ndm=2)
    mat = o3.nd_material.ElasticIsotropic(osi, e_mod=1.0, nu=1.0, rho=0.0)
    o3.nd_material.PlateFromPlaneStress(osi, pre_def_mat=mat, outof_plane_modulus=1.0)


def test_plate_rebar():
    osi = o3.OpenSeesInstance(ndm=2)
    mat = o3.uniaxial_material.Elastic(osi, 1.0)
    o3.nd_material.PlateRebar(osi, pre_def_mat=mat, sita=1.0)

