import o3seespy as o3  # for testing only
import pytest

@pytest.mark.skip()
def test_fluid_solid_porous_material():
    osi = o3.OpenSeesInstance(ndm=2)
    soil_mat = o3.nd_material.ElasticIsotropic(osi, e_mod=1.0, v=1.0, rho=0.0)
    o3.nd_material.FluidSolidPorousMaterial(osi, nd=1.0, soil_mat=soil_mat, combined_bulk_modul=1.0, pa=101.0)

