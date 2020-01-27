import o3seespy as o3  # for testing only
import pytest

@pytest.mark.skip()
def test_fluid_solid_porous_material():
    osi = o3.OpenseesInstance(dimensions=2)
    soil_mat = o3.nd_material.ElasticIsotropic(osi, big_e=1.0, v=1.0, rho=0.0)
    o3.nd_material.FluidSolidPorousMaterial(osi, nd=1.0, soil_mat=soil_mat, combined_bulk_modul=1.0, pa=101.0)

