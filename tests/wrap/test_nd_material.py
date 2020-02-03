import o3seespy as o3
import pytest


def test_can_set_PressureIndependMultiYield():
    osi = o3.OpenSeesInstance(ndm=2)

    # Define material
    nd = 2
    rho = 1400
    ref_shear_modul = 40e6
    ref_bulk_modul = 60e6
    cohesi = 60e3
    peak_shear_stra = 0.02
    o3.nd_material.PressureIndependMultiYield(osi, nd, rho, ref_shear_modul, ref_bulk_modul,
                                               cohesi, peak_shear_stra)


def test_can_set_PM4Sand():
    osi = o3.OpenSeesInstance(ndm=2)

    # Define material
    d_r = 0.5
    hpo = 0.4
    g0 = 555.
    p_atm = 101.3
    rho = 1.7
    o3.nd_material.PM4Sand(osi, d_r, g0, hpo, rho, p_atm)


def skip_test_can_set_stress_density_model():
    osi = o3.OpenSeesInstance(ndm=2)

    # Define material
    o3.nd_material.StressDensityModel(osi, den=1.8, e_init=0.73, area=250., n=0.6, nu=0.3, a1=0.58, b1=0.023,
                                       a2=230.0, b2=65., a3=79., b3=16., fd=4.0, mu_not=0.22, mu_cyc=0.0, sc=0.0055,
                                       big_m=0.607, p_atm=98.1)
