import o3seespy as o3  # for testing only
import pytest


def test_steel02():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.Steel02(osi, fy=1.0, e0=1.0, b=1.0, params=[15, 0.925, 0.15])



def test_hysteretic():
    osi = o3.OpenSeesInstance(ndm=2)
    p1 = [0.5, 0.5]
    p2 = [1.0, 1.0]
    p3 = [0, 1.5]
    n1 = [-0.5, -0.5]
    n2 = [-1.0, -1.0]
    n3 = [0, -1.5]
    o3.uniaxial_material.Hysteretic(osi, p1=p1, p2=p2, p3=p3, n1=n1, n2=n2, n3=n3, pinch_x=1, pinch_y=0, damage1=0, damage2=0)


def test_reinforcing_steel_ga_buck():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.ReinforcingSteelGABuck(osi, fy=1.0, fu=1.0, es=1.0, esh=1.0, eps_sh=1.0, eps_ult=1.0, lsr=1.0, beta=1.0, r=1.0, gamma=1.0)


def test_reinforcing_steel_dm_buck():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.ReinforcingSteelDMBuck(osi, fy=1.0, fu=1.0, es=1.0, esh=1.0, eps_sh=1.0, eps_ult=1.0, lsr_2=1, alpha=1.0)


def test_reinforcing_steel_cm_fatigue():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.ReinforcingSteelCMFatigue(osi, fy=1.0, fu=1.0, es=1.0, esh=1.0, eps_sh=1.0, eps_ult=1.0, cf=1.0, alpha_2=1, cd=1.0)


def test_reinforcing_steel_iso_hard():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.ReinforcingSteelIsoHard(osi, fy=1.0, fu=1.0, es=1.0, esh=1.0, eps_sh=1.0, eps_ult=1.0, a1=4.3, limit=1.0)


def test_reinforcing_steel_mp_curve_params():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.ReinforcingSteelMPCurveParams(osi, fy=1.0, fu=1.0, es=1.0, esh=1.0, eps_sh=1.0, eps_ult=1.0, r1=0.333, r2=18.0, r3=4.0)


@pytest.mark.skip()  # not connected in openseespy - works in tcl
def test_dodd_restrepo():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.DoddRestrepo(osi, fy=1.0, fsu=1.0, esh=1.0, esu=1.0, youngs=1.0, eshi=1.0, fshi=1.0, omega_fac=1.0)


def test_ramberg_osgood_steel():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.RambergOsgoodSteel(osi, fy=1.0, e0=1.0, a=1.0, n=1.0)


def test_steel_mpf():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.SteelMPF(osi, fyp=1.0, fyn=1.0, e0=1.0, bp=1.0, bn=1.0, params=[1.0, 1.0, 1.0], a1=0.0, a2=1.0, a3=0.0, a4=1.0)


def test_steel01thermal():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.Steel01Thermal(osi, fy=1.0, e0=1.0, b=1.0, a1=1.0, a2=1.0, a3=1.0, a4=1.0)

