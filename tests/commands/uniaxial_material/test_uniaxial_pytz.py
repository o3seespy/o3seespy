import o3seespy as o3  # for testing only


def test_py_simple1():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.PySimple1(osi, soil_type=1, pult=1.0, y50=1.0, cd=1.0, c=0.0)


def test_tz_simple1():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.TzSimple1(osi, soil_type=1, tult=1.0, z50=1.0, c=0.0)


def test_qz_simple1():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.QzSimple1(osi, qz_type=1, qult=1.0, z50=1.0, suction=0.0, c=0.0)


def test_py_liq1():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.PyLiq1(osi, soil_type=1, pult=1.0, y50=1.0, cd=1.0, c=1.0, p_res=1.0, ele1=1.0, ele2=1.0)


def test_tz_liq1():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.uniaxial_material.TzLiq1(osi, tz_type=1, tult=1.0, z50=1.0, c=1.0, ele1=1.0, ele2=1.0)

