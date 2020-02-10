import o3seespy as o3  # for testing only


def test_straight():
    osi = o3.OpenSeesInstance(ndm=2)
    start = [1.0, 1.0]
    end = [1.0, 1.0]
    rebar = o3.uniaxial_material.Steel01(osi, fy=60.0, e0=30000.0, b=0.02)
    o3.layer.Straight(osi, rebar, num_fiber=1, area_fiber=1.0, start=start, end=end)


def test_circ():
    osi = o3.OpenSeesInstance(ndm=2)
    center = [1.0, 1.0]
    rebar = o3.uniaxial_material.Steel01(osi, fy=60.0, e0=30000.0, b=0.02)
    o3.layer.Circ(osi, rebar, num_fiber=6, area_fiber=1.0, center=center, radius=1.0, ang=[0.0, 360.0-360/6])
