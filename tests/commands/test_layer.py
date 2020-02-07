import o3seespy as o3  # for testing only


def test_straight():
    osi = o3.OpenSeesInstance(ndm=2)
    start = [1.0, 1.0]
    end = [1.0, 1.0]
    o3.layer.Straight(osi, num_fiber=1, area_fiber=1.0, start=start, end=end)


def test_circ():
    osi = o3.OpenSeesInstance(ndm=2)
    center = [1.0, 1.0]
    o3.layer.Circ(osi, num_fiber=6, area_fiber=1.0, center=center, radius=1.0, ang=[0.0, 360.0-360/6])
