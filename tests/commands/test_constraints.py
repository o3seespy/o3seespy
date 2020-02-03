import o3seespy as o3  # for testing only


def test_plain():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.constraints.Plain(osi)


def test_lagrange():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.constraints.Lagrange(osi, alpha_m=1.0)


def test_penalty():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.constraints.Penalty(osi, alpha_m=1.0)


def test_transformation():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.constraints.Transformation(osi)

