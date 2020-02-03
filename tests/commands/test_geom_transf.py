import o3seespy as o3  # for testing only


def test_linear2d():
    osi = o3.OpenSeesInstance(ndm=2)
    d_i = [1.0, 1.0]
    d_j = [1.0, 1.0]
    o3.geom_transf.Linear2D(osi, d_i=d_i, d_j=d_j)

def test_linear3d():
    osi = o3.OpenSeesInstance(ndm=2)
    vecxz = [1.0, 1.0]
    d_i = [1.0, 1.0]
    d_j = [1.0, 1.0]
    o3.geom_transf.Linear3D(osi, vecxz=vecxz, d_i=d_i, d_j=d_j)


def test_p_delta2d():
    osi = o3.OpenSeesInstance(ndm=2)
    d_i = [1.0, 1.0]
    d_j = [1.0, 1.0]
    o3.geom_transf.PDelta2D(osi, d_i=d_i, d_j=d_j)

def test_p_delta3d():
    osi = o3.OpenSeesInstance(ndm=2)
    vecxz = [1.0, 1.0]
    d_i = [1.0, 1.0]
    d_j = [1.0, 1.0]
    o3.geom_transf.PDelta3D(osi, vecxz=vecxz, d_i=d_i, d_j=d_j)


def test_corotational2d():
    osi = o3.OpenSeesInstance(ndm=2)
    d_i = [1.0, 1.0]
    d_j = [1.0, 1.0]
    o3.geom_transf.Corotational2D(osi, d_i=d_i, d_j=d_j)

def test_corotational3d():
    osi = o3.OpenSeesInstance(ndm=2)
    vecxz = [1.0, 1.0]
    o3.geom_transf.Corotational3D(osi, vecxz=vecxz)

