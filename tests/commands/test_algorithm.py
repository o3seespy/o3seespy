import o3seespy as o3


def test_newton():
    osi = o3.OpenSeesInstance(ndm=2, ndf=2)
    o3.algorithm.Newton(osi)


def test_krylov_newton():
    osi = o3.OpenSeesInstance(ndm=2, ndf=2)
    o3.test_check.NormDispIncr(osi, tol=1.0e-3, max_iter=35, p_flag=0)
    o3.algorithm.KrylovNewton(osi)


if __name__ == '__main__':
    test_krylov_newton()
