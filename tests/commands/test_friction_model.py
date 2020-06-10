import o3seespy as o3  # for testing only


def test_coulomb():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.friction_model.Coulomb(osi, mu=1.0)


def test_vel_dependent():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.friction_model.VelDependent(osi, mu_slow=1.0, mu_fast=1.0, trans_rate=1.0)


def test_vel_normal_frc_dep():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.friction_model.VelNormalFrcDep(osi, a_slow=1.0, n_slow=1.0, a_fast=1.0, n_fast=1.0, alpha0=1.0, alpha1=1.0, alpha2=1.0, max_mu_fact=1.0)


def test_vel_pressure_dep():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.friction_model.VelPressureDep(osi, mu_slow=1.0, mu_fast0=1.0, big_a=1.0, delta_mu=1.0, alpha=1.0, trans_rate=1.0)


def test_vel_dep_multi_linear():
    osi = o3.OpenSeesInstance(ndm=2)
    vel_points = [0.0, 1.0]
    frn_points = [1.0, 1.0]
    o3.friction_model.VelDepMultiLinear(osi, vel_points=vel_points, frn_points=frn_points)


if __name__ == '__main__':
    test_vel_dep_multi_linear()
