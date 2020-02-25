import o3seespy as o3  # for testing only


def test_initial_state_analysis_wrapper():
    osi = o3.OpenSeesInstance(ndm=2)
    mat = o3.nd_material.ElasticIsotropic(osi, e_mod=1.0, nu=1.0, rho=0.0)
    o3.nd_material.InitialStateAnalysisWrapper(osi, n_d_mat=mat, n_dim=1)

