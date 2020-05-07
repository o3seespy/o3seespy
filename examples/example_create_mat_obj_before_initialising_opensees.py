import o3seespy as o3

using_prebuild = 1

if using_prebuild:  # Allows you to define the material outside of using opensees
    mat = o3.nd_material.ElasticIsotropic(None, e_mod=1.0, nu=1.0, rho=0.0)
    osi = o3.OpenSeesInstance(ndm=2)
    mat.build(osi)
else:
    osi = o3.OpenSeesInstance(ndm=2)
    o3.nd_material.ElasticIsotropic(osi, e_mod=1.0, nu=1.0, rho=0.0)
