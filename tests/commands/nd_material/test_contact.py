import o3seespy as o3  # for testing only


def test_contact_material2d():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.nd_material.ContactMaterial2D(osi, mu=1.0, g_mod=1.0, c=1.0, t=1.0)


def test_contact_material3d():
    osi = o3.OpenSeesInstance(ndm=3)
    o3.nd_material.ContactMaterial3D(osi, mu=1.0, g_mod=1.0, c=1.0, t=1.0)

