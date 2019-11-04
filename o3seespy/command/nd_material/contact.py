from o3seespy.command.nd_material.base_material import NDMaterialBase


class ContactMaterial2D(NDMaterialBase):

    def __init__(self, osi, mu, big_g, c, t):
        self.mu = float(mu)
        self.big_g = float(big_g)
        self.c = float(c)
        self.t = float(t)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.mu, self.big_g, self.c, self.t]
        self.to_process(osi)


class ContactMaterial3D(NDMaterialBase):

    def __init__(self, osi, mu, big_g, c, t):
        self.mu = float(mu)
        self.big_g = float(big_g)
        self.c = float(c)
        self.t = float(t)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.mu, self.big_g, self.c, self.t]
        self.to_process(osi)
