from o3seespy.command.element.base_element import ElementBase


class ElasticBeamColumn2D(ElementBase):
    op_type = 'elasticBeamColumn'

    def __init__(self, osi, ele_nodes, big_a, big_e, iz, transf, mass=None, c_mass=False):
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.big_a = float(big_a)
        self.big_e = float(big_e)
        self.iz = float(iz)
        self.transf = transf
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        self.c_mass = c_mass
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.big_a, self.big_e, self.iz, self.transf.tag]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        if getattr(self, 'c_mass'):
            self._parameters += ['-cMass']
        self.to_process(osi)

class ElasticBeamColumn3D(ElementBase):
    op_type = 'elasticBeamColumn'

    def __init__(self, osi, ele_nodes, big_a, big_e, big_g, big_j, iy, iz, transf, mass=None, c_mass=False):
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.big_a = float(big_a)
        self.big_e = float(big_e)
        self.big_g = float(big_g)
        self.big_j = float(big_j)
        self.iy = float(iy)
        self.iz = float(iz)
        self.transf = transf
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        self.c_mass = c_mass
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.big_a, self.big_e, self.big_g, self.big_j, self.iy, self.iz, self.transf.tag]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        if getattr(self, 'c_mass'):
            self._parameters += ['-cMass']
        self.to_process(osi)


class ModElasticBeam2Dmass(ElementBase):
    op_type = 'ModElasticBeam2d'

    def __init__(self, osi, ele_nodes, big_a, big_e, iz, k11, k33, k44, transf, mass_dens, c_mass=False):
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.big_a = float(big_a)
        self.big_e = float(big_e)
        self.iz = float(iz)
        self.k11 = float(k11)
        self.k33 = float(k33)
        self.k44 = float(k44)
        self.transf = transf
        self.mass_dens = float(mass_dens)
        self.c_mass = c_mass
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.big_a, self.big_e, self.iz, self.k11, self.k33, self.k44, self.transf.tag, '-mass', self.mass_dens]
        if getattr(self, 'c_mass'):
            self._parameters += ['-cMass']
        self.to_process(osi)


class ElasticTimoshenkoBeam2D(ElementBase):
    op_type = 'ElasticTimoshenkoBeam'

    def __init__(self, osi, ele_nodes, big_e, big_g, big_a, iz, avy, transf, mass=None, c_mass=False):
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.big_e = float(big_e)
        self.big_g = float(big_g)
        self.big_a = float(big_a)
        self.iz = float(iz)
        self.avy = float(avy)
        self.transf = transf
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        self.c_mass = c_mass
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.big_e, self.big_g, self.big_a, self.iz, self.avy, self.transf.tag]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        if getattr(self, 'c_mass'):
            self._parameters += ['-cMass']
        self.to_process(osi)

class ElasticTimoshenkoBeam3D(ElementBase):
    op_type = 'ElasticTimoshenkoBeam'

    def __init__(self, osi, ele_nodes, big_e, big_g, big_a, iz, jx, iy, iz_2, avy, avz, transf, mass=None, c_mass=False):
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.big_e = float(big_e)
        self.big_g = float(big_g)
        self.big_a = float(big_a)
        self.iz = float(iz)
        self.jx = float(jx)
        self.iy = float(iy)
        self.iz_2 = iz_2
        self.avy = float(avy)
        self.avz = float(avz)
        self.transf = transf
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        self.c_mass = c_mass
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.big_e, self.big_g, self.big_a, self.iz, self.jx, self.iy, self.iz_2, self.avy, self.avz, self.transf.tag]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        if getattr(self, 'c_mass'):
            self._parameters += ['-cMass']
        self.to_process(osi)



class DispBeamColumn(ElementBase):
    op_type = 'dispBeamColumn'

    def __init__(self, osi, i_node, j_node, transf, integration, c_mass=False, mass=None):
        self.i_node = i_node
        self.j_node = j_node
        self.transf = transf
        self.integration = integration
        self.c_mass = c_mass
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.i_node.tag, self.j_node.tag, self.transf.tag, self.integration.tag]
        if getattr(self, 'c_mass'):
            self._parameters += ['-cMass']
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        self.to_process(osi)


class ForceBeamColumn(ElementBase):
    op_type = 'forceBeamColumn'

    def __init__(self, osi, i_node, j_node, transf, integration, max_iter=None, tol=None, mass=None):
        self.i_node = i_node
        self.j_node = j_node
        self.transf = transf
        self.integration = integration
        self.max_iter = int(max_iter)
        self.tol = float(tol)
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.i_node.tag, self.j_node.tag, self.transf.tag, self.integration.tag, self.tol]
        if getattr(self, 'max_iter') is not None:
            self._parameters += ['-iter', self.max_iter]
        if getattr(self, 'tol') is not None:
            if getattr(self, 'max_iter') is None:
                raise ValueError('Cannot set: tol and not: max_iter')
            self._parameters += [self.tol]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        self.to_process(osi)


class NonlinearBeamColumn(ElementBase):
    op_type = 'nonlinearBeamColumn'

    def __init__(self, osi, i_node, j_node, num_intgr_pts, sec, transf, max_iter=None, tol=None, mass=None, int_type=None):
        self.i_node = i_node
        self.j_node = j_node
        self.num_intgr_pts = int(num_intgr_pts)
        self.sec = sec
        self.transf = transf
        self.max_iter = int(max_iter)
        self.tol = float(tol)
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        self.int_type = int_type
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.i_node.tag, self.j_node.tag, self.num_intgr_pts, self.sec.tag, self.transf.tag, self.tol]
        if getattr(self, 'max_iter') is not None:
            self._parameters += ['-iter', self.max_iter]
        if getattr(self, 'tol') is not None:
            if getattr(self, 'max_iter') is None:
                raise ValueError('Cannot set: tol and not: max_iter')
            self._parameters += [self.tol]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        if getattr(self, 'int_type') is not None:
            self._parameters += ['-integration', self.int_type]
        self.to_process(osi)


class DispBeamColumnInt(ElementBase):
    op_type = 'dispBeamColumnInt'

    def __init__(self, osi, ele_nodes, num_intgr_pts, sec, transf, c_rot, mass_dens=None):
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.num_intgr_pts = int(num_intgr_pts)
        self.sec = sec
        self.transf = transf
        self.c_rot = float(c_rot)
        if mass_dens is None:
            self.mass_dens = None
        else:
            self.mass_dens = float(mass_dens)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.num_intgr_pts, self.sec.tag, self.transf.tag, self.c_rot]
        if getattr(self, 'mass_dens') is not None:
            self._parameters += ['-mass', self.mass_dens]
        self.to_process(osi)


class MVLEM(ElementBase):
    op_type = 'MVLEM'

    def __init__(self, osi, dens, ele_nodes, m, c, thick=None, widths=None, rho=None, mat_concrete_tags=None, mat_steel_tags=None, mat_shear=None):
        self.dens = float(dens)
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.m = int(m)
        self.c = float(c)
        self.thick = thick
        self.widths = widths
        self.rho = rho
        self.mat_concrete_tags = [x.tag for x in mat_concrete_tags]
        self.mat_steel_tags = [x.tag for x in mat_steel_tags]
        self.mat_shear = mat_shear
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.dens, *self.ele_nodes, self.m, self.c]
        if getattr(self, 'thick') is not None:
            self._parameters += ['-thick', *self.thick]
        if getattr(self, 'widths') is not None:
            self._parameters += ['-width', *self.widths]
        if getattr(self, 'rho') is not None:
            self._parameters += ['-rho', *self.rho]
        if getattr(self, 'mat_concrete_tags') is not None:
            self._parameters += ['-matConcrete', *self.mat_concrete_tags]
        if getattr(self, 'mat_steel_tags') is not None:
            self._parameters += ['-matSteel', *self.mat_steel_tags]
        if getattr(self, 'mat_shear') is not None:
            self._parameters += ['-matShear', self.mat_shear.tag]
        self.to_process(osi)


class SFIMVLEM(ElementBase):
    op_type = 'SFI_MVLEM'

    def __init__(self, osi, ele_nodes, m, c, thick=None, widths=None, mat_tags=None):
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.m = int(m)
        self.c = float(c)
        self.thick = thick
        self.widths = widths
        self.mat_tags = mat_tags
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.m, self.c]
        if getattr(self, 'thick') is not None:
            self._parameters += ['-thick', *self.thick]
        if getattr(self, 'widths') is not None:
            self._parameters += ['-width', *self.widths]
        if getattr(self, 'mat_tags') is not None:
            self._parameters += ['-mat', *self.mat_tags]
        self.to_process(osi)
