from o3seespy.command.element.base_element import ElementBase


class ElasticBeamColumn(ElementBase):
    op_type = 'elasticBeamColumn'

    def __init__(self, osi, ele_nodes, big_a, big_e, iz, transf, mass=None, c_mass=False, big_g, big_j, iy):
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
        self.big_g = float(big_g)
        self.big_j = float(big_j)
        self.iy = float(iy)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.big_a, self.big_e, self.iz, self.transf.tag, self.big_g, self.big_j, self.iy]
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


class ElasticTimoshenkoBeam(ElementBase):
    op_type = 'ElasticTimoshenkoBeam'

    def __init__(self, osi, ele_nodes, big_e, big_g, big_a, iz, avy, transf, mass=None, c_mass=False, jx, iy, iz_2, avz, c_mas=False):
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
        self.jx = float(jx)
        self.iy = float(iy)
        self.iz_2 = iz_2
        self.avz = float(avz)
        self.c_mas = c_mas
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.big_e, self.big_g, self.big_a, self.iz, self.avy, self.transf.tag, self.jx, self.iy, self.iz_2, self.avz]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        if getattr(self, 'c_mass'):
            self._parameters += ['-cMass']
        if getattr(self, 'c_mas'):
            self._parameters += ['-cMas']
        self.to_process(osi)



class DispBeamColumn(ElementBase):
    op_type = 'dispBeamColumn'

    def __init__(self, osi, i_node, j_node, transf, integration, c_mass=False, mass=None):
        self.i_node = int(i_node)
        self.j_node = int(j_node)
        self.transf = transf
        self.integration = integration
        self.c_mass = c_mass
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.i_node, self.j_node, self.transf.tag, self.integration.tag]
        if getattr(self, 'c_mass'):
            self._parameters += ['-cMass']
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        self.to_process(osi)


class ForceBeamColumn(ElementBase):
    op_type = 'forceBeamColumn'

    def __init__(self, osi, i_node, j_node, transf, integration, tol=1e-12, mass=None, iter=None):
        self.i_node = int(i_node)
        self.j_node = int(j_node)
        self.transf = transf
        self.integration = integration
        self.tol = float(tol)
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        if iter is None:
            self.iter = None
        else:
            self.iter = float(iter)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.i_node, self.j_node, self.transf.tag, self.integration.tag, self.tol]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        if getattr(self, 'iter') is not None:
            self._parameters += ['-iter', self.iter]
        self.to_process(osi)


class NonlinearBeamColumnintegration(ElementBase):
    op_type = 'nonlinearBeamColumn'

    def __init__(self, osi, i_node, j_node, num_intgr_pts, sec, transf, tol=1e-12, mass=None, int_type, iter=None):
        self.i_node = int(i_node)
        self.j_node = int(j_node)
        self.num_intgr_pts = int(num_intgr_pts)
        self.sec = sec
        self.transf = transf
        self.tol = float(tol)
        if mass is None:
            self.mass = None
        else:
            self.mass = float(mass)
        self.int_type = int_type
        if iter is None:
            self.iter = None
        else:
            self.iter = float(iter)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.i_node, self.j_node, self.num_intgr_pts, self.sec.tag, self.transf.tag, self.tol, '-integration', self.int_type]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        if getattr(self, 'iter') is not None:
            self._parameters += ['-iter', self.iter]
        self.to_process(osi)


class DispBeamColumnIntmass(ElementBase):
    op_type = 'dispBeamColumnInt'

    def __init__(self, osi, ele_nodes, num_intgr_pts, sec, transf, c_rot, mass_dens):
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.num_intgr_pts = int(num_intgr_pts)
        self.sec = sec
        self.transf = transf
        self.c_rot = float(c_rot)
        self.mass_dens = float(mass_dens)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.num_intgr_pts, self.sec.tag, self.transf.tag, self.c_rot, '-mass', self.mass_dens]
        self.to_process(osi)


class MVLEMthick(ElementBase):
    op_type = 'MVLEM'

    def __init__(self, osi, dens, ele_nodes, m, c, thicknesses):
        self.dens = float(dens)
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.m = int(m)
        self.c = float(c)
        self.thicknesses = thicknesses
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.dens, *self.ele_nodes, self.m, self.c, '-thick', *self.thicknesses]
        self.to_process(osi)

class MVLEMwidth(ElementBase):
    op_type = 'MVLEM'

    def __init__(self, osi, dens, ele_nodes, m, c, widths):
        self.dens = float(dens)
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.m = int(m)
        self.c = float(c)
        self.widths = widths
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.dens, *self.ele_nodes, self.m, self.c, '-width', *self.widths]
        self.to_process(osi)

class MVLEMrho(ElementBase):
    op_type = 'MVLEM'

    def __init__(self, osi, dens, ele_nodes, m, c, reinforcing_ratios):
        self.dens = float(dens)
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.m = int(m)
        self.c = float(c)
        self.reinforcing_ratios = reinforcing_ratios
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.dens, *self.ele_nodes, self.m, self.c, '-rho', *self.reinforcing_ratios]
        self.to_process(osi)

class MVLEMmatConcrete(ElementBase):
    op_type = 'MVLEM'

    def __init__(self, osi, dens, ele_nodes, m, c, concrete_tags):
        self.dens = float(dens)
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.m = int(m)
        self.c = float(c)
        self.concrete_tags = concrete_tags
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.dens, *self.ele_nodes, self.m, self.c, '-matConcrete', *self.concrete_tags]
        self.to_process(osi)

class MVLEMmatSteel(ElementBase):
    op_type = 'MVLEM'

    def __init__(self, osi, dens, ele_nodes, m, c, steel_tags):
        self.dens = float(dens)
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.m = int(m)
        self.c = float(c)
        self.steel_tags = steel_tags
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.dens, *self.ele_nodes, self.m, self.c, '-matSteel', *self.steel_tags]
        self.to_process(osi)

class MVLEMmatShear(ElementBase):
    op_type = 'MVLEM'

    def __init__(self, osi, dens, ele_nodes, m, c, shear):
        self.dens = float(dens)
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.m = int(m)
        self.c = float(c)
        self.shear = shear
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.dens, *self.ele_nodes, self.m, self.c, '-matShear', self.shear.tag]
        self.to_process(osi)


class Thick(ElementBase):
    op_type = ''

    def __init__(self, osi, ele_nodes, m, c, thicknesses):
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.m = int(m)
        self.c = float(c)
        self.thicknesses = thicknesses
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.m, self.c, '-thick', *self.thicknesses]
        self.to_process(osi)

class Width(ElementBase):
    op_type = ''

    def __init__(self, osi, ele_nodes, m, c, widths):
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.m = int(m)
        self.c = float(c)
        self.widths = widths
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.m, self.c, '-width', *self.widths]
        self.to_process(osi)

class Mat(ElementBase):
    op_type = ''

    def __init__(self, osi, ele_nodes, m, c, material_tags):
        self.ele_nodes = [x.tag for x in ele_nodes]
        self.m = int(m)
        self.c = float(c)
        self.material_tags = material_tags
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.m, self.c, '-mat', *self.material_tags]
        self.to_process(osi)
