from o3seespy.command.element.base_element import ElementBase


class ElasticBeamColumn(ElementBase):

    def __init__(self, osi, ele_nodes, big_a, big_e, iz, transf, mass=None, cMass=False, big_g, big_j, iy):
        self.ele_nodes = ele_nodes
        self.big_a = float(big_a)
        self.big_e = float(big_e)
        self.iz = float(iz)
        self.transf = transf.tag
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
        if getattr(self, 'c_mass') is not None:
            self._parameters += ['-c_mass']
        self.to_process(osi)


class ModElasticBeam2dmass(ElementBase):

    def __init__(self, osi, ele_nodes, big_a, big_e, iz, k11, k33, k44, transf, mass_dens, cMass=False):
        self.ele_nodes = ele_nodes
        self.big_a = float(big_a)
        self.big_e = float(big_e)
        self.iz = float(iz)
        self.k11 = float(k11)
        self.k33 = float(k33)
        self.k44 = float(k44)
        self.transf = transf.tag
        self.mass_dens = float(mass_dens)
        self.c_mass = c_mass
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.big_a, self.big_e, self.iz, self.k11, self.k33, self.k44, self.transf.tag, '-mass', self.mass_dens]
        if getattr(self, 'c_mass') is not None:
            self._parameters += ['-c_mass']
        self.to_process(osi)


class ElasticTimoshenkoBeammass(ElementBase):

    def __init__(self, osi, ele_nodes, big_e, big_g, big_a, iz, avy, transf, mass_dens, cMass=False, jx, iy, iz_2, avz, mass=None, cMas=False):
        self.ele_nodes = ele_nodes
        self.big_e = float(big_e)
        self.big_g = float(big_g)
        self.big_a = float(big_a)
        self.iz = float(iz)
        self.avy = float(avy)
        self.transf = transf.tag
        self.mass_dens = mass_dens
        self.c_mass = c_mass
        self.jx = float(jx)
        self.iy = float(iy)
        self.iz_2 = iz_2
        self.avz = float(avz)
        self.mass = mass
        self.c_mas = c_mas
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.big_e, self.big_g, self.big_a, self.iz, self.avy, self.transf.tag, '-mass', self.mass_dens, self.jx, self.iy, self.iz_2, self.avz]
        if getattr(self, 'c_mass') is not None:
            self._parameters += ['-c_mass']
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        if getattr(self, 'c_mas') is not None:
            self._parameters += ['-c_mas']
        self.to_process(osi)



class DispBeamColumn(ElementBase):

    def __init__(self, osi, i_node, j_node, transf, integration, cMass=False, mass=0.0):
        self.i_node = int(i_node)
        self.j_node = int(j_node)
        self.transf = transf.tag
        self.integration = integration.tag
        self.c_mass = c_mass
        self.mass = float(mass)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.i_node, self.j_node, self.transf.tag, self.integration.tag]
        if getattr(self, 'c_mass') is not None:
            self._parameters += ['-c_mass']
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        self.to_process(osi)


class ForceBeamColumniter(ElementBase):

    def __init__(self, osi, i_node, j_node, transf, integration, max_iter=10, tol=1e-12, mass=0.0):
        self.i_node = int(i_node)
        self.j_node = int(j_node)
        self.transf = transf.tag
        self.integration = integration.tag
        self.max_iter = float(max_iter)
        self.tol = float(tol)
        self.mass = float(mass)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.i_node, self.j_node, self.transf.tag, self.integration.tag, '-iter', self.max_iter, self.tol]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        self.to_process(osi)


class NonlinearBeamColumniter(ElementBase):

    def __init__(self, osi, i_node, j_node, num_intgr_pts, sec, transf, max_iter=10, tol=1e-12, mass=0.0):
        self.i_node = int(i_node)
        self.j_node = int(j_node)
        self.num_intgr_pts = int(num_intgr_pts)
        self.sec = sec.tag
        self.transf = transf.tag
        self.max_iter = float(max_iter)
        self.tol = float(tol)
        self.mass = float(mass)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.i_node, self.j_node, self.num_intgr_pts, self.sec.tag, self.transf.tag, '-iter', self.max_iter, self.tol]
        if getattr(self, 'mass') is not None:
            self._parameters += ['-mass', self.mass]
        self.to_process(osi)

class NonlinearBeamColumnintegration(ElementBase):

    def __init__(self, osi, i_node, j_node, num_intgr_pts, sec, transf, int_type):
        self.i_node = int(i_node)
        self.j_node = int(j_node)
        self.num_intgr_pts = int(num_intgr_pts)
        self.sec = sec.tag
        self.transf = transf.tag
        self.int_type = int_type
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.i_node, self.j_node, self.num_intgr_pts, self.sec.tag, self.transf.tag, '-integration', self.int_type]
        self.to_process(osi)


class DispBeamColumnIntmass(ElementBase):

    def __init__(self, osi, ele_nodes, num_intgr_pts, sec, transf, c_rot, mass_dens):
        self.ele_nodes = ele_nodes
        self.num_intgr_pts = int(num_intgr_pts)
        self.sec = sec.tag
        self.transf = transf.tag
        self.c_rot = float(c_rot)
        self.mass_dens = float(mass_dens)
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.num_intgr_pts, self.sec.tag, self.transf.tag, self.c_rot, '-mass', self.mass_dens]
        self.to_process(osi)


class MVLEMthick(ElementBase):

    def __init__(self, osi, dens, ele_nodes, m, c, thicknesses):
        self.dens = float(dens)
        self.ele_nodes = ele_nodes
        self.m = int(m)
        self.c = float(c)
        self.thicknesses = thicknesses
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.dens, *self.ele_nodes, self.m, self.c, '-thick', *self.thicknesses]
        self.to_process(osi)

class MVLEMwidth(ElementBase):

    def __init__(self, osi, dens, ele_nodes, m, c, widths):
        self.dens = float(dens)
        self.ele_nodes = ele_nodes
        self.m = int(m)
        self.c = float(c)
        self.widths = widths
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.dens, *self.ele_nodes, self.m, self.c, '-width', *self.widths]
        self.to_process(osi)

class MVLEMrho(ElementBase):

    def __init__(self, osi, dens, ele_nodes, m, c, reinforcing_ratios):
        self.dens = float(dens)
        self.ele_nodes = ele_nodes
        self.m = int(m)
        self.c = float(c)
        self.reinforcing_ratios = reinforcing_ratios
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.dens, *self.ele_nodes, self.m, self.c, '-rho', *self.reinforcing_ratios]
        self.to_process(osi)

class MVLEMmatConcrete(ElementBase):

    def __init__(self, osi, dens, ele_nodes, m, c, concrete_tags):
        self.dens = float(dens)
        self.ele_nodes = ele_nodes
        self.m = int(m)
        self.c = float(c)
        self.concrete_tags = concrete_tags
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.dens, *self.ele_nodes, self.m, self.c, '-matConcrete', *self.concrete_tags]
        self.to_process(osi)

class MVLEMmatSteel(ElementBase):

    def __init__(self, osi, dens, ele_nodes, m, c, steel_tags):
        self.dens = float(dens)
        self.ele_nodes = ele_nodes
        self.m = int(m)
        self.c = float(c)
        self.steel_tags = steel_tags
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.dens, *self.ele_nodes, self.m, self.c, '-matSteel', *self.steel_tags]
        self.to_process(osi)

class MVLEMmatShear(ElementBase):

    def __init__(self, osi, dens, ele_nodes, m, c, shear):
        self.dens = float(dens)
        self.ele_nodes = ele_nodes
        self.m = int(m)
        self.c = float(c)
        self.shear = shear.tag
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, self.dens, *self.ele_nodes, self.m, self.c, '-matShear', self.shear.tag]
        self.to_process(osi)


class Thick(ElementBase):

    def __init__(self, osi, ele_nodes, m, c, thicknesses):
        self.ele_nodes = ele_nodes
        self.m = int(m)
        self.c = float(c)
        self.thicknesses = thicknesses
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.m, self.c, '-thick', *self.thicknesses]
        self.to_process(osi)

class Width(ElementBase):

    def __init__(self, osi, ele_nodes, m, c, widths):
        self.ele_nodes = ele_nodes
        self.m = int(m)
        self.c = float(c)
        self.widths = widths
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.m, self.c, '-width', *self.widths]
        self.to_process(osi)

class Mat(ElementBase):

    def __init__(self, osi, ele_nodes, m, c, material_tags):
        self.ele_nodes = ele_nodes
        self.m = int(m)
        self.c = float(c)
        self.material_tags = material_tags
        osi.n_ele += 1
        self._tag = osi.n_ele
        self._parameters = [self.op_type, self._tag, *self.ele_nodes, self.m, self.c, '-mat', *self.material_tags]
        self.to_process(osi)
