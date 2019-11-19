from o3seespy.command.nd_material import NDMaterialBase
from o3seespy.base_model import OpenseesObject


class SectionBase(OpenseesObject):
    op_base_type = "section"


class Elastic2D(SectionBase):
    op_type = 'Elastic'

    def __init__(self, osi, big_e, big_a, iz, big_g=0.0, alpha_y=0.0):
        self.big_e = float(big_e)
        self.big_a = float(big_a)
        self.iz = float(iz)
        self.big_g = float(big_g)
        self.alpha_y = float(alpha_y)
        osi.n_sect += 1
        self._tag = osi.n_sect
        self._parameters = [self.op_type, self._tag, self.big_e, self.big_a, self.iz, self.big_g, self.alpha_y]
        self.to_process(osi)


class Elastic3D(SectionBase):
    op_type = 'Elastic'

    def __init__(self, osi, big_e, big_a, iz, iy, big_g, big_j, alpha_y=0.0, alpha_z=0.0):
        self.big_e = float(big_e)
        self.big_a = float(big_a)
        self.iz = float(iz)
        self.iy = float(iy)
        self.big_g = float(big_g)
        self.big_j = float(big_j)
        self.alpha_y = float(alpha_y)
        self.alpha_z = float(alpha_z)
        osi.n_sect += 1
        self._tag = osi.n_sect
        self._parameters = [self.op_type, self._tag, self.big_e, self.big_a, self.iz, self.iy, self.big_g, self.big_j, self.alpha_y, self.alpha_z]
        self.to_process(osi)


class Fiber(SectionBase):
    op_type = 'Fiber'

    def __init__(self, osi, gj=None):
        if gj is None:
            self.gj = None
        else:
            self.gj = float(gj)
        osi.n_sect += 1
        self._tag = osi.n_sect
        self._parameters = [self.op_type, self._tag]
        if getattr(self, 'gj') is not None:
            self._parameters += ['-GJ', self.gj]
        self.to_process(osi)

class FiberThermal(SectionBase):
    op_type = 'FiberThermal'

    def __init__(self, osi, gj=None):
        self.gj = gj
        osi.n_sect += 1
        self._tag = osi.n_sect
        self._parameters = [self.op_type, self._tag]
        if getattr(self, 'gj') is not None:
            self._parameters += ['-GJ', self.gj]
        self.to_process(osi)


class NDFiber(SectionBase):
    op_type = 'NDFiber'

    def __init__(self, osi):
        osi.n_sect += 1
        self._tag = osi.n_sect
        self._parameters = [self.op_type, self._tag]
        self.to_process(osi)


class WFSection2D(SectionBase):
    op_type = 'WFSection2d'

    def __init__(self, osi, mat, d, tw, bf, tf, nfw, nff):
        self.mat = mat
        self.d = float(d)
        self.tw = float(tw)
        self.bf = float(bf)
        self.tf = float(tf)
        self.nfw = float(nfw)
        self.nff = float(nff)
        osi.n_sect += 1
        self._tag = osi.n_sect
        self._parameters = [self.op_type, self._tag, self.mat.tag, self.d, self.tw, self.bf, self.tf, self.nfw, self.nff]
        self.to_process(osi)


class RCSection2D(SectionBase):
    op_type = 'RCSection2d'

    def __init__(self, osi, core, cover, steel, d, b, cover_depth, atop, abot, aside, nfcore, nfcover, nfs):
        self.core = core
        self.cover = cover
        self.steel = steel
        self.d = float(d)
        self.b = float(b)
        self.cover_depth = float(cover_depth)
        self.atop = float(atop)
        self.abot = float(abot)
        self.aside = float(aside)
        self.nfcore = float(nfcore)
        self.nfcover = float(nfcover)
        self.nfs = float(nfs)
        osi.n_sect += 1
        self._tag = osi.n_sect
        self._parameters = [self.op_type, self._tag, self.core.tag, self.cover.tag, self.steel.tag, self.d, self.b, self.cover_depth, self.atop, self.abot, self.aside, self.nfcore, self.nfcover, self.nfs]
        self.to_process(osi)


class RCCircularSection(SectionBase):
    op_type = 'RCCircularSection'

    def __init__(self, osi, core, cover, steel, d, cover_depth, a_s, nrings_core, nrings_cover, newedges, nsteel):
        self.core = core
        self.cover = cover
        self.steel = steel
        self.d = float(d)
        self.cover_depth = float(cover_depth)
        self.a_s = float(a_s)
        self.nrings_core = int(nrings_core)
        self.nrings_cover = int(nrings_cover)
        self.newedges = int(newedges)
        self.nsteel = int(nsteel)
        osi.n_sect += 1
        self._tag = osi.n_sect
        self._parameters = [self.op_type, self._tag, self.core.tag, self.cover.tag, self.steel.tag, self.d, self.cover_depth, self.a_s, self.nrings_core, self.nrings_cover, self.newedges, self.nsteel]
        self.to_process(osi)


class Parallel(SectionBase):
    op_type = 'Parallel'

    def __init__(self, osi, tags):
        self.tags = tags
        osi.n_sect += 1
        self._tag = osi.n_sect
        self._parameters = [self.op_type, self._tag, *self.tags]
        self.to_process(osi)


class Aggregator(SectionBase):
    op_type = 'Aggregator'

    def __init__(self, osi, mats, section=None):
        self.mats = mats
        self.section = section
        osi.n_sect += 1
        self._tag = osi.n_sect
        self._parameters = [self.op_type, self._tag, *self.mats]
        if getattr(self, 'section') is not None:
            self._parameters += ['-section', self.section.tag]
        self.to_process(osi)


class Uniaxial(SectionBase):
    op_type = 'Uniaxial'

    def __init__(self, osi, mat, quantity):
        self.mat = mat
        self.quantity = quantity
        osi.n_sect += 1
        self._tag = osi.n_sect
        self._parameters = [self.op_type, self._tag, self.mat.tag, self.quantity]
        self.to_process(osi)


class ElasticMembranePlateSection(SectionBase):
    op_type = 'ElasticMembranePlateSection'

    def __init__(self, osi, big_e, nu, h, rho):
        self.big_e = float(big_e)
        self.nu = float(nu)
        self.h = float(h)
        self.rho = float(rho)
        osi.n_sect += 1
        self._tag = osi.n_sect
        self._parameters = [self.op_type, self._tag, self.big_e, self.nu, self.h, self.rho]
        self.to_process(osi)


class PlateFiber(SectionBase):
    op_type = 'PlateFiber'

    def __init__(self, osi, mat, h):
        self.mat = mat
        self.h = float(h)
        osi.n_sect += 1
        self._tag = osi.n_sect
        self._parameters = [self.op_type, self._tag, self.mat.tag, self.h]
        self.to_process(osi)


class Bidirectional(SectionBase):
    op_type = 'Bidirectional'

    def __init__(self, osi, big_e, fy, hiso, hkin, code1='Vy', code2='P'):
        self.big_e = float(big_e)
        self.fy = float(fy)
        self.hiso = float(hiso)
        self.hkin = float(hkin)
        self.code1 = code1
        self.code2 = code2
        osi.n_sect += 1
        self._tag = osi.n_sect
        self._parameters = [self.op_type, self._tag, self.big_e, self.fy, self.hiso, self.hkin, self.code1, self.code2]
        self.to_process(osi)


class Iso2spring(SectionBase):
    op_type = 'Iso2spring'

    def __init__(self, osi, tol, k1, fyo, k2o, kvo, hb, pe, po=0.0):
        self.tol = float(tol)
        self.k1 = float(k1)
        self.fyo = float(fyo)
        self.k2o = float(k2o)
        self.kvo = float(kvo)
        self.hb = float(hb)
        self.pe = float(pe)
        self.po = float(po)
        osi.n_sect += 1
        self._tag = osi.n_sect
        self._parameters = [self.op_type, self._tag, self.tol, self.k1, self.fyo, self.k2o, self.kvo, self.hb, self.pe, self.po]
        self.to_process(osi)


class LayeredShell(NDMaterialBase):
    op_type = 'LayeredShell'

    def __init__(self, osi, n_layers, mats):
        self.n_layers = int(n_layers)
        self.mats = mats
        osi.n_sect += 1
        self._tag = osi.n_sect
        self._parameters = [self.op_type, self._tag, self.n_layers, *self.mats]
        self.to_process(osi)
