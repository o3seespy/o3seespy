from o3seespy.base_model import OpenseesObject


class BeamIntegrationBase(OpenseesObject):
    op_base_type = "beamIntegration"


class HingeMidpoint(BeamIntegrationBase):
    op_type = 'HingeMidpoint'

    def __init__(self, osi, sect_i, lp_i, sect_j, lp_j, sect_e):
        self.sect_i = sect_i
        self.lp_i = float(lp_i)
        self.sect_j = sect_j
        self.lp_j = float(lp_j)
        self.sect_e = sect_e

        osi.n_integs += 1
        self._tag = osi.n_integs
        self._parameters = [self.op_type, self._tag, self.sect_i.tag, self.lp_i, self.sect_j.tag,
                            self.lp_j, self.sect_e.tag]
        self.to_process(osi)
