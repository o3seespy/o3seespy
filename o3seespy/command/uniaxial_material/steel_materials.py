from o3seespy.opensees_instance import OpenseesInstance
from o3seespy.command.uniaxial_material.base_material import UniaxialMaterialBase


class Steel01(UniaxialMaterialBase):
    op_type = "Steel01"

    def __init__(self, osi: OpenseesInstance, fy: float, e0: float, b: float, a1=None, a2=None, a3=None, a4=None):
        """
        Steel01 object

        Parameters
        ----------
        osi : opensees_pack.opensees_instance.OpenseesInstance object
            An instance of opensees
        fy
        e0
        b
        a1
        a2
        a3
        a4
        """
        self.fy = float(fy)
        self.e0 = float(e0)
        self.b = float(b)
        self.a_values = [a1, a2, a3, a4]
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.a4 = a4
        osi.n_mats += 1
        self._tag = osi.n_mats
        self._parameters = [self.op_type, self.tag, self.fy, self.e0, self.b]
        for a in self.a_values:
            if a is None:
                break
            self._parameters.append(a)
        self.to_process(osi)
