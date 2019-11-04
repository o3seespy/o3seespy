from o3seespy.command.uniaxial_material.base_material import UniaxialMaterialBase



class PySimple1(UniaxialMaterialBase):

    def __init__(self, osi, soil_type, pult, y50, cd, c=0.0):
        self.soil_type = int(soil_type)
        self.pult = float(pult)
        self.y50 = float(y50)
        self.cd = float(cd)
        self.c = float(c)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.soil_type, self.pult, self.y50, self.cd, self.c]
        self.to_process(osi)


class TzSimple1(UniaxialMaterialBase):

    def __init__(self, osi, soil_type, tult, z50, c=0.0):
        self.soil_type = int(soil_type)
        self.tult = float(tult)
        self.z50 = float(z50)
        self.c = float(c)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.soil_type, self.tult, self.z50, self.c]
        self.to_process(osi)


class QzSimple1(UniaxialMaterialBase):

    def __init__(self, osi, qz_type, qult, z50, suction=0.0, c=0.0):
        self.qz_type = int(qz_type)
        self.qult = float(qult)
        self.z50 = float(z50)
        self.suction = float(suction)
        self.c = float(c)
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.qz_type, self.qult, self.z50, self.suction, self.c]
        self.to_process(osi)


class PyLiq1(UniaxialMaterialBase):

    def __init__(self, osi, soil_type, pult, y50, cd, c, p_res, ele1, ele2, time_series=None):
        self.soil_type = int(soil_type)
        self.pult = float(pult)
        self.y50 = float(y50)
        self.cd = float(cd)
        self.c = float(c)
        self.p_res = float(p_res)
        self.ele1 = float(ele1)
        self.ele2 = float(ele2)
        self.time_series = time_series.tag
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.soil_type, self.pult, self.y50, self.cd, self.c, self.p_res, self.ele1, self.ele2]
        if getattr(self, 'time_series') is not None:
            self._parameters += ['-timeSeries', self.time_series]
        self.to_process(osi)


class TzLiq1(UniaxialMaterialBase):

    def __init__(self, osi, tz_type, tult, z50, c, ele1, ele2, time_series=None):
        self.tz_type = int(tz_type)
        self.tult = float(tult)
        self.z50 = float(z50)
        self.c = float(c)
        self.ele1 = float(ele1)
        self.ele2 = float(ele2)
        self.time_series = time_series.tag
        osi.n_mat += 1
        self._tag = osi.n_mat
        self._parameters = [self.op_type, self._tag, self.tz_type, self.tult, self.z50, self.c, self.ele1, self.ele2]
        if getattr(self, 'time_series') is not None:
            self._parameters += ['-timeSeries', self.time_series]
        self.to_process(osi)
