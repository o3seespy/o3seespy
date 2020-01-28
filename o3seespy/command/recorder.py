from o3seespy.base_model import OpenseesObject
import tempfile
import os


class RecorderBase(OpenseesObject):
    op_base_type = "recorder"


class NodeToFile(RecorderBase):
    op_type = "Node"

    def __init__(self, osi, fname, node, dofs, res_type, nsd=8, dt=None):
        self._parameters = [self.op_type, '-file', fname, '-precision', nsd, '-node', node.tag]
        if dt is not None:
            self._parameters += ['-dT', dt]
        self._parameters += ['-dof', *dofs, res_type]
        self.to_process(osi)


class NodesToFile(RecorderBase):
    op_type = "Node"

    def __init__(self, osi, fname, nodes, dofs, res_type, nsd=8):
        node_tags = [x.tag for x in nodes]
        self._parameters = [self.op_type, '-file', fname, '-precision', nsd, '-node', *node_tags, '-dof', *dofs, res_type]
        self.to_process(osi)


class NodeToArrayCache(RecorderBase):  # TODO: implement NodeToArray where data saved to memory and loaded as array without collect
    op_type = "Node"

    def __init__(self, osi, node, dofs, res_type, nsd=8, dt=None):
        self.tmpfname = tempfile.NamedTemporaryFile(delete=False).name
        self._parameters = [self.op_type, '-file', self.tmpfname, '-precision', nsd, '-node', node.tag]
        if dt is not None:
            self._parameters += ['-dT', dt]
        self._parameters += ['-dof', *dofs, res_type]
        self.to_process(osi)

    def collect(self):
        from numpy import loadtxt
        try:
            a = loadtxt(self.tmpfname, dtype=float)
        except ValueError as e:
            print('Warning: Need to run opy.wipe() before collecting arrays')
            raise ValueError(e)
        try:
            os.unlink(self.tmpfname)
        except PermissionError:
            print('Warning: Need to run opy.wipe() before collecting arrays')
        return a


class NodesToArrayCache(RecorderBase):  # TODO: implement NodeToArray where data saved to memory and loaded as array without collect
    op_type = "Node"

    def __init__(self, osi, nodes, dofs, res_type, nsd=8):
        node_tags = [x.tag for x in nodes]
        self.tmpfname = tempfile.NamedTemporaryFile(delete=False).name
        self._parameters = [self.op_type, '-file', self.tmpfname, '-precision', nsd, '-node', *node_tags, '-dof', *dofs, res_type]
        self.to_process(osi)

    def collect(self):
        from numpy import loadtxt
        try:
            a = loadtxt(self.tmpfname, dtype=float)
        except ValueError as e:
            print('Warning: Need to run opy.wipe() before collecting arrays')
            raise ValueError(e)
        try:
            os.unlink(self.tmpfname)
        except PermissionError:
            print('Warning: Need to run opy.wipe() before collecting arrays')
        return a


class ElementToFile(RecorderBase):
    op_type = "Element"

    def __init__(self, osi, fname, element, material=None, args=None, nsd=8):
        if args is None:
            args = []
        extra_pms = []
        if material is not None:
            extra_pms += ['material', material]
        self._parameters = [self.op_type, '-file', fname, '-precision', nsd, '-ele', element.tag, *extra_pms, *args]
        self.to_process(osi)


class ElementToArrayCache(RecorderBase):  # TODO: implement ElementToArray where data saved to memory and loaded as array without collect
    op_type = "Element"

    def __init__(self, osi, element, material=None, arg_vals=None, nsd=8, fname=None):
        if arg_vals is None:
            arg_vals = []
        extra_pms = []
        if material is not None:
            extra_pms += ['material', material]
        if fname is None:
            self.tmpfname = tempfile.NamedTemporaryFile(delete=False).name
        else:
            self.tmpfname = fname

        self._parameters = [self.op_type, '-file', self.tmpfname, '-precision', nsd, '-ele', element.tag, *extra_pms, *arg_vals]
        self.to_process(osi)

    def collect(self):
        from numpy import loadtxt
        try:
            a = loadtxt(self.tmpfname, dtype=float)
        except ValueError as e:
            print('Warning: Need to run opy.wipe() before collecting arrays')
            raise ValueError(e)
        try:
            os.unlink(self.tmpfname)
        except PermissionError:
            print('Warning: Need to run opy.wipe() before collecting arrays')
        return a
