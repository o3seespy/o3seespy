from o3seespy.base_model import OpenSeesObject
import tempfile
import os


class RecorderBase(OpenSeesObject):
    op_base_type = "recorder"
    
    
class RecorderToArrayCacheBase(RecorderBase):  # TODO: implement NodeToArray where data saved to memory and loaded as array without collect
    tmpfname = None
    
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


class NodeToFile(RecorderBase):
    op_type = "Node"

    def __init__(self, osi, fname, node, dofs, res_type, nsd=8, dt=None):
        self._parameters = [self.op_type, '-file', fname, '-precision', nsd, '-node', node.tag]
        if dt is not None:
            self._parameters.insert(5, '-dT')
            self._parameters.insert(6, dt)
        self._parameters += ['-dof', *dofs, res_type]
        self.to_process(osi)


class NodesToFile(RecorderBase):
    op_type = "Node"

    def __init__(self, osi, fname, nodes, dofs, res_type, nsd=8, dt=None):
        if nodes == 'all':
            node_tags = osi.to_process('getNodeTags', [])
        else:
            node_tags = [x.tag for x in nodes]
        self._parameters = [self.op_type, '-file', fname, '-precision', nsd, '-node', *node_tags, '-dof', *dofs, res_type]
        if dt is not None:
            self._parameters.insert(5, '-dT')
            self._parameters.insert(6, dt)
        self.to_process(osi)


class NodeToArrayCache(RecorderToArrayCacheBase):  # TODO: implement NodeToArray where data saved to memory and loaded as array without collect
    op_type = "Node"

    def __init__(self, osi, node, dofs, res_type, nsd=8, dt=None):
        self.tmpfname = tempfile.NamedTemporaryFile(delete=False).name
        self._parameters = [self.op_type, '-file', self.tmpfname, '-precision', nsd, '-node', node.tag]
        if dt is not None:
            self._parameters.insert(5, '-dT')
            self._parameters.insert(6, dt)
        self._parameters += ['-dof', *dofs, res_type]
        self.to_process(osi)


class NodesToArrayCache(RecorderToArrayCacheBase):  # TODO: implement NodeToArray where data saved to memory and loaded as array without collect
    op_type = "Node"

    def __init__(self, osi, nodes, dofs, res_type, nsd=8, dt=None):
        if nodes == 'all':
            node_tags = osi.to_process('getNodeTags', [])
        else:
            node_tags = [x.tag for x in nodes]
        self.tmpfname = tempfile.NamedTemporaryFile(delete=False).name
        self._parameters = [self.op_type, '-file', self.tmpfname, '-precision', nsd, '-node', *node_tags, '-dof', *dofs, res_type]
        if dt is not None:
            self._parameters.insert(5,'-dT')
            self._parameters.insert(6, dt)
        self.to_process(osi)


class ElementToFile(RecorderBase):
    op_type = "Element"

    def __init__(self, osi, fname, element, material=None, args=None, nsd=8, dt=None):
        if args is None:
            args = []
        extra_pms = []
        if material is not None:
            extra_pms += ['material', material]
        self._parameters = [self.op_type, '-file', fname, '-precision', nsd, '-ele', element.tag, *extra_pms, *args]
        if dt is not None:
            self._parameters.insert(5, '-dT')
            self._parameters.insert(6, dt)
        self.to_process(osi)


class ElementToArrayCache(RecorderToArrayCacheBase):  # TODO: implement ElementToArray where data saved to memory and loaded as array without collect
    op_type = "Element"

    def __init__(self, osi, element, material=None, arg_vals=None, nsd=8, fname=None, dt=None):
        if arg_vals is None:
            arg_vals = []
        extra_pms = []
        if material is not None:
            extra_pms += ['material', material]
        if fname is None:
            self.tmpfname = tempfile.NamedTemporaryFile(delete=False).name
        else:
            self.tmpfname = fname
        self.element = element
        self._parameters = [self.op_type, '-file', self.tmpfname, '-precision', nsd, '-ele', element.tag, *extra_pms, *arg_vals]
        if dt is not None:
            self._parameters.insert(5, '-dT')
            self._parameters.insert(6, dt)
        self.to_process(osi)

    # def collect(self):
    #     from numpy import loadtxt
    #     try:
    #         a = loadtxt(self.tmpfname, dtype=float)
    #     except ValueError as e:
    #         print('Warning: Need to run opy.wipe() before collecting arrays')
    #         raise ValueError(e)
    #     try:
    #         os.unlink(self.tmpfname)
    #     except PermissionError:
    #         print('Warning: Need to run opy.wipe() before collecting arrays')
    #     return a


class ElementsToArrayCache(RecorderToArrayCacheBase):
    op_type = "Element"

    def __init__(self, osi, elements, material=None, arg_vals=None, nsd=8, fname=None, dt=None):
        if arg_vals is None:
            arg_vals = []
        extra_pms = []
        if material is not None:
            extra_pms += ['material', material]
        if fname is None:
            self.tmpfname = tempfile.NamedTemporaryFile(delete=False).name
        else:
            self.tmpfname = fname
        self.ele_tags = [x.tag for x in elements]

        self._parameters = [self.op_type, '-file', self.tmpfname, '-precision', nsd, '-ele', *self.ele_tags, *extra_pms, *arg_vals]
        if dt is not None:
            self._parameters.insert(5, '-dT')
            self._parameters.insert(6, dt)
        self.to_process(osi)

    # def collect(self):
    #     from numpy import loadtxt
    #     try:
    #         a = loadtxt(self.tmpfname, dtype=float)
    #     except ValueError as e:
    #         print('Warning: Need to run opy.wipe() before collecting arrays')
    #         raise ValueError(e)
    #     try:
    #         os.unlink(self.tmpfname)
    #     except PermissionError:
    #         print('Warning: Need to run opy.wipe() before collecting arrays')
    #     return a