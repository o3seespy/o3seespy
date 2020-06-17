from o3seespy.base_model import OpenSeesObject
import tempfile
import os


class RecorderBase(OpenSeesObject):
    op_base_type = "recorder"
    
    
class RecorderToArrayCacheBase(RecorderBase):  # TODO: implement NodeToArray where data saved to memory and loaded as array without collect
    fname = None
    
    def collect(self, unlink=True):
        from numpy import loadtxt
        try:
            a = loadtxt(self.fname, dtype=float)
        except ValueError as e:
            print('Warning: Need to run opy.wipe() before collecting arrays')
            raise ValueError(e)
        if unlink:
            try:
                os.unlink(self.fname)
            except PermissionError:
                print('Warning: Need to run opy.wipe() before collecting arrays')
        return a


class NodeToFile(RecorderBase):
    op_type = "Node"

    def __init__(self, osi, fname, node, dofs, res_type, nsd=8, dt=None, time=False):
        """
        Records properties of a node and saves the results to a file

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        fname: str
            Full file name
        node: o3seespy.node.Node
        dofs: list
            A list of integers representing the degrees-of-freedom
        res_type: str
            Response type
        nsd: int
            Number of significant figures
        dt: float
            Time step
        time: bool
            If true the first column is the time
        """
        self.osi = osi
        self._parameters = [self.op_type, '-file', fname, '-precision', nsd, '-node', node.tag]
        if dt is not None:
            self._parameters.insert(5, '-dT')
            self._parameters.insert(6, dt)
        if time:
            self._parameters.insert(5, '-time')
        self._parameters += ['-dof', *dofs, res_type]
        self.to_process(osi)


class NodesToFile(RecorderBase):
    op_type = "Node"

    def __init__(self, osi, fname, nodes, dofs, res_type, nsd=8, dt=None, time=False):
        """
        Records properties of several nodes and saves the results to a file

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        fname: str
            Full file name
        node: list
            List of o3seespy.node.Node objects
        dofs: list
            A list of integers representing the degrees-of-freedom
        res_type: str
            Response type
        nsd: int
            Number of significant figures
        dt: float
            Time step
        time: bool
            If true the first column is the time
        """
        self.osi = osi
        if isinstance(nodes, str) and nodes == 'all':
            node_tags = osi.to_process('getNodeTags', [])
        else:
            node_tags = [x.tag for x in nodes]
        self._parameters = [self.op_type, '-file', fname, '-precision', nsd, '-node', *node_tags, '-dof', *dofs, res_type]
        if dt is not None:
            self._parameters.insert(5, '-dT')
            self._parameters.insert(6, dt)
        if time:
            self._parameters.insert(5, '-time')
        self.to_process(osi)


class NodeToArrayCache(RecorderToArrayCacheBase):  # TODO: implement NodeToArray where data saved to memory and loaded as array without collect
    op_type = "Node"

    def __init__(self, osi, node, dofs, res_type, nsd=8, dt=None, fname=None):
        """
        Records properties of a node and saves results to a numpy array

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        node: o3seespy.node.Node
        dofs: list
            A list of integers representing the degrees-of-freedom
        res_type: str
            Response type: `disp`, 'vel', 'accel', 'incrDisp', 'reaction', 'eigin i', 'rayleighForces'
        nsd: int
            Number of significant figures
        dt: float
            Time step
        fname: str
            Full file path where data should be stored and loaded from
        """
        self.osi = osi
        if fname is None:
            self.fname = tempfile.NamedTemporaryFile(delete=False).name
        else:
            self.fname = fname
        self._parameters = [self.op_type, '-file', self.fname, '-precision', nsd, '-node', node.tag]
        if dt is not None:
            self._parameters.insert(5, '-dT')
            self._parameters.insert(6, dt)
        self._parameters += ['-dof', *dofs, res_type]
        self.to_process(osi)


class NodesToArrayCache(RecorderToArrayCacheBase):  # TODO: implement NodeToArray where data saved to memory and loaded as array without collect
    op_type = "Node"

    def __init__(self, osi, nodes, dofs, res_type, nsd=8, dt=None, fname=None):
        """
        Records properties of several nodes and saves results to a numpy array

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        nodes: list
            A list of o3seespy.node.Node objects
        dofs: list
           A list of integers representing the degrees-of-freedom
        res_type: str
           Response type
        nsd: int
           Number of significant figures
        dt: float
           Time step
        """
        self.osi = osi
        if isinstance(nodes, str) and nodes == 'all':
            node_tags = osi.to_process('getNodeTags', [])
        else:
            node_tags = [x.tag for x in nodes]
        if fname is None:
            self.fname = tempfile.NamedTemporaryFile(delete=False).name
        else:
            self.fname = fname
        self._parameters = [self.op_type, '-file', self.fname, '-precision', nsd, '-node', *node_tags, '-dof', *dofs, res_type]
        if dt is not None:
            self._parameters.insert(5,'-dT')
            self._parameters.insert(6, dt)
        self.to_process(osi)


class TimeToArrayCache(RecorderBase):
    op_type = "Node"

    def __init__(self, osi, nsd=8, dt=None, fname=None):
        """
        Records the recorder time and saves to a numpy array

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        nsd: int
            Number of significant figures
        dt: float
            Time step
        """
        self.osi = osi
        if fname is None:
            self.fname = tempfile.NamedTemporaryFile(delete=False).name
        else:
            self.fname = fname
        self._parameters = [self.op_type, '-file', self.fname, '-precision', nsd, '-time', '-node', 1]
        if dt is not None:
            self._parameters.insert(5, '-dT')
            self._parameters.insert(6, dt)
        self._parameters += ['-dof', 1, 'accel']
        self.to_process(osi)

    def collect(self, unlink=True):
        from numpy import loadtxt
        try:
            a = loadtxt(self.fname, dtype=float)
        except ValueError as e:
            print('Warning: Need to run opy.wipe() before collecting arrays')
            raise ValueError(e)
        if unlink:
            try:
                os.unlink(self.fname)
            except PermissionError:
                print('Warning: Need to run opy.wipe() before collecting arrays')
        return a[:, 0]


class TimeToFile(RecorderBase):
    op_type = "Node"

    def __init__(self, osi, fname, nsd=8, dt=None):
        """
        Records the recorder time and saves to a numpy array

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        nsd: int
            Number of significant figures
        dt: float
            Time step
        """
        self.osi = osi
        self.fname = tempfile.NamedTemporaryFile(delete=False).name
        self._parameters = [self.op_type, '-file', fname, '-precision', nsd, '-time', '-node', 1]
        if dt is not None:
            self._parameters.insert(5, '-dT')
            self._parameters.insert(6, dt)
        self._parameters += ['-dof', 1, 'disp']
        self.to_process(osi)


class ElementToFile(RecorderBase):
    op_type = "Element"

    def __init__(self, osi, fname, ele, material=None, arg_vals=None, nsd=8, dt=None, time=False):
        """
        Records properties of an element and saves the results to a file

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        fname: str
            Full file name
        ele: o3seespy.element.BaseElement
            An o3seespy element
        material: -
        arg_vals: list
            Extra arguments passed to element recorder method
        nsd: int
            Number of significant figures
        dt: float
            Time step
        time: bool
            If true the first column is the time
        """
        self.osi = osi
        if arg_vals is None:
            arg_vals = []
        extra_pms = []
        if material is not None:
            extra_pms += ['material', material]
        self._parameters = [self.op_type, '-file', fname, '-precision', nsd, '-ele', ele.tag, *extra_pms, *arg_vals]
        if dt is not None:
            self._parameters.insert(5, '-dT')
            self._parameters.insert(6, dt)
        if time:
            self._parameters.insert(5, '-time')
        self.to_process(osi)


class ElementsToFile(RecorderBase):
    op_type = "Element"

    def __init__(self, osi, fname, eles, material=None, arg_vals=None, nsd=8, dt=None, time=False):
        """
        Records properties of an element and saves the results to a file

        Parameters
        ----------
        osi: o3seespy.OpenSeesInstance
        fname: str
            Full file name
        eles: list of o3seespy.element.BaseElement
            List of o3seespy elements
        material: -
        arg_vals: list
            Extra arguments passed to element recorder method
        nsd: int
            Number of significant figures
        dt: float
            Time step
        time: bool
            If true the first column is the time
        """
        self.osi = osi
        if arg_vals is None:
            arg_vals = []
        extra_pms = []
        if material is not None:
            extra_pms += ['material', material]
        self.ele_tags = [x.tag for x in eles]
        self._parameters = [self.op_type, '-file', fname, '-precision', nsd, '-ele', *self.ele_tags, *extra_pms, *arg_vals]
        if dt is not None:
            self._parameters.insert(5, '-dT')
            self._parameters.insert(6, dt)
        if time:
            self._parameters.insert(5, '-time')
        self.to_process(osi)


class ElementToArrayCache(RecorderToArrayCacheBase):  # TODO: implement ElementToArray where data saved to memory and loaded as array without collect
    op_type = "Element"

    def __init__(self, osi, ele, material=None, arg_vals=None, nsd=8, fname=None, dt=None):
        self.osi = osi
        if arg_vals is None:
            arg_vals = []
        self.arg_vals = [str(x) for x in arg_vals]
        extra_pms = []
        if material is not None:
            extra_pms += ['material', material]
        if fname is None:
            self.fname = tempfile.NamedTemporaryFile(delete=False).name
        else:
            self.fname = fname
        self.ele = ele
        self._parameters = [self.op_type, '-file', self.fname, '-precision', nsd, '-ele', ele.tag, *extra_pms, *self.arg_vals]
        if dt is not None:
            self._parameters.insert(5, '-dT')
            self._parameters.insert(6, dt)
        self.to_process(osi)


class ElementsToArrayCache(RecorderToArrayCacheBase):
    op_type = "Element"

    def __init__(self, osi, eles, material=None, arg_vals=None, nsd=8, fname=None, dt=None):
        self.osi = osi
        if arg_vals is None:
            arg_vals = []
        extra_pms = []
        if material is not None:
            extra_pms += ['material', material]
        if fname is None:
            self.fname = tempfile.NamedTemporaryFile(delete=False).name
        else:
            self.fname = fname
        self.ele_tags = [x.tag for x in eles]

        self._parameters = [self.op_type, '-file', self.fname, '-precision', nsd, '-ele', *self.ele_tags, *extra_pms, *arg_vals]
        if dt is not None:
            self._parameters.insert(5, '-dT')
            self._parameters.insert(6, dt)
        self.to_process(osi)



def load_recorder_options():
    folder_path = os.path.dirname(os.path.realpath(__file__))
    return open(folder_path + '/mat_recorder_options.csv')
