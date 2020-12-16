
def parameter(osi, p_args):
    """
    In DDM-based FE response sensitivity analysis, the sensitivity parameters can be material,geometry or discrete
    loading parameters. 

    Parameters
    ----------
    osi: o3seespy.OpenSeesInstance
    p_args: unk
            Depend on the object in the fe model encapsulating the desired parameters.
    """
    _parameters = [p_args]
    return osi.to_process("parameter", _parameters)


def add_to_parameter(osi, p_args):
    """
    In case that more objects (e.g., element, section) are mapped to an existing parameter,the  command can be used to
    relate these additional objects to the specific parameter.

    Parameters
    ----------
    osi: o3seespy.OpenSeesInstance
    p_args: unk
            Depend on the object in the fe model encapsulating the desired parameters.
    """
    _parameters = [p_args]
    return osi.to_process("addToParameter", _parameters)


def update_parameter(osi, new_value):
    """
    Once the parameters in FE model are defined, their value can be updated.

    Parameters
    ----------
    osi: o3seespy.OpenSeesInstance
    new_value: float
            The updated value to which the parameter needs to be set.

    Examples
    --------
    >>> import o3seespy as o3
    >>> # Example is currently not working
    >>> osi = o3.OpenSeesInstance(ndm=2)
    >>> o3.senscmds.update_parameter(osi, new_value=1.0)
    """
    new_value = float(new_value)
    _parameters = [new_value]
    return osi.to_process("updateParameter", _parameters)


def set_parameter(osi, new_value, start, end, args, eles: list=None):
    """
    set value for an element parameter

    Parameters
    ----------
    osi: o3seespy.OpenSeesInstance
    new_value: float
            The updated value to which the parameter needs to be set.
    start: int
            Start element object 
    end: int
            End element object 
    args: lists
            A list of strings for the element parameter
    eles: list, optional
            A list of element objects

    Examples
    --------
    >>> import o3seespy as o3
    >>> # Example is currently not working
    >>> osi = o3.OpenSeesInstance(ndm=2)
    >>> eles = [1, 1]
    >>> o3.senscmds.set_parameter(osi, new_value=1.0, eles=eles, start=1, end=1, args=1)
    """
    eles = [x.tag for x in eles]
    start = int(start)
    end = int(end)
    _parameters = ['-eleRange', start, end, *args]
    if eles is not None:
        _parameters += ['-ele', *eles]
    return osi.to_process("setParameter", _parameters)


def get_param_tags(osi):
    """
    Return a list of tags for all parameters

    Parameters
    ----------
    osi: o3seespy.OpenSeesInstance

    Examples
    --------
    >>> import o3seespy as o3
    >>> osi = o3.OpenSeesInstance(ndm=2)
    >>> o3.senscmds.get_param_tags(osi)
    """
    _parameters = []
    return osi.to_process("getParamTags", _parameters)


def get_param_value(osi):
    """
    Return the value of a parameter

    Parameters
    ----------
    osi: o3seespy.OpenSeesInstance

    Examples
    --------
    >>> import o3seespy as o3
    >>> # Example is currently not working
    >>> osi = o3.OpenSeesInstance(ndm=2)
    >>> o3.senscmds.get_param_value(osi)
    """
    _parameters = []
    return osi.to_process("getParamValue", _parameters)


def compute_gradients(osi):
    """
    This command is used to perform a sensitivity analysis. If the user wants

    Parameters
    ----------
    osi: o3seespy.OpenSeesInstance

    Examples
    --------
    >>> import o3seespy as o3
    >>> # Example is currently not working
    >>> osi = o3.OpenSeesInstance(ndm=2)
    >>> o3.senscmds.compute_gradients(osi)
    """
    _parameters = []
    return osi.to_process("computeGradients", _parameters)


def sensitivity_algorithm(osi):
    """
    This command is used to create a sensitivity algorithm.

    Parameters
    ----------
    osi: o3seespy.OpenSeesInstance

    Examples
    --------
    >>> import o3seespy as o3
    >>> # Example is currently not working
    >>> osi = o3.OpenSeesInstance(ndm=2)
    >>> o3.senscmds.sensitivity_algorithm(osi)
    """
    _parameters = []
    return osi.to_process("sensitivityAlgorithm", _parameters)


def get_sens_node_disp(osi, dof, param):
    """
    Returns the current displacement sensitivity to a parameter at a specified node.

    Parameters
    ----------
    osi: o3seespy.OpenSeesInstance
    dof: int
            Specific dof at the node (1 through ndf)
    param: obj
            Parameter object

    Examples
    --------
    >>> import o3seespy as o3
    >>> # Example is currently not working
    >>> osi = o3.OpenSeesInstance(ndm=2)
    >>> o3.senscmds.get_sens_node_disp(osi, dof=1, param=[])
    """
    dof = int(dof)
    _parameters = [dof, param.tag]
    return osi.to_process("sensNodeDisp", _parameters)


def get_sens_node_vel(osi, dof, param):
    """
    Returns the current velocity sensitivity to a parameter at a specified node.

    Parameters
    ----------
    osi: o3seespy.OpenSeesInstance
    dof: int
            Specific dof at the node (1 through ndf)
    param: obj
            Parameter object

    Examples
    --------
    >>> import o3seespy as o3
    >>> # Example is currently not working
    >>> osi = o3.OpenSeesInstance(ndm=2)
    >>> o3.senscmds.get_sens_node_vel(osi, dof=1, param=[])
    """
    dof = int(dof)
    _parameters = [dof, param.tag]
    return osi.to_process("sensNodeVel", _parameters)


def get_sens_node_accel(osi, dof, param):
    """
    Returns the current acceleration sensitivity to a parameter at a specified node.

    Parameters
    ----------
    osi: o3seespy.OpenSeesInstance
    dof: int
            Specific dof at the node (1 through ndf)
    param: obj
            Parameter object

    Examples
    --------
    >>> import o3seespy as o3
    >>> # Example is currently not working
    >>> osi = o3.OpenSeesInstance(ndm=2)
    >>> o3.senscmds.get_sens_node_accel(osi, dof=1, param=[])
    """
    dof = int(dof)
    _parameters = [dof, param.tag]
    return osi.to_process("sensNodeAccel", _parameters)


def get_sens_lambda(osi, param):
    """
    Returns the current load factor sensitivity to a parameter in a load pattern.

    Parameters
    ----------
    osi: o3seespy.OpenSeesInstance

    Examples
    --------
    >>> import o3seespy as o3
    >>> # Example is currently not working
    >>> osi = o3.OpenSeesInstance(ndm=2)
    >>> o3.senscmds.get_sens_lambda(osi, param=[])
    """
    _parameters = [param.tag]
    return osi.to_process("sensLambda", _parameters)


def get_sens_section_force(osi, sec_num, dof, param):
    """
    Returns the current section force sensitivity to a parameter at a specified element and section.

    Parameters
    ----------
    osi: o3seespy.OpenSeesInstance
    sec_num: int
            Section number
    dof: int
            Specific dof at the element (1 through element force ndf)
    param: obj
            Parameter object

    Examples
    --------
    >>> import o3seespy as o3
    >>> # Example is currently not working
    >>> osi = o3.OpenSeesInstance(ndm=2)
    >>> o3.senscmds.get_sens_section_force(osi, sec_num=1, dof=1, param=[])
    """
    sec_num = int(sec_num)
    dof = int(dof)
    _parameters = [sec_num, dof, param.tag]
    return osi.to_process("sensSectionForce", _parameters)


def get_sens_node_pressure(osi, param):
    """
    Returns the current pressure sensitivity to a parameter at a specified node.

    Parameters
    ----------
    osi: o3seespy.OpenSeesInstance
    param: obj
            Parameter object

    Examples
    --------
    >>> import o3seespy as o3
    >>> # Example is currently not working
    >>> osi = o3.OpenSeesInstance(ndm=2)
    >>> o3.senscmds.get_sens_node_pressure(osi, param=[])
    """
    _parameters = [param.tag]
    return osi.to_process("sensNodePressure", _parameters)
