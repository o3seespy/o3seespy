from inspect import signature
from collections import OrderedDict


def to_commands(op_base_type, parameters):
    para = []
    for i, e in enumerate(parameters):
        if isinstance(e, str):
            e = "'%s'" % e
        elif isinstance(e, float):
            e = '%.6g' % e
            if '.' not in e and 'e' not in e:
                e += '.0'
        para.append(str(e))
        if i > 40:  # avoid verbose print output
            break
    p_str = ', '.join(para)
    return 'opy.%s(%s)' % (op_base_type, p_str)


def to_py_file(osi, ofile='ofile.py', w_analyze=False):
    ofile = open(ofile, 'w')
    pstr = 'import openseespy.opensees as opy\n' + '\n'.join(osi.commands)
    if w_analyze:
        pstr += '\nopy.analyze(1, 0.1)\n'
    ofile.write(pstr)
    ofile.close()


def to_tcl_file(osi, ofile='ofile.tcl', w_analyze=False):
    ofile = open(ofile, 'w')
    pstr = '\n'.join(osi.commands)
    if w_analyze:
        pstr += '\nopy.analyze(1, 0.1)\n'
    tcl_str = py2tcl(pstr)
    ofile.write(tcl_str)
    ofile.close()


def get_o3_kwargs_from_obj(obj, o3_obj, custom=None, overrides=None):
    if custom is None:
        custom = {}
    if overrides is None:
        overrides = {}
    sig = signature(o3_obj)
    kwargs = OrderedDict()
    args = []
    sig_vals = sig.parameters.values()
    for p in sig_vals:
        if p.name in custom:
            pname = custom[p.name]
        else:
            pname = p.name
        if pname == 'osi':
            continue
        if pname in overrides:
            val = overrides[pname]
        else:
            try:
                val = getattr(obj, pname)
            except AttributeError as e:
                if p.default == p.empty:
                    raise AttributeError(e)
                else:
                    val = p.default
        if p.default == p.empty:
            args.append(val)
        else:
            if val is not None:
                kwargs[p.name] = val
    return args, kwargs


def has_o3_model_changed(cur_type, prev_type, cur_args, prev_args, cur_kwargs, prev_kwargs):
    import numpy as np
    changed = 0
    if cur_type != prev_type or len(cur_args) != len(prev_args) or len(cur_kwargs) != len(prev_kwargs):
        changed = 1
    else:
        for j, arg in enumerate(cur_args):
            if hasattr(arg, '__len__'):
                if len(arg) != len(prev_args[j]):
                    changed = 1
                    break
                for k, subarg in enumerate(arg):
                    if not np.isclose(subarg, prev_args[j][k]):
                        changed = 1
            elif not np.isclose(arg, prev_args[j]):
                changed = 1
                break
        for pm in cur_kwargs:
            if pm not in prev_kwargs:
                changed = 1
                break
            if hasattr(cur_kwargs[pm], '__len__'):
                if len(cur_kwargs[pm]) != len(prev_kwargs[pm]):
                    changed = 1
                    break
                for k, subarg in enumerate(cur_kwargs[pm]):
                    if not np.isclose(subarg, cur_kwargs[pm][k]):
                        changed = 1
                        break
            elif not np.isclose(cur_kwargs[pm], prev_kwargs[pm]):
                changed = 1
                break
    return changed


def py2tcl(pystr):
    """
    Converts openseespy script to tcl

    Returns
    -------

    """
    # new = '\n'.join(pystr.split()[1:])  # removes the import line
    new = pystr.replace('(', ' ')
    new = new.replace(')', ' ')
    new = new.replace('opy.', '')
    new = new.replace(',', '')
    new = new.replace("'", '')
    new = new.replace('"', '')
    return new
